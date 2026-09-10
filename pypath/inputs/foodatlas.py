#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
FoodAtlas data retrieval and processing module.
"""

import os
import json
import zipfile
import tempfile
import urllib.request
import pandas as pd

import pypath.share.session as session

_logger = session.Logger(name='foodatlas_input')
_log = _logger._log

FOODATLAS_URL = "https://foodatlasdownloadsstack-downloadsbucketb54b8c20-hfoegjgwag4w.s3.us-west-1.amazonaws.com/bundles/foodatlas-v4.9/foodatlas-v4.9.zip"

def get_foodatlas_data_dir():
    """
    Downloads and extracts the FoodAtlas archive to a temporary directory.
    Returns the path to the extracted directory containing the parquet files.
    """
    _log(f"Downloading FoodAtlas data from {FOODATLAS_URL}...")
    
    temp_dir = tempfile.mkdtemp(prefix="foodatlas_")
    zip_path = os.path.join(temp_dir, "foodatlas.zip")
    
    urllib.request.urlretrieve(FOODATLAS_URL, zip_path)
    
    _log(f"Extracting FoodAtlas archive to {temp_dir}...")
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)
        
    data_dir = temp_dir
    for root, dirs, files in os.walk(temp_dir):
        if "entities.parquet" in files:
            data_dir = root
            break
            
    return data_dir

def get_foodatlas_entities(file_path):
    """
    Reads entities.parquet and returns specific columns.
    Expected columns: foodatlas_id, entity_type, common_name, synonyms
    """
    columns_to_keep = ['foodatlas_id', 'entity_type', 'common_name', 'synonyms']
    
    df = pd.read_parquet(file_path, columns=columns_to_keep)
    
    return df

def get_foodatlas_relationships(file_path):
    """
    Reads relationships.parquet and returns specific columns.
    Expected columns: foodatlas_id, name
    """
    columns_to_keep = ['foodatlas_id', 'name']
    df = pd.read_parquet(file_path, columns=columns_to_keep)
    return df

def get_foodatlas_attestations(file_path):
    """
    Reads attestations.parquet and returns specific columns.
    Expected columns: attestation_id, evidence_id
    """
    columns_to_keep = ['attestation_id', 'evidence_id']
    df = pd.read_parquet(file_path, columns=columns_to_keep)
    return df

def get_foodatlas_evidence(file_path):
    """
    Reads evidence.parquet and returns specific columns.
    Expected columns: evidence_id, reference
    """
    columns_to_keep = ['evidence_id', 'reference']
    df = pd.read_parquet(file_path, columns=columns_to_keep)
    return df

def get_foodatlas_triplets(file_path):
    """
    Reads triplets.parquet and returns specific columns.
    Expected columns: head_id, relationship_id, tail_id, attestation_ids
    """
    columns_to_keep = ['head_id', 'relationship_id', 'tail_id', 'attestation_ids']
    df = pd.read_parquet(file_path, columns=columns_to_keep)
    return df

def get_foodatlas_interactions(allowed_relationships=None, data_dir=None):
    """
    Main function to read triplets.parquet, filter by allowed_relationships,
    merge with entities and relationships dataframes, and return a flattened format.
    """
    if data_dir is None:
        data_dir = get_foodatlas_data_dir()
        
    _log("Starting FoodAtlas data processing and merge...")
        
    triplets_path = os.path.join(data_dir, "triplets.parquet")
    entities_path = os.path.join(data_dir, "entities.parquet")
    rels_path = os.path.join(data_dir, "relationships.parquet")
    
    df_triplets = get_foodatlas_triplets(triplets_path)
    df_entities = get_foodatlas_entities(entities_path)
    df_rels = get_foodatlas_relationships(rels_path)
    
    if allowed_relationships:
        valid_rel_ids = df_rels[df_rels['name'].isin(allowed_relationships)]['foodatlas_id']
        df_triplets = df_triplets[df_triplets['relationship_id'].isin(valid_rel_ids)]
        
    df = df_triplets.merge(df_rels, left_on='relationship_id', right_on='foodatlas_id', how='left')
    df.rename(columns={'name': 'relationship_name'}, inplace=True)
    df.drop(columns=['foodatlas_id'], inplace=True) 
    
    df = df.merge(df_entities, left_on='head_id', right_on='foodatlas_id', how='left')
    df.rename(columns={
        'entity_type': 'head_entity_type',
        'common_name': 'head_common_name',
        'synonyms': 'head_synonyms'
    }, inplace=True)
    df.drop(columns=['foodatlas_id'], inplace=True)
    
    df = df.merge(df_entities, left_on='tail_id', right_on='foodatlas_id', how='left')
    df.rename(columns={
        'entity_type': 'tail_entity_type',
        'common_name': 'tail_common_name',
        'synonyms': 'tail_synonyms'
    }, inplace=True)
    df.drop(columns=['foodatlas_id'], inplace=True)
    
    attestations_path = os.path.join(data_dir, "attestations.parquet")
    evidence_path = os.path.join(data_dir, "evidence.parquet")
    
    df_attestations = get_foodatlas_attestations(attestations_path)
    df_evidence = get_foodatlas_evidence(evidence_path)
    
    df['attestation_ids'] = df['attestation_ids'].apply(
        lambda x: json.loads(x) if isinstance(x, str) else x
    )
    
    df = df.explode('attestation_ids')
    
    df = df.merge(df_attestations, left_on='attestation_ids', right_on='attestation_id', how='left')
    
    df = df.merge(df_evidence, on='evidence_id', how='left')
    
    groupby_cols = [
        'head_id', 'relationship_id', 'tail_id', 'relationship_name',
        'head_entity_type', 'head_common_name', 'head_synonyms',
        'tail_entity_type', 'tail_common_name', 'tail_synonyms'
    ]
    
    df = df.groupby(groupby_cols, dropna=False).agg({
        'reference': lambda x: list(x.dropna())
    }).reset_index()
    
    _log("FoodAtlas data processing completed successfully.")
    
    return df