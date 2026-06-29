#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

from __future__ import annotations

from typing import NamedTuple
import os 
import gzip
import requests

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.inputs.common as inputs_common
import pypath.share.cache as pypath_cache

_notavail = lambda x: None if x == 'Not Available' else x
_synonyms = lambda x: (
    tuple(sorted(y.strip() for y in x.split('|'))) if x else ()
)


class AdrecsAdr(NamedTuple):
    adr_class: str
    badd: str


class AdrecsChildParent(NamedTuple):
    child: AdrecsAdr
    parent: AdrecsAdr


class AdrecsDrugAdr(NamedTuple):
    drug_badd: str
    drug: str
    adr_badd: str
    adr: str


class AdrecsTerm(NamedTuple):
    adrecs_class: str
    badd: str
    name: str
    synonyms: tuple[str, ...]
    meddra: str


class AdrecsDrug(NamedTuple):
    badd: str
    drug: str
    synonyms: tuple[str, ...]
    drugbank: str
    pubchem_cid: str
    mesh: str
    kegg: str
    tdd: str


def adrecs_drug_identifiers() -> set[AdrecsDrug]:
    """
    Drug identifiers from the AdReCS database.

    Extracts IUPAC name, synonyms, DrugBank, MeSH, KEGG, and TDD IDs of drugs.
    http://www.bio-add.org/ADReCS/index.jsp

    Returns a set of AdrecsDrug namedtuples.
    """

    raw_data = _adrecs_base(
        url_key = 'drug_information', 
        record = AdrecsDrug,
        cell_range = 'A1:H2527',
        synonym_idx = [2],
    )

    result = set() 
    for drug in raw_data:
        if drug.badd:
            result.add(drug)
            
    return result

def adrecs_adr_ontology() -> set[AdrecsTerm]:
    """
    Adverse drug reaction (ADR) ontology from the AdReCS database.

    Returns a set of AdrecsTerm namedtuples.
    """

    raw_data = _adrecs_base(
        url_key = 'terminology', 
        record = AdrecsTerm,
        cell_range = 'A1:E13856',
        synonym_idx = [3],
    )

    result = set()
    for term in raw_data:
        if term.badd:
            result.add(term)
            
    return result

def adrecs_drug_adr() -> set[AdrecsDrugAdr]:
    """
    Drug-ADR pairs from the AdReCS database.
    """
    import pypath.share.cache as pypath_cache
    
    url = urls.urls['adrecs']['adrecs_drugs']
    result = set()
    
    cache_dir = pypath_cache.get_cachedir()

    local_file_patterns = [
        os.path.join(cache_dir, 'Drug_ADR_v3.3.txt.gz'),
        os.path.join(cache_dir, 'adrecs_adrecs_drugs_Drug_ADR.txt.gz')
    ]
    
    chosen_path = None

    for p in local_file_patterns:
        if os.path.exists(p):
            chosen_path = p
            break
            
    if chosen_path:

        try:
            with gzip.open(chosen_path, 'rt', encoding='utf-8') as f:
                _ = next(f)
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) < 4:
                        continue
                    record = AdrecsDrugAdr(*fields[:4])
                    result.add(record)

            if len(result) > 0:
                return result
        
        except Exception:
            result.clear()
           
    try:
        path = curl.Curl(url, silent=False, large=True)

        if path.outfile is None or not os.path.exists(path.outfile):
           
            return result
            
        with gzip.open(path.outfile, 'rt', encoding='utf-8') as f:

            _ = next(f) 

        for line in f:
            fields = line.strip().split('\t')

            if len(fields) < 4:
                continue

            record = AdrecsDrugAdr(*fields[:4])
            result.add(record)
        
    except Exception as e:

        return set()

    return result

def adrecs_hierarchy() -> set[AdrecsChildParent]:
    """
    Child-parent relationships between AdReCS ontology terms.
    """

    adr_ontology = adrecs_adr_ontology()

    child_adrs = {
        record.adrecs_class: record.badd
        for record in adr_ontology if record.adrecs_class
    }

    result = set()

    for field in adr_ontology:

        if not field.adrecs_class or '.' not in field.adrecs_class:
            continue

        parent_adrecs = field.adrecs_class.rsplit('.', 1)[0]
        parent_badd = child_adrs.get(parent_adrecs)

        if parent_badd is None:
            continue

        child_obj = AdrecsAdr(adr_class = field.adrecs_class, badd = field.badd)
        parent_obj = AdrecsAdr(adr_class = parent_adrecs, badd = child_adrs.get(parent_adrecs))
        relation = AdrecsChildParent(child = child_obj, parent = parent_obj)
        
        result.add(relation)

    return result


def _adrecs_base(
        url_key: str,
        record: type,
        cell_range: str,
        synonym_idx: list[int],) -> list[NamedTuple]:
    """
    Helper function which downloads and parses the excel files.
    """
    import tempfile
    
    url = urls.urls['adrecs'][url_key]
    result = []

    file_mapping = {
        'drug_information': 'Drug_information',
        'terminology': 'ADR_ontology'
    }
    
    local_filename = file_mapping.get(url_key)
    excel_path = None

    if local_filename:

        cache_dir = pypath_cache.get_cachedir()
        potential_path = os.path.join(cache_dir, local_filename)

        if os.path.exists(potential_path):
                       
            try:
                
                contents = inputs_common.read_xls(potential_path, cell_range = cell_range)
                
                if contents and len(contents) > 1:
                    excel_path = potential_path

            except Exception:
                excel_path = None

    if not excel_path:
        try:
            headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}
            r = requests.get(url, headers=headers, timeout=90)
            
            if r.status_code != 200:        
                return result
                
            with tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx') as temp_file:
                
                temp_file.write(r.content)
                excel_path = temp_file.name

        except Exception:
            return result

    try:
        contents = inputs_common.read_xls(excel_path, cell_range = cell_range)
    
    except Exception:
    
        return result

    for line in contents[1:]:
        line = [_notavail(x) for x in line]

        for isyn in synonym_idx:
            line[isyn] = _synonyms(line[isyn])

        try:
            result.append(record(*line))
        except TypeError:
            continue

    return result
