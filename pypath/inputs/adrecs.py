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

from typing import Generator, NamedTuple

import collections

# In accordance with the guidelines, the pandas dependency has been removed.

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.inputs.common as inputs_common


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
    synonyms: tuple[str]
    meddra: str


class AdrecsDrug(NamedTuple):
    badd: str
    drug: str
    synonyms: tuple[str]
    drugbank: str
    pubchem_cid: str
    mesh: str
    kegg: str
    tdd: str

# MAIN INTEGRATION FUNCTIONS

def adrecs_drug_identifiers() -> dict[str, set[AdrecsDrug]]:
    """
    Drug identifiers from the AdReCS database.

    Extracts IUPAC name, synonyms, DrugBank, MeSH, KEGG, and TDD IDs of drugs.
    http://www.bio-add.org/ADReCS/index.jsp

    Returns:
       A dictionary mapping primary BADD IDs to a set of AdrecsDrug namedtuples.
    """

    raw_data = _adrecs_base(
        url_key = 'drug_information', 
        record = AdrecsDrug,
        cell_range = 'A1:H2527',
        synonym_idx = [2],
    )

    result = collections.defaultdict(set) 
    for drug in raw_data:
        if drug.badd:
            result[drug.badd].add(drug)
            
    return dict(result)

def adrecs_adr_ontology() -> dict[str, set[AdrecsTerm]]:
    """
    Adverse drug reaction (ADR) ontology from the AdReCS database.

    Extracts ADR term classification, name, synonyms, and MedDRA identifiers.

    Returns:
        A dictionary mapping ADR BADD IDs to a set of AdrecsTerm namedtuples.
    """

    raw_data = _adrecs_base(
        url_key = 'terminology', 
        record = AdrecsTerm,
        cell_range = 'A1:E13856',
        synonym_idx = [3],
    )

    result = collections.defaultdict(set)
    for term in raw_data:
        if term.badd:
            result[term.badd].add(term)
            
    return dict(result)

def adrecs_drug_adr() -> dict[str, set[AdrecsDrugAdr]]:
    """
    Drug-ADR pairs from the AdReCS database.

    Parses the mapping between drugs and their respective adverse reactions.

    Returns:
        A dictionary mapping drug BADD IDs to a set of AdrecsDrugAdr namedtuples.
    """

    url = urls.urls['adrecs']['adrecs_drugs']
    c = curl.Curl(url, large = True, silent = False)
    _ = next(c.result) 

    result = collections.defaultdict(set)
    
    for line in c.result:
        fields = line.strip().split('\t')


        if len(fields) < 4:
            continue
            
        record = AdrecsDrugAdr(*fields)
        result[record.drug_badd].add(record)
        
    return dict(result)

def adrecs_hierarchy() -> dict[str, set[AdrecsChildParent]]:
    """
    Child-parent relationships between AdReCS ontology terms.

    Reconstructs the hierarchical tree structure of adverse reactions.

    Returns:
         A dictionary mapping child BADD IDs to a set of AdrecsChildParent namedtuples.
    """

    url = urls.urls['adrecs']['terminology']
    path = curl.Curl(url, silent = True, large = True)
    contents = inputs_common.read_xls(path.outfile, cell_range = 'A1:E13856')
    
    adr_ontology = []
    for line in contents[1:]:
        line = [_notavail(x) for x in line]
        line[3] = _synonyms(line[3])
        adr_ontology.append(AdrecsTerm(*line))

    child_adrs = {
        record.adrecs_class: record.badd
        for record in adr_ontology
    }

    result = collections.defaultdict(set)

    for field in adr_ontology:
        if '.' not in field.adrecs_class:
            continue

        parent_adrecs = field.adrecs_class.rsplit('.', 1)[0]
        
        child_obj = AdrecsAdr(adr_class = field.adrecs_class, badd = field.badd)
        parent_obj = AdrecsAdr(adr_class = parent_adrecs, badd = child_adrs.get(parent_adrecs))
        
        relation = AdrecsChildParent(child = child_obj, parent = parent_obj)
        
        result[field.badd].add(relation)

    return dict(result)

#HELPER FUNCTIONS

def _adrecs_base(
        url_key: str,
        record: type,
        cell_range: str,
        synonym_idx: list[int],) -> list[tuple]:

    """
    Helper function which downloads and parses the excel files. 
    
    Critical Fixation: NameError cause 'record_name' block completely removed.
    """

    url = urls.urls['adrecs'][url_key]
    path = curl.Curl(url, silent = False, large = True)
    contents = inputs_common.read_xls(path.outfile, cell_range = cell_range)
    result = []

    for line in contents[1:]:
        line = [_notavail(x) for x in line]

        for isyn in synonym_idx:
            line[isyn] = _synonyms(line[isyn])

        result.append(record(*line))

    return result
