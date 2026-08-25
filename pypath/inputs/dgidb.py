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

import csv
import collections
import json

import bs4

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.share.session as session

_logger = session.Logger(name='dgidb_input')
_log = _logger._log

DGIDB_GRAPHQL_QUERY = """
    query getInteractions($chunk_size: Int!, $after_cursor: String) {
      interactions(first: $chunk_size, after: $after_cursor) {
        pageInfo {
          endCursor
          hasNextPage
        }
        edges {
          node {
            interactionScore
            interactionTypes {
              type
            }
            sources {
              sourceDbName
            }
            drug {
              name
              conceptId
            }
            gene {
              name
              conceptId
            }
            publications {
              pmid
            }
          }
        }
      }
    }
"""

def dgidb_graphql_request(variables: dict) -> dict:
    """
    Sends a request to the DGIdb GraphQL API and returns the result as JSON.
    """
    url = 'https://dgidb.org/api/graphql' 
    req_headers = {
        'Accept-Encoding': 'gzip, deflate, br',
        'Content-Type': 'application/json',
        'Connection': 'keep-alive'
    }
    
    query_param = {
        'query': DGIDB_GRAPHQL_QUERY,
        'variables': variables
    }

    binary_data = json.dumps(query_param).encode('utf-8')
    c = curl.Curl(url=url, silent=True, req_headers=req_headers, binary_data=binary_data)
    
    result = json.loads(c.result)
    return result.get('data', {})


def dgidb_interactions(chunk_size: int = 1000) -> list[tuple]:
    """
    Retrieves drug-gene interactions from DGIdb via GraphQL API.

    Returns:
        A list with tuples. Tuples are dgidb interactions
    """

    result = set()

    DgidbInteraction = collections.namedtuple(
        'DgidbInteraction',
        [
            'genesymbol',
            'gene_id',
            'resource',
            'type',
            'drug_name',
            'drug_id',
            'score',
            'pmid'
        ],
    )

    variables = {
        'chunk_size': chunk_size,
        'after_cursor': "", 
    }
    
    page_num = 1

    while True:
        _log(f"Fetching DGIdb data from API... Page #{page_num}")
        
        response = dgidb_graphql_request(variables)
        interactions_data = response.get('interactions') or {}
        edges = interactions_data.get('edges') or []
        page_info = interactions_data.get('pageInfo') or {}

        if not edges:
            break

        for edge in edges:
            node = edge.get('node') or {}
            
            drug = node.get('drug') or {}
            gene = node.get('gene') or {}
            pubs = node.get('publications') or []
            sources = node.get('sources') or []
            int_types = node.get('interactionTypes') or []
            
            drug_name = drug.get('name')
            drug_concept_id = drug.get('conceptId')
            gene_name = gene.get('name')
            gene_concept_id = gene.get('conceptId')
            score = node.get('interactionScore')
            
            pmids = [str(p.get('pmid')) for p in pubs if p.get('pmid')]
            pmids_str = ",".join(pmids) if pmids else None
            
            source_names = [s.get('sourceDbName') for s in sources if s.get('sourceDbName')]
            resource_str = ",".join(source_names) if source_names else None
            
            type_names = [t.get('type') for t in int_types if t.get('type')]
            type_str = ",".join(type_names) if type_names else None
            
            if drug_concept_id and gene_concept_id:
                dgidb_interaction = DgidbInteraction(
                    genesymbol = gene_name,
                    gene_id = gene_concept_id, 
                    resource = resource_str,
                    type = type_str,
                    drug_name = drug_name, 
                    drug_id = drug_concept_id,
                    score = score, 
                    pmid = pmids_str, 
                )
                result.add(dgidb_interaction)
        
        if not page_info.get('hasNextPage'):
            break
            
        variables['after_cursor'] = page_info.get('endCursor')
        page_num += 1

    return list(result)


def dgidb_annotations():
    """
    Downloads druggable protein annotations from DGIdb.
    """

    DgidbAnnotation = collections.namedtuple(
        'DgidbAnnotation',
        ['category'],
    )


    url = urls.urls['dgidb']['categories']
    c = curl.Curl(url = url, silent = False, large = True)
    data = csv.DictReader(c.result, delimiter = '\t')

    result = collections.defaultdict(set)

    for rec in data:

        uniprots = mapping.map_name(
            rec['name'], 
            'genesymbol',
            'uniprot',
        )

        for uniprot in uniprots:
            result[uniprot].add(
                DgidbAnnotation(
                    category = rec['name-2'] 
                )
            )

    return dict(result)


def get_dgidb_old():
    """
    Deprecated. Will be removed soon.

    Downloads and processes the list of all human druggable proteins.
    Returns a list of GeneSymbols.
    """

    genesymbols = []
    url = urls.urls['dgidb']['main_url']
    c = curl.Curl(url, silent = False)
    html = c.result
    soup = bs4.BeautifulSoup(html, 'html.parser')
    cats = [
        o.attrs['value']
        for o in soup.find('select', {'id': 'gene_categories'})
        .find_all('option')
    ]

    for cat in cats:
        url = urls.urls['dgidb']['url'] % cat
        c = curl.Curl(url)
        html = c.result
        soup = bs4.BeautifulSoup(html, 'html.parser')
        trs = soup.find('tbody').find_all('tr')
        genesymbols.extend([tr.find('td').text.strip() for tr in trs])

    return mapping.map_names(genesymbols, 'genesymbol', 'uniprot')