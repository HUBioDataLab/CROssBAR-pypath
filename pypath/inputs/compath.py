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

from typing import Generator, Literal

import collections
import pandas as pd

import pypath.resources.urls as urls
import pypath.share.curl as curl


def compath_mappings(
        source_db: Literal['kegg', 'wikipathways', 'reactome'] | None = None,
        target_db: Literal['kegg', 'wikipathways', 'reactome'] | None = None,
        return_df: bool = False,
    ) -> Generator[tuple] | pd.DataFrame:
    """
    Cross-database pathway to pathway mappings from Compath.

    Compath contains proposed and accepted mappings by the users/curators
    between pairs of pathways across databases. The source and target
    databases specify the direction of the mapping.

    Args:
        source_db:
            Name of the source database.
        target_db:
            Name of the target database.
        return_df:
            Return a pandas data frame.

    Returns:
        Tuples of pathway-to-pathway mappings.
    """

    result = _compath_mappings(source_db, target_db)

    return pd.DataFrame(result) if return_df else result


def _compath_mappings(
        source_db: Literal['kegg', 'wikipathways', 'reactome'] | None = None,
        target_db: Literal['kegg', 'wikipathways', 'reactome'] | None = None,
    ) -> Generator[tuple]:

    fields = (
        'pathway1',
        'pathway_id_1',
        'source_db',
        'relation',
        'pathway2',
        'pathway_id_2',
        'target_db',
    )
    record = collections.namedtuple('CompathPathwayToPathway', fields)
    for file_url in urls.urls['compath']['github_urls']:

        c = curl.Curl(file_url, large = True)

        if c.result is None:
            continue

        lines = list(c.result)

        for line in lines[1:]:  # basligi atla

            parts = line.strip().split(',')

            if len(parts) < 7:
                continue

            src_res, src_id, src_name, relation, tgt_res, tgt_id, tgt_name = parts[:7]

            src_res = 'kegg' if src_res == 'kegg.pathway' else src_res
            tgt_res = 'kegg' if tgt_res == 'kegg.pathway' else tgt_res

            if (
                (source_db is None or src_res == source_db) and
                (target_db is None or tgt_res == target_db)
            ):

                yield record(
                    pathway1 = src_name,
                    pathway_id_1 = src_id,
                    source_db = src_res,
                    relation = relation,
                    pathway2 = tgt_name,
                    pathway_id_2 = tgt_id,
                    target_db = tgt_res,
                )