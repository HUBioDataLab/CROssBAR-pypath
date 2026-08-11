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

import csv
from typing import Union

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.internals.intera as intera
import pypath.utils.taxonomy as taxonomy

def _split_semicolon(raw, expected_length=None):
    """
    Splits a semicolon-delimited CORUM field.

    A trailing semicolon is treated as an export artifact only when
    dropping the resulting spurious empty entry makes the result
    match `expected_length`; otherwise the trailing empty entry is
    kept, since it may be a real "no value" marker for the last
    subunit. When `expected_length` is None (used for the UniProt ID
    field itself, which defines the expected length for the other
    fields), a lone trailing empty entry is always dropped, since an
    empty UniProt ID carries no information and gets filtered out
    anyway.
    """

    parts = raw.split(';')

    if expected_length is None:

        if parts and parts[-1] == '':

            parts = parts[:-1]

    elif len(parts) == expected_length + 1 and parts[-1] == '':

        parts = parts[:-1]

    return parts

def corum_complexes(organism: Union[int, str] = "all"):
    """

    Retrieves the "Complete complexes" dataset. The
    file uses a new snake_case schema (`complex_id`, `subunits_uniprot_id`,
    etc.).

    Note 1 : records are kept separate, keyed by CORUM's own `complex_id`
    (unique per row); we don't aggregate rows that happen to share a
    name and component set, so that records CORUM lists separately
    stay separate here too.

    Note 2 : when the top-level `organism` field can't be resolved to an
    NCBI Taxonomy ID (e.g. "MINK", ambiguous between American mink,
    452646, and European mink, 9666), we fall back to
    `subunits_organism`: if every subunit agrees on one resolvable
    organism, we use that; otherwise the complex is skipped rather
    than guessed.

    Note 3 : mixed-species complexes (subunits from more than one
    organism) can't have a single species, so CORUM labels their
    `organism` as "Mammalia" instead. This resolves fine on its own to
    NCBI Taxonomy ID 40674, so no special handling is needed - it's
    included with `organism="all"` and correctly excluded when filtering
    for one species (e.g. `organism=9606`).

    Args:
        organism: NCBI Taxonomy ID of the organism to keep; complexes
            of other organisms are discarded. `None` or `"all"`
            disables filtering (all organisms kept).

    Returns:
        A dict of `intera.Complex` objects keyed by CORUM's `complex_id`.
    """

    if organism in (None, 'all'):
        organism = None
    else:
        organism = taxonomy.ensure_ncbi_tax_id(organism)

    url = urls.urls['corum']['url']
    c = curl.Curl(url, large = True, silent = False)

    tab = csv.DictReader(c.result, delimiter = '\t')

    complexes = {}

    for rec in tab:

        tax_id = taxonomy.ensure_ncbi_tax_id(rec['organism'])

        if tax_id is None:

            # Note 2 is applied here.
            subunit_organisms = {
                o.strip()
                for o in _split_semicolon(rec['subunits_organism'])
                if o.strip()
            }

            if len(subunit_organisms) == 1:

                tax_id = taxonomy.ensure_ncbi_tax_id(subunit_organisms.pop())

        if tax_id is None:

            raise ValueError(
                f'Could not resolve organism for CORUM complex '
                f'`{rec["complex_id"]}` (organism = `{rec["organism"]}`, '
                f'subunits_organism = `{rec["subunits_organism"]}`).'
            )

        # organism resolved successfully, just not the one requested
        if organism and tax_id != organism:

            continue

        uniprots_raw = _split_semicolon(rec['subunits_uniprot_id'])
        genesymbols_raw = _split_semicolon(
            rec['subunits_gene_name'], len(uniprots_raw)
        )

        # subunit fields are positionally aligned (same index = same
        # subunit); a length mismatch means the alignment can't be
        # trusted, so we fail loudly rather than guess
        if len(uniprots_raw) == len(genesymbols_raw):

            subunit_pairs = zip(uniprots_raw, genesymbols_raw)

        else:

            raise ValueError(
                f'Mismatched subunit field lengths for CORUM complex '
                f'`{rec["complex_id"]}`: subunits_uniprot_id has '
                f'{len(uniprots_raw)} entries, subunits_gene_name has '
                f'{len(genesymbols_raw)}.'
            )

        # filter as pairs, not independently, so positions stay in sync
        subunit_pairs = [(u, g) for u, g in subunit_pairs if u]

        if not subunit_pairs:

            continue

        gene_names = {u: g for u, g in subunit_pairs}
        stoich_raw = _split_semicolon(
    rec['subunits_stoechiometrie'], len(uniprots_raw)
)
 
        # unlike uniprot/genesymbol, a stoichiometry mismatch doesn't
        # invalidate the whole record - we just can't trust the
        # per-subunit values, so every subunit falls back to the
        # "one copy, count unknown" default individually
        if len(stoich_raw) == len(uniprots_raw):
 
            stoich_pairs = zip(uniprots_raw, stoich_raw)
 
        else:
 
            stoich_pairs = zip(uniprots_raw, [''] * len(uniprots_raw))
 
        stoichiometry = {}
 
        for u, s in stoich_pairs:
 
            if not u:
 
                continue
 
            try:
 
                stoichiometry[u] = int(s)
 
            except (TypeError, ValueError):
 
                # empty or non-numeric entry (CORUM leaves this field
                # blank often); default to one copy
                stoichiometry[u] = 1


        pubmeds = {p for p in _split_semicolon(rec['pmid']) if p}
 
        evi_raw = _split_semicolon(rec['functions_evi'])
        pmid_raw = _split_semicolon(rec['functions_pmid'], len(evi_raw))
        goid_raw = _split_semicolon(rec['functions_go_id'], len(evi_raw))

        # functions_evi / functions_pmid / functions_go_id are a
        # triple, positionally aligned (each index = one GO
        # annotation).We store them as a set of (go_id, evidence,
        # pmid) tuples rather than three parallel lists: this makes
        # the record self-aligned (no index bookkeeping) instead of
        # three lists that can silently drift out of sync.
        if len(evi_raw) == len(pmid_raw) == len(goid_raw):

            go_annotations = {
                (g, e, p) for e, p, g in zip(evi_raw, pmid_raw, goid_raw)
                if g
            }

        else:

            raise ValueError(
                f'Mismatched GO annotation field lengths for CORUM '
                f'complex `{rec["complex_id"]}`: functions_evi has '
                f'{len(evi_raw)}, functions_pmid has {len(pmid_raw)}, '
                f'functions_go_id has {len(goid_raw)}.'
            )

        cplex = intera.Complex(
            name = rec['complex_name'],
            components = stoichiometry,
            sources = 'CORUM',
            references = pubmeds,
            ncbi_tax_id = tax_id,
            ids = rec['complex_id'],
            attrs = {
                'synonyms': rec['synonyms'],
                'cell_line': rec['cell_line'],
                'comment': rec['comment_complex'],
                'gene_names': gene_names,
                'go_annotations': go_annotations,
            },
        )

        complexes[rec['complex_id']] = cplex

    return complexes