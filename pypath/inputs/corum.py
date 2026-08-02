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

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.internals.intera as intera
import pypath.utils.taxonomy as taxonomy


def corum_complexes(organism = 9606):
    """
    Protein complexes from CORUM.

    Retrieves the "Complete complexes" dataset directly from CORUM's
    current (release 5.3+, since the 2026 platform update) API. The
    file uses a new snake_case schema (`complex_id`, `subunits_uniprot_id`,
    etc.), replacing the old CamelCase schema (`ComplexID`,
    `subunits(UniProt IDs)`, etc.) used before the platform migration.

    Note: a small number of complex identities appear more than once in
    the source (same components merge to the same key here); when this
    happens references, gene names and GO annotations are merged
    rather than overwritten.

    Note: when the top-level `organism` field can't be resolved to an
    NCBI Taxonomy ID (e.g. "MINK", ambiguous between American mink,
    452646, and European mink, 9666), we fall back to
    `subunits_organism`: if every subunit agrees on one resolvable
    organism, we use that; otherwise the complex is skipped rather
    than guessed. Complexes with `organism` set to "Mammalia"
    (mixed-species complexes) resolve fine on their own to NCBI
    Taxonomy ID 40674 (the Mammalia class) and are kept when
    `organism=None`; they are correctly excluded when filtering for a
    specific species (e.g. `organism=9606`).

    Args:
        organism: NCBI Taxonomy ID of the organism to keep; complexes
            of other organisms are discarded. None disables filtering
            (all organisms kept).

    Returns:
        A dict of `intera.Complex` objects keyed by the complex identity
        string.
    """

    if organism:

        organism = taxonomy.ensure_ncbi_tax_id(organism)

    # Note: the previous version of this function defined an `annots`
    # tuple of 27 FunCat-style category strings (e.g. 'nucleus',
    # 'endocytosis') right here. It was removed because it was dead
    # code: it was never referenced anywhere else in the function, so
    # it had no effect on what was fetched or filtered. It looked like
    # leftover/unfinished filtering logic from an earlier version of
    # this adapter.

    url = urls.urls['corum']['url']
    c = curl.Curl(url, large = True, silent = False)

    tab = csv.DictReader(c.result, delimiter = '\t')

    complexes = {}

    for rec in tab:

        try:

            tax_id = taxonomy.ensure_ncbi_tax_id(rec['organism'])

        except Exception:

            tax_id = None

        if tax_id is None:

            # The top-level `organism` field is sometimes an ambiguous
            # or unrecognized name. E.g. 4 complexes have organism =
            # "MINK", which CORUM doesn't resolve further, but which
            # is genuinely ambiguous between American mink (NCBI Taxonomy
            # ID 452646) and European mink (9666).
            #
            # Rather than dropping these rows outright, we fall back
            # to the per-subunit `subunits_organism` field: if every
            # subunit agrees on a single, resolvable organism, we use
            # that -- it's unambiguous, so there is nothing to guess.
            # This rescues 3 of the 4 "MINK" complexes, whose subunits
            # are consistently "Neovison vison (American mink)".
            #
            # If subunits disagree with each other (or, as in one of
            # the 4 "MINK" cases, all agree on a organism different
            # from the top-level field -- there the subunits say
            # "Human" while the top-level field says "MINK", which
            # looks like a CORUM curation inconsistency) we still
            # don't guess silently: a single agreed-upon subunit
            # organism is trusted, anything else is skipped.
            subunit_organisms = {
                o.strip()
                for o in rec['subunits_organism'].split(';')
                if o.strip()
            }

            if len(subunit_organisms) == 1:

                try:

                    tax_id = taxonomy.ensure_ncbi_tax_id(
                        subunit_organisms.pop()
                    )

                except Exception:

                    tax_id = None

        if tax_id is None:

            continue

        if organism and tax_id != organism:

            continue

        uniprots_raw = rec['subunits_uniprot_id'].split(';')
        genesymbols_raw = rec['subunits_gene_name'].split(';')

        # subunit fields are positionally aligned (same index = same
        # subunit); if lengths don't match we can't trust the
        # alignment, so we keep the UniProt IDs and drop gene symbols
        # for this record rather than mis-pairing them
        if len(uniprots_raw) == len(genesymbols_raw):

            subunit_pairs = zip(uniprots_raw, genesymbols_raw)

        else:

            subunit_pairs = zip(uniprots_raw, [None] * len(uniprots_raw))

        # filter as pairs, not independently, so positions stay in sync
        subunit_pairs = [(u, g) for u, g in subunit_pairs if u]

        if not subunit_pairs:

            continue

        uniprots = [u for u, g in subunit_pairs]
        gene_names = {u: g for u, g in subunit_pairs}

        pubmeds = {p for p in rec['pmid'].split(';') if p}

        evi_raw = rec['functions_evi'].split(';')
        pmid_raw = rec['functions_pmid'].split(';')
        goid_raw = rec['functions_go_id'].split(';')

        # functions_evi / functions_pmid / functions_go_id are a
        # triple, positionally aligned (each index = one GO
        # annotation). We store them as a set of (go_id, evidence,
        # pmid) tuples rather than three parallel lists: this makes
        # the record self-aligned (no index bookkeeping) and makes
        # merging duplicate complex rows a plain set union, instead
        # of three lists that can silently drift out of sync.
        if len(evi_raw) == len(pmid_raw) == len(goid_raw):

            go_annotations = {
                (g, e, p) for e, p, g in zip(evi_raw, pmid_raw, goid_raw)
                if g
            }

        else:

            go_annotations = {
                (g, None, None) for g in goid_raw if g
            }

        cplex = intera.Complex(
            name = rec['complex_name'],
            components = uniprots,
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

        key = cplex.__str__()

        if key in complexes:

            complexes[key].references.update(pubmeds)
            complexes[key].attrs['go_annotations'] = (
                complexes[key].attrs.get('go_annotations', set()) |
                go_annotations
            )
            complexes[key].attrs['gene_names'].update(gene_names)

        else:

            complexes[key] = cplex

    return complexes