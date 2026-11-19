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

from typing import Generator, Iterable
import collections
import datetime

import pypath.share.session as session
import pypath.share.curl as curl
import pypath.resources.urls as urls

_logger = session.Logger(name = 'ncbi_taxonomy_input')

ROOT_TAX_ID = 1

UNRANKED = {'no rank', 'clade'}

NAME_CLASSES = {
    'scientific name': 'scientific_name',
    'common name': 'common_name',
    'genbank common name': 'genbank_common_name',
}

NcbiTaxon = collections.namedtuple(
    'NcbiTaxon',
    [
        'tax_id',
        'scientific_name',
        'common_name',
        'genbank_common_name',
        'rank',
        'parent_tax_id',
    ],
)
NcbiTaxon.__new__.__defaults__ = (None,) * 6


def _taxdump(files_needed: list[str]) -> tuple[dict, str]:
    """
    Downloads the NCBI Taxonomy taxdump archive and opens the requested
    members of it.

    Args:
        files_needed:
            Names of the `.dmp` files to open from the archive.

    Returns:
        A dict with file names as keys and binary file objects as values,
        and the release date of the taxdump, i.e. the modification time of
        the first requested file, in YYYY-MM-DD format.
    """

    url = urls.urls['ncbi_taxonomy']['taxdump']

    c = curl.Curl(
        url,
        large = True,
        silent = False,
        compr = 'tgz',
        files_needed = files_needed,
    )

    missing = set(files_needed) - set(c.result or {})

    if missing:

        msg = (
            'Files missing from NCBI Taxonomy taxdump archive: %s. '
            'The cache file `%s` might be corrupted.' % (
                ', '.join(sorted(missing)),
                c.cache_file_name,
            )
        )
        _logger._log(msg)
        raise RuntimeError(msg)

    mtime = c.tarfile.getmember(files_needed[0]).mtime
    release_date = datetime.datetime.fromtimestamp(
        mtime,
        tz = datetime.timezone.utc,
    ).strftime('%Y-%m-%d')

    return c.result, release_date


def _iter_dmp(fileobj) -> Generator[list[str], None, None]:
    """
    Iterates the records of a taxdump `.dmp` file. Fields are separated
    by `\\t|\\t` and each line ends with `\\t|`.
    """

    for line in fileobj:

        line = line.decode('utf-8').rstrip('\r\n')

        if line.endswith('\t|'):

            line = line[:-2]

        yield line.split('\t|\t')


def ncbi_taxonomy(
        ranks: set[str] | None = None,
        include_no_rank: bool = True,
        keep_root: bool = False,
        tax_ids: Iterable[int] | None = None,
        release_date: bool = False,
    ) -> dict[int, NcbiTaxon] | tuple[dict[int, NcbiTaxon], str]:
    """
    Taxa from the NCBI Taxonomy database, with names, rank and parent.

    Args:
        ranks:
            Keep only taxa of these ranks (e.g. ``{'species', 'genus'}``).
            If None, taxa of all ranks are kept.
        include_no_rank:
            Keep the taxa without a named rank (``no rank`` and ``clade``).
            The top level taxa ``cellular organisms`` and ``Viruses`` have
            the ranks ``cellular root`` and ``acellular root``, hence they
            are kept even if this is False.
        keep_root:
            Keep the root even if the filters above exclude it. It is then
            the closest kept ancestor of the taxa without any other kept
            ancestor.
        tax_ids:
            Keep only these taxa and all of their ancestors up to the root,
            regardless of their rank. If given, ``ranks`` and
            ``include_no_rank`` are not used. IDs missing from NCBI
            Taxonomy (e.g. merged or deleted taxa) are skipped.
        release_date:
            Return also the release date of the taxdump.

    If any taxon is excluded by the filters above, the ``parent_tax_id``
    of its descendants points to their closest ancestor which is kept,
    hence the lineage up to the root remains connected. The
    ``parent_tax_id`` is None for the root and for taxa without any kept
    ancestor.

    Returns:
        A dict with NCBI Taxonomy IDs as keys and taxon records as values.
        If ``release_date`` is True, a tuple of this dict and the release
        date of the taxdump in YYYY-MM-DD format.
    """

    files, taxdump_date = _taxdump(['nodes.dmp', 'names.dmp'])

    _logger._log('Processing NCBI Taxonomy nodes.')

    nodes = {}

    for rec in _iter_dmp(files['nodes.dmp']):

        nodes[int(rec[0])] = (int(rec[1]), rec[2])

    if tax_ids is not None:

        if ranks is not None or not include_no_rank:

            _logger._log(
                '`tax_ids` is given, `ranks` and `include_no_rank` '
                'are not used.'
            )

        tax_ids = {int(tax_id) for tax_id in tax_ids}
        missing = tax_ids - nodes.keys()

        if missing:

            _logger._log(
                '%u of the %u requested taxa are not in NCBI Taxonomy, '
                'e.g. %s.' % (
                    len(missing),
                    len(tax_ids),
                    ', '.join(str(i) for i in sorted(missing)[:10]),
                )
            )

        kept = set()

        for tax_id in tax_ids - missing:

            # the lineage until the first taxon which is already kept
            while tax_id not in kept:

                kept.add(tax_id)
                parent = nodes[tax_id][0]

                if parent == tax_id:

                    break

                tax_id = parent

    else:

        def keep(rank):

            return (
                (ranks is None or rank in ranks) and
                (include_no_rank or rank not in UNRANKED)
            )

        kept = {tax_id for tax_id, (_, rank) in nodes.items() if keep(rank)}

    if keep_root:

        kept.add(ROOT_TAX_ID)

    _logger._log('Processing NCBI Taxonomy names.')

    names = collections.defaultdict(dict)

    for rec in _iter_dmp(files['names.dmp']):

        tax_id = int(rec[0])
        field = NAME_CLASSES.get(rec[3])

        if field and tax_id in kept:

            # the first one of the multiple common names is kept
            names[tax_id].setdefault(field, rec[1])

    # closest kept ancestor of each node, filled lazily
    kept_ancestor = {}

    def closest_kept_ancestor(tax_id):

        path = []
        parent = nodes[tax_id][0]
        ancestor = None

        while parent != tax_id:

            if parent in kept:

                ancestor = parent
                break

            if parent in kept_ancestor:

                ancestor = kept_ancestor[parent]
                break

            path.append(parent)
            tax_id, parent = parent, nodes[parent][0]

        for unkept in path:

            kept_ancestor[unkept] = ancestor

        return ancestor

    result = {}

    for tax_id, (parent, rank) in nodes.items():

        if tax_id not in kept:

            continue

        if parent == tax_id:

            parent = None

        elif parent not in kept:

            parent = closest_kept_ancestor(tax_id)

        result[tax_id] = NcbiTaxon(
            tax_id = tax_id,
            rank = rank,
            parent_tax_id = parent,
            **names.get(tax_id, {}),
        )

    if release_date:

        return result, taxdump_date

    return result


def ncbi_taxonomy_merged() -> dict[int, int]:
    """
    Taxa merged into other taxa in NCBI Taxonomy.

    Returns:
        A dict with the old NCBI Taxonomy IDs as keys and the IDs of the
        taxa they have been merged into as values.
    """

    files, _ = _taxdump(['merged.dmp'])

    return {
        int(rec[0]): int(rec[1])
        for rec in _iter_dmp(files['merged.dmp'])
    }

