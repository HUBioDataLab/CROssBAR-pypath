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

from typing import Generator, Optional

import re
import csv
import json
import collections

import pypath.resources.urls as urls
import pypath.inputs.common as inputs_common
import pypath.share.curl as curl
import pypath.share.session as session

_logger = session.Logger(name = 'isdb_input')
_log = _logger._log

# the most recent release at the time of writing, used only if the list of
# the available releases can not be retrieved from GitHub
ISDB_FALLBACK_VERSION = '2026_05_09'

# the releases are dated files in the `versions` directory of the ISDB
# repository; the search indices in the same directory are not releases
ISDB_VERSION_FILE = re.compile(r'ISDB_(\d{4}_\d{2}_\d{2})\.tsv\.gz$')

ISDB_COLUMNS = (
    'Serial.Number',
    'Taxonomy.ID.A',
    'Organism.A',
    'UniProt.ID.A',
    'Protein.Name.A',
    'Taxonomy.ID.B',
    'Organism.B',
    'UniProt.ID.B',
    'Protein.Name.B',
    'Interaction.Type',
    'Ontology.ID',
    'Reference',
    'Database',
)

# the processed record yielded by the `_interactions` functions, defined
# once here rather than per record so all records share one type
IsdbInteraction = collections.namedtuple(
    'IsdbInteraction',
    (
        'uniprot_a',
        'uniprot_b',
        'isoform_a',
        'isoform_b',
        'taxon_a',
        'taxon_b',
        'interaction_type',
        'ontology_id',
        'pubmeds',
        'databases',
    ),
)


def isdb_latest_version() -> str:
    """
    The most recent ISDB release available on GitHub.

    The releases are distributed as dated files, their names carry the
    release date in `YYYY_MM_DD` format, hence the most recent one is the
    last in alphabetical order. If the release list can not be retrieved,
    falls back to `ISDB_FALLBACK_VERSION`.

    Returns
        (str): The date of the most recent release, e.g. `2026_05_09`.
    """

    try:

        c = curl.Curl(
            urls.urls['isdb']['versions'],
            silent = True,
            large = False,
            cache = False,
        )

        versions = sorted(
            match.group(1)
            for item in json.loads(c.result)
            if (match := ISDB_VERSION_FILE.match(item['name']))
        )

        if versions:

            _log('Latest ISDB release on GitHub: `%s`.' % versions[-1])

            return versions[-1]

        _log('No ISDB release found in the GitHub repository.')

    except Exception:

        _log('Failed to retrieve the list of ISDB releases:')
        _logger._log_traceback()

    _log('Falling back to ISDB release `%s`.' % ISDB_FALLBACK_VERSION)

    return ISDB_FALLBACK_VERSION


def isdb_raw(version: Optional[str] = None) -> Generator[tuple, None, None]:
    """
    All records of one ISDB release, without filtering.

    The distributed file is a gzipped, quoted TSV with a byte order mark,
    all its fields are carried over to the named tuples as they are.

    Args
        version: Date of an ISDB release, e.g. `2026_05_09`. By default the
            most recent release is used.

    Yields
        (tuple): Named tuples, each representing one ISDB record with all
            its fields as strings.
    """

    IsdbRecord = collections.namedtuple(
        'IsdbRecord',
        (
            'serial_number',
            'taxonomy_id_a',
            'organism_a',
            'uniprot_id_a',
            'protein_name_a',
            'taxonomy_id_b',
            'organism_b',
            'uniprot_id_b',
            'protein_name_b',
            'interaction_type',
            'ontology_id',
            'reference',
            'database',
        ),
    )

    version = version or isdb_latest_version()
    url = urls.urls['isdb']['url'] % version

    # the file is written by R: character fields are quoted, hence a plain
    # split by tabs would carry the quotes over into the data; the encoding
    # strips the byte order mark from the first column name
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        encoding = 'utf-8-sig',
    )

    records = csv.reader(c.result, delimiter = '\t', quotechar = '"')
    columns = tuple(next(records))

    if columns != ISDB_COLUMNS:

        msg = (
            'Unexpected columns in the ISDB release `%s`: expected `%s`, '
            'got `%s`.' % (
                version,
                ', '.join(ISDB_COLUMNS),
                ', '.join(columns),
            )
        )
        _log(msg)
        raise RuntimeError(msg)

    malformed = 0

    for record in records:

        if len(record) != len(ISDB_COLUMNS):

            malformed += 1
            continue

        yield IsdbRecord(*(field.strip() for field in record))

    if malformed:

        _log(
            'Skipped %u malformed records in the ISDB release `%s`.' % (
                malformed,
                version,
            )
        )


def _isdb_pubmeds(reference: str) -> set:
    """
    PubMed IDs from the `Reference` field of an ISDB record.

    The field is a `|` separated mix of `pubmed:`, `PMID:` (sometimes with
    a comma separated list), bare PMIDs, `doi:`, `imex:` and `mint:` items,
    a few records contain line breaks. Only the numeric PubMed IDs are kept,
    values such as `pubmed:DIP-17212S` or `PMID:nan` are dropped.
    """

    pubmeds = set()

    for item in re.split(r'[|\r\n]+', reference):

        item = item.strip()

        if item.startswith('pubmed:'):

            candidates = (item[7:],)

        elif item.startswith('PMID:'):

            candidates = item[5:].split(',')

        else:

            candidates = (item,)

        pubmeds.update(c.strip() for c in candidates if c.strip().isdigit())

    return pubmeds


def _isdb_databases(database: str) -> set:
    """
    Source databases from the `Database` field of an ISDB record, which
    is `|` separated when a record is supported by more than one resource.
    """

    return {db.strip() for db in database.split('|') if db.strip()}


def _isdb_interaction(record: tuple) -> tuple:
    """
    Processed record from a raw ISDB record, with only the fields needed
    for an interaction and in the form the other PPI input modules use.

    UniProt isoform suffixes are split off into `isoform_a` and
    `isoform_b`, taxonomy IDs are integers, references and databases are
    sets. A partner without UniProt ID has `None` in `uniprot_a` or
    `uniprot_b`, its taxon still identifies it.
    """

    uniprot_a, isoform_a = (
        inputs_common._try_isoform(record.uniprot_id_a)
            if record.uniprot_id_a else
        (None, None)
    )
    uniprot_b, isoform_b = (
        inputs_common._try_isoform(record.uniprot_id_b)
            if record.uniprot_id_b else
        (None, None)
    )

    return IsdbInteraction(
        uniprot_a = uniprot_a,
        uniprot_b = uniprot_b,
        isoform_a = isoform_a,
        isoform_b = isoform_b,
        taxon_a = int(record.taxonomy_id_a),
        taxon_b = int(record.taxonomy_id_b),
        interaction_type = record.interaction_type,
        ontology_id = record.ontology_id,
        pubmeds = _isdb_pubmeds(record.reference),
        databases = _isdb_databases(record.database),
    )


def _isdb_interactions(
        version: Optional[str] = None,
        organism: Optional[int] = None,
    ) -> Generator[tuple, None, None]:
    """
    Processed ISDB records, optionally restricted to one organism.

    Args
        version: Date of an ISDB release, e.g. `2026_05_09`. By default the
            most recent release is used.
        organism: NCBI Taxonomy ID. If provided, only the records where
            at least one of the partners belongs to this organism are
            yielded; by default all the records are, as ISDB is an
            interspecies resource.

    Yields
        (tuple): Named tuples, see `IsdbInteraction` for the fields.
    """

    # the taxonomy IDs are compared as strings on the raw records, so the
    # processed record is built only for the ones kept
    organism = None if organism is None else str(int(organism))

    for record in isdb_raw(version = version):

        if (
            organism is None or
            organism in (record.taxonomy_id_a, record.taxonomy_id_b)
        ):

            yield _isdb_interaction(record)


def isdb_ppi_interactions(
        version: Optional[str] = None,
        organism: Optional[int] = None,
    ) -> Generator[tuple, None, None]:
    """
    Protein level interactions from ISDB.

    Records where both partners are identified by a UniProt ID. The presence
    of the UniProt IDs splits the records into three disjoint sets, see also
    `isdb_protein_organism_interactions` and
    `isdb_organism_organism_interactions`.

    Args
        version: Date of an ISDB release, e.g. `2026_05_09`. By default the
            most recent release is used.
        organism: NCBI Taxonomy ID. If provided, only the records where
            at least one of the partners belongs to this organism are
            yielded; by default all the records are, as ISDB is an
            interspecies resource.

    Yields
        (tuple): `IsdbInteraction` named tuples, each representing one
            interaction.
    """

    for i in _isdb_interactions(version = version, organism = organism):

        if i.uniprot_a and i.uniprot_b:

            yield i


def isdb_protein_organism_interactions(
        version: Optional[str] = None,
        organism: Optional[int] = None,
    ) -> Generator[tuple, None, None]:
    """
    Interactions between a protein and an organism from ISDB.

    Records where only one of the partners is identified by a UniProt ID:
    one end is a protein, the other one is known only at the level of the
    taxon. Most of them come from host-pathogen resources (phi-base, PHISTO,
    VirHostNet), the rest from general interaction resources where the
    partner is not a UniProt protein. Which of them represent a host-pathogen
    relationship is not decided here, the records are passed on as they are.

    Args
        version: Date of an ISDB release, e.g. `2026_05_09`. By default the
            most recent release is used.
        organism: NCBI Taxonomy ID. If provided, only the records where
            at least one of the partners belongs to this organism are
            yielded; by default all the records are, as ISDB is an
            interspecies resource.

    Yields
        (tuple): `IsdbInteraction` named tuples, each representing one
            interaction.
    """

    for i in _isdb_interactions(version = version, organism = organism):

        if bool(i.uniprot_a) != bool(i.uniprot_b):

            yield i


def isdb_organism_organism_interactions(
        version: Optional[str] = None,
        organism: Optional[int] = None,
    ) -> Generator[tuple, None, None]:
    """
    Organism level interactions from ISDB.

    Records where neither partner is identified by a UniProt ID: ecological
    and epidemiological relationships such as `hasHost`, `preysOn` or
    `pathogenOf` between two taxa.

    Args
        version: Date of an ISDB release, e.g. `2026_05_09`. By default the
            most recent release is used.
        organism: NCBI Taxonomy ID. If provided, only the records where
            at least one of the partners belongs to this organism are
            yielded; by default all the records are, as ISDB is an
            interspecies resource.

    Yields
        (tuple): `IsdbInteraction` named tuples, each representing one
            interaction.
    """

    for i in _isdb_interactions(version = version, organism = organism):

        if not i.uniprot_a and not i.uniprot_b:

            yield i
