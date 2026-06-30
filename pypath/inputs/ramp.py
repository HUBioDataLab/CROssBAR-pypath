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

"""
Access the RaMP metabolomic pathway and metabolite database.
"""

from typing import IO, Literal

import gzip
import json
import os
import pathlib
import pprint
import re
import sqlite3

import pandas as pd

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session
import pypath.share.common as common
import pypath.formats.sqldump as sqldump

_log = session.Logger(name = 'ramp_input')._log

_FALLBACK_VERSION = '3.0.7'


def _latest_ramp_version() -> str:
    """
    Retrieve the latest RaMP SQLite version number from GitHub API.
    """

    try:
        api_url = urls.urls['ramp']['github_api']
        c = curl.Curl(api_url, silent=True, large=False)

        if c.result is None:
            raise ValueError('Empty response from GitHub API')

        contents = json.loads(c.result)
        versions = sorted(
            m.group(1)
            for item in contents
            if (m := re.search(r'RaMP_SQLite_v([\d.]+)\.sqlite\.gz', item['name']))
        )

        if versions:
            _log(f'Latest RaMP SQLite version from GitHub: {versions[-1]}')
            return versions[-1]

    except Exception as e:
        _log(f'Could not determine latest RaMP version from GitHub: {e}')

    return _FALLBACK_VERSION


def ramp_sqlite(version: str | None = None) -> sqlite3.Connection:
    """
    Download and open the RaMP database in SQLite format from GitHub.

    Downloads from https://github.com/ncats/RaMP-DB/tree/main/db and
    caches in ~/.cache/pypath/. This is the preferred data source when
    the NIH REST API (rampdb.nih.gov) is unavailable.

    Args:
        version:
            Specific RaMP version to download (e.g. '3.0.7').
            If None, the latest version is looked up from GitHub.

    Returns:
        An open sqlite3 connection to the RaMP database.
    """

    version = version or _latest_ramp_version()

    cache_dir = pathlib.Path.home() / '.cache' / 'pypath'
    cache_dir.mkdir(parents=True, exist_ok=True)

    sqlite_path = cache_dir / f'RaMP_SQLite_v{version}.sqlite'

    # Treat an empty or very small file as corrupted from a previous failed attempt
    if sqlite_path.exists() and sqlite_path.stat().st_size < 1024:
        _log(f'Removing incomplete RaMP SQLite cache at {sqlite_path}')
        sqlite_path.unlink()

    if not sqlite_path.exists():

        sqlite_url = urls.urls['ramp']['github_sqlite'] % version
        _log(f'Downloading RaMP SQLite v{version} from {sqlite_url}')

        c = curl.Curl(sqlite_url, large=True, silent=False)

        if not hasattr(c, 'cache_file_name') or not os.path.exists(c.cache_file_name):
            raise RuntimeError(f'Failed to download RaMP SQLite v{version}')

        # Validate gzip magic bytes — catches Git LFS pointer downloads
        with open(c.cache_file_name, 'rb') as f:
            magic = f.read(2)

        if magic != b'\x1f\x8b':
            raise RuntimeError(
                f'Downloaded file is not a valid gzip archive (magic={magic!r}). '
                f'The URL may have returned a Git LFS pointer instead of the actual file.'
            )

        # Decompress into a temp file; only rename on success to avoid partial sqlite
        tmp_path = sqlite_path.with_suffix('.tmp')
        tmp_path.unlink(missing_ok=True)

        _log('Decompressing RaMP SQLite...')

        try:
            with gzip.open(c.cache_file_name, 'rb') as gz_in:
                with open(tmp_path, 'wb') as out:
                    while chunk := gz_in.read(1024 * 1024):
                        out.write(chunk)
            tmp_path.rename(sqlite_path)
        except Exception:
            tmp_path.unlink(missing_ok=True)
            raise

        _log(f'RaMP SQLite v{version} ready at {sqlite_path}')

    return sqlite3.connect(str(sqlite_path))


def _ramp_sqldump() -> IO:
    """
    Download the RaMP metabolomic pathway and metabolite database.
    """

    url = urls.urls['ramp']['url']
    c = curl.Curl(url, large = True, silent = False, compr = 'gz')

    return c._gzfile_mode_r


def ramp_raw(
        tables: list[str] = None,
        sqlite: bool = False,
        **kwargs
    ) -> dict[str, pd.DataFrame, sqlite3.Connection]:
    """
    Retrieve RaMP database contents from raw SQL dump.

    Args:
        tables:
            One or more tables to retrieve. If None, all tables are retrieved.
        sqlite:
            Return an SQLite database instead of a pandas DataFrame.
        kwargs:
            Options for the SQLite database: this way you can point to a new
            or existing database, while by default, an in-memory, temporary
            database is used.

    Returns:
        Either a dictionary with the table names as keys and  pandas dataframes
        as values, or an SQLite database connection.
    """

    fp = _ramp_sqldump()

    return sqldump.tables(
        sqldump = fp,
        tables = tables,
        return_df = True,
        return_sqlite = sqlite,
        con_param = kwargs,
        source_id = (fp.name, f'{os.path.getmtime(fp.name):.0f}'),
    )


def ramp_list_tables() -> dict[str, list[str]]:
    """
    List the tables of the RaMP database from SQL dump.
    """

    return sqldump.list_tables(_ramp_sqldump())


def ramp_show_tables() -> None:
    """
    Show the tables of the RaMP database from SQL dump.
    """

    pprint.pprint(ramp_list_tables())


def ramp_mapping(
        id_type_a: str,
        id_type_b: str,
        return_df: bool = False,
        curies: bool = False,
    ) -> dict[str, set[str]] | pd.DataFrame:
    """
    Retrieve the mapping between two identifiers.

    Args:
        id_type_a:
            The identifier type of the first identifier.
        id_type_b:
            The identifier type of the second identifier.
        return_df:
            Return a pandas DataFrame instead of a dictionary.
        curies:
            Do not remove CURIEs from the identifiers.

    Returns:
        A dictionary with the mapping between the two identifiers.
    """

    query = (
        'SELECT DISTINCT a.sourceId as id_type_a, b.sourceId as id_type_b '
        'FROM '
        '   (SELECT sourceId, rampId '
        '    FROM source '
        f'   WHERE geneOrCompound = "compound" AND IDtype = "{id_type_a}") a '
        'JOIN '
        '   (SELECT sourceId, rampId '
        '    FROM source '
        f'   WHERE geneOrCompound = "compound" AND IDtype = "{id_type_b}") b '
        'ON a.rampId = b.rampId;'
    )

    try:
        con = ramp_sqlite()
        _log('RaMP mapping: using SQLite backend')
    except Exception as e:
        _log(f'RaMP SQLite unavailable ({e}), falling back to SQL dump')
        con = ramp_raw(tables = 'source', sqlite = True)

    df = pd.read_sql_query(query, con)

    if not curies:

        df[df.columns] = df[df.columns].apply(
            lambda y: [x.split(':', maxsplit = 1)[-1] for x in y],
        )

    return (
        df
            if return_df else
        df.groupby('id_type_a')['id_type_b'].apply(set).to_dict()
    )


def ramp_id_types(
        entity_type: Literal['gene', 'compound'] | None = None,
    ) -> set[str]:
    """
    List the identifier types of the RaMP database.

    Uses the SQLite backend (GitHub) as primary source, falls back
    to the SQL dump if SQLite is unavailable.
    """

    query = (
        'SELECT DISTINCT(IDtype) AS id_type FROM source' +
        (f' WHERE geneOrCompound = "{entity_type}"' if entity_type else '')
    )

    try:
        con = ramp_sqlite()
        df = pd.read_sql_query(query, con)
        con.close()
        return set(df['id_type'])
    except Exception as e:
        _log(f'RaMP SQLite unavailable ({e}), falling back to SQL dump for ID types')

    sqldump_query = (
        'SELECT DISTINCT(s.IDtype) as id_type FROM source s' +
        (f' WHERE geneOrCompound = "{entity_type}";' if entity_type else ';')
    )
    con = ramp_raw(tables = 'source', sqlite = True)
    df = pd.read_sql_query(sqldump_query, con)

    return set(df['id_type'])


def ramp_id_types_2(
        entity_type: Literal['gene', 'compound'] | None = None,
    ) -> set[str]:
    """
    List the identifier types of the RaMP database via REST API.

    Tries the NIH REST API first (fast, no download required). Falls back
    to ramp_id_types() (SQLite/SQL dump) if the API is unavailable.
    """

    entity_types = {
        'compound': 'Metabolites',
        'gene': 'Genes/Proteins',
    }

    url = urls.urls['ramp']['api'] % 'id-types'
    c = curl.Curl(url, silent = True, large = False)

    try:
        return {
            id_type.strip()
            for i in json.loads(c.result)['data']
            if not entity_type or i['analyteType'] == entity_types[entity_type]
            for id_type in i['idTypes'].split(',')
        }
    except Exception as e:
        _log(f'RaMP REST API unavailable ({e}), falling back to SQLite/SQL dump')
        return ramp_id_types(entity_type)
