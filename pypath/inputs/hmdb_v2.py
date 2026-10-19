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

"""
Metabolites, proteins and their associations from the Human Metabolome
Database (HMDB), as named tuples like in the other input modules.

Independent of the older `pypath.inputs.hmdb` package, which parses the
same XML files into nested dicts and data frames.
"""

from typing import Generator, Literal, Optional

import re
import zipfile
import collections

import lxml.etree as etree

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session

_logger = session.Logger(name = 'hmdb_v2_input')
_log = _logger._log

HMDB_XMLNS = '{http://www.hmdb.ca}'

HMDB_RECORD_TAGS = {
    'metabolites': 'metabolite',
    'proteins': 'protein',
}

# identifiers of the metabolites in other resources: field names in the
# records and tag names in the metabolites XML
HMDB_METABOLITE_ID_TAGS = {
    'chebi': 'chebi_id',
    'pubchem': 'pubchem_compound_id',
    'drugbank': 'drugbank_id',
    'kegg': 'kegg_id',
    'vmh': 'vmh_id',
    'bigg': 'bigg_id',
    'chemspider': 'chemspider_id',
    'foodb': 'foodb_id',
    'pdb': 'pdb_id',
    'biocyc': 'biocyc_id',
    'wikipedia': 'wikipedia_id',
    'metlin': 'metlin_id',
    'knapsack': 'knapsack_id',
    'phenol_explorer': 'phenol_explorer_compound_id',
    'fbonto': 'fbonto_id',
}

# the two kinds of concentrations are in separate containers, and the
# subjects are described by differently named tags: container, age, sex
# and condition
HMDB_CONCENTRATION_TAGS = {
    'normal': (
        'normal_concentrations',
        'subject_age',
        'subject_sex',
        'subject_condition',
    ),
    'abnormal': (
        'abnormal_concentrations',
        'patient_age',
        'patient_sex',
        'patient_information',
    ),
}

_NUM = r'(\d+(?:\.\d*)?(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?)'
_PM = r'\s*(?:\+/-|\+-|±)\s*'
_RANGE = _NUM + r'\s*[-–]\s*' + _NUM

# the forms of the concentration values, each with the fields its groups
# fill in, tried in this order
HMDB_CONCENTRATION_FORMS = (
    (re.compile(_NUM), ('value',)),
    (re.compile(_NUM + _PM + _NUM), ('value', 'value_sd')),
    (re.compile(_RANGE), ('value_min', 'value_max')),
    # mean (SD): all values in this form come from seven papers, the two
    # behind 92% of them (PMID 24704061, 26010610) give mean (SD) in the
    # tables the values were taken from, as does PMID 23434268
    (re.compile(_NUM + r'\s*\(\s*' + _NUM + r'\s*\)'), ('value', 'value_sd')),
    # value (low-high): no single convention, the values in this form come
    # from more than 500 sources; in the largest ones the value is a mean
    # or a median, the interval is the full range (PMID 24023812, 28959601),
    # the interquartile range (PMID 27294788) or percentiles (PMID 8087979)
    (
        re.compile(_NUM + r'\s*\(\s*' + _RANGE + r'\s*\)'),
        ('value', 'value_min', 'value_max'),
    ),
    (
        re.compile(_NUM + _PM + _NUM + r'\s*\(\s*' + _RANGE + r'\s*\)'),
        ('value', 'value_sd', 'value_min', 'value_max'),
    ),
    (re.compile(r'<=?\s*' + _NUM), ('value_max',)),
    (re.compile(r'>=?\s*' + _NUM), ('value_min',)),
)

# the records yielded by the functions below, each defined once here rather
# than per record so all records of a function share one type
HmdbMetabolite = collections.namedtuple(
    'HmdbMetabolite',
    (
        'hmdb_id',
        'secondary_ids',
        'name',
        'synonyms',
        'status',
        'formula',
        'avg_mw',
        'mono_mw',
        'smiles',
        'inchi',
        'inchikey',
        'iupac_name',
        'cas',
        'kingdom',
        'super_class',
        'class_',
        'sub_class',
        'direct_parent',
        'cellular_locations',
        'biospecimen_locations',
        'tissue_locations',
    ) + tuple(HMDB_METABOLITE_ID_TAGS),
)

HmdbMetaboliteProtein = collections.namedtuple(
    'HmdbMetaboliteProtein',
    (
        'hmdb_id',
        'hmdb_protein_id',
        'uniprot_id',
        'gene_symbol',
        'protein_name',
        'protein_type',
        'pubmed_ids',
    ),
)

HmdbMetaboliteDisease = collections.namedtuple(
    'HmdbMetaboliteDisease',
    (
        'hmdb_id',
        'disease_name',
        'omim_id',
        'pubmed_ids',
        'reference_texts',
    ),
)

HmdbMetabolitePathway = collections.namedtuple(
    'HmdbMetabolitePathway',
    (
        'hmdb_id',
        'pathway_name',
        'smpdb_id',
        'kegg_map_id',
    ),
)

HmdbConcentration = collections.namedtuple(
    'HmdbConcentration',
    (
        'hmdb_id',
        'kind',
        'biospecimen',
        'value_raw',
        'value',
        'value_sd',
        'value_min',
        'value_max',
        'units',
        'age',
        'sex',
        'condition',
        'comment',
        'pubmed_ids',
    ),
)

HmdbProtein = collections.namedtuple(
    'HmdbProtein',
    (
        'hmdb_protein_id',
        'uniprot_id',
        'gene_symbol',
        'hgnc_id',
        'protein_type',
        'general_function',
        'specific_function',
    ),
)


def _hmdb_records(
        dataset: Literal['metabolites', 'proteins'],
    ) -> Generator[etree._Element, None, None]:
    """
    The records of one HMDB XML file, one element at a time.

    HMDB serves its files behind a bot challenge, so the download fails
    when it is not done in a web browser. In that case the zip file has
    to be downloaded manually and saved in the pypath cache, the error
    message tells the exact path. The same message is given if the file
    in the cache is not a zip file, e.g. a web page saved by mistake.

    Each element is cleared after it has been processed, together with the
    ones before it, as the metabolites XML is too large to keep in memory.
    """

    xml_name = f'hmdb_{dataset}.xml'
    url = urls.urls['hmdb'][dataset]
    cache_path = curl.Curl(
        url,
        setup = False,
        call = False,
        process = False,
    ).cache_file_name
    manual = 'download `%s` in a web browser and save it as `%s`' % (
        url,
        cache_path,
    )

    try:

        c = curl.Curl(
            url,
            large = True,
            silent = False,
            slow = True,
            default_mode = 'rb',
            files_needed = [xml_name],
        )

    except zipfile.BadZipFile:

        msg = 'The HMDB file in the cache is not a zip file: %s.' % manual
        _log(msg)
        raise RuntimeError(msg)

    if c.download_failed or not c.result or xml_name not in c.result:

        msg = (
            'Could not open `%s`. HMDB blocks automated downloads: %s.' % (
                xml_name,
                manual,
            )
        )
        _log(msg)
        raise RuntimeError(msg)

    tag = f'{HMDB_XMLNS}{HMDB_RECORD_TAGS[dataset]}'

    for _, record in etree.iterparse(c.result[xml_name], tag = tag):

        # a tag of the same name nested deeper is part of a record, it
        # must not be cleared before its record is processed
        if record.getparent().getparent() is not None:

            continue

        yield record

        record.clear()

        while record.getprevious() is not None:

            del record.getparent()[0]


def _hmdb_path(path: str) -> str:
    """
    An element path with the HMDB namespace added to each tag.
    """

    return '/'.join(f'{HMDB_XMLNS}{tag}' for tag in path.split('/'))


def _hmdb_text(elem: etree._Element, path: str) -> Optional[str]:
    """
    Text of the first element at `path`, `None` if it is missing or empty.
    """

    text = elem.findtext(_hmdb_path(path))

    return (text.strip() or None) if text else None


def _hmdb_texts(elem: etree._Element, path: str) -> set:
    """
    Texts of all elements at `path`, without the empty ones.
    """

    return {
        e.text.strip()
        for e in elem.iterfind(_hmdb_path(path))
        if e.text and e.text.strip()
    }


def _hmdb_float(elem: etree._Element, path: str) -> Optional[float]:
    """
    Number from the element at `path`, `None` if it is missing or not
    a number.
    """

    try:

        return float(_hmdb_text(elem, path))

    except (TypeError, ValueError):

        return None


def _hmdb_references(elem: etree._Element) -> tuple[set, set]:
    """
    PubMed IDs and reference texts from the `references` of an element.

    Only the numeric PubMed IDs are kept. References without PubMed ID
    still contribute their text.
    """

    pubmeds = set()
    texts = set()

    for reference in elem.iterfind(_hmdb_path('references/reference')):

        pubmed = _hmdb_text(reference, 'pubmed_id')
        text = _hmdb_text(reference, 'reference_text')

        if pubmed and pubmed.isdigit():

            pubmeds.add(pubmed)

        if text:

            texts.add(text)

    return pubmeds, texts


def _hmdb_concentration_value(value_raw: Optional[str]) -> dict:
    """
    Numbers from a concentration value.

    The values are free text such as `9.8`, `9.8 +/- 1.2`, `9.8 (1.2)`,
    `4.0-15.0`, `9.8 (4.0-15.0)` or `< 0.5`. Values in none of the forms in
    `HMDB_CONCENTRATION_FORMS` give no numbers, their text is still
    available as `value_raw`.
    """

    values = dict.fromkeys(('value', 'value_sd', 'value_min', 'value_max'))

    if value_raw:

        for form, fields in HMDB_CONCENTRATION_FORMS:

            if match := form.fullmatch(value_raw):

                values.update(zip(fields, map(float, match.groups())))
                break

    return values


def _hmdb_metabolite(record: etree._Element) -> tuple:
    """
    Metabolite record from an element of the metabolites XML.
    """

    locations = 'biological_properties/%s_locations/%s'

    return HmdbMetabolite(
        hmdb_id = _hmdb_text(record, 'accession'),
        secondary_ids = _hmdb_texts(record, 'secondary_accessions/accession'),
        name = _hmdb_text(record, 'name'),
        synonyms = _hmdb_texts(record, 'synonyms/synonym'),
        status = _hmdb_text(record, 'status'),
        formula = _hmdb_text(record, 'chemical_formula'),
        avg_mw = _hmdb_float(record, 'average_molecular_weight'),
        # sic, the tag name is misspelled in HMDB
        mono_mw = _hmdb_float(record, 'monisotopic_molecular_weight'),
        smiles = _hmdb_text(record, 'smiles'),
        inchi = _hmdb_text(record, 'inchi'),
        inchikey = _hmdb_text(record, 'inchikey'),
        iupac_name = _hmdb_text(record, 'iupac_name'),
        cas = _hmdb_text(record, 'cas_registry_number'),
        kingdom = _hmdb_text(record, 'taxonomy/kingdom'),
        super_class = _hmdb_text(record, 'taxonomy/super_class'),
        class_ = _hmdb_text(record, 'taxonomy/class'),
        sub_class = _hmdb_text(record, 'taxonomy/sub_class'),
        direct_parent = _hmdb_text(record, 'taxonomy/direct_parent'),
        cellular_locations = _hmdb_texts(
            record,
            locations % ('cellular', 'cellular'),
        ),
        biospecimen_locations = _hmdb_texts(
            record,
            locations % ('biospecimen', 'biospecimen'),
        ),
        tissue_locations = _hmdb_texts(
            record,
            locations % ('tissue', 'tissue'),
        ),
        **{
            field: _hmdb_text(record, tag)
            for field, tag in HMDB_METABOLITE_ID_TAGS.items()
        },
    )


def hmdb_metabolites() -> Generator[tuple, None, None]:
    """
    Metabolites from HMDB.

    Identity, structure, chemical taxonomy and locations of the metabolites,
    with their identifiers in other resources.

    Yields
        (tuple): `HmdbMetabolite` named tuples, one for each metabolite.
            Fields with more than one value are sets, the molecular weights
            are floats, missing values are `None`.
    """

    for record in _hmdb_records('metabolites'):

        yield _hmdb_metabolite(record)


def hmdb_metabolite_pathway() -> Generator[tuple, None, None]:
    """
    Pathways of the metabolites from HMDB.

    Yields
        (tuple): `HmdbMetabolitePathway` named tuples, one for each
            metabolite-pathway pair.
    """

    for record in _hmdb_records('metabolites'):

        hmdb_id = _hmdb_text(record, 'accession')

        for pathway in record.iterfind(
            _hmdb_path('biological_properties/pathways/pathway')
        ):

            yield HmdbMetabolitePathway(
                hmdb_id = hmdb_id,
                pathway_name = _hmdb_text(pathway, 'name'),
                smpdb_id = _hmdb_text(pathway, 'smpdb_id'),
                kegg_map_id = _hmdb_text(pathway, 'kegg_map_id'),
            )


def hmdb_metabolite_disease() -> Generator[tuple, None, None]:
    """
    Diseases associated with the metabolites in HMDB.

    Yields
        (tuple): `HmdbMetaboliteDisease` named tuples, one for each
            metabolite-disease pair. The references are in two sets: the
            PubMed IDs, and the texts of all references, including the
            ones without PubMed ID.
    """

    for record in _hmdb_records('metabolites'):

        hmdb_id = _hmdb_text(record, 'accession')

        for disease in record.iterfind(_hmdb_path('diseases/disease')):

            pubmeds, texts = _hmdb_references(disease)

            yield HmdbMetaboliteDisease(
                hmdb_id = hmdb_id,
                disease_name = _hmdb_text(disease, 'name'),
                omim_id = _hmdb_text(disease, 'omim_id'),
                pubmed_ids = pubmeds,
                reference_texts = texts,
            )


def hmdb_concentrations(
        kind: Literal['all', 'normal', 'abnormal'] = 'all',
    ) -> Generator[tuple, None, None]:
    """
    Concentrations of the metabolites in biospecimens from HMDB.

    HMDB distinguishes normal concentrations, measured in healthy subjects,
    and abnormal ones, measured in patients. The original value is kept in
    `value_raw`, the numbers are taken from it into `value`, `value_sd`,
    `value_min` and `value_max` where its form is recognized.

    The values are collected from many sources, and their meaning follows
    the source: the units differ between records, so the numbers are
    comparable only within the same `units`; `value` may be a mean or a
    median; the spread after `+/-` may be a standard deviation or a
    standard error; and `value_min` and `value_max` are the bounds of the
    interval as given, which may be the lowest and highest value measured,
    the interquartile range or other percentiles. Where HMDB records the
    statistics used, it is in `comment`, e.g. `Median (range)`.

    Args
        kind: Either `normal`, `abnormal`, or `all` for both.

    Yields
        (tuple): `HmdbConcentration` named tuples, one for each
            concentration. For the normal ones `condition` is the condition
            of the subject, for the abnormal ones the patient information.
    """

    if kind == 'all':

        kinds = tuple(HMDB_CONCENTRATION_TAGS)

    elif kind in HMDB_CONCENTRATION_TAGS:

        kinds = (kind,)

    else:

        raise ValueError(
            'Unknown concentration kind: `%s`. Use one of `all`, `%s`.' % (
                kind,
                '`, `'.join(HMDB_CONCENTRATION_TAGS),
            )
        )

    for record in _hmdb_records('metabolites'):

        hmdb_id = _hmdb_text(record, 'accession')

        for _kind in kinds:

            container, age, sex, condition = HMDB_CONCENTRATION_TAGS[_kind]

            for conc in record.iterfind(
                _hmdb_path(f'{container}/concentration')
            ):

                value_raw = _hmdb_text(conc, 'concentration_value')
                pubmeds, _ = _hmdb_references(conc)

                yield HmdbConcentration(
                    hmdb_id = hmdb_id,
                    kind = _kind,
                    biospecimen = _hmdb_text(conc, 'biospecimen'),
                    value_raw = value_raw,
                    units = _hmdb_text(conc, 'concentration_units'),
                    age = _hmdb_text(conc, age),
                    sex = _hmdb_text(conc, sex),
                    condition = _hmdb_text(conc, condition),
                    comment = _hmdb_text(conc, 'comment'),
                    pubmed_ids = pubmeds,
                    **_hmdb_concentration_value(value_raw),
                )


def hmdb_proteins() -> Generator[tuple, None, None]:
    """
    Proteins from HMDB.

    Yields
        (tuple): `HmdbProtein` named tuples, one for each protein.
    """

    for record in _hmdb_records('proteins'):

        yield HmdbProtein(
            hmdb_protein_id = _hmdb_text(record, 'accession'),
            uniprot_id = _hmdb_text(record, 'uniprot_id'),
            gene_symbol = _hmdb_text(record, 'gene_name'),
            hgnc_id = _hmdb_text(record, 'hgnc_id'),
            protein_type = _hmdb_text(record, 'protein_type'),
            general_function = _hmdb_text(record, 'general_function'),
            specific_function = _hmdb_text(record, 'specific_function'),
        )


def hmdb_metabolite_protein() -> Generator[tuple, None, None]:
    """
    Proteins associated with the metabolites in HMDB.

    The associations are listed in the metabolites XML, while the
    references supporting them are in the proteins XML, hence both files
    are processed: first the proteins, to collect the references.

    Yields
        (tuple): `HmdbMetaboliteProtein` named tuples, one for each
            metabolite-protein pair. `pubmed_ids` is empty if the proteins
            XML has no reference for the pair.
    """

    pubmeds = collections.defaultdict(set)

    for record in _hmdb_records('proteins'):

        protein_id = _hmdb_text(record, 'accession')

        for ref in record.iterfind(
            _hmdb_path('metabolite_references/metabolite_reference')
        ):

            hmdb_id = _hmdb_text(ref, 'metabolite/accession')
            pubmed = _hmdb_text(ref, 'reference/pubmed_id')

            if hmdb_id and pubmed and pubmed.isdigit():

                pubmeds[(hmdb_id, protein_id)].add(pubmed)

    for record in _hmdb_records('metabolites'):

        hmdb_id = _hmdb_text(record, 'accession')

        for protein in record.iterfind(
            _hmdb_path('protein_associations/protein')
        ):

            protein_id = _hmdb_text(protein, 'protein_accession')

            yield HmdbMetaboliteProtein(
                hmdb_id = hmdb_id,
                hmdb_protein_id = protein_id,
                uniprot_id = _hmdb_text(protein, 'uniprot_id'),
                gene_symbol = _hmdb_text(protein, 'gene_name'),
                protein_name = _hmdb_text(protein, 'name'),
                protein_type = _hmdb_text(protein, 'protein_type'),
                pubmed_ids = set(pubmeds.get((hmdb_id, protein_id), ())),
            )


def _hmdb_values(value) -> set:
    """
    The values of a record field as a set of strings.
    """

    if value is None:

        return set()

    if isinstance(value, set):

        return value

    return {str(value)}


def hmdb_mapping(id_a: str, id_b: str) -> dict[str, set[str]]:
    """
    Translation table between two fields of the HMDB records.

    Both fields must be in the same kind of record: either fields of
    `HmdbMetabolite`, e.g. `hmdb_id`, `chebi`, `kegg` or `inchikey`, or
    fields of `HmdbProtein`, e.g. `hmdb_protein_id`, `uniprot_id` or
    `gene_symbol`. For fields with more than one value, such as
    `secondary_ids` or `synonyms`, each value becomes a separate key.

    Args
        id_a: The field to translate from.
        id_b: The field to translate to.

    Returns
        (dict): Values of `id_a` as keys, sets of the values of `id_b` as
            values. Records without value in either field are left out.
    """

    for records, record_type in (
        (hmdb_metabolites, HmdbMetabolite),
        (hmdb_proteins, HmdbProtein),
    ):

        if id_a in record_type._fields and id_b in record_type._fields:

            break

    else:

        raise ValueError(
            'No HMDB record has both `%s` and `%s`. Use two fields of '
            '`HmdbMetabolite` or two fields of `HmdbProtein`.' % (id_a, id_b)
        )

    mapping = collections.defaultdict(set)

    for record in records():

        values_b = _hmdb_values(getattr(record, id_b))

        if values_b:

            for value_a in _hmdb_values(getattr(record, id_a)):

                mapping[value_a].update(values_b)

    return dict(mapping)
