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

from typing import Generator, Literal, Optional, Union

import os
import zipfile
import collections

import lxml.etree as etree

import pypath.resources.urls as urls
import pypath.share.common as common
import pypath.share.curl as curl
import pypath.share.session as session
import pypath.share.settings as settings

_logger = session.Logger(name = 'hmdb_v2_input')
_log = _logger._log

HMDB_XMLNS = '{http://www.hmdb.ca}'

HMDB_RECORD_TAGS = {
    'metabolites': 'metabolite',
    'proteins': 'protein',
}

# the status of each metabolite: the first two have been found in human
# samples, `expected` ones are expected but not yet found in humans,
# `predicted` ones are only predicted
HMDB_STATUSES = ('quantified', 'detected', 'expected', 'predicted')

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

# fields of the metabolite-protein records describing the protein: field
# names in the records and tag names, the same in both XMLs
HMDB_PAIR_PROTEIN_TAGS = {
    'uniprot_id': 'uniprot_id',
    'gene_symbol': 'gene_name',
    'protein_name': 'name',
    'protein_type': 'protein_type',
}

# the records yielded by the functions below, each defined once here rather
# than per record so all records of a function share one type
HmdbMetabolite = collections.namedtuple(
    'HmdbMetabolite',
    (
        'hmdb_id',
        'secondary_ids',
        'name',
        'synonyms',
        'description',
        'status',
        'formula',
        'avg_mw',
        'mono_mw',
        'state',
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
        'value',
        'units',
        'age',
        'sex',
        'condition',
        'comment',
        'pubmed_ids',
    ),
)

HmdbMetaboliteOntology = collections.namedtuple(
    'HmdbMetaboliteOntology',
    (
        'hmdb_id',
        'root',
        'term',
        'level',
        'parent_term',
    ),
)

HmdbExperimentalProperty = collections.namedtuple(
    'HmdbExperimentalProperty',
    (
        'hmdb_id',
        'kind',
        'value',
        'source',
    ),
)

HmdbProtein = collections.namedtuple(
    'HmdbProtein',
    (
        'hmdb_protein_id',
        'name',
        'uniprot_id',
        'gene_symbol',
        'hgnc_id',
        'protein_type',
        'general_function',
        'specific_function',
    ),
)

HmdbProteinPathway = collections.namedtuple(
    'HmdbProteinPathway',
    (
        'hmdb_protein_id',
        'pathway_name',
        'smpdb_id',
        'kegg_map_id',
    ),
)


def _hmdb_records(
        dataset: Literal['metabolites', 'proteins'],
    ) -> Generator[etree._Element, None, None]:
    """
    The records of one HMDB XML file, one element at a time.

    HMDB serves its files behind a bot challenge, so the download fails
    when it is not done in a web browser. The zip file then has to be
    downloaded manually, and either its path set in the
    `hmdb_metabolites_path` or `hmdb_proteins_path` setting, or the file
    saved in the pypath cache; the error message tells the exact path.
    The same message is given if the file is not a zip file, e.g. a web
    page saved by mistake.

    Each element is cleared after it has been processed, together with the
    ones before it, as the metabolites XML is too large to keep in memory.
    """

    xml_name = f'hmdb_{dataset}.xml'
    url = urls.urls['hmdb'][dataset]
    setting = f'hmdb_{dataset}_path'
    path = settings.get(setting)
    cache_path = curl.Curl(
        url,
        setup = False,
        call = False,
        process = False,
    ).cache_file_name
    manual = (
        'download `%s` in a web browser, then either set its path in the '
        '`%s` setting or save it as `%s`' % (url, setting, cache_path)
    )

    # a path that does not exist would be taken for a URL by `curl.Curl`
    if path and not os.path.isfile(path):

        msg = 'The file in the `%s` setting does not exist: `%s`.' % (
            setting,
            path,
        )
        _log(msg)
        raise RuntimeError(msg)

    try:

        c = curl.Curl(
            path or url,
            large = True,
            silent = False,
            slow = True,
            default_mode = 'rb',
            files_needed = [xml_name],
        )

    except zipfile.BadZipFile:

        msg = 'The HMDB file `%s` is not a zip file: %s.' % (
            path or cache_path,
            manual,
        )
        _log(msg)
        raise RuntimeError(msg)

    if c.download_failed or not c.result or xml_name not in c.result:

        msg = (
            'No `%s` in the file in the `%s` setting, `%s`: %s.' % (
                xml_name,
                setting,
                path,
                manual,
            )
                if path else
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


def _hmdb_status(status: Optional[Union[str, set]]) -> set:
    """
    The metabolite statuses to keep, as a set; empty if all are kept.

    Raises on a status not in `HMDB_STATUSES`, as a misspelled one would
    silently filter out every record.
    """

    status = common.to_set(status)
    unknown = status - set(HMDB_STATUSES)

    if unknown:

        raise ValueError(
            'Unknown metabolite status: `%s`. Use any of `%s`.' % (
                '`, `'.join(sorted(unknown)),
                '`, `'.join(HMDB_STATUSES),
            )
        )

    return status


def _hmdb_metabolite_records(
        status: set,
    ) -> Generator[etree._Element, None, None]:
    """
    The records of the metabolites XML, only the ones with a status in
    `status` if it is not empty.
    """

    for record in _hmdb_records('metabolites'):

        if not status or _hmdb_text(record, 'status') in status:

            yield record


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
        description = _hmdb_text(record, 'description'),
        status = _hmdb_text(record, 'status'),
        formula = _hmdb_text(record, 'chemical_formula'),
        avg_mw = _hmdb_float(record, 'average_molecular_weight'),
        # sic, the tag name is misspelled in HMDB
        mono_mw = _hmdb_float(record, 'monisotopic_molecular_weight'),
        state = _hmdb_text(record, 'state'),
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


def hmdb_metabolites(
        status: Optional[Union[str, set]] = None,
    ) -> Generator[tuple, None, None]:
    """
    Metabolites from HMDB.

    Identity, description, structure, physical state, chemical taxonomy and
    locations of the metabolites, with their identifiers in other resources.

    Args
        status: Only the metabolites with these statuses, either one or a
            set of `quantified`, `detected`, `expected` and `predicted`,
            e.g. `{'quantified', 'detected'}` for the ones found in human
            samples. By default all metabolites are included.

    Yields
        (tuple): `HmdbMetabolite` named tuples, one for each metabolite.
            Fields with more than one value are sets, the molecular weights
            are floats, missing values are `None`.
    """

    status = _hmdb_status(status)

    for record in _hmdb_metabolite_records(status):

        yield _hmdb_metabolite(record)


def hmdb_metabolite_pathway(
        status: Optional[Union[str, set]] = None,
    ) -> Generator[tuple, None, None]:
    """
    Pathways of the metabolites from HMDB.

    Args
        status: Only the metabolites with these statuses, either one or a
            set of `quantified`, `detected`, `expected` and `predicted`,
            e.g. `{'quantified', 'detected'}` for the ones found in human
            samples. By default all metabolites are included.

    Yields
        (tuple): `HmdbMetabolitePathway` named tuples, one for each
            metabolite-pathway pair.
    """

    status = _hmdb_status(status)

    for record in _hmdb_metabolite_records(status):

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


def hmdb_metabolite_disease(
        status: Optional[Union[str, set]] = None,
    ) -> Generator[tuple, None, None]:
    """
    Diseases associated with the metabolites in HMDB.

    Args
        status: Only the metabolites with these statuses, either one or a
            set of `quantified`, `detected`, `expected` and `predicted`,
            e.g. `{'quantified', 'detected'}` for the ones found in human
            samples. By default all metabolites are included.

    Yields
        (tuple): `HmdbMetaboliteDisease` named tuples, one for each
            metabolite-disease pair. The references are in two sets: the
            PubMed IDs, and the texts of all references, including the
            ones without PubMed ID.
    """

    status = _hmdb_status(status)

    for record in _hmdb_metabolite_records(status):

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
        status: Optional[Union[str, set]] = None,
    ) -> Generator[tuple, None, None]:
    """
    Concentrations of the metabolites in biospecimens from HMDB.

    HMDB distinguishes normal concentrations, measured in healthy subjects,
    and abnormal ones, measured in patients.

    `value` is the concentration as HMDB gives it, as text: a single
    number, or one with a spread or an interval, e.g. `9.8`,
    `7.7 +/- 1.9`, `50.65(14.14)`, `8.3 (2.4-28.4)`, `10-17` or `<500`.
    It is not split into numbers, as their meaning follows the source the
    value was collected from:

    - in `x (y)`, `y` is the standard deviation: all values in this form
      come from seven papers, and the tables behind 92% of them (PMID
      24704061, 26010610) give mean (SD), as does PMID 23434268;
    - in `x (a-b)`, `x` is a mean or a median, and the interval is the
      full range (e.g. PMID 24023812, 28959601), the interquartile range
      (PMID 27294788) or percentiles (PMID 8087979), depending on the
      source;
    - the spread after `+/-` may be a standard deviation or a standard
      error.

    Where HMDB records the statistics used, it is in `comment`, e.g.
    `Median (range)`. The units differ between records, see `units`.

    Args
        kind: Either `normal`, `abnormal`, or `all` for both.
        status: Only the metabolites with these statuses, either one or a
            set of `quantified`, `detected`, `expected` and `predicted`,
            e.g. `{'quantified', 'detected'}` for the ones found in human
            samples. By default all metabolites are included.

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

    status = _hmdb_status(status)

    for record in _hmdb_metabolite_records(status):

        hmdb_id = _hmdb_text(record, 'accession')

        for _kind in kinds:

            container, age, sex, condition = HMDB_CONCENTRATION_TAGS[_kind]

            for conc in record.iterfind(
                _hmdb_path(f'{container}/concentration')
            ):

                pubmeds, _ = _hmdb_references(conc)

                yield HmdbConcentration(
                    hmdb_id = hmdb_id,
                    kind = _kind,
                    biospecimen = _hmdb_text(conc, 'biospecimen'),
                    value = _hmdb_text(conc, 'concentration_value'),
                    units = _hmdb_text(conc, 'concentration_units'),
                    age = _hmdb_text(conc, age),
                    sex = _hmdb_text(conc, sex),
                    condition = _hmdb_text(conc, condition),
                    comment = _hmdb_text(conc, 'comment'),
                    pubmed_ids = pubmeds,
                )


def _hmdb_ontology(
        hmdb_id: str,
        elem: etree._Element,
        root: Optional[str] = None,
        parent_term: Optional[str] = None,
    ) -> Generator[tuple, None, None]:
    """
    Ontology records of a metabolite from a term and all terms below it.
    """

    term = _hmdb_text(elem, 'term')
    root = root or term
    level = _hmdb_text(elem, 'level')

    yield HmdbMetaboliteOntology(
        hmdb_id = hmdb_id,
        root = root,
        term = term,
        level = int(level) if level and level.isdigit() else None,
        parent_term = parent_term,
    )

    for descendant in elem.iterfind(_hmdb_path('descendants/descendant')):

        yield from _hmdb_ontology(hmdb_id, descendant, root, term)


def hmdb_metabolite_ontology(
        status: Optional[Union[str, set]] = None,
    ) -> Generator[tuple, None, None]:
    """
    ChemFOnt classification of the metabolites from HMDB.

    The `ontology` of a metabolite is a tree with up to four roots:
    `Disposition` (source, e.g. endogenous, food, or the organism it comes
    from, biological location and route of exposure), `Role` (e.g. drug,
    nutrient, toxin), `Process` and `Physiological effect`. Each term of
    the tree is one record, the roots and the intermediate terms included,
    so the higher levels can be used for filtering. In HMDB 5.0 a term has
    the same parent under every metabolite, hence a term is identified by
    its name.

    The definitions and synonyms of the terms are not included, as they
    describe the terms, not the metabolites.

    Args
        status: Only the metabolites with these statuses, either one or a
            set of `quantified`, `detected`, `expected` and `predicted`,
            e.g. `{'quantified', 'detected'}` for the ones found in human
            samples. By default all metabolites are included.

    Yields
        (tuple): `HmdbMetaboliteOntology` named tuples, one for each term
            of each metabolite. `root` is the root term the term belongs
            to, `level` is 1 for the roots, `parent_term` is `None` for the
            roots.
    """

    status = _hmdb_status(status)

    for record in _hmdb_metabolite_records(status):

        hmdb_id = _hmdb_text(record, 'accession')

        for root in record.iterfind(_hmdb_path('ontology/root')):

            yield from _hmdb_ontology(hmdb_id, root)


def hmdb_experimental_properties(
        status: Optional[Union[str, set]] = None,
    ) -> Generator[tuple, None, None]:
    """
    Experimental physico-chemical properties of the metabolites from HMDB.

    Four kinds of properties: `logp`, `water_solubility`, `melting_point`
    and `boiling_point`. As for the concentrations, `value` is the text as
    HMDB gives it, with its unit and conditions, e.g. `200 g/kg`,
    `187.0 - 189.5 °C` or `139.80 °C. @  760.00 mm Hg`; `logp` has no unit.
    Not all values are measured: some are marked `(est)` in `value`, or
    come from the source `Extrapolated`. The predicted properties of the
    metabolites XML are not included.

    Args
        status: Only the metabolites with these statuses, either one or a
            set of `quantified`, `detected`, `expected` and `predicted`,
            e.g. `{'quantified', 'detected'}` for the ones found in human
            samples. By default all metabolites are included.

    Yields
        (tuple): `HmdbExperimentalProperty` named tuples, one for each
            property of each metabolite. `source` is a reference, a
            database or a URL, `None` if HMDB gives none.
    """

    status = _hmdb_status(status)

    for record in _hmdb_metabolite_records(status):

        hmdb_id = _hmdb_text(record, 'accession')

        for prop in record.iterfind(
            _hmdb_path('experimental_properties/property')
        ):

            yield HmdbExperimentalProperty(
                hmdb_id = hmdb_id,
                kind = _hmdb_text(prop, 'kind'),
                value = _hmdb_text(prop, 'value'),
                source = _hmdb_text(prop, 'source'),
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
            name = _hmdb_text(record, 'name'),
            uniprot_id = _hmdb_text(record, 'uniprot_id'),
            gene_symbol = _hmdb_text(record, 'gene_name'),
            hgnc_id = _hmdb_text(record, 'hgnc_id'),
            protein_type = _hmdb_text(record, 'protein_type'),
            general_function = _hmdb_text(record, 'general_function'),
            specific_function = _hmdb_text(record, 'specific_function'),
        )


def hmdb_protein_pathway() -> Generator[tuple, None, None]:
    """
    Pathways of the proteins from HMDB.

    Yields
        (tuple): `HmdbProteinPathway` named tuples, one for each
            protein-pathway pair.
    """

    for record in _hmdb_records('proteins'):

        protein_id = _hmdb_text(record, 'accession')

        for pathway in record.iterfind(_hmdb_path('pathways/pathway')):

            yield HmdbProteinPathway(
                hmdb_protein_id = protein_id,
                pathway_name = _hmdb_text(pathway, 'name'),
                smpdb_id = _hmdb_text(pathway, 'smpdb_id'),
                kegg_map_id = _hmdb_text(pathway, 'kegg_map_id'),
            )


def _hmdb_pair_protein(elem: etree._Element) -> dict:
    """
    The protein fields of a metabolite-protein record, from a protein in
    either XML.
    """

    return {
        field: _hmdb_text(elem, tag)
        for field, tag in HMDB_PAIR_PROTEIN_TAGS.items()
    }


def hmdb_metabolite_protein(
        status: Optional[Union[str, set]] = None,
    ) -> Generator[tuple, None, None]:
    """
    Proteins associated with the metabolites in HMDB.

    Both XML files list the associations: the metabolites XML under each
    metabolite (`protein_associations`), the proteins XML under each
    protein (`metabolite_associations`). The proteins XML has about twice
    as many, most of the additional ones are of predicted metabolites. The
    references supporting the associations are only in the proteins XML.

    Hence both files are processed: first the proteins, to collect their
    associations, references and fields, then the metabolites. For each
    metabolite its pairs from the metabolites XML are yielded first, then
    the ones only in the proteins XML, with the protein fields taken from
    the protein record. The status of a metabolite is only known from the
    metabolites XML, so a pair of a metabolite missing from it would be
    left out; in HMDB 5.0 there is none.

    Args
        status: Only the metabolites with these statuses, either one or a
            set of `quantified`, `detected`, `expected` and `predicted`,
            e.g. `{'quantified', 'detected'}` for the ones found in human
            samples. By default all metabolites are included.

    Yields
        (tuple): `HmdbMetaboliteProtein` named tuples, one for each
            metabolite-protein pair. `pubmed_ids` is empty if the proteins
            XML has no reference for the pair.
    """

    # checked before the proteins XML is read, not after
    status = _hmdb_status(status)
    pubmeds = collections.defaultdict(set)
    proteins = {}
    associations = collections.defaultdict(list)

    for record in _hmdb_records('proteins'):

        protein_id = _hmdb_text(record, 'accession')
        proteins[protein_id] = _hmdb_pair_protein(record)

        for metabolite in record.iterfind(
            _hmdb_path('metabolite_associations/metabolite')
        ):

            hmdb_id = _hmdb_text(metabolite, 'accession')

            if hmdb_id:

                associations[hmdb_id].append(protein_id)

        for ref in record.iterfind(
            _hmdb_path('metabolite_references/metabolite_reference')
        ):

            hmdb_id = _hmdb_text(ref, 'metabolite/accession')
            pubmed = _hmdb_text(ref, 'reference/pubmed_id')

            if hmdb_id and pubmed and pubmed.isdigit():

                pubmeds[(hmdb_id, protein_id)].add(pubmed)

    for record in _hmdb_metabolite_records(status):

        hmdb_id = _hmdb_text(record, 'accession')
        pairs = [
            (
                _hmdb_text(protein, 'protein_accession'),
                _hmdb_pair_protein(protein),
            )
            for protein in record.iterfind(
                _hmdb_path('protein_associations/protein')
            )
        ]
        listed = {protein_id for protein_id, _ in pairs}
        pairs.extend(
            (protein_id, proteins[protein_id])
            for protein_id in associations.pop(hmdb_id, ())
            if protein_id not in listed
        )

        for protein_id, protein in pairs:

            yield HmdbMetaboliteProtein(
                hmdb_id = hmdb_id,
                hmdb_protein_id = protein_id,
                pubmed_ids = set(pubmeds.get((hmdb_id, protein_id), ())),
                **protein,
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
