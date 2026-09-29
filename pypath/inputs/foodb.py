"""
FooDB (https://foodb.ca) scraper - tek dosyada.
import csv
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Optional

import scrapy
from pydantic import BaseModel, field_validator
from scrapy.crawler import CrawlerProcess
from tqdm import tqdm

OUTPUT_DIR = Path("outputs")
USE_CACHE = True


def cell_text(row, td_index):
    parts = row.xpath(f"./td[{td_index}]//text()").getall()
    return " ".join(p.strip() for p in parts if p.strip())

class FoodListModel(BaseModel):
    food_id: str
    name: str
    scientific_name: Optional[str] = None

    @field_validator("scientific_name", mode="before")
    @classmethod
    def blank_to_none(cls, v):
        return None if not v or v == "Not Available" else v

    @field_validator("food_id")
    @classmethod
    def validate_food_id(cls, v: str) -> str:
        if not re.match(r"^FOOD\d+$", v):
            raise ValueError(f"Invalid food_id format: {v}")
        return v

    @field_validator("name")
    @classmethod
    def validate_name(cls, v: str) -> str:
        if not v:
            raise ValueError("Name cannot be empty.")
        return v


class FoodDetailModel(BaseModel):
    primary_id: str
    name: str
    scientific_name: Optional[str] = None
    compounds: Optional[list[str]] = None

    @field_validator("scientific_name", mode="before")
    @classmethod
    def blank_to_none(cls, v):
        return None if not v or v == "Not Available" else v

    @field_validator("compounds", mode="before")
    @classmethod
    def empty_list_to_none(cls, v):
        return None if not v else v


class CompoundListModel(BaseModel):
    foodb_id: str
    name: str
    cas_number: Optional[str] = None
    foods: Optional[list[str]] = None

    @field_validator("cas_number", mode="before")
    @classmethod
    def blank_string_to_none(cls, v):
        return None if not v else v

    @field_validator("foods", mode="before")
    @classmethod
    def empty_list_to_none(cls, v):
        return None if not v else v

    @field_validator("foodb_id")
    @classmethod
    def validate_foodb_id(cls, v: str) -> str:
        if not re.match(r"^FDB\d+$", v):
            raise ValueError(f"Invalid foodb_id format: {v}")
        return v

    @field_validator("name")
    @classmethod
    def validate_name(cls, v: str) -> str:
        if not v:
            raise ValueError("Name cannot be empty.")
        return v


class CompoundDetailModel(BaseModel):
    primary_id: str
    foodb_name: str
    description: Optional[str] = None
    cas_number: Optional[str] = None
    synonyms: Optional[list[str]] = None
    chemical_formula: Optional[str] = None
    iupac_name: Optional[str] = None
    inchi_identifier: Optional[str] = None
    inchi_key: Optional[str] = None
    isomeric_smiles: Optional[str] = None
    average_molecular_weight: Optional[str] = None
    monoisotopic_molecular_weight: Optional[str] = None
    chembl_id: Optional[str] = None
    kegg_compound_id: Optional[str] = None
    pubchem_compound_id: Optional[str] = None
    pubchem_substance_id: Optional[str] = None
    chebi_id: Optional[str] = None
    drugbank_id: Optional[str] = None
    hmdb_id: Optional[str] = None
    associated_foods: Optional[list[str]] = None

    @field_validator(
        "description", "cas_number", "chemical_formula", "iupac_name",
        "inchi_identifier", "inchi_key", "isomeric_smiles",
        "average_molecular_weight", "monoisotopic_molecular_weight",
        "chembl_id", "kegg_compound_id", "pubchem_compound_id",
        "pubchem_substance_id", "chebi_id", "drugbank_id", "hmdb_id",
        mode="before",
    )
    @classmethod
    def not_available_to_none(cls, v):
        return None if not v or v == "Not Available" else v

    @field_validator("synonyms", "associated_foods", mode="before")
    @classmethod
    def empty_list_to_none(cls, v):
        return None if not v else v

    @field_validator("foodb_name")
    @classmethod
    def validate_name(cls, v: str) -> str:
        if not v:
            raise ValueError("FooDB name cannot be empty.")
        return v

class FoodsSpider(scrapy.Spider):
    name = "foods"
    allowed_domains = ["foodb.ca"]
    start_urls = ["https://foodb.ca/foods"]

    custom_settings = {
        "FEEDS": {
            str(OUTPUT_DIR / "foods_list.csv"): {
                "format": "csv",
                "encoding": "utf8",
                "item_classes": [FoodListModel],
            },
        }
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pbar = None

    def parse(self, response):
        last_page_href = response.css(
            ".pagination.pagination-sm .last.next a::attr(href)"
        ).get()
        if not last_page_href:
            raise ValueError("Could not find the last page link.")

        match = re.search(r"page=(\d+)", last_page_href)
        if not match:
            raise ValueError("Could not extract the last page number from the link.")

        last_page_number = int(match.group(1))
        self.pbar = tqdm(total=last_page_number, desc="Scraping food pages", unit="page")

        for page in range(1, last_page_number + 1):
            url = (
                "https://foodb.ca/foods"
                if page == 1
                else f"https://foodb.ca/foods?page={page}"
            )
            yield scrapy.Request(url, callback=self.parse_foods)

    def parse_foods(self, response):
        rows = response.css(".table-standard.table-condensed tbody tr")
        for row in rows:
            food_id = cell_text(row, 1)
            name = cell_text(row, 2)
            scientific_name = cell_text(row, 3)
            if not food_id or not name:
                self.logger.warning(
                    "Bos food satiri atlandi (url=%s, food_id=%r, name=%r)",
                    response.url, food_id, name,
                )
                continue

            yield FoodListModel(
                food_id=food_id,
                name=name,
                scientific_name=scientific_name,
            )

        if self.pbar:
            self.pbar.update(1)

    def closed(self, reason):
        if self.pbar:
            self.pbar.close()

class CompoundsSpider(scrapy.Spider):
    name = "compounds"
    allowed_domains = ["foodb.ca"]
    start_urls = ["https://foodb.ca/compounds"]

    custom_settings = {
        "FEEDS": {
            str(OUTPUT_DIR / "compounds_list.csv"): {
                "format": "csv",
                "encoding": "utf8",
                "item_classes": [CompoundListModel],
            },
            str(OUTPUT_DIR / "compounds_detail.jsonl"): {
                "format": "jsonlines",
                "encoding": "utf8",
                "item_classes": [CompoundDetailModel],
            },
        }
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pbar = None

    def parse(self, response):
        last_page_href = response.css(
            ".pagination.pagination-sm .last.next a::attr(href)"
        ).get()
        if not last_page_href:
            raise ValueError("Could not find the last page link.")

        match = re.search(r"page=(\d+)", last_page_href)
        if not match:
            raise ValueError("Could not extract the last page number from the link.")

        last_page_number = int(match.group(1))
        self.pbar = tqdm(total=last_page_number, desc="Scraping compound pages", unit="page")

        for page in range(1, last_page_number + 1):
            url = (
                "https://foodb.ca/compounds"
                if page == 1
                else f"https://foodb.ca/compounds?page={page}"
            )
            yield scrapy.Request(url, callback=self.parse_compounds)

    def parse_compounds(self, response):
        rows = response.css(".table-standard.table-condensed tbody tr")
        for row in rows:
            foodb_id = cell_text(row, 1)
            name = cell_text(row, 2)
            cas_number = cell_text(row, 4)
            foods = row.xpath(
                "./td[6]/div[@class='expandable-list']/ul/li/a/text()"
            ).getall()
            if not foods:
                foods = row.xpath("./td[6]//a/text()").getall()

            if not foodb_id or not name:
                self.logger.warning(
                    "Bos compound satiri atlandi (url=%s, foodb_id=%r, name=%r)",
                    response.url, foodb_id, name,
                )
                continue

            yield CompoundListModel(
                foodb_id=foodb_id,
                name=name,
                cas_number=cas_number,
                foods=foods,
            )

            detail_url = f"https://foodb.ca/compounds/{foodb_id}.xml"
            yield scrapy.Request(
                detail_url,
                callback=self.parse_compound_detail,
                cb_kwargs={"foodb_id": foodb_id},
            )

        if self.pbar:
            self.pbar.update(1)

    def parse_compound_detail(self, response, foodb_id):
        xml_selector = scrapy.Selector(text=response.text, type="xml")
        xml_selector.remove_namespaces()

        def first(*xpaths):
            for xp in xpaths:
                val = xml_selector.xpath(xp).get()
                if val:
                    val = val.strip()
                    if val:
                        return val
            return None

        primary_id = first("/compound/accession/text()") or foodb_id
        foodb_name = first("/compound/name/text()") or ""
        description = first("/compound/description/text()")
        cas_number = first("/compound/cas_registry_number/text()")
        synonyms = [
            s.strip()
            for s in xml_selector.xpath("/compound/synonyms/synonym/text()").getall()
            if s.strip()
        ]
        chemical_formula = first("/compound/chemical_formula/text()")
        iupac_name = first("/compound/iupac_name/text()")
        inchi_identifier = first("/compound/inchi/text()")
        inchi_key = first("/compound/inchikey/text()")
        isomeric_smiles = first("/compound/smiles/text()")
        average_molecular_weight = first("/compound/average_molecular_weight/text()")
        monoisotopic_molecular_weight = first(
            "/compound/monisotopic_moleculate_weight/text()"
        )
        chembl_id = first("/compound/chembl_id/text()")
        kegg_compound_id = first("/compound/kegg_id/text()")
        pubchem_compound_id = first("/compound/pubchem_compound_id/text()")
        pubchem_substance_id = first("/compound/pubchem_substance_id/text()")
        chebi_id = first("/compound/chebi_id/text()")
        drugbank_id = first("/compound/drugbank_id/text()")
        hmdb_id = first("/compound/hmdb_id/text()")
        associated_foods = [
            f.strip()
            for f in xml_selector.xpath("/compound/foods/food/name/text()").getall()
            if f.strip()
        ]

        yield CompoundDetailModel(
            primary_id=primary_id,
            foodb_name=foodb_name,
            description=description,
            cas_number=cas_number,
            synonyms=synonyms,
            chemical_formula=chemical_formula,
            iupac_name=iupac_name,
            inchi_identifier=inchi_identifier,
            inchi_key=inchi_key,
            isomeric_smiles=isomeric_smiles,
            average_molecular_weight=average_molecular_weight,
            monoisotopic_molecular_weight=monoisotopic_molecular_weight,
            chembl_id=chembl_id,
            kegg_compound_id=kegg_compound_id,
            pubchem_compound_id=pubchem_compound_id,
            pubchem_substance_id=pubchem_substance_id,
            chebi_id=chebi_id,
            drugbank_id=drugbank_id,
            hmdb_id=hmdb_id,
            associated_foods=associated_foods,
        )

    def closed(self, reason):
        if self.pbar:
            self.pbar.close()

def build_food_compounds():
    foods_list_path = OUTPUT_DIR / "foods_list.csv"
    compounds_detail_path = OUTPUT_DIR / "compounds_detail.jsonl"
    output_path = OUTPUT_DIR / "foods_detail.jsonl"

    if not foods_list_path.exists() or not compounds_detail_path.exists():
        print(
            "UYARI: foods_list.csv veya compounds_detail.jsonl bulunamadi, "
            "foods_detail.jsonl olusturulamadi."
        )
        return

    food_to_compounds = defaultdict(list)
    with compounds_detail_path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            compound_name = record.get("foodb_name")
            associated_foods = record.get("associated_foods") or []
            if not compound_name:
                continue
            for food_name in associated_foods:
                food_to_compounds[food_name].append(compound_name)

    written = 0
    with foods_list_path.open(encoding="utf-8", newline="") as in_fh, \
         output_path.open("w", encoding="utf-8") as out_fh:

        reader = csv.DictReader(in_fh)
        for row in reader:
            food_id = row.get("food_id", "").strip()
            name = row.get("name", "").strip()
            scientific_name = row.get("scientific_name") or None
            compounds = food_to_compounds.get(name)

            detail = FoodDetailModel(
                primary_id=food_id,
                name=name,
                scientific_name=scientific_name,
                compounds=compounds,
            )
            out_fh.write(detail.model_dump_json() + "\n")
            written += 1

    print(f"{written} food kaydi yazildi: {output_path}")


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(exist_ok=True)

    process = CrawlerProcess(
        settings={
            "ROBOTSTXT_OBEY": True,
            "CONCURRENT_REQUESTS": 32,
            "CONCURRENT_REQUESTS_PER_DOMAIN": 16,
            "AUTOTHROTTLE_ENABLED": True,
            "AUTOTHROTTLE_START_DELAY": 0.1,
            "AUTOTHROTTLE_MAX_DELAY": 10,
            "AUTOTHROTTLE_TARGET_CONCURRENCY": 8.0,
            "RETRY_ENABLED": True,
            "RETRY_TIMES": 3,
            "USER_AGENT": (
                "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
            ),
            "FEED_EXPORT_ENCODING": "utf-8",
            "LOG_LEVEL": "INFO",
            "HTTPCACHE_ENABLED": USE_CACHE,
            "HTTPCACHE_DIR": str(OUTPUT_DIR / ".scrapy_cache"),
            "HTTPCACHE_EXPIRATION_SECS": 0,
            "HTTPCACHE_STORAGE": "scrapy.extensions.httpcache.FilesystemCacheStorage",
        }
    )
    process.crawl(FoodsSpider)
    process.crawl(CompoundsSpider)
    process.start()
    build_food_compounds()
