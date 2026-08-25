from __future__ import annotations
from collections import namedtuple

import requests
import urllib3
import pandas as pd
import bs4
import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

HEADERS = {
    "Accept": "application/json, text/javascript, */*; q=0.01",
    "Content-Type": "application/x-www-form-urlencoded; charset=UTF-8",
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)",
    "X-Requested-With": "XMLHttpRequest",
}

def _fetch_datatables_records(url: str) -> list[dict]:
    records_per_page = 50
    draw = 1
    all_data = []

    while True:
        start_index = (draw - 1) * records_per_page
        payload = {
            "draw": draw,
            "start": start_index,
            "length": records_per_page,
            "search[value]": "",
            "search[regex]": "false",
        }

        max_retries = 3
        page_success = False
        
        for attempt in range(max_retries):
            try:
                response = requests.post(url, headers=HEADERS, data=payload, verify=False, timeout=60)
                
                if response.status_code == 200:
                    data = response.json()
                    if data and "data" in data and len(data["data"]) > 0:
                        all_data.extend(data["data"])
                        total_records = data.get("recordsFiltered", 0)
                        
                        if start_index + records_per_page >= total_records:
                            return all_data 
                        
                        page_success = True
                        break 
                    else:
                        return all_data 
                else:
                    print(f"\n[API ERROR] Server rejected: {response.status_code}")
                    break 
                    
            except Exception as e:
                import time
                print(f"\n[WARNING] {url} - Server did not respond ({e}). Attempt: {attempt+1}/{max_retries}...")
                time.sleep(5)
        
        if not page_success:
            break
            
        draw += 1

    return all_data

def _get_interaction_descriptions() -> dict:
    raw_interactions = _fetch_datatables_records(urls.urls['ddinter_v2']['interaction_source'])
    return {str(inter.get("id")): inter.get("interaction", "") for inter in raw_interactions}

def _ensure_hashable(data):
    if isinstance(data, (dict, list, set)):
        return tuple(data)
    return data

def _scrape_legacy_identifiers(drug_id: str) -> dict:
    mapping_dict = {}
    url = urls.urls['ddinter']['mapping'] % drug_id
    
    try:
        c = curl.Curl(url, silent=True, large=True)
        
        if hasattr(c, 'fileobj') and c.fileobj:
            html_content = c.fileobj
        elif hasattr(c, 'result') and c.result:
            html_content = c.result
        else:
            response = requests.get(url, verify=False, timeout=10)
            if response.status_code == 200:
                html_content = response.text
            else:
                return mapping_dict 

        soup = bs4.BeautifulSoup(html_content, 'html.parser')
        refs = soup.find_all('a')
        
        mapping_targets = ['drugbank', 'chembl', 'pubchem']
        
        for link in refs:
            href = link.get('href', '')
            for target in mapping_targets:
                if target in href:
                    mapping_dict[target] = href.split('/')[-1]
                    break
                    
    except Exception as e:
        pass
        
    return mapping_dict

def ddinter_mappings_v2(return_df: bool = False) -> list[tuple] | pd.DataFrame:
    raw_drugs = _fetch_datatables_records(urls.urls['ddinter_v2']['drugs'])
    
    fields = ('ddinter', 'drugbank', 'chembl', 'pubchem')
    record = namedtuple('DdinterV2Identifiers', fields, defaults=(None,) * len(fields))
    result = set()

    for drug in raw_drugs:
        ddinter_id = drug.get("internalID")
        
        if ddinter_id:
            scraped_ids = _scrape_legacy_identifiers(ddinter_id)
            
            result.add(record(
                ddinter=ddinter_id,
                drugbank=scraped_ids.get('drugbank'),
                chembl=scraped_ids.get('chembl'),
                pubchem=scraped_ids.get('pubchem')
            ))

    return pd.DataFrame(result) if return_df else list(result)

def ddinter_interactions_v2(return_df: bool = False) -> list[tuple] | pd.DataFrame:
    raw_drugs = _fetch_datatables_records(urls.urls['ddinter_v2']['drugs'])
    drug_map = {d.get("internalID"): d for d in raw_drugs if d.get("internalID")}
    desc_map = _get_interaction_descriptions()
    
    result = set()
    record = namedtuple(
        'DdinterV2Interaction',
        ('drug1_id', 'drug1_name', 'drug2_id', 'drug2_name', 'level', 'actions'),
        defaults=None,
    )
    
    for drug_a_id, drug_a_info in drug_map.items():
        url = urls.urls['ddinter_v2']['drug_interactions'] % drug_a_id
        drug_a_pairs = _fetch_datatables_records(url)
        
        drug1_id = drug_a_id
        drug1_name = drug_a_info.get("name")
        
        for pair in drug_a_pairs:
            inter_id = str(pair.get("interaction_id"))
            
            mechanisms = [k for k in ('absorption', 'metabolism', 'excretion', 'distribution', 'synergistic_effect', 'antagonistic_effect') if pair.get(k) == '1']
            description = desc_map.get(inter_id, "")
            
            actions = tuple(mechanisms) + (description,) if description else tuple(mechanisms)
            
            interaction = record(
                drug1_id=drug1_id,
                drug1_name=drug1_name,
                drug2_id=_ensure_hashable(pair.get("drug_id")),
                drug2_name=_ensure_hashable(pair.get("drug_name")),
                level=_ensure_hashable(pair.get("level")),
                actions=_ensure_hashable(actions),
            )
            result.add(interaction)

    return pd.DataFrame(result) if return_df else list(result)

def ddinter_n_drugs_v2() -> int:
    raw_drugs = _fetch_datatables_records(urls.urls['ddinter_v2']['drugs'])
    return len(raw_drugs)

def ddinter_identifiers_v2(drug: str) -> list[tuple]:
    fields = ('ddinter', 'drugbank', 'chembl', 'pubchem')
    record = namedtuple('DdinterV2Identifiers', fields, defaults=(None,) * len(fields))
    scraped_ids = _scrape_legacy_identifiers(drug)
    
    result = record(
        ddinter=drug,
        drugbank=scraped_ids.get('drugbank'),
        chembl=scraped_ids.get('chembl'),
        pubchem=scraped_ids.get('pubchem')
    )
    return [result]

def ddinter_drug_interactions_v2(drug: str, return_df: bool = False) -> list[tuple] | pd.DataFrame:
    desc_map = _get_interaction_descriptions()
    
    raw_drugs = _fetch_datatables_records(urls.urls['ddinter_v2']['drugs'])
    drug1_name = next((d.get("name") for d in raw_drugs if d.get("internalID") == drug), None)
    
    record = namedtuple(
        'DdinterV2Interaction',
        ('drug1_id', 'drug1_name', 'drug2_id', 'drug2_name', 'level', 'actions'),
        defaults=None,
    )
    result = set()
    url = urls.urls['ddinter_v2']['drug_interactions'] % drug
    drug_pairs = _fetch_datatables_records(url)
    
    if not drug_pairs:
        return pd.DataFrame(result) if return_df else list(result)
        
    for pair in drug_pairs:
        inter_id = str(pair.get("interaction_id"))
        mechanisms = [k for k in ('absorption', 'metabolism', 'excretion', 'distribution', 'synergistic_effect', 'antagonistic_effect') if pair.get(k) == '1']
        description = desc_map.get(inter_id, "")
        
        actions = tuple(mechanisms) + (description,) if description else tuple(mechanisms)
        
        interaction = record(
            drug1_id=drug,
            drug1_name=drug1_name,  
            drug2_id=_ensure_hashable(pair.get("drug_id")),
            drug2_name=_ensure_hashable(pair.get("drug_name")),
            level=_ensure_hashable(pair.get("level")),
            actions=_ensure_hashable(actions),
        )
        result.add(interaction)

    return pd.DataFrame(result) if return_df else list(result)

def ddinter_disease_interactions_v2(return_df: bool = False) -> list[tuple] | pd.DataFrame:
    raw_drugs = _fetch_datatables_records(urls.urls['ddinter_v2']['drugs'])
    drug_map = {d.get("internalID"): d.get("name") for d in raw_drugs if d.get("internalID")}
    
    record = namedtuple(
        'DdinterV2DiseaseInteraction',
        ('drug_id', 'drug_name', 'disease_name', 'level', 'actions'), 
        defaults=None,
    )
    result = set()

    for d_id, d_name in drug_map.items():
        url = urls.urls['ddinter_v2']['disease_interactions'] % d_id
        records = _fetch_datatables_records(url)
        
        for row in records:
            disease_name = row.get('diseaseName', '').strip()
            
            description = row.get('text', '').strip()
            actions = (description,) if description else ('',)
                
            interaction = record(
                drug_id=_ensure_hashable(d_id),
                drug_name=_ensure_hashable(d_name),
                disease_name=_ensure_hashable(disease_name),
                level=_ensure_hashable(int(row.get('level') or 0)),
                actions=_ensure_hashable(actions)
            )
            result.add(interaction)
            
    return pd.DataFrame(result) if return_df else list(result)

def ddinter_food_interactions_v2(return_df: bool = False) -> list[tuple] | pd.DataFrame:
    raw_drugs = _fetch_datatables_records(urls.urls['ddinter_v2']['drugs'])
    drug_map = {d.get("internalID"): d.get("name") for d in raw_drugs if d.get("internalID")}
    
    record = namedtuple(
        'DdinterV2FoodInteraction',
        ('drug_id', 'drug_name', 'food_name', 'level', 'actions'),
        defaults=None,
    )
    result = set()

    for d_id, d_name in drug_map.items():
        url = urls.urls['ddinter_v2']['food_interactions'] % d_id
        records = _fetch_datatables_records(url)
        
        for row in records:
            description = row.get('newInteraction', '').strip()
            
            actions = (description,) if description else ('',)
                
            interaction = record(
                drug_id=_ensure_hashable(d_id),
                drug_name=_ensure_hashable(d_name),
                food_name=_ensure_hashable(row.get('foodName', '')),
                level=_ensure_hashable(int(row.get('level') or 0)),
                actions=_ensure_hashable(actions)
            )
            result.add(interaction)
            
    return pd.DataFrame(result) if return_df else list(result)