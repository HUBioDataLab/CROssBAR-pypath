import json
import urllib.request
import ssl
import collections
import pypath.share.session as session

_logger = session.Logger(name='disbiome_input')
_log = _logger._log

DisbiomeExperiment = collections.namedtuple(
    'DisbiomeExperiment',
    [
        'experiment_id', 'disease_id', 'organism_id', 'location_id',
        'host_id', 'publication_id', 'methodological_detail_id',
        'organism_name', 'organism_ncbi_id',
        'disease_name', 'meddra_id', 'meddra_level',
        'qualitative_outcome', 'subject_value', 'control_value',
        'ratio', 'response_name', 'response_unit',
        'method_name', 'sample_name', 'control_name', 'host_type',
        'location_name', 'pubmed_id'
    ]
)

def fetch_json_from_url(url):
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req, context=ctx, timeout=30) as response:
            return json.loads(response.read().decode('utf-8'))
    except Exception as e:
        _log(f"Error fetching {url}: {e}")
        return []

def get_disbiome_locations():
    data = fetch_json_from_url("https://disbiome.ugent.be:8080/location")
    return {item.get('location_id') or item.get('id'): item.get('name') for item in data if item}

def get_disbiome_publications():
    data = fetch_json_from_url("https://disbiome.ugent.be:8080/publication")
    pub_dict = {}
    for item in data:
        pub_id = item.get('publication_id') or item.get('id')
        pubmed_url = item.get('pubmed_url')
        pubmed_id = None
        if pubmed_url:
            parts = [p for p in pubmed_url.strip().split('/') if p]
            if parts:
                pubmed_id = parts[-1]
        pub_dict[pub_id] = pubmed_id
    return pub_dict

def get_disbiome_experiments():
    _log("Starting Disbiome experiments data retrieval...")
    
    locations_map = get_disbiome_locations()
    publications_map = get_disbiome_publications()
    
    url = "https://disbiome.ugent.be:8080/experiment"
    data = fetch_json_from_url(url)
    
    if not data:
        return []
        
    experiments = set()
    
    for item in data:
        try:
            exp = DisbiomeExperiment(
                experiment_id=item.get('experiment_id'), 
                disease_id=item.get('disease_id'),
                organism_id=item.get('organism_id'),
                location_id=item.get('location_id'),
                host_id=item.get('host_id'),
                publication_id=item.get('publication_id'),
                methodological_detail_id=item.get('methodological_detail_id'),
                organism_name=item.get('organism_name'),
                organism_ncbi_id=item.get('organism_ncbi_id'),
                disease_name=item.get('disease_name'),
                meddra_id=item.get('meddra_id'),
                meddra_level=item.get('meddra_level'),
                qualitative_outcome=item.get('qualitative_outcome'),
                subject_value=item.get('subject_value'),
                control_value=item.get('control_value'),
                ratio=item.get('ratio'),
                response_name=item.get('response_name'),
                response_unit=item.get('response_unit'),
                method_name=item.get('method_name'),
                sample_name=item.get('sample_name'),
                control_name=item.get('control_name'),
                host_type=item.get('host_type'),
                location_name=locations_map.get(item.get('location_id')),
                pubmed_id=publications_map.get(item.get('publication_id'))
            )
            experiments.add(exp)
        except Exception as e:
            _log(f"Error parsing item: {e}")
            
    _log(f"Successfully loaded {len(experiments)} Disbiome experiments.")
    return list(experiments)