import requests
import pandas as pd
from typing import Union

from tqdm import tqdm

import pypath.share.session as session

class ExpressionAtlasExperimentManager:
    def __init__(self, 
                 base_url: str, 
                 experimental_factors: Union[list[str], None] = ["organism part", "cell type", "cell line", "compound", "infect", "disease"],
                 ) -> None:
        
        self._logger = session.Logger(name='inputs.expression_atlas.experiement_manager')
        self._log = self._logger._log
        
        self.base_url = base_url
        
        self.EXPERIMENTAL_FACTORS_TO_INCLUDE = experimental_factors
        

        self.EXPERIMENT_LINKS = {
            'RNASEQ_MRNA_BASELINE': ('https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv','https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/ExperimentDesignFile.Baseline/experiment-design'),
            'PROTEOMICS_BASELINE': ('https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/ExperimentDownloadSupplier.Proteomics/tsv', 'https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/ExperimentDesignFile.Baseline/experiment-design'),
            'PROTEOMICS_BASELINE_DIA': ('https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/ExperimentDownloadSupplier.Proteomics/tsv', 'https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/ExperimentDesignFile.Baseline/experiment-design'),
            'RNASEQ_MRNA_DIFFERENTIAL': ('https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/DifferentialSecondaryDataFiles.RnaSeq/analytics', 'https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/ExperimentDesignFile.RnaSeq/experiment-design'),
            'MICROARRAY_1COLOUR_MRNA_DIFFERENTIAL': ('https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/DifferentialSecondaryDataFiles.Microarray/analytics', 'https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/ExperimentDesignFile.Microarray/experiment-design'),
            'MICROARRAY_2COLOUR_MRNA_DIFFERENTIAL': ('https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/DifferentialSecondaryDataFiles.Microarray/analytics', 'https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/ExperimentDesignFile.Microarray/experiment-design'),
            'MICROARRAY_1COLOUR_MICRORNA_DIFFERENTIAL': ('https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/DifferentialSecondaryDataFiles.Microarray/analytics','https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/ExperimentDesignFile.Microarray/experiment-design' ),
            'PROTEOMICS_DIFFERENTIAL': ('https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/ExperimentDownloadSupplier.BulkDifferential/tsv','https://www.ebi.ac.uk/gxa/experiments-content/{}/resources/ExperimentDesignFile.RnaSeq/experiment-design')
        }
    
    def fetch_experiment_data(self) -> pd.DataFrame:
        """Fetch experiment data from the API."""
        self._log(f"Fetching experiment data from the API: {self.base_url}")
        response = requests.get(self.base_url)
        response.raise_for_status()
        self._log(f"Experiment data fetched successfully.")
        return pd.DataFrame(response.json()["experiments"])
    
    def get_matching_factors(self, factors) -> Union[str, None]:
        """Return matching factors if any, else None."""
        if isinstance(factors, list) and len(factors) == 1:
            matching_factors = [factor for factor in factors if factor in self.EXPERIMENTAL_FACTORS_TO_INCLUDE]
            if matching_factors:
                return matching_factors[0]
        
        return None
    
    def filter_experiments(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filter experiments based on specific terms and return a grouped DataFrame with matching factors."""

        df['matching_factors'] = df["experimentalFactors"].apply(self.get_matching_factors)
        df = df[df['matching_factors'].notna()]
        df = df[~df['rawExperimentType'].str.contains("PROTEOMICS", case=False)]
        return df.groupby(['rawExperimentType', 'matching_factors'])['experimentAccession'].apply(list).reset_index()
        
    def generate_links_table(self, df: pd.DataFrame) -> pd.DataFrame:
        """Generate a DataFrame with experiment links."""
        results = []
        for _, row in tqdm(df.iterrows(), total=df.shape[0], desc="Generating links"):
            experiment_type = row['rawExperimentType']
            matching_factors = row['matching_factors']
            for accession in row['experimentAccession']:
                data_link, design_link = (link.format(accession) for link in self.EXPERIMENT_LINKS[experiment_type])
                results.append({
                    'experimentAccession': accession,
                    'data_link': data_link,
                    'design_link': design_link,
                    'matching_factors': matching_factors
                })
                
        return pd.DataFrame(results)
    
    def separate_experiment_types(self, df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Separate results into baseline and differential dataframes."""
        df_baseline_links = df[df['data_link'].str.contains("Baseline")]
        df_differential_links = df[df['data_link'].str.contains("Differential")]
        return df_baseline_links.reset_index(drop=True), df_differential_links.reset_index(drop=True)
    
    def __call__(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        df = self.filter_experiments(self.fetch_experiment_data())
        links_table = self.generate_links_table(df)
        return self.separate_experiment_types(links_table)

