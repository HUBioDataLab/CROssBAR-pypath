import requests
import aiohttp
import asyncio

import pandas as pd

from tqdm import tqdm
from tqdm.asyncio import tqdm_asyncio
from typing import Union





class ExpressionAtlasExperimentManager:
    def __init__(self, 
                 base_url: str, 
                 experimental_factors: Union[list[str], None] = ["organism part", "cell type", "cell line", "compound", "infect", "disease"],
                 ) -> None:
        
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
        response = requests.get(self.base_url)
        response.raise_for_status()
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
        
    def check_link(self, url: str) -> bool:
        """Check if a URL is accessible."""
        try:
            response = requests.head(url, allow_redirects=True)
            return response.status_code == 200
        except requests.RequestException:
            return False
        
    def validate_links(self, df: pd.DataFrame) -> pd.DataFrame:
        """Validate experiment links and return results DataFrame."""
        results = []
        for _, row in tqdm(df.iterrows(), total=df.shape[0], desc="Checking links"):
            experiment_type = row['rawExperimentType']
            matching_factors = row['matching_factors']
            for accession in row['experimentAccession']:
                data_link, design_link = (link.format(accession) for link in self.EXPERIMENT_LINKS[experiment_type])
                results.append({
                    'experimentAccession': accession,
                    'data_link': data_link,
                    'is_data_link_accessible': self.check_link(data_link),
                    'design_link': design_link,
                    'is_design_link_accessible': self.check_link(design_link),
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
        validated_df = self.validate_links(df)
        return self.separate_experiment_types(validated_df)

    


class AsyncExpressionAtlasExperimentManager(ExpressionAtlasExperimentManager):
    def __init__(self, 
                 base_url: str, 
                 experimental_factors: Union[list[str], None] = ["organism part", "cell type", "cell line", "compound", "infect", "disease"],
                 ) -> None:
        
        super().__init__(base_url, experimental_factors)


    
    async def check_link(self, url: str) -> bool:
        """Asynchronously check if a URL is accessible."""
        try:
            async with aiohttp.ClientSession() as session:
                async with session.head(url, allow_redirects=True) as response:
                    return response.status == 200
        except aiohttp.ClientError:
            return False
        
    
    async def validate_links(self, df: pd.DataFrame) -> pd.DataFrame:
        """Validate experiment links asynchronously and return results DataFrame."""
        tasks = []
        results = []

        async def validate(row):
            experiment_type = row['rawExperimentType']
            matching_factors = row['matching_factors']
            for accession in row['experimentAccession']:
                data_link, design_link = (link.format(accession) for link in self.EXPERIMENT_LINKS[experiment_type])
                is_data_link_accessible = await self.check_link(data_link)
                is_design_link_accessible = await self.check_link(design_link)
                results.append({
                    'experimentAccession': accession,
                    'data_link': data_link,
                    'is_data_link_accessible': is_data_link_accessible,
                    'design_link': design_link,
                    'is_design_link_accessible': is_design_link_accessible,
                    'matching_factors': matching_factors
                })

        for _, row in df.iterrows():
            tasks.append(validate(row))

        # await asyncio.gather(*tasks)
        await tqdm_asyncio.gather(*tasks, total=len(tasks), desc="Validating experiment links")
        return pd.DataFrame(results)
    
    async def __call__(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        df = self.filter_experiments(self.fetch_experiment_data())
        validated_df = await self.validate_links(df)
        return self.separate_experiment_types(validated_df)
