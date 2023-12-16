import GEOparse
import logging

import pandas as pd

from ibd.core.platforms.utils import get_platform


class Dataset:
    id: str
    raw_dataset: pd.DataFrame
    data: pd.DataFrame
    ensembl_release: int

    def __init__(self, id, ensembl_release):
        self.id = id
        self.ensembl_release = ensembl_release

    def load_raw_data(self, cache_dir='./data/geo_cache'):
        logging.info(f'Loading dataset {self.id}')
    
        gse = GEOparse.get_GEO(geo=self.id, destdir=cache_dir)
        
        self.raw_dataset = gse

    def process(self):
        platform = get_platform(self.platform_id)

        self.data = platform.process(self)

    def persist(self):
        logging.info(f'Saving dataset {self.id}')

        self.data.to_csv(f'./db/{self.id}.csv', index=False)

    @property
    def platform_id(self):
        return self.raw_dataset.metadata['platform_id'][0]