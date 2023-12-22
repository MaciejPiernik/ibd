import GEOparse
import logging

import pandas as pd

from ibd.core.platforms.Platform import Platform

GEOparse.logger.set_verbosity("ERROR")


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
        platform = Platform(self.platform_id)

        self.data = platform.process(self)

    def persist(self):
        logging.info(f'Saving dataset {self.id}')

        self.data.to_parquet(f'./db/{self.id}.parquet')

    @property
    def platform_id(self):
        return self.raw_dataset.metadata['platform_id'][0]