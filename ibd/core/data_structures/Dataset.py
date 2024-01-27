import GEOparse
import logging

import pandas as pd

import ibd.core.data_processing.metadata_processors

from ibd.core.data_processing.gpt_metadata_processor import GPTMetadataProcessor
from ibd.core.platforms.Platform import Platform


GEOparse.logger.set_verbosity("ERROR")


class Dataset:
    id: str
    raw_dataset: pd.DataFrame
    metadata: pd.DataFrame
    data: pd.DataFrame


    @property
    def platform_id(self):
        return self.raw_dataset.metadata['platform_id'][0]


    def __init__(self, id):
        self.id = id


    def load_raw_data(self, cache_dir='./data/geo_cache'):
        logging.info(f'Loading dataset {self.id}')
    
        gse = GEOparse.get_GEO(geo=self.id, destdir=cache_dir)
        
        self.raw_dataset = gse


    def process(self):
        platform = Platform(self.platform_id)

        self.data = platform.process(self)
        self.metadata = self.process_metadata()


    def process_metadata(self):
        logging.info(f'Processing metadata for dataset {self.id}')
        # processor = self.get_metadata_processor()
        processor = GPTMetadataProcessor()

        return processor.process(self.raw_dataset.phenotype_data)


    def persist(self):
        logging.info(f'Saving dataset {self.id}')

        self.data.to_parquet(f'./db/{self.id}.parquet')
        self.metadata.to_parquet(f'./db/{self.id}_metadata_gpt.parquet')


    def get_metadata_processor(self):
        class_ = getattr(ibd.core.data_processing.metadata_processors, f'{self.id}_MetadataProcessor')

        return class_()