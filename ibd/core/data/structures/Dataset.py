from io import StringIO
import os
import GEOparse
import logging

import pandas as pd

import ibd.core.data.processing.metadata_processors

from ibd.core.data.processing.gpt_metadata_processor import GPTMetadataProcessor
from ibd.core.platforms.Platform import Platform


GEOparse.logger.set_verbosity("ERROR")


class Dataset:
    id: str
    raw_dataset: pd.DataFrame
    metadata: pd.DataFrame
    data: pd.DataFrame
    batch_size: int


    @property
    def platform_id(self):
        return self.raw_dataset.metadata['platform_id'][0]


    def __init__(self, id, batch_size=None):
        self.id = id
        self.batch_size = batch_size


    def load_raw_data(self, cache_dir='./data/geo_cache'):
        logging.info(f'Loading dataset {self.id}')
    
        gse = GEOparse.get_GEO(geo=self.id, destdir=cache_dir)
        
        self.raw_dataset = gse


    def process(self, recompute_data=False, recompute_metadata=False):
        logging.info(f'Processing dataset {self.id}')

        self.load_raw_data()

        platform = Platform(self.platform_id)

        if self.data_in_db() and not recompute_data:
            logging.info(f'Dataset {self.id} data already in database, skipping')
        else:
            self.data = platform.process(self)
            self.persist()

        if self.metadata_in_db() and not recompute_metadata:
            logging.info(f'Dataset {self.id} metadata already in database, skipping')
        else:
            self.metadata = self.process_metadata()
            self.persist_metadata()


    def process_metadata(self):
        logging.info(f'Processing metadata for dataset {self.id}')
        processor = self.get_metadata_processor()

        if self.batch_size is None:
            metadata = processor.process(self.raw_dataset.phenotype_data)
        else:
            metadata_batches = [self.raw_dataset.phenotype_data[i:i + self.batch_size] for i in range(0, len(self.raw_dataset.phenotype_data), self.batch_size)]

            metadata = pd.DataFrame()
            for batch in metadata_batches:
                batch_metadata_str = processor.process(batch)
                batch_metadata_buf = StringIO(batch_metadata_str)
                batch_metadata_df = pd.read_csv(batch_metadata_buf)

                metadata = pd.concat([metadata, batch_metadata_df])

        return metadata


    def persist(self):
        logging.info(f'Saving dataset {self.id}')

        self.data.to_parquet(f'./db/{self.id}.parquet')


    def persist_metadata(self):
        logging.info(f'Saving metadata for dataset {self.id}')

        self.metadata.to_csv(f'./db/{self.id}_metadata.csv')


    def get_metadata_processor(self):
        class_ = getattr(ibd.core.data.processing.metadata_processors, f'{self.id}_MetadataProcessor')

        return class_()
    

    def data_in_db(self):
        return os.path.exists(f'./db/{self.id}.parquet')
    

    def metadata_in_db(self):
        return os.path.exists(f'./db/{self.id}_metadata.csv')