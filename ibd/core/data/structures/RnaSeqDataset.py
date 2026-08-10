import logging
import os
import re
import urllib.request

import pandas as pd

import ibd.core.data.processing.rnaseq_processors


GEO_SERIES = 'https://ftp.ncbi.nlm.nih.gov/geo/series'


class RnaSeqDataset:
    """A GEO RNA-seq series, processed into the same form as the microarray datasets.

    Mirrors Dataset: the raw data is fetched from GEO, a per-dataset processor
    turns it into a gene-level expression matrix and harmonised metadata, and
    both are written to ./db_rnaseq/. RNA-seq series deposit their expression
    data as supplementary files rather than in the series matrix, so the
    processors read those directly.
    """

    id: str
    raw_metadata: str
    metadata: pd.DataFrame
    data: pd.DataFrame


    def __init__(self, id):
        self.id = id


    @property
    def url(self):
        prefix = re.sub(r'\d{1,3}$', 'nnn', self.id)

        return f'{GEO_SERIES}/{prefix}/{self.id}'


    def load_raw_data(self, cache_dir='./data/geo_cache_rnaseq'):
        """Cache the sample annotation from the head of the series matrix."""
        logging.info(f'Loading dataset {self.id}')

        os.makedirs(cache_dir, exist_ok=True)
        path = f'{cache_dir}/{self.id}.txt'

        if not os.path.exists(path):
            with urllib.request.urlopen(f'{self.url}/matrix/{self.id}_series_matrix.txt.gz') as response:
                import gzip
                text = gzip.decompress(response.read()).decode('utf-8', errors='replace')
            open(path, 'w', encoding='utf-8').write(text.split('!series_matrix_table_begin')[0])

        self.raw_metadata = open(path, encoding='utf-8').read()


    def process(self, recompute_data=False, recompute_metadata=False):
        logging.info(f'Processing dataset {self.id}')

        self.load_raw_data()

        if self.data_in_db() and not recompute_data and self.metadata_in_db() and not recompute_metadata:
            logging.info(f'Dataset {self.id} already in database, skipping')
            return

        self.data, self.metadata = self.get_processor().process(self)
        self.persist()
        self.persist_metadata()


    def supplementary(self, filename, cache_dir='./db_rnaseq_raw'):
        """Path to a supplementary file, downloading it on first use."""
        os.makedirs(cache_dir, exist_ok=True)
        path = f'{cache_dir}/{filename}'

        if not os.path.exists(path):
            logging.info(f'Downloading {filename}')
            urllib.request.urlretrieve(f'{self.url}/suppl/{filename}', path)

        return path


    def sample_field(self, field):
        """One row of the series matrix header, with one entry per sample."""
        for line in self.raw_metadata.splitlines():
            if line.startswith(f'!Sample_{field}'):
                return [value.strip('"') for value in line.split('\t')[1:]]


    def sample_characteristic(self, name):
        """One characteristics field of the series matrix, with one entry per sample."""
        for line in self.raw_metadata.splitlines():
            if line.startswith('!Sample_characteristics_ch1') and f'"{name}:' in line:
                values = [value.strip('"') for value in line.split('\t')[1:]]
                return [value.split(':', 1)[1].strip() for value in values]


    def persist(self):
        logging.info(f'Saving dataset {self.id}')

        self.data.to_parquet(f'./db_rnaseq/{self.id}.parquet')


    def persist_metadata(self):
        logging.info(f'Saving metadata for dataset {self.id}')

        self.metadata.to_csv(f'./db_rnaseq/{self.id}_metadata.csv')


    def get_processor(self):
        class_ = getattr(ibd.core.data.processing.rnaseq_processors, f'{self.id}_RnaSeqProcessor')

        return class_()


    def has_processor(self):
        return hasattr(ibd.core.data.processing.rnaseq_processors, f'{self.id}_RnaSeqProcessor')


    def data_in_db(self):
        return os.path.exists(f'./db_rnaseq/{self.id}.parquet')


    def metadata_in_db(self):
        return os.path.exists(f'./db_rnaseq/{self.id}_metadata.csv')
