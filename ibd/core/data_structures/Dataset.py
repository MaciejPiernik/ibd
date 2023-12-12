import pandas as pd

import logging

from ibd.core.platforms.platforms import get_platform

class Dataset:
    id: str
    raw_dataset: pd.DataFrame
    data: pd.DataFrame
    platform_id: str

    def __init__(self, id, raw_dataset):
        self.id = id
        self.raw_dataset = raw_dataset
        self.platform_id = raw_dataset.metadata['platform_id'][0]

    def process(self):
        platform = get_platform(self.platform_id)

        self.data = platform.process(self)

    def persist(self):
        logging.info(f'Saving dataset {self.id}')

        self.data.to_csv(f'./db/{self.id}.csv', index=False)