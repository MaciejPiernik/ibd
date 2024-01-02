import logging

import pandas as pd

from ibd.core.platforms.GPL1708 import GPL1708
from ibd.core.platforms.GPL570 import GPL570
from ibd.core.platforms.GPL6244 import GPL6244


class Platform:
    def __init__(self, platform_id: str):
        self.id = platform_id
        self._instance = self.get_platform(platform_id)


    def process(self, dataset):
        logging.info(f'Extracting expression data {dataset.id}')

        expr_matrix = dataset.raw_dataset.pivot_samples('VALUE')
        expr_matrix_t = expr_matrix.T
        expr_df = pd.DataFrame(expr_matrix_t)

        locations = dataset.raw_dataset.gpls[next(iter(dataset.raw_dataset.gpls))].table.set_index('ID')

        id_lists = self._instance.get_entrez_id_lists(locations)

        ids = id_lists.apply(pd.Series, 1).stack()
        ids = ids.reset_index(level=1, drop=True)
        ids = ids.rename('ENTREZ_ID')

        expr_df = expr_df.transpose().join(ids, how='inner').set_index(ids.name).transpose()

        expr_df = expr_df.groupby(expr_df.columns, axis=1).mean()

        return expr_df
    

    def get_platform(self, platform_id):
        if platform_id in ['GPL570', 'GPL17996']:
            return GPL570()
        elif platform_id == 'GPL1708':
            return GPL1708()
        elif platform_id == 'GPL6244':
            return GPL6244()
        else:
            raise NotImplementedError(f'{platform_id} is not implemented yet')