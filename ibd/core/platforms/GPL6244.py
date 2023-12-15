import logging

import pandas as pd
from pyensembl import EnsemblRelease


class GPL6244():
    def process(self, dataset):
        logging.info(f'Processing dataset {dataset.id}')

        expr_matrix = dataset.raw_dataset.pivot_samples('VALUE')
        expr_matrix_t = expr_matrix.T
        expr_df = pd.DataFrame(expr_matrix_t)
        transcript_ids = dataset.raw_dataset.gpls[next(iter(dataset.raw_dataset.gpls))].table.set_index('ID').ENSEMBL_ID.dropna()
        expr_df = pd.concat([expr_df.transpose(), transcript_ids], axis=1, join='inner').set_index('ENSEMBL_ID').transpose()

        return expr_df