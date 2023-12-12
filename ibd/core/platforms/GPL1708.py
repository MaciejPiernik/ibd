import logging

import pandas as pd


class GPL1708():
    def process(self, dataset):
        logging.info(f'Processing dataset {dataset.id}')

        expr_matrix = dataset.raw_dataset.pivot_samples('VALUE')
        expr_matrix_t = expr_matrix.T
        expr_df = pd.DataFrame(expr_matrix_t)
        gene_symbols = dataset.raw_dataset.gpls[next(iter(dataset.raw_dataset.gpls))].table.set_index('ID').GENE_SYMBOL.dropna()
        expr_df = pd.concat([expr_df.transpose(), gene_symbols], axis=1, join='inner').set_index('GENE_SYMBOL').transpose()

        return expr_df