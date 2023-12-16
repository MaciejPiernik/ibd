import logging

import pandas as pd

from ibd.core.ensembl.mapping import locations_to_genes


class GPL1708():
    def process(self, dataset):
        logging.info(f'Processing dataset {dataset.id}')

        expr_matrix = dataset.raw_dataset.pivot_samples('VALUE')
        expr_matrix_t = expr_matrix.T
        expr_df = pd.DataFrame(expr_matrix_t)
        locations = dataset.raw_dataset.gpls[next(iter(dataset.raw_dataset.gpls))].table.set_index('ID').CHROMOSOMAL_LOCATION

        genes = locations_to_genes(locations, dataset.ensembl_release)

        expr_df = pd.concat([expr_df.transpose(), genes], axis=1, join='inner').set_index('ENSEMBL_ID').transpose()

        return expr_df