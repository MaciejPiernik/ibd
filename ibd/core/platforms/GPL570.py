import logging

import pandas as pd

from ibd.core.ensembl.mapping import locations_to_genes


class GPL570():
    def process(self, dataset):
        logging.info(f'Processing dataset {dataset.id}')

        raise NotImplementedError('GPL570 is not implemented yet')