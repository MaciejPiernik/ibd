import logging

import pandas as pd

from tqdm import tqdm
from pyensembl import EnsemblRelease


class GPL1708():
    def process(self, dataset):
        logging.info(f'Processing dataset {dataset.id}')

        expr_matrix = dataset.raw_dataset.pivot_samples('VALUE')
        expr_matrix_t = expr_matrix.T
        expr_df = pd.DataFrame(expr_matrix_t)
        locations = dataset.raw_dataset.gpls[next(iter(dataset.raw_dataset.gpls))].table.set_index('ID').CHROMOSOMAL_LOCATION

        data = EnsemblRelease(54)

        genes = []

        for loc in tqdm(locations, total=len(locations)):
            if pd.isnull(loc):
                genes.append(None)
                continue

            loc_split = loc.split(':')
            chr = loc_split[0][3:]
            locs = loc_split[1].split('-')
            start = int(locs[0])
            end = int(locs[1])

            gene = data.genes_at_locus(contig=chr, position=start, end=end)

            if len(gene) > 0:
                genes.append(gene[0].gene_id)
            else:
                genes.append(None)

        genes = pd.Series(genes, index=locations.index).dropna().rename('ENSEMBL_ID')

        expr_df = pd.concat([expr_df.transpose(), genes], axis=1, join='inner').set_index('ENSEMBL_ID').transpose()

        return expr_df