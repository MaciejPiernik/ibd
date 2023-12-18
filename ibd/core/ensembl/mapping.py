import pandas as pd

from pyensembl import EnsemblRelease
from tqdm import tqdm


def locations_to_genes(locations, ensembl_release):
    ensembl = EnsemblRelease(ensembl_release)

    genes = []

    for loc in tqdm(locations, total=len(locations)):
        if pd.isnull(loc):
            genes.append(None)
            continue
        elif loc == 'control':
            genes.append(loc)
            continue

        loc_split = loc.split(':')
        chr = loc_split[0][3:]
        locs = loc_split[1].split('-')
        start = int(locs[0])
        end = int(locs[1])

        gene = ensembl.genes_at_locus(contig=chr, position=start, end=end)

        if len(gene) > 0:
            genes.append(gene[0].gene_id)
        else:
            genes.append(None)

    genes = pd.Series(genes, index=locations.index).dropna().rename('ENSEMBL_ID')

    return genes