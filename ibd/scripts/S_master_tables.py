"""
Build the master gene-level and pathway-level result tables.

Every downstream table and figure needs the same handful of facts about a gene
or a pathway, currently scattered across five microarray analyses, three
RNA-seq analyses, the rank-aggregation output, the knee filter, the clustering
and the layer classification. This assembles all of it into one row per gene
and one row per pathway, so that anything downstream is a lookup rather than
another merge.

Columns are named <quantity>_<analysis>, with the analyses abbreviated:

  infl    inflamed UC vs controls          cd      inflamed CD vs controls
  uninf   uninflamed UC vs controls        uccd    UC vs CD, directly
  direct  inflamed vs uninflamed UC

and an rna_ prefix for the three contrasts that have an RNA-seq counterpart.
Significance is the paper's criterion throughout (q < 0.05, direction
consistency >= 0.8, presence in at least half the contributing cohorts).

The gene table also carries the activity-layer classification and its RNA-seq
grading from S_gene_layers.py, which must be run first.

Run:
    python ibd/scripts/S_master_tables.py
"""

import argparse
import math
import os
import re

import numpy as np
import pandas as pd

MICRO = [('infl', 'uc_vs_hc'), ('uninf', 'uninf_vs_hc'), ('direct', 'infl_vs_uninfl'),
         ('cd', 'cd_vs_hc'), ('uccd', 'uc_vs_cd')]
RNA = [('rna_infl', 'rnaseq_uc_vs_hc'), ('rna_uninf', 'rnaseq_uninf_vs_hc'),
       ('rna_direct', 'rnaseq_infl_vs_uninfl')]
HGNC_URL = ('https://www.genenames.org/cgi-bin/download/custom?'
            'col=gd_app_sym&col=gd_locus_group&status=Approved&'
            'hgnc_datea=&hgnc_dateb=&order_by=gd_app_sym_sort&'
            'format=text&submit=submit')


def significance(table, value_column):
    """The paper's three-part criterion, and the presence threshold behind it."""
    minimum = math.ceil(table['datasets'].max() / 2)
    significant = (
        (table['q_value'] < 0.05)
        & (table['direction_ratio'] >= 0.8)
        & (table['datasets'] >= minimum)
    )
    return significant, minimum


def gene_block(results, path, suffix, full):
    """One analysis' gene-level results, with RRA where the fuller set is wanted."""
    table = pd.read_csv(f'{results}/{path}/gene_meta_all.csv')
    table['significant'], minimum = significance(table, 'hedges_g')

    keep = ['gene', 'hedges_g', 'se', 'q_value', 'datasets', 'significant']
    if full:
        keep += ['I2', 'direction_ratio']
        rra = pd.read_csv(f'{results}/{path}/gene_rra.csv')[['gene', 'rra_rho', 'rra_q']]
        table = table.merge(rra, on='gene', how='left')
        keep += ['rra_rho', 'rra_q']

    print(f'  {path:22s} {int(table["datasets"].max()):2d} cohorts, presence >= {minimum}, '
          f'{int(table["significant"].sum()):5d} significant')

    return table[keep].rename(columns={c: f'{c}_{suffix}' for c in keep[1:]})


def pathway_block(results, path, suffix, full):
    """One analysis' pathway-level results, with knee and cluster status."""
    table = pd.read_csv(f'{results}/{path}/pathway_meta_all.csv')
    table['significant'], minimum = significance(table, 'combined_nes')

    keep = ['term', 'combined_nes', 'se', 'q_value', 'datasets', 'significant']
    if full:
        keep += ['I2', 'direction_ratio']
        knee = f'{results}/{path}/pathway_meta_significant_knee.csv'
        if os.path.exists(knee):
            table['knee'] = table['term'].isin(pd.read_csv(knee)['term'])
            keep.append('knee')
        clusters = f'{results}/{path}/hierarchical_pathway_clusters.csv'
        if os.path.exists(clusters):
            member = {}
            for row in pd.read_csv(clusters).itertuples():
                for term in str(row.Pathways).split(';'):
                    member[term.strip()] = row.Representative
            table['cluster'] = table['term'].map(member)
            keep.append('cluster')

    print(f'  {path:22s} {int(table["datasets"].max()):2d} cohorts, presence >= {minimum}, '
          f'{int(table["significant"].sum()):5d} significant')

    return table[keep].rename(columns={c: f'{c}_{suffix}' for c in keep[1:]})


def tidy_flags(table):
    """Outer joins leave NaN in the boolean columns; absent means not significant."""
    for column in table.columns:
        if column.startswith(('significant_', 'replicated_', 'knee')):
            table[column] = table[column].fillna(False).astype(bool)
    return table


def add_replication(table, key):
    """Whether each RNA-seq counterpart confirms the microarray call."""
    for suffix in ['infl', 'uninf', 'direct']:
        micro, rna = f'{key}_{suffix}', f'{key}_rna_{suffix}'
        same_sign = np.sign(table[micro]) == np.sign(table[rna])
        table[f'replicated_{suffix}'] = (
            table[f'significant_{suffix}'].fillna(False)
            & table[f'significant_rna_{suffix}'].fillna(False)
            & same_sign)
    return table


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', default='results/revision_v2')
    parser.add_argument('--output-dir', default='results/revision_v2')
    parser.add_argument('--hgnc', default=None,
                        help='Path to HGNC TSV; if omitted, downloads from genenames.org')
    args = parser.parse_args()

    # ── genes ───────────────────────────────────────────────────────────
    print('Gene-level, microarray:')
    genes = None
    for suffix, path in MICRO:
        part = gene_block(args.results, path, suffix, full=True)
        genes = part if genes is None else genes.merge(part, on='gene', how='outer')

    print('Gene-level, RNA-seq:')
    for suffix, path in RNA:
        genes = genes.merge(gene_block(args.results, path, suffix, full=False),
                            on='gene', how='outer')

    genes = add_replication(tidy_flags(genes), 'hedges_g')

    layers = pd.read_csv(f'{args.results}/gene_layers.csv')[
        ['gene', 'constitutive', 'inflammation_dependent', 'unresolved', 'direction',
         'support_constitutive', 'support_inflammation_dependent']]
    genes = genes.merge(layers, on='gene', how='left')
    for flag in ['constitutive', 'inflammation_dependent', 'unresolved']:
        genes[flag] = genes[flag].fillna(False).astype(bool)

    coding = pd.read_csv(args.hgnc or HGNC_URL, sep='\t')
    coding = set(coding.loc[coding['Locus group'] == 'protein-coding gene', 'Approved symbol'])
    genes['protein_coding'] = genes['gene'].isin(coding)

    genes = tidy_flags(genes).sort_values('hedges_g_infl', key=lambda c: c.abs(),
                                          ascending=False)
    os.makedirs(args.output_dir, exist_ok=True)
    genes.to_csv(f'{args.output_dir}/gene_master.csv', index=False)
    print(f'\n  {len(genes)} genes x {len(genes.columns)} columns '
          f'-> {args.output_dir}/gene_master.csv')

    # ── pathways ────────────────────────────────────────────────────────
    print('\nPathway-level, microarray:')
    pathways = None
    for suffix, path in MICRO:
        part = pathway_block(args.results, path, suffix, full=True)
        pathways = part if pathways is None else pathways.merge(part, on='term', how='outer')

    print('Pathway-level, RNA-seq:')
    for suffix, path in RNA:
        pathways = pathways.merge(pathway_block(args.results, path, suffix, full=False),
                                  on='term', how='outer')

    pathways = add_replication(tidy_flags(pathways), 'combined_nes')

    lead = pd.read_csv(f'{args.results}/uc_vs_hc/pathway_meta_all.csv')[
        ['term', 'lead_genes_core']]
    pathways = pathways.merge(lead, on='term', how='left')
    pathways['term_clean'] = pathways['term'].str.replace(r' R-HSA-\d+$', '', regex=True)

    pathways = tidy_flags(pathways).sort_values('combined_nes_infl', key=lambda c: c.abs(),
                                                ascending=False)
    pathways.to_csv(f'{args.output_dir}/pathway_master.csv', index=False)
    print(f'\n  {len(pathways)} pathways x {len(pathways.columns)} columns '
          f'-> {args.output_dir}/pathway_master.csv')


if __name__ == '__main__':
    main()
