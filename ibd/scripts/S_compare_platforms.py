"""
Compare the microarray meta-analysis with the independent RNA-seq validation.

For the inflamed and the uninflamed comparison separately, reports agreement of
gene-level effect sizes and of pathway-level enrichment scores between the two
platforms, and the effective number of cohorts contributing to each
meta-analysis.

The effective number of cohorts is 1 / sum(w^2) over the normalised
random-effects weights. It states how many cohorts a pooled estimate really
rests on: with k cohorts of equal weight it equals k, and when one cohort
dominates it approaches 1.

Run:
    python ibd/scripts/compare_platforms.py
"""

import os

import numpy as np
import pandas as pd

MICRO = 'results/revision_v2'
RNA = 'results/revision_v2'
OUT = 'results/revision_v2/comparison'

COMPARISONS = [
    ('inflamed', 'uc_vs_hc', 'rnaseq_uc_vs_hc'),
    ('uninflamed', 'uninf_vs_hc', 'rnaseq_uninf_vs_hc'),
]


def effective_cohorts(directory):
    """Weight share per cohort and the effective number of cohorts."""
    summary = pd.read_csv(f'{directory}/dataset_summary.csv')
    tau2 = pd.read_csv(f'{directory}/gene_meta_all.csv')['tau2'].median()

    weight = 1 / (tau2 + 1 / summary['n_uc'] + 1 / summary['n_ctrl'])
    share = weight / weight.sum()

    return pd.DataFrame({'dataset': summary['dataset'], 'weight_share': share}), 1 / (share ** 2).sum()


def gene_agreement(micro_dir, rna_dir):
    micro = pd.read_csv(f'{micro_dir}/gene_meta_all.csv')[['gene', 'hedges_g', 'q_value']]
    rna = pd.read_csv(f'{rna_dir}/gene_meta_all.csv')[['gene', 'hedges_g', 'q_value']]
    micro_significant = set(pd.read_csv(f'{micro_dir}/gene_meta_significant.csv')['gene'])
    rna_significant = set(pd.read_csv(f'{rna_dir}/gene_meta_significant.csv')['gene'])

    merged = micro.merge(rna, on='gene', suffixes=('_micro', '_rna')).dropna()
    merged['significant_micro'] = merged['gene'].isin(micro_significant)
    merged['significant_rna'] = merged['gene'].isin(rna_significant)
    merged['same_direction'] = np.sign(merged['hedges_g_micro']) == np.sign(merged['hedges_g_rna'])

    hits = merged[merged['significant_micro']]

    return merged, {
        'shared_genes': len(merged),
        'r_all': merged['hedges_g_micro'].corr(merged['hedges_g_rna']),
        'microarray_significant': len(hits),
        'r_significant': hits['hedges_g_micro'].corr(hits['hedges_g_rna']),
        'same_direction': hits['same_direction'].mean(),
        'replicated': (hits['same_direction'] & hits['significant_rna']).mean(),
    }


def pathway_agreement(micro_dir, rna_dir):
    micro = pd.read_csv(f'{micro_dir}/pathway_meta_all.csv')[['term', 'combined_nes', 'q_value']]
    rna = pd.read_csv(f'{rna_dir}/pathway_meta_all.csv')[['term', 'combined_nes', 'q_value']]

    merged = micro.merge(rna, on='term', suffixes=('_micro', '_rna')).dropna()
    merged['same_direction'] = (
        np.sign(merged['combined_nes_micro']) == np.sign(merged['combined_nes_rna']))
    hits = merged[merged['q_value_micro'] < 0.05]

    return merged, {
        'shared_pathways': len(merged),
        'r_all': merged['combined_nes_micro'].corr(merged['combined_nes_rna']),
        'microarray_significant': len(hits),
        'same_direction': hits['same_direction'].mean(),
        'replicated': (hits['same_direction'] & (hits['q_value_rna'] < 0.05)).mean(),
    }


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    summaries = []

    for label, micro_name, rna_name in COMPARISONS:
        micro_dir, rna_dir = f'{MICRO}/{micro_name}', f'{RNA}/{rna_name}'
        print(f'\n===== {label.upper()} UC vs controls =====')

        for platform, directory in [('microarray', micro_dir), ('RNA-seq', rna_dir)]:
            shares, effective = effective_cohorts(directory)
            print(f'{platform}: {len(shares)} cohorts, effective {effective:.2f}, '
                  f'largest {shares["weight_share"].max():.0%}')
            shares.to_csv(f'{OUT}/weights_{label}_{platform}.csv', index=False)

        genes, gene_stats = gene_agreement(micro_dir, rna_dir)
        print(f'genes:    {gene_stats["shared_genes"]} shared, r = {gene_stats["r_all"]:.2f} | '
              f'of {gene_stats["microarray_significant"]} microarray-significant: '
              f'r = {gene_stats["r_significant"]:.2f}, '
              f'{gene_stats["same_direction"]:.1%} same direction, '
              f'{gene_stats["replicated"]:.1%} replicated')

        pathways, pathway_stats = pathway_agreement(micro_dir, rna_dir)
        print(f'pathways: {pathway_stats["shared_pathways"]} shared, r = {pathway_stats["r_all"]:.2f} | '
              f'of {pathway_stats["microarray_significant"]} microarray-significant: '
              f'{pathway_stats["same_direction"]:.1%} same direction, '
              f'{pathway_stats["replicated"]:.1%} replicated')

        genes.to_csv(f'{OUT}/genes_{label}.csv', index=False)
        pathways.to_csv(f'{OUT}/pathways_{label}.csv', index=False)
        summaries.append({'comparison': label, **gene_stats,
                          **{f'pathway_{k}': v for k, v in pathway_stats.items()}})

    pd.DataFrame(summaries).to_csv(f'{OUT}/summary.csv', index=False)
    print(f'\nWritten to {OUT}/')
