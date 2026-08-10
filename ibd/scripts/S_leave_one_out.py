"""
Leave-one-dataset-out sensitivity analysis.

Reviewer 1 asked whether the inverse-variance weighting lets the largest cohorts
drive the pooled estimates, leaving the smaller ones with little influence, and
asked for the analysis to be repeated with the largest datasets left out.

Two things are reported for every analysis:

  - how much of the pooled weight each cohort actually carries, under the
    random-effects weights used here and under fixed-effect weights for
    comparison. A large between-study variance pulls random-effects weights
    towards equality, so the concern applies far less than it would to a
    fixed-effect model.

  - the result of repeating the meta-analysis with each cohort omitted in turn:
    how many genes and pathways stay significant, how much of the full list is
    retained, and how closely the effect sizes agree with the full analysis.

Only the pooling step is repeated. Per-dataset statistics and enrichment scores
are unchanged, so the folds reuse the pipeline's own meta-analysis functions.

Run:
    python ibd/scripts/S_leave_one_out.py
"""

import argparse
import math
import os

import numpy as np
import pandas as pd

import experiments.Experiment as exp
from experiments.Experiment import build_config

ANALYSES = ['uc_vs_hc', 'uninf_vs_hc', 'infl_vs_uninfl', 'cd_vs_hc', 'uc_vs_cd']


def weight_shares(per_dataset, tau2):
    """Share of the pooled weight carried by each cohort, per weighting scheme."""
    variance = per_dataset.pivot(index='gene', columns='dataset', values='hedges_g_var')

    random = 1 / variance.add(tau2, axis=0)
    fixed = 1 / variance

    return pd.DataFrame({
        'random_effects': (random.div(random.sum(axis=1), axis=0)).mean(),
        'fixed_effect': (fixed.div(fixed.sum(axis=1), axis=0)).mean(),
    })


def run(name, results_dir):
    directory = f'{results_dir}/{name}'
    per_dataset = pd.read_csv(f'{directory}/per_dataset_gene_stats.csv')
    cohorts = sorted(per_dataset['dataset'].unique())

    full_meta = pd.read_csv(f'{directory}/gene_meta_all.csv')
    full_significant = set(pd.read_csv(f'{directory}/gene_meta_significant.csv')['gene'])
    full_effect = full_meta.set_index('gene')['hedges_g']

    print(f'\n===== {name}: {len(cohorts)} cohorts, {len(full_significant)} significant genes =====')

    shares = weight_shares(per_dataset, full_meta.set_index('gene')['tau2'])
    print('weight carried by each cohort (random-effects / fixed-effect):')
    for cohort, row in shares.sort_values('random_effects', ascending=False).iterrows():
        print(f'  {cohort:11s} {row.random_effects:5.1%} / {row.fixed_effect:5.1%}')

    rows = []
    for dropped in cohorts:
        kept = per_dataset[per_dataset['dataset'] != dropped]
        config = build_config(data_dir='./db', output_dir='/tmp/loo',
                              min_datasets=math.ceil((len(cohorts) - 1) / 2))
        fold = exp._filter_genes(exp._meta_analyse_genes(kept, config), config)
        fold_significant = set(fold['gene'])

        common = full_effect.index.intersection(fold.set_index('gene').index)
        correlation = full_effect.loc[common].corr(fold.set_index('gene')['hedges_g'].loc[common])
        retained = len(full_significant & fold_significant)

        rows.append({
            'analysis': name,
            'omitted': dropped,
            'weight_share': shares.loc[dropped, 'random_effects'],
            'significant': len(fold_significant),
            'retained': retained / len(full_significant) if full_significant else np.nan,
            'correlation': correlation,
        })
        print(f'  -{dropped:11s} significant {len(fold_significant):5d}  '
              f'retained {rows[-1]["retained"]:5.1%}  r = {correlation:.3f}')

    return pd.DataFrame(rows), shares.assign(analysis=name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', default='results/revision_v2')
    parser.add_argument('--output-dir', default='results/revision_v2')
    parser.add_argument('--analyses', default=','.join(ANALYSES))
    args = parser.parse_args()

    folds, shares = [], []
    for name in args.analyses.split(','):
        fold, share = run(name, args.results)
        folds.append(fold)
        shares.append(share)

    folds = pd.concat(folds, ignore_index=True)
    shares = pd.concat(shares)

    os.makedirs(args.output_dir, exist_ok=True)
    folds.to_csv(f'{args.output_dir}/leave_one_out.csv', index=False)
    shares.to_csv(f'{args.output_dir}/cohort_weights.csv')

    print('\nSummary across analyses:')
    for name, group in folds.groupby('analysis', sort=False):
        print(f'  {name:16s} retained {group.retained.min():5.1%}-{group.retained.max():5.1%}, '
              f'r {group.correlation.min():.3f}-{group.correlation.max():.3f}, '
              f'largest cohort weight {group.weight_share.max():.1%}')
    print(f'\nWritten to {args.output_dir}/')


if __name__ == '__main__':
    main()
