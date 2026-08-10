"""
Classify genes as constitutive or inflammation-dependent using the direct
inflamed-versus-uninflamed contrast.

The original classification inferred inflammation dependence indirectly, from a
gene's effect against controls when inflamed and when uninflamed. Those two
analyses have very different power, so a gene can look inflammation-dependent
merely because the uninflamed comparison was unable to detect it. The direct
contrast asks the same question within datasets and is better powered for it.

  inflammation-dependent  differs between inflamed and uninflamed UC
  constitutive            altered in uninflamed UC, no further change when inflamed
  inflammation-amplified  altered in uninflamed UC and further changed when inflamed
  undetermined            altered in inflamed UC, but neither of the above is shown

Run:
    python ibd/scripts/S_classify_layers.py
"""

import argparse
import math
import os

import numpy as np
import pandas as pd

COLUMNS = ['gene', 'hedges_g', 'se', 'q_value', 'datasets', 'direction_ratio']


def load(path, suffix):
    """Meta-analysis results plus a significance flag on the paper's criteria."""
    table = pd.read_csv(f'{path}/gene_meta_all.csv')[COLUMNS]
    minimum = math.ceil(table['datasets'].max() / 2)
    table['significant'] = (
        (table['q_value'] < 0.05)
        & (table['direction_ratio'] >= 0.8)
        & (table['datasets'] >= minimum)
    )
    print(f'  {path}: {int(table["datasets"].max())} cohorts, presence >= {minimum}, '
          f'{int(table["significant"].sum())} significant')

    return table.rename(columns={c: f'{c}_{suffix}' for c in COLUMNS[1:] + ['significant']})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', default='results/revision_v2')
    parser.add_argument('--output-dir', default='results/revision_v2')
    args = parser.parse_args()

    print('Loading:')
    inflamed = load(f'{args.results}/uc_vs_hc', 'infl')
    uninflamed = load(f'{args.results}/uninf_vs_hc', 'uninf')
    direct = load(f'{args.results}/infl_vs_uninfl', 'direct')

    genes = inflamed.merge(uninflamed, on='gene').merge(direct, on='gene')

    genes['layer'] = 'not significant'
    genes.loc[genes['significant_infl'] | genes['significant_uninf'], 'layer'] = 'undetermined'
    genes.loc[genes['significant_uninf'] & ~genes['significant_direct'], 'layer'] = 'constitutive'
    genes.loc[genes['significant_uninf'] & genes['significant_direct'], 'layer'] = 'inflammation-amplified'
    genes.loc[~genes['significant_uninf'] & genes['significant_direct'], 'layer'] = 'inflammation-dependent'

    order = ['constitutive', 'inflammation-amplified', 'inflammation-dependent',
             'undetermined', 'not significant']
    print(f'\nClassification over {len(genes)} genes measured in all three analyses:')
    for layer in order:
        subset = genes[genes['layer'] == layer]
        if len(subset):
            print(f'  {layer:24s} {len(subset):6d}   median |g| inflamed {subset["hedges_g_infl"].abs().median():.2f}')

    # How the direct contrast revises the indirect classification.
    print('\nGenes previously called inflammation-dependent, by the direct contrast:')
    previously = genes[genes['significant_infl'] & ~genes['significant_uninf']]
    confirmed = previously['significant_direct']
    print(f'  confirmed as inflammation-dependent: {int(confirmed.sum())} of {len(previously)} '
          f'({100 * confirmed.mean():.1f}%)')
    print(f'  undetermined: {int((~confirmed).sum())} — significant against controls when inflamed, '
          f'but not shown to differ from uninflamed UC')

    os.makedirs(args.output_dir, exist_ok=True)
    genes.to_csv(f'{args.output_dir}/layer_classification.csv', index=False)
    print(f'\nWritten to {args.output_dir}/layer_classification.csv')


if __name__ == '__main__':
    main()
