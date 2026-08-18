"""
Classify every gene into an activity layer.

Each gene is assigned by a triage between the three UC-specific contrasts:
inflamed vs uninflamed (direct), inflamed vs controls, and uninflamed vs
controls. A gene has to be significant in inflamed vs controls to be classified
at all; which of the other two contrasts it is also significant in then says
whether the change is already present without inflammation, and whether
inflammation adds anything to it.

  inflammation-amplified  significant in all three contrasts: present in
                          uninflamed tissue, and further changed when inflamed
  inflammation-dependent  not significant in uninflamed vs controls: only
                          present when inflammation is active
  constitutive            not significant in inflamed vs uninflamed: present in
                          uninflamed tissue, and no further change when inflamed
  undetermined            significant in at least one contrast, but in a
                          combination that places it in none of the layers
  not significant         significant in none of the three contrasts

Significance is the criterion used throughout the paper: q < 0.05, direction
consistency >= 0.8, and presence in at least half of the datasets contributing
to that contrast.

The pooled effect additionally has to carry the same sign in every contrast in
which the gene is significant. This is a condition between contrasts; agreement
between the datasets within a contrast is already covered by the direction
consistency in the significance criterion above. The remaining contrast is, for
each layer, the one the gene was not shown to differ in, so its sign carries no
information and is not required to agree.

Run:
    python ibd/scripts/S_gene_layers.py
"""

import argparse
import math
import os

import numpy as np
import pandas as pd

COLUMNS = ['gene', 'hedges_g', 'q_value', 'datasets', 'direction_ratio']
LAYERS = ['constitutive', 'inflammation-amplified', 'inflammation-dependent',
          'undetermined', 'not significant']
CLASSIFIED = LAYERS[:3]
OUTPUT = ['gene', 'layer', 'direction', 'hedges_g_infl', 'q_value_infl',
          'hedges_g_uninf', 'q_value_uninf', 'hedges_g_direct', 'q_value_direct']


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

    # Everything is judged relative to the inflamed contrast, so a gene that is
    # not significant there cannot be placed in any layer.
    altered = genes['significant_infl']
    sign = np.sign(genes['hedges_g_infl'])
    agrees_uninf = sign == np.sign(genes['hedges_g_uninf'])
    agrees_direct = sign == np.sign(genes['hedges_g_direct'])

    genes['layer'] = 'not significant'
    genes.loc[
        altered | genes['significant_uninf'] | genes['significant_direct'],
        'layer'] = 'undetermined'
    genes.loc[
        altered & genes['significant_uninf'] & ~genes['significant_direct'] & agrees_uninf,
        'layer'] = 'constitutive'
    genes.loc[
        altered & genes['significant_uninf'] & genes['significant_direct']
        & agrees_uninf & agrees_direct,
        'layer'] = 'inflammation-amplified'
    genes.loc[
        altered & ~genes['significant_uninf'] & genes['significant_direct'] & agrees_direct,
        'layer'] = 'inflammation-dependent'

    genes['direction'] = pd.Series(np.where(sign > 0, 'up', 'down'),
                                   index=genes.index).where(genes['layer'].isin(CLASSIFIED))

    print(f'\nClassification over {len(genes)} genes measured in all three contrasts:')
    for layer in LAYERS:
        subset = genes[genes['layer'] == layer]
        counts = subset['direction'].value_counts()
        split = f'{counts.get("up", 0)} up, {counts.get("down", 0)} down' if layer in CLASSIFIED else ''
        print(f'  {layer:24s} {len(subset):6d}   {split}')

    genes['layer'] = pd.Categorical(genes['layer'], categories=LAYERS, ordered=True)
    genes = genes.sort_values(
        ['layer', 'hedges_g_infl'], key=lambda c: c.abs() if c.name == 'hedges_g_infl' else c,
        ascending=[True, False])

    os.makedirs(args.output_dir, exist_ok=True)
    genes[OUTPUT].to_csv(f'{args.output_dir}/gene_layers.csv', index=False)
    print(f'\nWritten to {args.output_dir}/gene_layers.csv')


if __name__ == '__main__':
    main()
