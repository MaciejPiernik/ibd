"""
Sensitivity of the constitutive gene set to the effect-size ratio window.

The original classification identified constitutive genes in two steps: a Wald
test finding no significant difference between the inflamed and uninflamed
effect sizes, then a window on their ratio requiring the two to be numerically
close. The Wald test alone establishes only absence of evidence for a
difference, which underpowered data pass easily, and the window was added to
require approximate equality.

Reviewer 1 asked why that window was set at 0.70-1.30 and whether other widths
give overlapping gene lists. This varies the width and reports the resulting
sets. Because the criterion is a symmetric interval on a single statistic, the
sets are necessarily nested: the width sets stringency, not composition.

Run:
    python ibd/scripts/S_ratio_window_sensitivity.py
"""

import argparse
import os

import pandas as pd

WINDOWS = [(0.50, 1.50), (0.60, 1.40), (0.70, 1.30), (0.80, 1.20), (0.85, 1.15), (0.90, 1.10)]
REFERENCE = (0.70, 1.30)


def constitutive(genes, low, high):
    return set(genes[(genes['category'] == 'constitutive') & genes['ratio'].between(low, high)]['gene'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--comparison', default='results/revision_v2/comparison_stats.csv')
    parser.add_argument('--output-dir', default='results/revision_v2')
    args = parser.parse_args()

    genes = pd.read_csv(args.comparison)
    genes['ratio'] = genes['hedges_g_uninfl'] / genes['hedges_g_infl']

    candidates = genes[genes['category'] == 'constitutive']
    print(f'Genes significant in both analyses with no significant difference: {len(candidates)}')
    print(f'Ratio: median {candidates["ratio"].median():.2f}, '
          f'IQR {candidates["ratio"].quantile(0.25):.2f}-{candidates["ratio"].quantile(0.75):.2f}')

    reference = constitutive(genes, *REFERENCE)
    rows = []
    for low, high in WINDOWS:
        selected = constitutive(genes, low, high)
        rows.append({
            'window': f'{low:.2f}-{high:.2f}',
            'genes': len(selected),
            'shared_with_reference': len(selected & reference),
            'nested': selected <= reference or reference <= selected,
        })

    summary = pd.DataFrame(rows)
    print()
    print(summary.to_string(index=False))
    print(f'\nAll sets nested: {bool(summary["nested"].all())}')

    os.makedirs(args.output_dir, exist_ok=True)
    summary.to_csv(f'{args.output_dir}/ratio_window_sensitivity.csv', index=False)
    print(f'Written to {args.output_dir}/ratio_window_sensitivity.csv')


if __name__ == '__main__':
    main()
