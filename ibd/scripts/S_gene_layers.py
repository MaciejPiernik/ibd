"""
Classify genes into the two activity layers, graded by RNA-seq support.

Each gene is judged on the three UC contrasts:

  A  inflamed UC vs controls
  B  uninflamed UC vs controls
  C  inflamed vs uninflamed UC (the direct contrast)

and carries two independent properties, each a positive claim:

  constitutive            significant in A and in B, same direction: the change
                          is already present before inflammation is
  inflammation-dependent  significant in A and in C, same direction:
                          inflammation measurably changes the gene

These are not exclusive. A gene present at baseline that inflammation pushes
further carries both, and is reported as a subset of each rather than as a
third category. A gene significant in A but carrying neither property is
unresolved: the data place it in no layer, which is stated rather than assumed
to mean absence.

C is not an independent axis --- it tracks (g_A - g_B) closely --- but inflamed
and uninflamed biopsies are paired within datasets, so it measures that axis
more precisely than the difference of two separate analyses, and it lets
inflammation dependence be asserted positively instead of inferred from a
non-significant B.

The RNA-seq cohorts grade each claim rather than filtering it. Requiring
significance on both platforms would discard a signal that replicates in
direction but not in significance, since the uninflamed RNA-seq arm has three
cohorts. Each property is therefore graded on the contrast that defines it:
B for constitutive, C for inflammation-dependent.

Run:
    python ibd/scripts/S_gene_layers.py
"""

import argparse
import math
import os

import numpy as np
import pandas as pd

CONTRASTS = [
    ('infl', 'uc_vs_hc', 'rnaseq_uc_vs_hc'),
    ('uninf', 'uninf_vs_hc', 'rnaseq_uninf_vs_hc'),
    ('direct', 'infl_vs_uninfl', 'rnaseq_infl_vs_uninfl'),
]
# Every contrast the grading depends on is carried through, so a support call can
# be audited from this file alone.
OUTPUT = ['gene', 'constitutive', 'inflammation_dependent', 'unresolved', 'direction',
          'support_constitutive', 'support_inflammation_dependent',
          'hedges_g_infl', 'q_value_infl', 'hedges_g_rna_infl', 'q_value_rna_infl',
          'hedges_g_uninf', 'q_value_uninf', 'hedges_g_rna_uninf', 'q_value_rna_uninf',
          'hedges_g_direct', 'q_value_direct', 'hedges_g_rna_direct', 'q_value_rna_direct']


def load(results, path, suffix, quiet=False):
    """Meta-analysis results plus a significance flag on the paper's criteria."""
    table = pd.read_csv(f'{results}/{path}/gene_meta_all.csv')[
        ['gene', 'hedges_g', 'q_value', 'datasets', 'direction_ratio']]
    minimum = math.ceil(table['datasets'].max() / 2)
    table['significant'] = (
        (table['q_value'] < 0.05)
        & (table['direction_ratio'] >= 0.8)
        & (table['datasets'] >= minimum)
    )
    if not quiet:
        print(f'  {path:22s} {int(table["datasets"].max()):2d} cohorts, presence >= {minimum}, '
              f'{int(table["significant"].sum()):5d} significant')

    keep = ['gene', 'hedges_g', 'q_value', 'significant']
    return table[keep].rename(columns={c: f'{c}_{suffix}' for c in keep[1:]})


def grade(genes, holds, contrasts):
    """How far the RNA-seq cohorts support a claim.

    A layer is defined by a conjunction of contrasts, so the claim is only
    fully supported when the RNA-seq cohorts confirm every one of them. Grading
    on the layer-defining contrast alone would call a gene strongly supported
    while the other half of its definition failed to replicate.
    """
    agree, replicated, tested = [], [], []
    for contrast in contrasts:
        micro = genes[f'hedges_g_{contrast}']
        rna = genes[f'hedges_g_rna_{contrast}']
        same_sign = np.sign(micro) == np.sign(rna)
        agree.append(same_sign)
        replicated.append(same_sign & genes[f'significant_rna_{contrast}'].fillna(False))
        tested.append(rna.notna())

    all_agree = pd.concat(agree, axis=1).all(axis=1)
    all_tested = pd.concat(tested, axis=1).all(axis=1)
    n_replicated = pd.concat(replicated, axis=1).sum(axis=1)

    support = pd.Series('', index=genes.index)
    support[holds] = 'unsupported'
    support[holds & all_agree] = 'direction-consistent'
    support[holds & all_agree & (n_replicated > 0)] = 'partial'
    support[holds & all_agree & (n_replicated == len(contrasts))] = 'strong'
    support[holds & ~all_tested] = 'not tested'

    return support


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', default='results/revision_v2')
    parser.add_argument('--output-dir', default='results/revision_v2')
    args = parser.parse_args()

    print('Microarray (the discovery cohorts):')
    genes = None
    for name, micro, _ in CONTRASTS:
        part = load(args.results, micro, name)
        genes = part if genes is None else genes.merge(part, on='gene')

    print('RNA-seq (grading only, never filtering):')
    for name, _, rna in CONTRASTS:
        # Left join: a gene missing from the RNA-seq corpus is untested, not refuted.
        genes = genes.merge(load(args.results, rna, f'rna_{name}'), on='gene', how='left')

    sign = np.sign(genes['hedges_g_infl'])
    genes['constitutive'] = (
        genes['significant_infl'] & genes['significant_uninf']
        & (sign == np.sign(genes['hedges_g_uninf'])))
    genes['inflammation_dependent'] = (
        genes['significant_infl'] & genes['significant_direct']
        & (sign == np.sign(genes['hedges_g_direct'])))
    genes['unresolved'] = (
        genes['significant_infl']
        & ~genes['constitutive'] & ~genes['inflammation_dependent'])

    genes['direction'] = pd.Series(np.where(sign > 0, 'up', 'down'), index=genes.index).where(
        genes['constitutive'] | genes['inflammation_dependent'] | genes['unresolved'])
    # Each layer is graded on the whole conjunction that defines it.
    genes['support_constitutive'] = grade(
        genes, genes['constitutive'], ['infl', 'uninf'])
    genes['support_inflammation_dependent'] = grade(
        genes, genes['inflammation_dependent'], ['infl', 'direct'])

    both = genes['constitutive'] & genes['inflammation_dependent']
    print(f'\nClassification over {len(genes)} genes measured in all three UC contrasts:')
    print(f'  constitutive            {int(genes["constitutive"].sum()):6d}   '
          f'of which {int(both.sum())} are also inflammation-dependent')
    print(f'  inflammation-dependent  {int(genes["inflammation_dependent"].sum()):6d}   '
          f'of which {int(both.sum())} are also constitutive')
    print(f'  unresolved              {int(genes["unresolved"].sum()):6d}   '
          f'altered when inflamed, placed in neither layer')

    print('\nRNA-seq support:')
    for label, column in [('constitutive', 'support_constitutive'),
                          ('inflammation-dependent', 'support_inflammation_dependent')]:
        counts = genes.loc[genes[column] != '', column].value_counts()
        detail = ', '.join(f'{counts.get(k, 0)} {k}'
                           for k in ['strong', 'partial', 'direction-consistent',
                                     'unsupported', 'not tested'])
        print(f'  {label:24s} {detail}')

    genes = genes.sort_values(
        ['constitutive', 'inflammation_dependent', 'hedges_g_infl'],
        key=lambda c: c.abs() if c.name == 'hedges_g_infl' else c,
        ascending=[False, False, False])

    os.makedirs(args.output_dir, exist_ok=True)
    genes[OUTPUT].to_csv(f'{args.output_dir}/gene_layers.csv', index=False)
    print(f'\nWritten to {args.output_dir}/gene_layers.csv')


if __name__ == '__main__':
    main()
