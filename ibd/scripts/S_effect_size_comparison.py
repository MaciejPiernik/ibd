"""
Compare inflamed-vs-HC and uninflamed-vs-HC meta-analyses.

For each gene present in both analyses, computes a Wald test for the
difference in effect sizes:

    z = (g_infl - g_uninfl) / sqrt(SE_infl^2 + SE_uninfl^2)

Outputs:
  1. comparison_stats.csv    - full gene-level comparison table
  2. comparison_summary.tex  - LaTeX summary stats + top constitutive/
                               inflammation-dependent gene tables
  3. Printed summary to stdout

Usage:
    python ibd/scripts/0X_effect_size_comparison.py \
        --inflamed   results/paper/uc_vs_hc/gene_meta_all.csv \
        --uninflamed results/paper/uninf_vs_hc/gene_meta_all.csv \
        --output-dir results/paper/
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import norm
def _bh_correction(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    order = np.argsort(pvals)
    qvals = np.empty(n)
    qvals[order] = np.minimum(1, np.minimum.accumulate(
        (pvals[order] * n / np.arange(1, n + 1))[::-1]
    )[::-1])
    return qvals


def _escape_latex(s: str) -> str:
    return s.replace('_', r'\_').replace('&', r'\&').replace('%', r'\%')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inflamed', required=True, help='gene_meta_all.csv from inflamed vs HC')
    parser.add_argument('--uninflamed', required=True, help='gene_meta_all.csv from uninflamed vs HC')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    parser.add_argument('--top-n', type=int, default=20, help='Genes per table')
    parser.add_argument('--encoding', default='latin1')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    infl = pd.read_csv(args.inflamed, encoding=args.encoding)
    uninfl = pd.read_csv(args.uninflamed, encoding=args.encoding)

    # Merge on gene
    merged = infl.merge(uninfl, on='gene', suffixes=('_infl', '_uninfl'))
    print(f"Genes in inflamed analysis:   {len(infl)}")
    print(f"Genes in uninflamed analysis: {len(uninfl)}")
    print(f"Genes in both:                {len(merged)}")

    # Wald test for difference
    merged['delta_g'] = merged['hedges_g_infl'] - merged['hedges_g_uninfl']
    merged['delta_se'] = np.sqrt(merged['se_infl']**2 + merged['se_uninfl']**2)
    merged['delta_z'] = merged['delta_g'] / merged['delta_se']
    merged['delta_p'] = 2 * norm.sf(np.abs(merged['delta_z']))

    # BH correction
    mask = np.isfinite(merged['delta_p'])
    merged.loc[~mask, 'delta_q'] = np.nan
    merged.loc[mask, 'delta_q'] = _bh_correction(merged.loc[mask, 'delta_p'].values)

    # Classification
    # Significance in each analysis
    # Inflamed: q<0.05, direction consistency >=0.8, present in >=50% datasets (6/12)
    # Uninflamed: q<0.05, present in >=50% datasets (3/5)
    merged['sig_infl'] = (
        (merged['q_value_infl'] < 0.05)
        & (merged['direction_ratio_infl'] >= 0.8)
        & (merged['datasets_infl'] >= 6)
    )
    merged['sig_uninfl'] = (
        (merged['q_value_uninfl'] < 0.05)
        & (merged['datasets_uninfl'] >= 3)
    )

    # Categories
    merged['category'] = 'not significant'
    merged.loc[
        merged['sig_infl'] & ~merged['sig_uninfl'],
        'category'
    ] = 'inflammation-dependent'
    merged.loc[
        merged['sig_infl'] & merged['sig_uninfl'] & (merged['delta_q'] < 0.05),
        'category'
    ] = 'inflammation-amplified'
    merged.loc[
        merged['sig_infl'] & merged['sig_uninfl'] & (merged['delta_q'] >= 0.05),
        'category'
    ] = 'constitutive'
    merged.loc[
        ~merged['sig_infl'] & merged['sig_uninfl'],
        'category'
    ] = 'uninflamed-only'

    # Filter non-coding for display
    coding_mask = ~merged['gene'].str.match(r'^(MIR\d|LOC\d|LINC\d|SNORD|SCARNA)')

    # ── Summary stats ──
    print(f"\n{'='*60}")
    print("COMPARISON SUMMARY")
    print(f"{'='*60}")
    for cat in ['constitutive', 'inflammation-amplified', 'inflammation-dependent',
                 'uninflamed-only', 'not significant']:
        n = (merged['category'] == cat).sum()
        print(f"  {cat:30s}: {n:,}")

    # Direction agreement among genes significant in both
    both_sig = merged[merged['sig_infl'] & merged['sig_uninfl']]
    same_dir = (both_sig['direction_infl'] == both_sig['direction_uninfl']).sum() if len(both_sig) > 0 else 0
    print(f"\nGenes significant in both: {len(both_sig)}")
    if len(both_sig) > 0:
        print(f"  Same direction: {same_dir} ({same_dir/len(both_sig)*100:.1f}%)")
    else:
        print(f"  Same direction: N/A")

    # Correlation of effect sizes
    from scipy.stats import spearmanr, pearsonr
    both_finite = merged[
        np.isfinite(merged['hedges_g_infl']) & np.isfinite(merged['hedges_g_uninfl'])
    ]
    r_pearson, _ = pearsonr(both_finite['hedges_g_infl'], both_finite['hedges_g_uninfl'])
    r_spearman, _ = spearmanr(both_finite['hedges_g_infl'], both_finite['hedges_g_uninfl'])
    print(f"\nEffect size correlation (all shared genes):")
    print(f"  Pearson r:  {r_pearson:.3f}")
    print(f"  Spearman r: {r_spearman:.3f}")

    # ── Save full table ──
    out_cols = [
        'gene',
        'hedges_g_infl', 'se_infl', 'q_value_infl', 'I2_infl',
        'datasets_infl', 'direction_infl',
        'hedges_g_uninfl', 'se_uninfl', 'q_value_uninfl', 'I2_uninfl',
        'datasets_uninfl', 'direction_uninfl',
        'delta_g', 'delta_se', 'delta_z', 'delta_p', 'delta_q',
        'sig_infl', 'sig_uninfl', 'category',
    ]
    csv_path = os.path.join(args.output_dir, 'comparison_stats.csv')
    merged[out_cols].to_csv(csv_path, index=False)
    print(f"\nFull comparison table: {csv_path}")

    # ── LaTeX tables ──
    n = args.top_n
    lines = []

    # Summary as comments
    lines.append('% ── Inflamed vs Uninflamed comparison summary ──')
    lines.append(f'% Genes in both analyses: {len(merged)}')
    for cat in ['constitutive', 'inflammation-amplified', 'inflammation-dependent',
                 'uninflamed-only', 'not significant']:
        lines.append(f'%   {cat}: {(merged["category"]==cat).sum()}')
    lines.append(f'% Pearson r (effect sizes): {r_pearson:.3f}')
    lines.append(f'% Spearman r (effect sizes): {r_spearman:.3f}')
    lines.append(f'% Genes sig in both, same direction: {same_dir}/{len(both_sig)} '
                 f'({same_dir/max(len(both_sig),1)*100:.1f}%)')
    lines.append('')

    # Top constitutive genes (significant in both, delta_q >= 0.05, by |g_infl|)
    const = merged[coding_mask & (merged['category'] == 'constitutive')].copy()
    const['abs_g'] = const['hedges_g_infl'].abs()
    const_top = const.nlargest(n, 'abs_g')

    lines.append(r'\begin{table*}[t]')
    lines.append(r'\centering')
    lines.append(r'\scriptsize')
    lines.append(r'\renewcommand{\arraystretch}{1.15}')
    lines.append(r'\begin{tabular}{lrrrrr}')
    lines.append(r'\hline')
    lines.append(r"Gene & $g_{\mathrm{infl}}$ & $g_{\mathrm{uninfl}}$ "
                 r"& $\Delta g$ & $q_{\Delta}$ & Direction \\")
    lines.append(r'\hline')
    for _, row in const_top.iterrows():
        gene = _escape_latex(row['gene'])
        lines.append(
            f'  {gene} & {row["hedges_g_infl"]:+.2f} & {row["hedges_g_uninfl"]:+.2f} '
            f'& {row["delta_g"]:+.2f} & {row["delta_q"]:.2f} & {row["direction_infl"]} \\\\'
        )
    lines.append(r'\hline')
    lines.append(r'\end{tabular}')
    lines.append(r'\caption{Top ' + str(n) + r' constitutive genes: significantly dysregulated '
                 r'in both inflamed and uninflamed UC relative to healthy controls, '
                 r'with no significant difference in effect size between conditions '
                 r'($q_{\Delta} \geq 0.05$, Wald test). Ranked by $|g_{\mathrm{infl}}|$.}')
    lines.append(r'\label{tab:constitutive}')
    lines.append(r'\end{table*}')
    lines.append('')

    # Top inflammation-dependent genes (sig in inflamed only, by |g_infl|)
    dep = merged[coding_mask & (merged['category'] == 'inflammation-dependent')].copy()
    dep['abs_g'] = dep['hedges_g_infl'].abs()
    dep_top = dep.nlargest(n, 'abs_g')

    lines.append(r'\begin{table*}[t]')
    lines.append(r'\centering')
    lines.append(r'\scriptsize')
    lines.append(r'\renewcommand{\arraystretch}{1.15}')
    lines.append(r'\begin{tabular}{lrrrrr}')
    lines.append(r'\hline')
    lines.append(r"Gene & $g_{\mathrm{infl}}$ & $g_{\mathrm{uninfl}}$ "
                 r"& $\Delta g$ & $q_{\Delta}$ & Direction \\")
    lines.append(r'\hline')
    for _, row in dep_top.iterrows():
        gene = _escape_latex(row['gene'])
        dq = row['delta_q']
        dq_str = f'{dq:.2e}' if dq < 0.01 else f'{dq:.2f}'
        lines.append(
            f'  {gene} & {row["hedges_g_infl"]:+.2f} & {row["hedges_g_uninfl"]:+.2f} '
            f'& {row["delta_g"]:+.2f} & {dq_str} & {row["direction_infl"]} \\\\'
        )
    lines.append(r'\hline')
    lines.append(r'\end{tabular}')
    lines.append(r'\caption{Top ' + str(n) + r' inflammation-dependent genes: significantly '
                 r'dysregulated in inflamed UC but not in uninflamed UC relative to healthy '
                 r'controls. Ranked by $|g_{\mathrm{infl}}|$.}')
    lines.append(r'\label{tab:inflammation_dependent}')
    lines.append(r'\end{table*}')

    tex_path = os.path.join(args.output_dir, 'comparison_tables.tex')
    with open(tex_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"LaTeX tables: {tex_path}")


if __name__ == '__main__':
    main()