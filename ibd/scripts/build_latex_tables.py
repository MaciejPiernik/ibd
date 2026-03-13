"""
Generate gene-level LaTeX tables for the UC meta-analysis paper.

Outputs:
  1. Summary statistics (LaTeX comments for Results text)
  2. Top upregulated genes table
  3. Top downregulated genes table

Usage:
    python ibd/scripts/build_latex_tables.py \
        --meta results/paper/uc_vs_hc/gene_meta_all.csv \
        --rra results/paper/uc_vs_hc/gene_rra.csv \
        --pathways results/paper/uc_vs_hc/pathway_meta_significant.csv \
        --output-dir results/paper
"""

import argparse
import os

import pandas as pd
from scipy.stats import spearmanr

MIN_DATASETS = 6
DIRECTION_THRESHOLD = 0.8
ALPHA = 0.05


def get_protein_coding_genes(path: str | None = None) -> set[str]:
    """Return a set of HGNC-approved protein-coding gene symbols.

    Parameters
    ----------
    path : str or None
        Local path to an HGNC custom-download TSV (must contain columns
        ``Approved symbol`` and ``Locus group``).  If *None*, downloads
        the current list from genenames.org.
    """
    HGNC_URL = (
        'https://www.genenames.org/cgi-bin/download/custom?'
        'col=gd_app_sym&col=gd_locus_group&status=Approved&'
        'hgnc_datea=&hgnc_dateb=&order_by=gd_app_sym_sort&'
        'format=text&submit=submit'
    )
    src = path if path else HGNC_URL
    df = pd.read_csv(src, sep='\t')
    return set(df.loc[df['Locus group'] == 'protein-coding gene', 'Approved symbol'])


def _escape_latex(s: str) -> str:
    return s.replace('_', r'\_').replace('&', r'\&').replace('%', r'\%')


def _shorten_theme(s: str, maxlen: int = 1000) -> str:
    s = s.split(' R-HSA-')[0] if ' R-HSA-' in s else s
    if len(s) > maxlen:
        s = s[:maxlen - 1] + '…'
    return s


def build_gene_to_pathway(pw_knee: pd.DataFrame) -> dict:
    """Map each gene to its highest-|NES| pathway via all leading-edge genes."""
    gene_best = {}  # gene -> (abs_nes, pathway_name)
    for _, prow in pw_knee.iterrows():
        term = prow['term'].split(' R-HSA-')[0]
        abs_nes = abs(prow['combined_nes'])
        for col in ['lead_genes_all', 'lead_genes_core']:
            gs = str(prow.get(col, ''))
            if gs and gs != 'nan':
                for g in gs.split(';'):
                    g = g.strip()
                    if g:
                        prev = gene_best.get(g, (0, None))
                        if abs_nes > prev[0]:
                            gene_best[g] = (abs_nes, term)
    return {g: name for g, (_, name) in gene_best.items() if name}


def make_top_genes_table(df: pd.DataFrame, direction: str, n: int = 20) -> str:
    mask = ~df['gene'].str.match(r'^(MIR\d|LOC\d|LINC\d|SNORD|SCARNA)')
    sub = df[mask & (df['direction'] == direction)].nlargest(n, 'abs_g')

    label = 'upregulated' if direction == 'up' else 'downregulated'

    lines = []
    lines.append(r'\begin{table*}[t]')
    lines.append(r'\centering')
    lines.append(r'\scriptsize')
    lines.append(r'\renewcommand{\arraystretch}{1.15}')
    lines.append(r'\begin{tabular}{lrl}')
    lines.append(r'\hline')
    lines.append(r"Gene & $g$ & Pathway \\")
    lines.append(r'\hline')

    for _, row in sub.iterrows():
        gene = _escape_latex(row['gene'])
        g = row['hedges_g']
        theme = row.get('pathway', None)
        if pd.isna(theme) or theme is None:
            theme = '---'
        else:
            theme = _shorten_theme(theme)
        theme = _escape_latex(theme)
        lines.append(f'  {gene} & {g:+.2f} & {theme} \\\\')

    lines.append(r'\hline')
    lines.append(r'\end{tabular}')
    lines.append(
        r'\caption{Top ' + str(n) + ' ' + label +
        r" genes by Hedges' $g$ from random-effects meta-analysis. "
        r"Pathway indicates the highest-scoring enriched Reactome pathway containing "
        r"the gene in its leading edge; ``---'' denotes genes not present in any enriched pathway.}")
    lines.append(r'\label{tab:top_' + direction + '}')
    lines.append(r'\end{table*}')
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta', default='results/uc_vs_hc/gene_meta_all.csv', help='Gene-level meta-analysis results')
    parser.add_argument('--rra', default='results/uc_vs_hc/gene_rra.csv', help='RRA results')
    parser.add_argument('--pathways', default='results/uc_vs_hc/pathway_meta_all.csv', help='Pathway results')
    parser.add_argument('--output-dir', default='results/uc_vs_hc', help='Directory to save LaTeX tables')
    parser.add_argument('--top-n', type=int, default=20)
    parser.add_argument('--encoding', default='latin1')
    parser.add_argument('--hgnc', default=None,
                        help='Path to HGNC TSV; if omitted, downloads from genenames.org')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    meta_all = pd.read_csv(args.meta, encoding=args.encoding)
    rra = pd.read_csv(args.rra, encoding=args.encoding)
    pw_knee = pd.read_csv(args.pathways, encoding=args.encoding)

    # Keep only protein-coding genes
    pc = get_protein_coding_genes(args.hgnc)

    # Filter significant
    meta_sig = meta_all[
        (meta_all['datasets'] >= MIN_DATASETS)
        & (meta_all['direction_ratio'] >= DIRECTION_THRESHOLD)
        & (meta_all['q_value'] < ALPHA)
    ].copy()
    meta_sig['abs_g'] = meta_sig['hedges_g'].abs()

    # Map genes to best pathway
    gene_to_pw = build_gene_to_pathway(pw_knee)
    meta_sig['pathway'] = meta_sig['gene'].map(gene_to_pw)

    # RRA concordance
    merged = meta_sig.merge(rra[['gene', 'rra_rho', 'rra_q']], on='gene', how='left')
    rra_r, _ = spearmanr(
        merged['abs_g'].rank(ascending=False),
        merged['rra_rho'].rank(ascending=True),
    )

    n_up = (meta_sig['direction'] == 'up').sum()
    n_down = (meta_sig['direction'] == 'down').sum()
    n_mapped = meta_sig['pathway'].notna().sum()

    summary = '\n'.join([
        '% ── Summary statistics for Results section ──',
        f'% Total genes in meta-analysis: {len(meta_all)}',
        f'% Significant genes: {len(meta_sig)}',
        f'%   Upregulated: {n_up}',
        f'%   Downregulated: {n_down}',
        f'%   Median datasets per gene: {meta_sig["datasets"].median():.0f}',
        f'%   Median |Hedges\' g|: {meta_sig["abs_g"].median():.2f}',
        f'%   Median I^2: {meta_sig["I2"].median():.2f}',
        f'% RRA concordance (Spearman rho): {rra_r:.3f}',
        f'% Genes mapped to a pathway: {n_mapped}',
        f'% Genes not in any pathway: {len(meta_sig) - n_mapped}',
    ])

    merged = merged[merged['gene'].isin(pc)]

    table_up = make_top_genes_table(merged, 'up', n=args.top_n)
    table_down = make_top_genes_table(merged, 'down', n=args.top_n)

    out_path = os.path.join(args.output_dir, 'gene_tables.tex')
    with open(out_path, 'w') as f:
        f.write(summary + '\n\n')
        f.write(table_up + '\n\n')
        f.write(table_down + '\n')

    print(summary.replace('% ', ''))
    print(f'\nLaTeX tables written to {out_path}')


if __name__ == '__main__':
    main()