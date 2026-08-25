"""
Scatter plot of the constitutive layer, graded by RNA-seq support.

Each point is a protein-coding constitutive gene from S_gene_layers.py: altered
against controls both when inflamed and when uninflamed, in the same direction.
Non-coding symbols are excluded from the plot entirely.

Distance below the y = x diagonal is how much larger the effect is once
inflammation is present, so the genes furthest below it are those that also
carry the inflammation-dependent property. Colour grades how well the
independent RNA-seq cohorts support the uninflamed effect --- the claim that
makes a gene constitutive, and the weakest link in the analysis.

Usage:
    python ibd/scripts/06_fig_scatter_comparison.py \
        --input results/revision_v2/gene_layers.csv \
        --output figures/scatter_comparison.pdf
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from adjustText import adjust_text

# Graded by how well the RNA-seq cohorts support the uninflamed effect, which is
# the contrast that makes a gene constitutive.
SUPPORT = [
    ("strong", "#2166ac", "Both contrasts replicated"),
    ("partial", "#67a9cf", "One contrast replicated"),
    ("direction-consistent", "#d1e5f0", "Same direction, neither significant"),
    ("unsupported", "#d6604d", "Not supported"),
    ("not tested", "#bdbdbd", "Not tested in RNA-seq"),
]


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="results/revision_v2/gene_layers.csv",
                        help="gene_layers.csv from S_gene_layers.py")
    parser.add_argument("--output", default="figures/scatter_comparison.pdf")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--labels", type=int, default=10,
                        help="how many of the largest effects to name, on top of the "
                             "genes the RNA-seq cohorts settle either way")
    parser.add_argument("--hgnc", default=None,
                        help="Path to HGNC TSV; if omitted, downloads from genenames.org")
    args = parser.parse_args()

    genes = pd.read_csv(args.input)
    genes = genes[genes.gene.isin(get_protein_coding_genes(args.hgnc))]
    shown = genes[genes.constitutive]

    fig, ax = plt.subplots(figsize=(7, 7))

    for name, colour, label in SUPPORT:
        sub = shown[shown.support_constitutive == name]
        if sub.empty:
            continue
        ax.scatter(sub.hedges_g_infl, sub.hedges_g_uninf,
                   s=55, alpha=0.85, c=colour, edgecolors="k", linewidths=0.4,
                   zorder=3, label=f"{label} ($n = {len(sub)}$)")

    # Labelling all of them is unreadable, so name the genes a reader can act on:
    # those the RNA-seq cohorts settle either way, plus the largest effects.
    settled = shown[shown.support_constitutive.isin(["strong", "unsupported"])]
    largest = shown.reindex(shown.hedges_g_uninf.abs().sort_values(ascending=False).index)
    labelled = pd.concat([settled, largest.head(args.labels)]).drop_duplicates("gene")
    texts = [ax.text(r.hedges_g_infl, r.hedges_g_uninf, r.gene, fontsize=8, color="k",
                     fontweight="medium", zorder=5)
             for r in labelled.itertuples()]
    adjust_text(
        texts, ax=ax,
        x=shown.hedges_g_infl.values, y=shown.hedges_g_uninf.values,
        arrowprops=dict(arrowstyle="-", color="grey", lw=0.5, alpha=0.8),
        expand=(3.3, 3.5), force_text=(0.4, 0.6),
    )

    lim = max(shown.hedges_g_infl.abs().max(),
              shown.hedges_g_uninf.abs().max()) * 1.15
    ax.plot([-lim, lim], [-lim, lim], "k--", alpha=0.4, lw=1, label="$y = x$")
    ax.axhline(0, color="grey", lw=0.5, alpha=0.3)
    ax.axvline(0, color="grey", lw=0.5, alpha=0.3)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Hedges' $g$ (inflamed UC vs HC)", fontsize=12)
    ax.set_ylabel("Hedges' $g$ (uninflamed UC vs HC)", fontsize=12)
    ax.set_aspect("equal")
    ax.legend(fontsize=9, loc="upper left")

    plt.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved → {args.output}")


if __name__ == "__main__":
    main()
