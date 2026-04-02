"""
Volcano plot for the inflamed UC vs control gene-level meta-analysis.

Produces a publication-quality figure showing Hedges' g vs -log10(q) for all
tested genes, with significant genes color-coded by direction and the most
extreme genes labelled.

Usage:
    python ibd/scripts/volcano_plot.py
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

# ── style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "figure.dpi": 300,
})

UP_COLOR = "#c0392b"
DOWN_COLOR = "#2471a3"
NS_COLOR = "#bdc3c7"
LABEL_UP_COLOR = "#922b21"
LABEL_DOWN_COLOR = "#1a5276"

# Genes to label — top extremes from the paper tables + key discussion genes
LABEL_GENES_UP = [
    "LPCAT1", "TIMP1", "CD55", "SLC6A14", "CXCL1", "BACE2", "CFB",
    "DUOX2", "KYNU", "IFITM3", "PI3", "LCN2",
]
LABEL_GENES_DOWN = [
    "PDE6A", "DPP10", "UGT2A3", "SLC26A2", "ACSF2", "AQP8",
    "CLDN8", "SLC22A5", "PPARGC1A",
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
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="results/paper/uc_vs_hc",
                     help="Directory containing gene_meta_all.csv")
    ap.add_argument("--output", default="figures/volcano.pdf")
    args = ap.parse_args()

    df = pd.read_csv(f"{args.data_dir}/gene_meta_all.csv")
    df["neg_log10_q"] = -np.log10(df["q_value"].clip(lower=1e-300))
    df = df[df["gene"].isin(get_protein_coding_genes())]

    label_genes_top = df.sort_values("hedges_g").head(10)["gene"].tolist() + \
                      df.sort_values("hedges_g").tail(10)["gene"].tolist()
    
    # label_genes_top = [g for g in label_genes_top if g not in LABEL_GENES_UP and g not in LABEL_GENES_DOWN]

    # Determine minimum dataset coverage (≥50% of total datasets in analysis)
    min_datasets = df["datasets"].max() / 2

    sig = ((df["q_value"] < 0.05)
           & (df["direction_ratio"] >= 0.8)
           & (df["datasets"] >= min_datasets))
    up = sig & (df["direction"] == "up")
    down = sig & (df["direction"] == "down")
    ns = ~(up | down)

    fig, ax = plt.subplots(figsize=(6, 4.5))

    ax.scatter(df.loc[ns, "hedges_g"], df.loc[ns, "neg_log10_q"],
               s=3, c=NS_COLOR, alpha=0.4, linewidths=0, rasterized=True)
    ax.scatter(df.loc[up, "hedges_g"], df.loc[up, "neg_log10_q"],
               s=5, c=UP_COLOR, alpha=0.5, linewidths=0, rasterized=True)
    ax.scatter(df.loc[down, "hedges_g"], df.loc[down, "neg_log10_q"],
               s=5, c=DOWN_COLOR, alpha=0.5, linewidths=0, rasterized=True)

    # ── labels ─────────────────────────────────────────────────────────────
    texts = []
    label_genes = LABEL_GENES_UP + LABEL_GENES_DOWN
    for gene in label_genes:
        row = df[df["gene"] == gene]
        if row.empty:
            continue
        x = row["hedges_g"].values[0]
        y = row["neg_log10_q"].values[0]
        color = LABEL_UP_COLOR if gene in LABEL_GENES_UP else LABEL_DOWN_COLOR
        ax.scatter(x, y, s=18, c=color, zorder=5, linewidths=0.4,
                   edgecolors="white")
        texts.append(ax.text(x, y, gene, fontsize=6.5, color=color,
                             fontstyle="italic", fontweight="medium"))

    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="0.5", lw=0.4), force_text=(0.4, 0.6), expand=(1.2, 1.4))

    # ── threshold lines ───────────────────────────────────────────────────
    ax.axhline(-np.log10(0.05), color="0.6", ls="--", lw=0.5)

    # ── counts annotation ─────────────────────────────────────────────────
    n_up = up.sum()
    n_down = down.sum()
    ax.text(0.98, 0.97,
            f"{n_up:,} up\n{n_down:,} down",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox=dict(facecolor="white", edgecolor="0.7", pad=3, lw=0.5))

    ax.set_xlabel("Hedges' $g$ (meta-analytic effect size)")
    ax.set_ylabel("$-\\log_{10}\\,q$")
    ax.set_title("Inflamed UC vs controls — gene-level meta-analysis",
                 fontsize=10, fontweight="bold", pad=8)

    # keep axes symmetric and tidy
    xlim = max(abs(df["hedges_g"]))
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(-0.5, df["neg_log10_q"].quantile(0.999) * 1.05)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(args.output, bbox_inches="tight")
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
