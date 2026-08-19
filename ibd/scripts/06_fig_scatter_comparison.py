"""
Scatter plot of the two activity layers present in uninflamed mucosa.

Each point is a gene assigned by the activity-layer triage (S_gene_layers.py).
Both layers are altered in uninflamed UC; they differ in whether the direct
inflamed-versus-uninflamed contrast finds a further change. Vertical distance
below the y = x diagonal is how much inflammation adds.

The constitutive layer is defined by the absence of a significant further
change, not by numerical equality of the two effects, so these genes scatter
below the diagonal rather than sitting on it: the direct contrast is drawn from
six cohorts and is underpowered to resolve differences of this size.

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

LAYERS = [
    ("inflammation-amplified", "#d6604d", "Inflammation-amplified"),
    ("constitutive", "#2166ac", "Constitutive"),
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="results/revision_v2/gene_layers.csv",
                        help="gene_layers.csv from S_gene_layers.py")
    parser.add_argument("--output", default="figures/scatter_comparison.pdf")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    genes = pd.read_csv(args.input)
    shown = genes[genes.layer.isin([name for name, _, _ in LAYERS])]

    fig, ax = plt.subplots(figsize=(7, 7))

    for name, colour, label in LAYERS:
        sub = shown[shown.layer == name]
        ax.scatter(sub.hedges_g_infl, sub.hedges_g_uninf,
                   s=55, alpha=0.85, c=colour, edgecolors="k", linewidths=0.4,
                   zorder=3, label=f"{label} ($n = {len(sub)}$)")

    # The constitutive layer is the one the text discusses gene by gene.
    texts = [ax.text(r.hedges_g_infl, r.hedges_g_uninf, r.gene, fontsize=8, color="k", fontweight="medium", zorder=5)
             for r in shown[shown.layer == "constitutive"].itertuples()]
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
