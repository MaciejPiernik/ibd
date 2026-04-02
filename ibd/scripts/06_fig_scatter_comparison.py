"""
Scatter plot of constitutive genes: inflamed vs uninflamed effect sizes.

Genes colored by ratio deviation from 1 (proximity to y=x diagonal).
Labels highlight high-confidence constitutive genes (ratio ~1, large effect).

Usage:
    python ibd/scripts/06_fig_scatter_comparison.py \
        --input results/paper/comparison_stats.csv \
        --output figures/scatter_comparison.pdf
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="comparison_stats.csv")
    parser.add_argument("--output", default="figures/scatter_comparison.pdf")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    fin = df[
        np.isfinite(df.hedges_g_infl) & np.isfinite(df.hedges_g_uninfl)
    ].copy()
    const = fin[fin["category"] == "constitutive"].copy()
    const["ratio"] = const["hedges_g_uninfl"] / const["hedges_g_infl"]
    ratio_dev = (const["ratio"] - 1).abs()

    fig, ax = plt.subplots(figsize=(7, 7))

    sc = ax.scatter(
        const.hedges_g_infl, const.hedges_g_uninfl,
        c=ratio_dev, cmap="RdYlBu_r", s=55, alpha=0.85, zorder=3,
        edgecolors="k", linewidths=0.4, vmin=0, vmax=0.8,
    )

    # ── Labels: near-diagonal with large effect ──
    to_label = const[const["ratio"].between(0.85, 1.15)]

    texts = []
    for _, r in to_label.iterrows():
        texts.append(
            ax.text(
                r.hedges_g_infl, r.hedges_g_uninfl, r.gene,
                fontsize=8, color="k", fontweight="medium", zorder=5,
            )
        )

    adjust_text(
        texts, ax=ax,
        x=const.hedges_g_infl.values,
        y=const.hedges_g_uninfl.values,
        arrowprops=dict(arrowstyle="-", color="grey", lw=0.5, alpha=0.8),
        expand=(5, 5),
        force_text=(1, 1),
        # max_move=None,
        # prevent_crossings=True,
        # zorder=100
    )

    lim = max(const.hedges_g_infl.abs().max(),
              const.hedges_g_uninfl.abs().max()) * 1.15
    ax.plot([-lim, lim], [-lim, lim], "k--", alpha=0.4, lw=1, label="$y = x$")
    ax.axhline(0, color="grey", lw=0.5, alpha=0.3)
    ax.axvline(0, color="grey", lw=0.5, alpha=0.3)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Hedges' $g$ (inflamed UC vs HC)", fontsize=12)
    ax.set_ylabel("Hedges' $g$ (uninflamed UC vs HC)", fontsize=12)
    ax.set_aspect("equal")
    ax.legend(fontsize=9, loc="upper left")

    cb = fig.colorbar(sc, ax=ax, shrink=0.75, pad=0.02)
    cb.set_label(
        "$|g_{\\mathrm{uninfl}} \\,/\\, g_{\\mathrm{infl}} - 1|$",
        fontsize=11,
    )

    plt.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved → {args.output}")


if __name__ == "__main__":
    main()