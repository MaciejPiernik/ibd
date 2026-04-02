#!/usr/bin/env python3
"""
generate_supplementary.py

Master script to generate ALL supplementary materials for the UC meta-analysis paper.

Outputs:
  output/supplementary.pdf                      — Supplementary document with figures and tables
  output/Table_S1_genes_uc_vs_hc.csv            — Gene-level results: UC inflamed vs HC
  output/Table_S2_genes_uninf_vs_hc.csv         — Gene-level results: UC uninflamed vs HC
  output/Table_S3_genes_cd_vs_hc.csv            — Gene-level results: CD inflamed vs HC
  output/Table_S4_genes_uc_vs_cd.csv            — Gene-level results: UC vs CD direct
  output/Table_S5_pathways_uc_vs_hc.csv         — Significant pathways: UC inflamed vs HC
  output/Table_S6_pathways_uninf_vs_hc.csv      — Significant pathways: UC uninflamed vs HC
  output/Table_S7_pathways_cd_vs_hc.csv         — Significant pathways: CD inflamed vs HC
  output/Table_S8_pathways_uc_vs_cd.csv         — Significant pathways: UC vs CD direct

Usage:
    python generate_supplementary.py --data-dir /path/to/data
"""

import argparse, os, sys, subprocess, colorsys, textwrap
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings("ignore")

# ── Styling ───────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--data-dir", default="results/paper/",
                   help="Root directory containing uc_vs_hc/, uninf_vs_hc/, cd_vs_hc/, uc_vs_cd/")
    p.add_argument("--output-dir", default="output")
    return p.parse_args()


# ═══════════════════════════════════════════════════════════════════════════════
#  HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def load_csv(path):
    """Load CSV, trying to handle space-prefixed filenames from macOS zips."""
    p = Path(path)
    if p.exists():
        return pd.read_csv(p, encoding="latin1")
    # Try with leading space
    alt = p.parent / f" {p.name}"
    if alt.exists():
        return pd.read_csv(alt, encoding="latin1")
    raise FileNotFoundError(f"Cannot find {path} or {alt}")


def clean_pathway_name(term):
    """Remove R-HSA-XXXXX suffix from pathway names."""
    return term.split(" R-HSA")[0]


def escape_latex(s):
    """Escape special LaTeX characters."""
    s = str(s)
    for ch, repl in [("&", r"\&"), ("_", r"\_"), ("%", r"\%"),
                      ("#", r"\#"), ("$", r"\$"), ("{", r"\{"),
                      ("}", r"\}"), ("~", r"\textasciitilde{}"),
                      ("^", r"\textasciicircum{}")]:
        s = s.replace(ch, repl)
    return s


def desaturate(rgba, factor=0.30):
    r, g, b, a = rgba
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    l_new = l + (0.85 - l) * (1 - factor)
    s_new = s * factor
    r2, g2, b2 = colorsys.hls_to_rgb(h, l_new, s_new)
    return (r2, g2, b2, a)


# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE S1: I² DIAGNOSTICS
# ═══════════════════════════════════════════════════════════════════════════════

def generate_figure_s1(data_dir, fig_dir):
    print("Generating Figure S1: I² diagnostics...")
    meta = load_csv(f"{data_dir}/uc_vs_hc/gene_meta_significant.csv")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    # Panel A: Histogram of I²
    ax = axes[0]
    i2_vals = meta.I2.values
    ax.hist(i2_vals, bins=50, color="#4393c3", edgecolor="white",
            linewidth=0.3, alpha=0.85)
    ax.axvline(np.median(i2_vals), color="#b2182b", linestyle="--",
               linewidth=1.2, label=f"Median = {np.median(i2_vals):.2f}")
    ax.set_xlabel("I²")
    ax.set_ylabel("Number of genes")
    ax.set_title("A.  Heterogeneity distribution", fontweight="bold", loc="left")
    ax.legend(frameon=False)

    # Panel B: |g| vs I²
    ax = axes[1]
    abs_g = np.abs(meta.hedges_g.values)
    i2 = meta.I2.values
    neg_log_q = -np.log10(meta.q_value.values + 1e-300)
    neg_log_q_cap = np.clip(neg_log_q, 0, 30)

    sc = ax.scatter(abs_g, i2, c=neg_log_q_cap, cmap="RdYlBu_r", s=4,
                    alpha=0.4, edgecolors="none", rasterized=True)
    cbar = plt.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("−log₁₀(q)")
    ax.set_xlabel("|Hedges' g|")
    ax.set_ylabel("I²")
    ax.set_title("B.  Effect size vs heterogeneity", fontweight="bold", loc="left")

    rho, _ = stats.spearmanr(abs_g, i2)
    ax.text(0.97, 0.05, f"Spearman ρ = {rho:.2f}",
            transform=ax.transAxes, fontsize=8, ha="right",
            bbox=dict(facecolor="white", edgecolor="#cccccc", alpha=0.9,
                      boxstyle="round,pad=0.3"))

    plt.tight_layout()
    out = f"{fig_dir}/figure_s1.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  → {out}")


# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE S2: FOREST PLOTS
# ═══════════════════════════════════════════════════════════════════════════════

def forest_plot_multi(genes, per_df, meta_df, outpath, total_width=12.6):
    """Multi-panel forest plot with fixed total figure width."""
    available = [g for g in genes if g in per_df.gene.values and g in meta_df.gene.values]
    if not available:
        print(f"  WARNING: no genes found for forest plot, skipping {outpath}")
        return
    n_genes = len(available)
    panel_width = total_width / n_genes
    fig, axes = plt.subplots(1, n_genes, figsize=(total_width, 5.5), sharey=False)
    if n_genes == 1:
        axes = [axes]

    for ax, gene_name in zip(axes, available):
        gdata = per_df[per_df.gene == gene_name].copy()
        gmeta = meta_df[meta_df.gene == gene_name].iloc[0]
        gdata = gdata.sort_values("hedges_g", ascending=True).reset_index(drop=True)

        y_positions = np.arange(len(gdata))
        ci_lo = gdata.hedges_g - 1.96 * np.sqrt(gdata.hedges_g_var)
        ci_hi = gdata.hedges_g + 1.96 * np.sqrt(gdata.hedges_g_var)

        weights = 1.0 / gdata.hedges_g_var.clip(lower=1e-6)
        marker_sizes = 30 + 150 * (weights / weights.max())

        for i, (y, g, lo, hi, ms) in enumerate(
                zip(y_positions, gdata.hedges_g, ci_lo, ci_hi, marker_sizes)):
            ax.plot([lo, hi], [y, y], color="#4a4a4a", linewidth=0.8, zorder=1)
            ax.scatter(g, y, s=ms, color="#2166ac", edgecolors="white",
                       linewidth=0.5, zorder=2)

        y_combined = len(gdata) + 0.8
        meta_lo, meta_hi, meta_g = gmeta.ci_lo, gmeta.ci_hi, gmeta.hedges_g
        diamond_x = [meta_lo, meta_g, meta_hi, meta_g]
        diamond_y = [y_combined, y_combined + 0.35, y_combined, y_combined - 0.35]
        ax.fill(diamond_x, diamond_y, color="#b2182b", alpha=0.85, zorder=3)

        ax.axvline(0, color="grey", linestyle="--", linewidth=0.6, alpha=0.6, zorder=0)

        labels = list(gdata.dataset)
        yticks = list(y_positions) + [y_combined]
        yticklabels = labels + ["Combined"]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize=7.5)
        ax.set_xlabel("Hedges' g")
        ax.set_title(f"{gene_name}", fontweight="bold", fontsize=11)

        txt = f"g = {meta_g:.2f}  [{meta_lo:.2f}, {meta_hi:.2f}]\n"
        txt += f"I² = {gmeta.I2:.0%}   τ² = {gmeta.tau2:.2f}\n"
        txt += f"q = {gmeta.q_value:.1e}   k = {int(gmeta.datasets)}"
        ax.text(0.97, 0.02, txt, transform=ax.transAxes, fontsize=7,
                verticalalignment="bottom", horizontalalignment="right",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#f7f7f7",
                          edgecolor="#cccccc", alpha=0.9))

        for i, (y, row) in enumerate(zip(y_positions, gdata.itertuples())):
            ax.text(1.0, y, f"  {int(row.n_uc)}/{int(row.n_ctrl)}",
                    transform=ax.get_yaxis_transform(), fontsize=6.5,
                    verticalalignment="center", color="#666666")

        ax.set_ylim(-0.8, y_combined + 1.0)

    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def generate_figure_s2(data_dir, fig_dir):
    print("Generating Figure S2: Forest plots...")
    per_gene = load_csv(f"{data_dir}/uc_vs_hc/per_dataset_gene_stats_significant.csv")
    meta_gene = load_csv(f"{data_dir}/uc_vs_hc/gene_meta_significant.csv")

    # Panel A: metabolic/mitochondrial
    forest_plot_multi(
        ["PPARGC1A", "OPA1", "ACSF2"],
        per_gene, meta_gene,
        f"{fig_dir}/figure_s2a.pdf"
    )

    # Panel B: top upregulated + diverse
    forest_plot_multi(
        ["LPCAT1", "DUOX2", "SLC6A14"],
        per_gene, meta_gene,
        f"{fig_dir}/figure_s2b.pdf"
    )

    # Panel C: immune + barrier
    forest_plot_multi(
        ["CXCL1", "AQP8"],
        per_gene, meta_gene,
        f"{fig_dir}/figure_s2c.pdf"
    )

    print(f"  → {fig_dir}/figure_s2[a-c].pdf")


# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE S3: STOUFFER SENSITIVITY
# ═══════════════════════════════════════════════════════════════════════════════

def generate_figure_s3(data_dir, fig_dir):
    print("Generating Figure S3: Stouffer sensitivity...")
    from statsmodels.stats.multitest import multipletests

    per_gsea = load_csv(f"{data_dir}/uc_vs_hc/per_dataset_gsea.csv")
    pw_all = load_csv(f"{data_dir}/uc_vs_hc/pathway_meta_all.csv")
    per_gene = load_csv(f"{data_dir}/uc_vs_hc/per_dataset_gene_stats_significant.csv")

    # Get sample sizes per dataset
    ds_n = {}
    for ds in per_gene.dataset.unique():
        sub = per_gene[per_gene.dataset == ds].iloc[0]
        ds_n[ds] = sub.n_uc + sub.n_ctrl

    # Stouffer's method per pathway
    pathway_terms = per_gsea.Term.unique()
    stouffer_results = []

    for term in pathway_terms:
        pdata = per_gsea[per_gsea.Term == term].copy()
        if len(pdata) < 2:
            continue
        n_pos = (pdata.NES > 0).sum()
        direction = "up" if n_pos >= len(pdata) / 2 else "down"

        z_scores, weights = [], []
        for _, row in pdata.iterrows():
            p = max(min(row["NOM p-val"], 1 - 1e-10), 1e-10)
            nes = row.NES
            if direction == "up":
                p_one = p / 2 if nes > 0 else 1 - p / 2
            else:
                p_one = p / 2 if nes < 0 else 1 - p / 2
            z = stats.norm.ppf(1 - p_one)
            w = np.sqrt(ds_n.get(row.dataset, 50))
            z_scores.append(z)
            weights.append(w)

        z_scores, weights = np.array(z_scores), np.array(weights)
        z_combined = np.sum(weights * z_scores) / np.sqrt(np.sum(weights**2))
        p_combined = 2 * (1 - stats.norm.cdf(abs(z_combined)))
        # Sign the z-score so downregulated pathways have negative z
        z_signed = z_combined if direction == "up" else -z_combined

        stouffer_results.append({
            "term": term,
            "z_stouffer": z_signed,
            "p_stouffer": p_combined,
            "direction": direction,
            "n_datasets": len(pdata),
        })

    stouffer_df = pd.DataFrame(stouffer_results)
    _, stouffer_df["q_stouffer"], _, _ = multipletests(
        stouffer_df.p_stouffer, method="fdr_bh"
    )

    merged = stouffer_df.merge(
        pw_all[["term", "combined_nes", "q_value"]],
        on="term", how="inner"
    )
    merged["sig_reml"] = merged.q_value < 0.05
    merged["sig_stouffer"] = merged.q_stouffer < 0.05

    reml_sig = merged[merged.sig_reml]
    concordance = reml_sig.sig_stouffer.mean()
    r_z, _ = stats.pearsonr(merged.combined_nes, merged.z_stouffer)

    print(f"  Stouffer concordance: {concordance:.1%}, r = {r_z:.2f}")

    # Plot
    fig, ax = plt.subplots(figsize=(5.5, 5))
    nonsig = merged[~merged.sig_reml & ~merged.sig_stouffer]
    ax.scatter(nonsig.combined_nes, nonsig.z_stouffer, s=6, alpha=0.25,
               color="#cccccc", edgecolors="none", rasterized=True,
               label="Neither significant")

    both_df = merged[merged.sig_reml & merged.sig_stouffer]
    ax.scatter(both_df.combined_nes, both_df.z_stouffer, s=10, alpha=0.5,
               color="#2166ac", edgecolors="none", rasterized=True,
               label="Both significant")

    reml_only_df = merged[merged.sig_reml & ~merged.sig_stouffer]
    if len(reml_only_df) > 0:
        ax.scatter(reml_only_df.combined_nes, reml_only_df.z_stouffer, s=15,
                   alpha=0.7, color="#d6604d", edgecolors="none", marker="s",
                   label=f"REML only (n={len(reml_only_df)})")

    stouffer_only_df = merged[~merged.sig_reml & merged.sig_stouffer]
    if len(stouffer_only_df) > 0:
        ax.scatter(stouffer_only_df.combined_nes, stouffer_only_df.z_stouffer,
                   s=15, alpha=0.7, color="#4dac26", edgecolors="none",
                   marker="^", label=f"Stouffer only (n={len(stouffer_only_df)})")

    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.set_xlabel("REML combined NES")
    ax.set_ylabel("Stouffer's weighted z")
    ax.set_title(f"Pathway meta-analysis: REML vs Stouffer\n"
                 f"(Pearson r = {r_z:.2f})", fontsize=10)
    ax.legend(fontsize=7, loc="upper left", frameon=True, fancybox=True)

    plt.tight_layout()
    out = f"{fig_dir}/figure_s3.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  → {out}")

    return concordance, r_z


# ═══════════════════════════════════════════════════════════════════════════════
#  CONSTITUTIVE ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

def compute_constitutive_analysis(data_dir):
    """Compute constitutive vs inflammation-dependent gene classification."""
    print("Computing constitutive gene analysis...")

    infl = load_csv(f"{data_dir}/uc_vs_hc/gene_meta_all.csv")
    uninfl = load_csv(f"{data_dir}/uninf_vs_hc/gene_meta_all.csv")
    infl_sig = load_csv(f"{data_dir}/uc_vs_hc/gene_meta_significant.csv")
    uninfl_sig = load_csv(f"{data_dir}/uninf_vs_hc/gene_meta_significant.csv")

    # Merge ALL shared genes
    from statsmodels.stats.multitest import multipletests
    m = infl[["gene", "hedges_g", "se", "q_value", "direction"]].merge(
        uninfl[["gene", "hedges_g", "se", "q_value", "direction"]],
        on="gene", suffixes=("_infl", "_uninfl"), how="inner"
    )

    # Wald test on ALL shared genes (BH correction across full set)
    m["delta_g"] = m.hedges_g_infl - m.hedges_g_uninfl
    m["delta_se"] = np.sqrt(m.se_infl**2 + m.se_uninfl**2)
    m["delta_z"] = m.delta_g / m.delta_se
    m["delta_p"] = 2 * (1 - stats.norm.cdf(np.abs(m.delta_z)))
    _, m["delta_q"], _, _ = multipletests(m.delta_p, method="fdr_bh")

    # Now filter to genes significant in both analyses, same direction
    sig_infl_genes = set(infl_sig.gene)
    sig_uninfl_genes = set(uninfl_sig.gene)
    both_sig = sig_infl_genes & sig_uninfl_genes
    same_dir = m[m.gene.isin(both_sig)].copy()
    same_dir = same_dir[same_dir.direction_infl == same_dir.direction_uninfl]

    # Classify (delta_q already from genome-wide BH correction)
    same_dir["ratio"] = same_dir.hedges_g_uninfl / same_dir.hedges_g_infl

    constitutive = same_dir[
        (same_dir.delta_q >= 0.05) &
        (same_dir.ratio >= 0.70) & (same_dir.ratio <= 1.30)
    ].copy()

    amplified = same_dir[same_dir.delta_q < 0.05].copy()

    print(f"  Genes significant in both: {len(same_dir)}")
    print(f"  Constitutive (ratio 0.7-1.3): {len(constitutive)}")
    print(f"  Inflammation-amplified: {len(amplified)}")

    return constitutive, amplified, same_dir


# ═══════════════════════════════════════════════════════════════════════════════
#  SCATTER PLOT (used as main figure, generated here for completeness)
# ═══════════════════════════════════════════════════════════════════════════════

def generate_scatter_comparison(same_dir, fig_dir):
    """Generate the constitutive scatter comparison plot."""
    print("Generating scatter comparison plot...")

    const = same_dir[same_dir.delta_q >= 0.05].copy()
    const["ratio"] = const.hedges_g_uninfl / const.hedges_g_infl
    ratio_dev = (const["ratio"] - 1).abs()

    fig, ax = plt.subplots(figsize=(7, 7))
    sc = ax.scatter(
        const.hedges_g_infl, const.hedges_g_uninfl,
        c=ratio_dev, cmap="RdYlBu_r", s=55, alpha=0.85, zorder=3,
        edgecolors="k", linewidths=0.4, vmin=0, vmax=0.8,
    )

    # Labels for near-diagonal genes
    to_label = const[const["ratio"].between(0.7, 1.3)]
    try:
        from adjustText import adjust_text
        texts = []
        for _, r in to_label.iterrows():
            texts.append(ax.text(
                r.hedges_g_infl, r.hedges_g_uninfl, r.gene,
                fontsize=8, color="k", fontweight="medium", zorder=5,
            ))
        adjust_text(
            texts, ax=ax,
            x=const.hedges_g_infl.values,
            y=const.hedges_g_uninfl.values,
            arrowprops=dict(arrowstyle="-", color="grey", lw=0.5, alpha=0.8),
            expand=(5, 5), force_text=(1, 1),
        )
    except ImportError:
        pass

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
    cb.set_label("$|g_{\\mathrm{uninfl}} \\,/\\, g_{\\mathrm{infl}} - 1|$", fontsize=11)

    out = f"{fig_dir}/scatter_comparison.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  → {out}")


# ═══════════════════════════════════════════════════════════════════════════════
#  PATHWAY HEATMAP (butterfly layout)
# ═══════════════════════════════════════════════════════════════════════════════

def generate_pathway_heatmap(data_dir, fig_dir):
    """Generate three-way pathway comparison heatmap."""
    print("Generating pathway heatmap (butterfly layout)...")

    SIG_THR = 0.05

    uc_cl = load_csv(f"{data_dir}/uc_vs_hc/hierarchical_pathway_clusters.csv")
    un_cl = load_csv(f"{data_dir}/uninf_vs_hc/hierarchical_pathway_clusters.csv")
    cd_cl = load_csv(f"{data_dir}/cd_vs_hc/hierarchical_pathway_clusters.csv")

    all_reps = set()
    for cl in [uc_cl, un_cl, cd_cl]:
        all_reps.update(cl["Representative"].tolist())

    uc_pw = load_csv(f"{data_dir}/uc_vs_hc/pathway_meta_all.csv")
    un_pw = load_csv(f"{data_dir}/uninf_vs_hc/pathway_meta_all.csv")
    cd_pw = load_csv(f"{data_dir}/cd_vs_hc/pathway_meta_all.csv")

    uc_nes = dict(zip(uc_pw["term"], uc_pw["combined_nes"]))
    un_nes = dict(zip(un_pw["term"], un_pw["combined_nes"]))
    cd_nes = dict(zip(cd_pw["term"], cd_pw["combined_nes"]))
    uc_qval = dict(zip(uc_pw["term"], uc_pw["q_value"]))
    un_qval = dict(zip(un_pw["term"], un_pw["q_value"]))
    cd_qval = dict(zip(cd_pw["term"], cd_pw["q_value"]))

    data = []
    for term in all_reps:
        name = clean_pathway_name(term)
        if len(name) > 55:
            name = name[:52] + "..."
        uc_val = uc_nes.get(term, np.nan)
        un_val = un_nes.get(term, np.nan)
        cd_val = cd_nes.get(term, np.nan)
        sort_nes = uc_val if not np.isnan(uc_val) else (
            cd_val if not np.isnan(cd_val) else un_val)

        uc_sig = uc_qval.get(term, np.nan)
        un_sig = un_qval.get(term, np.nan)
        cd_sig = cd_qval.get(term, np.nan)
        uc_sig = (uc_sig < SIG_THR) if not np.isnan(uc_sig) else np.nan
        un_sig = (un_sig < SIG_THR) if not np.isnan(un_sig) else np.nan
        cd_sig = (cd_sig < SIG_THR) if not np.isnan(cd_sig) else np.nan

        data.append({"name": name,
                      "uc": uc_val, "un": un_val, "cd": cd_val,
                      "uc_sig": uc_sig, "un_sig": un_sig, "cd_sig": cd_sig,
                      "sort_nes": sort_nes})

    df = pd.DataFrame(data)
    df_up = df[df["sort_nes"] > 0].sort_values("sort_nes", ascending=False).reset_index(drop=True)
    df_down = df[df["sort_nes"] <= 0].sort_values("sort_nes", ascending=True).reset_index(drop=True)

    mat_up, sig_up = df_up[["uc","un","cd"]].values, df_up[["uc_sig","un_sig","cd_sig"]].values
    mat_down, sig_down = df_down[["uc","un","cd"]].values, df_down[["uc_sig","un_sig","cd_sig"]].values
    n_up, n_down = len(df_up), len(df_down)
    n_max = max(n_up, n_down)

    def pad(mat, sig, names, n_actual, n_target):
        if n_actual < n_target:
            pad_n = n_target - n_actual
            mat = np.vstack([mat, np.full((pad_n, 3), np.nan)])
            sig = np.vstack([sig, np.full((pad_n, 3), np.nan)])
            names = names + [""] * pad_n
        return mat, sig, names

    mat_up, sig_up, names_up = pad(mat_up, sig_up, df_up["name"].tolist(), n_up, n_max)
    mat_down, sig_down, names_down = pad(mat_down, sig_down, df_down["name"].tolist(), n_down, n_max)

    all_vals = np.concatenate([mat_up.ravel(), mat_down.ravel()])
    vmax = np.nanmax(np.abs(all_vals[np.isfinite(all_vals)]))
    cmap = plt.cm.RdBu_r
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    col_labels = ["UC\ninfl.", "UC\nuninfl.", "CD\ninfl."]
    cell_w, cell_h, gap = 1.0, 1.0, 0.4
    x_up = [0, cell_w, 2 * cell_w]
    x_down = [3 * cell_w + gap, 4 * cell_w + gap, 5 * cell_w + gap]

    fig_w, fig_h = 13.5, max(5.5, 0.36 * n_max + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    def draw_cell(x, y, val, significant):
        if np.isnan(val):
            rect = plt.Rectangle((x, y), cell_w, cell_h, facecolor="#ededed",
                                  edgecolor="white", linewidth=0.8)
            ax.add_patch(rect)
            return
        base_color = cmap(norm(val))
        if significant is True or significant == 1.0:
            face, hatch = base_color, None
        else:
            face, hatch = desaturate(base_color, factor=0.1), "////"

        rect = plt.Rectangle((x, y), cell_w, cell_h, facecolor=face,
                              edgecolor="white", linewidth=0.8, hatch=hatch)
        if hatch:
            rect.set_edgecolor("white")
            ax.add_patch(rect)
            hatch_rect = plt.Rectangle((x, y), cell_w, cell_h, facecolor="none",
                                        edgecolor="#888888", linewidth=0, hatch=hatch, alpha=0.45)
            ax.add_patch(hatch_rect)
        else:
            ax.add_patch(rect)

        text = f"{val:+.1f}"
        brightness = 0.299 * face[0] + 0.587 * face[1] + 0.114 * face[2]
        txt_color = "white" if brightness < 0.55 else "#222222"
        ax.text(x + cell_w / 2, y + cell_h / 2, text, ha="center", va="center",
                fontsize=6.5, color=txt_color, fontweight="medium")

    for i in range(n_max):
        y = (n_max - 1 - i) * cell_h
        for j in range(3):
            if i < n_up:
                draw_cell(x_up[j], y, mat_up[i, j], sig_up[i, j])
            else:
                ax.add_patch(plt.Rectangle((x_up[j], y), cell_w, cell_h,
                             facecolor="#fafafa", edgecolor="#fafafa", linewidth=0))
            if i < n_down:
                draw_cell(x_down[j], y, mat_down[i, j], sig_down[i, j])
            else:
                ax.add_patch(plt.Rectangle((x_down[j], y), cell_w, cell_h,
                             facecolor="#fafafa", edgecolor="#fafafa", linewidth=0))

    for i, name in enumerate(names_up):
        if name:
            ax.text(x_up[0] - 0.15, (n_max - 1 - i) * cell_h + cell_h / 2,
                    name, ha="right", va="center", fontsize=6.5)
    for i, name in enumerate(names_down):
        if name:
            ax.text(x_down[2] + cell_w + 0.15, (n_max - 1 - i) * cell_h + cell_h / 2,
                    name, ha="left", va="center", fontsize=6.5)

    header_y = n_max * cell_h + 0.15
    for j, label in enumerate(col_labels):
        ax.text(x_up[j] + cell_w / 2, header_y, label, ha="center", va="bottom",
                fontsize=7, fontweight="bold")
        ax.text(x_down[j] + cell_w / 2, header_y, label, ha="center", va="bottom",
                fontsize=7, fontweight="bold")

    title_y = n_max * cell_h + 1.5
    ax.text((x_up[0] + x_up[2] + cell_w) / 2, title_y, "Upregulated",
            ha="center", va="bottom", fontsize=9, fontweight="bold", color="#b22222")
    ax.text((x_down[0] + x_down[2] + cell_w) / 2, title_y, "Downregulated",
            ha="center", va="bottom", fontsize=9, fontweight="bold", color="#1a5276")

    cbar_ax = fig.add_axes([0.4123, 0.1, 0.2, 0.018])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Normalized Enrichment Score (NES)", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)

    ax.set_xlim(x_up[0] - 0.05, x_down[2] + cell_w + 0.05)
    ax.set_ylim(-1.0, n_max * cell_h + 2.2)
    ax.set_aspect("equal")
    ax.axis("off")

    out = f"{fig_dir}/pathway_heatmap.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  → {out} ({n_up} up + {n_down} down = {n_up + n_down} pathways)")


# ═══════════════════════════════════════════════════════════════════════════════
#  UC vs CD CONCORDANCE PLOT
# ═══════════════════════════════════════════════════════════════════════════════

def generate_uc_cd_concordance(data_dir, fig_dir):
    """Generate concordance scatter plot for direct vs indirect UC-CD comparison."""
    print("Generating UC vs CD concordance plot...")

    uc = load_csv(f"{data_dir}/uc_vs_hc/gene_meta_all.csv")[["gene", "hedges_g"]].rename(
        columns={"hedges_g": "g_uc"})
    cd = load_csv(f"{data_dir}/cd_vs_hc/gene_meta_all.csv")[["gene", "hedges_g"]].rename(
        columns={"hedges_g": "g_cd"})
    direct = load_csv(f"{data_dir}/uc_vs_cd/gene_meta_all.csv")[["gene", "hedges_g"]].rename(
        columns={"hedges_g": "g_direct"})

    merged = uc.merge(cd, on="gene").merge(direct, on="gene")
    merged["delta_indirect"] = merged["g_uc"] - merged["g_cd"]
    r = merged["delta_indirect"].corr(merged["g_direct"])
    print(f"  Pearson r = {r:.3f}, n = {len(merged)}")

    label_genes = {
        "PPARGC1A": "metabolic", "PPARGC1B": "metabolic",
        "CPT1A": "metabolic", "CPT2": "metabolic",
        "HADHA": "metabolic", "HADHB": "metabolic",
        "ACADS": "metabolic", "ACO2": "metabolic",
        "CS": "metabolic", "SDHB": "metabolic",
        "CXCL1": "immune", "LCN2": "immune",
        "S100A8": "immune", "IL1B": "immune", "TNF": "immune",
        "TTC9": "other", "BMAL2": "other", "FGFR2": "other",
    }
    colors_map = {"metabolic": "#2166ac", "immune": "#b2182b", "other": "#4d4d4d"}

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(merged.delta_indirect, merged.g_direct, s=3, alpha=0.12,
               color="#999999", rasterized=True, zorder=1)

    for gene, cat in label_genes.items():
        row = merged[merged.gene == gene]
        if row.empty:
            continue
        x, y = row.delta_indirect.values[0], row.g_direct.values[0]
        ax.scatter(x, y, s=30, color=colors_map[cat], edgecolors="white",
                   linewidths=0.3, zorder=3)

    try:
        from adjustText import adjust_text
        texts = []
        for gene, cat in label_genes.items():
            row = merged[merged.gene == gene]
            if row.empty:
                continue
            x, y = row.delta_indirect.values[0], row.g_direct.values[0]
            texts.append(ax.text(x, y, gene, fontsize=5.5, color=colors_map[cat],
                                 fontweight="bold", zorder=4))
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.4),
                    expand=(1.4, 1.6), force_text=(0.4, 0.5))
    except ImportError:
        pass

    lim = max(abs(merged.delta_indirect).max(), abs(merged.g_direct).max()) * 1.05
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.6, alpha=0.4, zorder=0)
    ax.axhline(0, color="#cccccc", lw=0.5, zorder=0)
    ax.axvline(0, color="#cccccc", lw=0.5, zorder=0)
    ax.set_xlabel(r"Indirect ($g_{\mathrm{UC\,vs\,HC}}$ − $g_{\mathrm{CD\,vs\,HC}}$)", fontsize=9)
    ax.set_ylabel(r"Direct ($g_{\mathrm{UC\,vs\,CD}}$)", fontsize=9)
    ax.set_title(f"Gene-level concordance (Pearson r = {r:.2f})", fontsize=10)
    ax.tick_params(labelsize=8)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#2166ac", markersize=6, label="Metabolic"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#b2182b", markersize=6, label="Immune"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#4d4d4d", markersize=6, label="Other"),
    ]
    ax.legend(handles=legend_elements, fontsize=7, loc="upper left", framealpha=0.8)

    fig.tight_layout()
    out = f"{fig_dir}/uc_cd_concordance.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  → {out}")


# ═══════════════════════════════════════════════════════════════════════════════
#  SUPPLEMENTARY TABLES (CSV)
# ═══════════════════════════════════════════════════════════════════════════════

def generate_csv_tables(data_dir, output_dir):
    """Generate supplementary CSV tables for gene and pathway results."""
    print("Generating supplementary CSV tables...")

    cols_gene = ["gene", "datasets", "hedges_g", "se", "ci_lo", "ci_hi",
                 "p_value", "q_value", "tau2", "I2", "direction",
                 "direction_ratio"]

    for label, subdir, fname in [
        ("S1", "uc_vs_hc", "gene_meta_significant.csv"),
        ("S2", "uninf_vs_hc", "gene_meta_significant.csv"),
        ("S3", "cd_vs_hc", "gene_meta_significant.csv"),
        ("S4", "uc_vs_cd", "gene_meta_significant.csv"),
    ]:
        df = load_csv(f"{data_dir}/{subdir}/{fname}")
        use_cols = [c for c in cols_gene if c in df.columns]
        out = f"{output_dir}/Table_{label}_genes_{subdir}.csv"
        df[use_cols].to_csv(out, index=False, float_format="%.6g")
        print(f"  → {out} ({len(df)} genes)")

    # Pathway tables S5--S8: separate significant pathway CSV per comparison
    cols_pw = ["term", "datasets", "combined_nes", "se", "ci_lo", "ci_hi",
               "p_value", "q_value", "tau2", "I2", "direction",
               "direction_ratio", "lead_genes_core"]

    for label, subdir, desc in [
        ("S5", "uc_vs_hc", "UC inflamed vs healthy controls"),
        ("S6", "uninf_vs_hc", "UC uninflamed vs healthy controls"),
        ("S7", "cd_vs_hc", "CD inflamed vs healthy controls"),
        ("S8", "uc_vs_cd", "UC vs CD direct comparison"),
    ]:
        df = load_csv(f"{data_dir}/{subdir}/pathway_meta_significant.csv")
        use_cols = [c for c in cols_pw if c in df.columns]
        out = f"{output_dir}/Table_{label}_pathways_{subdir}.csv"
        df[use_cols].to_csv(out, index=False, float_format="%.6g")
        print(f"  → {out} ({len(df)} significant pathways — {desc})")


# ═══════════════════════════════════════════════════════════════════════════════
#  LATEX SUPPLEMENTARY DOCUMENT
# ═══════════════════════════════════════════════════════════════════════════════

def cluster_table_latex(data_dir, subdir, comparison_label):
    """Generate LaTeX for a pathway cluster table."""
    cl = load_csv(f"{data_dir}/{subdir}/hierarchical_pathway_clusters.csv")
    up_cl = cl[cl.Mean_NES > 0].sort_values("Mean_NES", ascending=False)
    dn_cl = cl[cl.Mean_NES < 0].sort_values("Mean_NES")

    lines = []
    lines.append(r"\begin{tabular}{p{7cm}rr}")
    lines.append(r"\toprule")
    lines.append(r"Cluster representative & Pathways & Mean NES \\")
    lines.append(r"\midrule")

    for _, row in up_cl.iterrows():
        name = escape_latex(clean_pathway_name(row["Representative"]))
        lines.append(f"  {name} & {row['Size']:.0f} & ${row['Mean_NES']:+.2f}$ \\\\")

    if len(dn_cl) > 0:
        lines.append(r"\midrule")
        for _, row in dn_cl.iterrows():
            name = escape_latex(clean_pathway_name(row["Representative"]))
            lines.append(f"  {name} & {row['Size']:.0f} & ${row['Mean_NES']:+.2f}$ \\\\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    return "\n".join(lines)


def constitutive_table_latex(df, direction):
    """Generate LaTeX for constitutive or amplified gene table."""
    sub = df.copy()
    if direction == "constitutive_up":
        sub = sub[sub.hedges_g_uninfl > 0].sort_values("hedges_g_uninfl", ascending=False)
    elif direction == "constitutive_down":
        sub = sub[sub.hedges_g_uninfl < 0].sort_values("hedges_g_uninfl")
    else:
        sub = sub.sort_values("delta_g", key=abs, ascending=False)

    lines = []
    lines.append(r"\begin{tabular}{lrrrrr}")
    lines.append(r"\toprule")
    lines.append(r"Gene & $g_{\mathrm{infl}}$ & $g_{\mathrm{uninfl}}$ & "
                 r"Ratio & $\Delta g$ & $q_{\Delta}$ \\")
    lines.append(r"\midrule")

    for _, r in sub.iterrows():
        gene = escape_latex(r.gene)
        ratio = r.hedges_g_uninfl / r.hedges_g_infl if r.hedges_g_infl != 0 else np.nan
        lines.append(
            f"  {gene} & ${r.hedges_g_infl:+.2f}$ & ${r.hedges_g_uninfl:+.2f}$ & "
            f"${ratio:.2f}$ & ${r.delta_g:+.2f}$ & ${r.delta_q:.3f}$ \\\\"
        )

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    return "\n".join(lines)


def build_supplementary_latex(fig_dir, output_dir, data_dir,
                              constitutive_df, amplified_df):
    """Build the complete supplementary LaTeX document."""
    print("Building supplementary LaTeX document...")

    # Pathway cluster tables for all 4 comparisons
    uc_cluster_table = cluster_table_latex(data_dir, "uc_vs_hc", "UC inflamed vs HC")
    uninfl_cluster_table = cluster_table_latex(data_dir, "uninf_vs_hc", "Uninflamed UC vs HC")
    cd_cluster_table = cluster_table_latex(data_dir, "cd_vs_hc", "CD vs HC")
    uccd_cluster_table = cluster_table_latex(data_dir, "uc_vs_cd", "UC vs CD")

    const_up = constitutive_df[constitutive_df.hedges_g_uninfl > 0]
    const_down = constitutive_df[constitutive_df.hedges_g_uninfl < 0]

    const_up_table = constitutive_table_latex(const_up, "constitutive_up")
    const_down_table = constitutive_table_latex(const_down, "constitutive_down")
    amplified_table = constitutive_table_latex(amplified_df, "amplified")

    # Count significant pathways per comparison for the text
    pw_counts = {}
    for subdir, key in [("uc_vs_hc", "uc"), ("uninf_vs_hc", "uninfl"),
                         ("cd_vs_hc", "cd"), ("uc_vs_cd", "uccd")]:
        try:
            df = load_csv(f"{data_dir}/{subdir}/pathway_meta_significant.csv")
            pw_counts[key] = len(df)
        except Exception:
            pw_counts[key] = "?"

    doc = r"""\documentclass[11pt,a4paper]{article}
\usepackage[margin=2cm]{geometry}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{hyperref}
\usepackage{xcolor}
\usepackage{float}

\hypersetup{colorlinks=true, linkcolor=blue!60!black, urlcolor=blue!60!black}

\renewcommand{\arraystretch}{1.15}
\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\thetable}{S\arabic{table}}

\title{\textbf{Supplementary Materials}\\[0.3em]
\large Constitutive and inflammation-dependent transcriptomic signatures\\
in ulcerative colitis: a multi-cohort meta-analysis of mucosal biopsies}
\author{Piernik et al.}
\date{}

\begin{document}
\maketitle
\tableofcontents
\clearpage

% ══════════════════════════════════════════════════════════════════
\section{Supplementary Figures}
% ══════════════════════════════════════════════════════════════════

% ── Figure S1 ──
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{""" + os.path.abspath(f"{fig_dir}/figure_s1.pdf") + r"""}
\caption{\textbf{Heterogeneity diagnostics for the inflamed UC vs control gene-level meta-analysis.}
\textbf{(A)}~Distribution of $I^2$ across 7{,}727 significantly dysregulated genes.
The median $I^2$ of 0.76 reflects variation in effect-size magnitude across the 12 datasets
and 9 microarray platforms, rather than inconsistency in direction (controlled by the
$\geq 80\%$ direction consistency filter).
\textbf{(B)}~Relationship between absolute effect size ($|g|$) and heterogeneity ($I^2$).
Colour encodes statistical significance ($-\log_{10} q$).
Genes with large effects tend to have high $I^2$ (Spearman $\rho = 0.74$),
while genes ranking highest by $q$-value have low $I^2$, reflecting complementary
properties of the two ranking criteria.}
\label{fig:s1}
\end{figure}
\clearpage

% ── Figure S2 ──
\begin{figure}[H]
\centering
\includegraphics[width=0.8\textwidth]{""" + os.path.abspath(f"{fig_dir}/figure_s2a.pdf") + r"""}

\vspace{0.8em}
\includegraphics[width=0.8\textwidth]{""" + os.path.abspath(f"{fig_dir}/figure_s2b.pdf") + r"""}

\vspace{0.8em}
\includegraphics[width=0.8\textwidth]{""" + os.path.abspath(f"{fig_dir}/figure_s2c.pdf") + r"""}
\caption{\textbf{Forest plots for representative genes from the inflamed UC vs control meta-analysis.}
Each panel shows per-dataset effect sizes (Hedges'~$g$, blue circles scaled by inverse variance)
with 95\% confidence intervals, and the random-effects combined estimate (red diamond).
Numbers at right indicate sample sizes ($n_{\mathrm{UC}}/n_{\mathrm{Ctrl}}$).
\textbf{Top row:} PPARGC1A, OPA1, and ACSF2 (downregulated metabolic genes).
\textbf{Middle row:} LPCAT1, DUOX2, and SLC6A14 (top upregulated genes).
\textbf{Bottom row:} CXCL1 (top upregulated immune) and AQP8 (downregulated barrier gene).}
\label{fig:s2}
\end{figure}
\clearpage

% ── Figure S3 ──
\begin{figure}[H]
\centering
\includegraphics[width=0.7\textwidth]{""" + os.path.abspath(f"{fig_dir}/figure_s3.pdf") + r"""}
\caption{\textbf{Sensitivity analysis: REML random-effects vs Stouffer's weighted $z$-method
for pathway-level meta-analysis.}
Each point represents one Reactome pathway.
The $x$-axis shows the REML combined normalized enrichment score (NES); the $y$-axis shows
Stouffer's weighted $z$-statistic (weights = $\sqrt{n}$).
Blue points are significant by both methods ($q < 0.05$); red squares are REML-only;
green triangles are Stouffer-only; grey points are non-significant by both.
91.3\% of REML-significant pathways were confirmed by the Stouffer method,
with 100\% direction concordance and Pearson $r = 0.89$, supporting the validity
of the NES variance approximation used in the REML framework.}
\label{fig:s3}
\end{figure}
\clearpage

% ══════════════════════════════════════════════════════════════════
\section{Supplementary Tables}
% ══════════════════════════════════════════════════════════════════

% ── Tables S1--S4 ──
\subsection*{Tables S1--S4: Gene-level meta-analysis results}

Full gene-level meta-analysis results for all four comparisons are provided as
separate CSV files:

\begin{itemize}
\item \textbf{Table~S1}: UC inflamed vs healthy controls (7{,}727 significant genes)
\item \textbf{Table~S2}: UC uninflamed vs healthy controls (165 significant genes)
\item \textbf{Table~S3}: CD inflamed vs healthy controls (2{,}494 significant genes)
\item \textbf{Table~S4}: UC vs CD direct comparison (1{,}831 significant genes)
\end{itemize}

\noindent Each file contains: gene symbol, number of datasets, Hedges'~$g$, standard error,
95\% confidence interval, $p$-value, $q$-value (Benjamini--Hochberg), $\tau^2$, $I^2$,
direction of effect, and direction consistency ratio.

% ── Tables S5--S8 ──
\subsection*{Tables S5--S8: Pathway-level meta-analysis results}

Significant pathway-level meta-analysis results ($q < 0.05$) for all four comparisons
are provided as separate CSV files:

\begin{itemize}
\item \textbf{Table~S5}: UC inflamed vs healthy controls (""" + str(pw_counts["uc"]) + r""" significant pathways)
\item \textbf{Table~S6}: UC uninflamed vs healthy controls (""" + str(pw_counts["uninfl"]) + r""" significant pathways)
\item \textbf{Table~S7}: CD inflamed vs healthy controls (""" + str(pw_counts["cd"]) + r""" significant pathways)
\item \textbf{Table~S8}: UC vs CD direct comparison (""" + str(pw_counts["uccd"]) + r""" significant pathways)
\end{itemize}

\noindent Each file contains: Reactome pathway term, number of datasets, combined NES, standard error,
95\% CI, $p$-value, $q$-value, $\tau^2$, $I^2$, direction, direction ratio, and core leading-edge genes.

\clearpage

% ── Table S9 ──
\subsection*{Table S9: Pathway clusters --- UC inflamed vs healthy controls}

\setcounter{table}{8}
\begin{table}[H]
\centering
\scriptsize
""" + uc_cluster_table + r"""
\caption{Pathway clusters from the inflamed UC vs healthy controls analysis.
Each row is a community of co-enriched Reactome pathways identified by Louvain clustering
on shared leading-edge genes. Positive NES = upregulated in inflamed UC relative to controls;
negative NES = downregulated.}
\label{tab:s9}
\end{table}
\clearpage

% ── Table S10 ──
\subsection*{Table S10: Pathway clusters --- uninflamed UC vs healthy controls}

\begin{table}[H]
\centering
\scriptsize
""" + uninfl_cluster_table + r"""
\caption{Pathway clusters from the uninflamed UC vs healthy controls analysis.
Each row is a community of co-enriched Reactome pathways identified by Louvain clustering
on shared leading-edge genes. Positive NES = upregulated in uninflamed UC relative to controls;
negative NES = downregulated.}
\label{tab:s10}
\end{table}
\clearpage

% ── Table S11 ──
\subsection*{Table S11: Pathway clusters --- CD inflamed vs healthy controls}

\begin{table}[H]
\centering
\scriptsize
""" + cd_cluster_table + r"""
\caption{Pathway clusters from the CD inflamed vs healthy controls analysis.
Each row is a community of co-enriched Reactome pathways identified by Louvain clustering
on shared leading-edge genes. Positive NES = upregulated in CD relative to controls;
negative NES = downregulated.}
\label{tab:s11}
\end{table}
\clearpage

% ── Table S12 ──
\subsection*{Table S12: Pathway clusters --- UC vs CD direct comparison}

\begin{table}[H]
\centering
\scriptsize
""" + uccd_cluster_table + r"""
\caption{Pathway clusters from the UC vs CD direct comparison analysis.
Each row is a community of co-enriched Reactome pathways identified by Louvain clustering
on shared leading-edge genes. Positive NES = upregulated in UC relative to CD;
negative NES = downregulated in UC relative to CD.}
\label{tab:s12}
\end{table}
\clearpage

% ── Table S13 ──
\subsection*{Table S13: Constitutive genes}

Genes significant in both the inflamed and uninflamed UC vs control meta-analyses,
with no significant difference between conditions (Wald $q_{\Delta} \geq 0.05$)
and effect-size ratio within $0.70$--$1.30$.

\begin{table}[H]
\centering
\scriptsize
\caption{Constitutive genes (""" + str(len(constitutive_df)) + r""" genes: """ + str(len(const_up)) + r""" upregulated, """ + str(len(const_down)) + r""" downregulated).
Ratio $= g_{\mathrm{uninfl}} / g_{\mathrm{infl}}$.}
\label{tab:s13}
""" + const_up_table + r"""

\vspace{0.5em}

""" + const_down_table + r"""
\end{table}
\clearpage

% ── Table S14 ──
\subsection*{Table S14: Inflammation-amplified genes}

Genes significant in both the inflamed and uninflamed UC vs control meta-analyses
whose effect sizes differed significantly between conditions (Wald $q_{\Delta} < 0.05$).

\begin{table}[H]
\centering
\scriptsize
\caption{Inflammation-amplified genes (""" + str(len(amplified_df)) + r""" genes).
$\Delta g = g_{\mathrm{infl}} - g_{\mathrm{uninfl}}$.}
\label{tab:s14}
""" + amplified_table + r"""
\end{table}

\end{document}
"""

    tex_path = f"{output_dir}/supplementary.tex"
    with open(tex_path, "w") as f:
        f.write(doc)
    print(f"  → {tex_path}")
    return tex_path


def compile_latex(tex_path, output_dir):
    """Compile LaTeX to PDF."""
    print("Compiling LaTeX to PDF...")
    for _ in range(2):  # two passes for TOC
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode",
             "-output-directory", output_dir, tex_path],
            capture_output=True, text=True
        )
    pdf_path = tex_path.replace(".tex", ".pdf")
    if os.path.exists(pdf_path):
        print(f"  → {pdf_path}")
    else:
        print(f"  ERROR: PDF compilation failed")
        print(result.stdout[-2000:] if result.stdout else "")
        print(result.stderr[-2000:] if result.stderr else "")


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    args = parse_args()
    data_dir = args.data_dir
    output_dir = args.output_dir
    fig_dir = f"{output_dir}/figures"
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    # ── Figures ──
    generate_figure_s1(data_dir, fig_dir)
    generate_figure_s2(data_dir, fig_dir)
    generate_figure_s3(data_dir, fig_dir)

    # ── Main-text figures (also regenerated for completeness) ──
    constitutive_df, amplified_df, same_dir = compute_constitutive_analysis(data_dir)
    generate_scatter_comparison(same_dir, fig_dir)
    generate_pathway_heatmap(data_dir, fig_dir)
    generate_uc_cd_concordance(data_dir, fig_dir)

    # ── CSV supplementary tables ──
    generate_csv_tables(data_dir, output_dir)

    # ── LaTeX supplementary document ──
    tex_path = build_supplementary_latex(fig_dir, output_dir, data_dir,
                                          constitutive_df, amplified_df)
    compile_latex(tex_path, output_dir)

    print("\n" + "=" * 70)
    print("DONE. Generated files:")
    for f in sorted(Path(output_dir).rglob("*")):
        if f.is_file() and not f.name.endswith((".aux", ".log", ".out", ".toc")):
            print(f"  {f}")


if __name__ == "__main__":
    main()