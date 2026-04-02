#!/usr/bin/env python3
"""
08_fig_uc_cd_concordance.py

Generate figures and tables for the direct UC vs CD comparison.

Inputs (from uc_vs_cd/ analysis directory):
  - uc_vs_cd/dataset_summary.csv
  - uc_vs_cd/gene_meta_all.csv
  - uc_vs_cd/gene_meta_significant.csv
  - uc_vs_cd/pathway_meta_significant_knee_with_clusters.csv
  - uc_vs_cd/hierarchical_pathway_clusters.csv

Also requires indirect comparison results:
  - uc_vs_hc/gene_meta_all.csv
  - cd_vs_hc/gene_meta_all.csv

Outputs:
  - figures/uc_cd_concordance.pdf  — scatter of indirect Δg vs direct g
  - tables/uc_vs_cd_clusters.tex   — pathway cluster table
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────
BASE = Path("results/paper/")
DIRECT = BASE / "uc_vs_cd"
UC_HC = BASE / "uc_vs_hc"
CD_HC = BASE / "cd_vs_hc"
FIGDIR = Path("figures")
TABDIR = BASE / "tables"
FIGDIR.mkdir(exist_ok=True)
TABDIR.mkdir(exist_ok=True)


# ══════════════════════════════════════════════════════════════════
#  FIGURE: Concordance between indirect (g_UC − g_CD) and direct g
# ══════════════════════════════════════════════════════════════════

uc = pd.read_csv(UC_HC / "gene_meta_all.csv")[["gene", "hedges_g"]].rename(
    columns={"hedges_g": "g_uc"}
)
cd = pd.read_csv(CD_HC / "gene_meta_all.csv")[["gene", "hedges_g"]].rename(
    columns={"hedges_g": "g_cd"}
)
direct = pd.read_csv(DIRECT / "gene_meta_all.csv")[["gene", "hedges_g"]].rename(
    columns={"hedges_g": "g_direct"}
)

merged = uc.merge(cd, on="gene").merge(direct, on="gene")
merged["delta_indirect"] = merged["g_uc"] - merged["g_cd"]

r = merged["delta_indirect"].corr(merged["g_direct"])
print(f"Pearson r (indirect Δg vs direct g): {r:.3f}")
print(f"Genes in concordance plot: {len(merged)}")

# Key genes to label
label_genes = {
    # Metabolic
    "PPARGC1A": "metabolic", "PPARGC1B": "metabolic",
    "CPT1A": "metabolic", "CPT2": "metabolic",
    "HADHA": "metabolic", "HADHB": "metabolic",
    "ACADS": "metabolic", "ACO2": "metabolic",
    "CS": "metabolic", "SDHB": "metabolic",
    "NDUFS1": "metabolic", "COX5A": "metabolic",
    "ATP5F1B": "metabolic", "ETFDH": "metabolic",
    # Immune
    "CXCL1": "immune", "LCN2": "immune",
    "S100A8": "immune", "IL1B": "immune",
    "TNF": "immune",
    # Translation / top hits
    "TTC9": "other", "BMAL2": "other",
    "RMDN2": "other", "FGFR2": "other",
}

colors_map = {"metabolic": "#2471a3", "immune": "#c0392b", "other": "#4d4d4d"}

fig, ax = plt.subplots(figsize=(5.5, 5.5))

# Background scatter
ax.scatter(
    merged["delta_indirect"], merged["g_direct"],
    s=3, alpha=0.12, color="#999999", rasterized=True, zorder=1
)

# Highlighted genes
for gene, cat in label_genes.items():
    row = merged[merged.gene == gene]
    if row.empty:
        continue
    x, y = row.delta_indirect.values[0], row.g_direct.values[0]
    ax.scatter(x, y, s=30, color=colors_map[cat], edgecolors="white",
               linewidths=0.3, zorder=3)

# Labels with manual nudge
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
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.4),
                expand=(1.7, 1.7), force_text=(2, 2))
except ImportError:
    # Fallback: no adjustText
    for gene, cat in label_genes.items():
        row = merged[merged.gene == gene]
        if row.empty:
            continue
        x, y = row.delta_indirect.values[0], row.g_direct.values[0]
        ax.annotate(gene, (x, y), fontsize=5, color=colors_map[cat],
                    textcoords="offset points", xytext=(4, 4))

# Diagonal
lim = max(abs(merged["delta_indirect"]).max(), abs(merged["g_direct"]).max())
ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.6, alpha=0.4, zorder=0)
ax.axhline(0, color="#bdc3c7", lw=0.5, zorder=0)
ax.axvline(0, color="#bdc3c7", lw=0.5, zorder=0)

ax.set_xlabel(r"Indirect difference ($g_{\mathrm{UC\,vs\,HC}}$ − $g_{\mathrm{CD\,vs\,HC}}$)", fontsize=9)
ax.set_ylabel(r"Direct comparison ($g_{\mathrm{UC\,vs\,CD}}$)", fontsize=9)
ax.set_title(f"Gene-level concordance (Pearson r = {r:.2f})", fontsize=10)
ax.tick_params(labelsize=8)

plot_lim = 1.8
ax.set_xlim(-plot_lim, plot_lim)
ax.set_ylim(-plot_lim, plot_lim)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#2471a3", markersize=6, label="Metabolic"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#c0392b", markersize=6, label="Immune"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#4d4d4d", markersize=6, label="Other"),
]
ax.legend(handles=legend_elements, fontsize=7, loc="upper left", framealpha=0.8)

fig.tight_layout()
fig.savefig(FIGDIR / "uc_cd_concordance.pdf", dpi=300, bbox_inches="tight")
print(f"Saved {FIGDIR / 'uc_cd_concordance.pdf'}")
plt.close()


# ══════════════════════════════════════════════════════════════════
#  TABLE: Pathway clusters from direct UC vs CD comparison
# ══════════════════════════════════════════════════════════════════

cl = pd.read_csv(DIRECT / "hierarchical_pathway_clusters.csv")
pk = pd.read_csv(DIRECT / "pathway_meta_significant_knee_with_clusters.csv")

lines = []
lines.append(r"\begin{table}[t]")
lines.append(r"\centering")
lines.append(r"\scriptsize")
lines.append(r"\renewcommand{\arraystretch}{1.15}")
lines.append(r"\begin{tabular}{p{6.4cm}rr}")
lines.append(r"\hline")
lines.append(r"Cluster representative & Pathways & Mean NES \\")
lines.append(r"\hline")

# Upregulated clusters
up_cl = cl[cl.Mean_NES > 0].sort_values("Mean_NES", ascending=False)
for _, row in up_cl.iterrows():
    name = row["Representative"].split(" R-HSA")[0]
    # Escape ampersands and underscores for LaTeX
    name = name.replace("&", r"\&").replace("_", r"\_")
    lines.append(f"  {name} & {row['Size']:.0f} & ${row['Mean_NES']:+.2f}$ \\\\")

lines.append(r"\hline")

# Downregulated clusters
dn_cl = cl[cl.Mean_NES < 0].sort_values("Mean_NES")
for _, row in dn_cl.iterrows():
    name = row["Representative"].split(" R-HSA")[0]
    name = name.replace("&", r"\&").replace("_", r"\_")
    lines.append(f"  {name} & {row['Size']:.0f} & ${row['Mean_NES']:+.2f}$ \\\\")

lines.append(r"\hline")
lines.append(r"\end{tabular}")
lines.append(
    r"\caption{Pathway clusters from the direct UC\,vs\,CD comparison. "
    r"Positive NES indicates higher enrichment in UC; negative indicates higher in CD. "
    r"The analysis comprised five datasets containing both UC and CD samples (293~total samples).}"
)
lines.append(r"\label{tab:pathway_clusters_uc_cd}")
lines.append(r"\end{table}")

table_tex = "\n".join(lines)
(TABDIR / "uc_vs_cd_clusters.tex").write_text(table_tex)
print(f"Saved {TABDIR / 'uc_vs_cd_clusters.tex'}")


# ══════════════════════════════════════════════════════════════════
#  Summary statistics for text
# ══════════════════════════════════════════════════════════════════

ds = pd.read_csv(DIRECT / "dataset_summary.csv")
gs = pd.read_csv(DIRECT / "gene_meta_significant.csv")

print("\n=== Summary for Results text ===")
print(f"Datasets: {len(ds)}")
print(f"Total samples: {ds.n_total.sum()}, UC: {ds.n_uc.sum()}, CD: {ds.n_cd.sum()}")
print(f"Sig genes: {len(gs)} (up in UC: {(gs.hedges_g>0).sum()}, down: {(gs.hedges_g<0).sum()})")
print(f"Median |g|: {gs.hedges_g.abs().median():.2f}")
print(f"Median I²: {gs.I2.median():.3f}")
print(f"Concordance r: {r:.2f}")
print(f"Knee-point pathways: {len(pk)}, Clusters: {len(cl)}")
print(f"  Upregulated clusters: {len(up_cl)}, Downregulated: {len(dn_cl)}")

# Key pathway NES in direct
pa = pd.read_csv(DIRECT / "pathway_meta_all.csv")
key_pw = [
    "Respiratory Electron Transport R-HSA-611105",
    "Citric Acid (TCA) Cycle And Respiratory Electron Transport R-HSA-1428517",
    "Mitochondrial Fatty Acid Beta-Oxidation R-HSA-77289",
    "Glucuronidation R-HSA-156588",
    "Eukaryotic Translation Elongation R-HSA-156842",
    "Peptide Chain Elongation R-HSA-156902",
    "Interferon Alpha/Beta Signaling R-HSA-909733",
    "Complement Cascade R-HSA-166658",
    "PD-1 Signaling R-HSA-389948",
]
print("\nKey pathways:")
for pw in key_pw:
    row = pa[pa.term == pw]
    if not row.empty:
        name = pw.split(" R-HSA")[0]
        print(f"  {name:<55s}  NES={row.combined_nes.values[0]:+.2f}  q={row.q_value.values[0]:.4f}")

# PPARGC1A/B
print("\nKey genes:")
for g in ["PPARGC1A", "PPARGC1B", "ESRRA"]:
    row = gs[gs.gene == g]
    if not row.empty:
        print(f"  {g}: g={row.hedges_g.values[0]:+.3f}, I²={row.I2.values[0]:.3f}, q={row.q_value.values[0]:.2e}")
