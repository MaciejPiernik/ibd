#!/usr/bin/env python3
"""
Generate figures for the UC meta-analysis paper:
1. Forest plots for key genes
2. I² diagnostic plots (histogram + scatter)
3. Stouffer's sensitivity check for NES pathway meta-analysis
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import warnings
warnings.filterwarnings('ignore')

# ── Styling ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 9,
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

OUTDIR = '/home/claude'

# ── Load data ────────────────────────────────────────────────────────────────
per_gene = pd.read_csv('/mnt/user-data/uploads/per_dataset_gene_stats_significant.csv')
meta_gene = pd.read_csv('/mnt/user-data/uploads/gene_meta_significant.csv')
per_gsea = pd.read_csv('/mnt/user-data/uploads/per_dataset_gsea.csv')
pw_sig = pd.read_csv('/mnt/user-data/uploads/pathway_meta_significant.csv')
pw_all_infl = pd.read_csv('/mnt/user-data/uploads/pathway_meta_all_inflamed.csv')
pw_all_uninfl = pd.read_csv('/mnt/user-data/uploads/pathway_meta_all_uninflamed.csv')

# ═════════════════════════════════════════════════════════════════════════════
# 1. FOREST PLOTS
# ═════════════════════════════════════════════════════════════════════════════

def forest_plot_multi(genes, per_df, meta_df, filename):
    """Create a multi-panel forest plot for several genes."""
    n_genes = len(genes)
    fig, axes = plt.subplots(1, n_genes, figsize=(4.2 * n_genes, 5.5), sharey=False)
    if n_genes == 1:
        axes = [axes]

    for ax, gene_name in zip(axes, genes):
        gdata = per_df[per_df.gene == gene_name].copy()
        gmeta = meta_df[meta_df.gene == gene_name].iloc[0]

        # Sort datasets by effect size for readability
        gdata = gdata.sort_values('hedges_g', ascending=True).reset_index(drop=True)

        y_positions = np.arange(len(gdata))
        ci_lo = gdata.hedges_g - 1.96 * np.sqrt(gdata.hedges_g_var)
        ci_hi = gdata.hedges_g + 1.96 * np.sqrt(gdata.hedges_g_var)

        # Per-dataset estimates
        # Scale marker by inverse variance (larger = more precise)
        weights = 1.0 / gdata.hedges_g_var
        marker_sizes = 30 + 150 * (weights / weights.max())

        for i, (y, g, lo, hi, ms) in enumerate(zip(y_positions, gdata.hedges_g, ci_lo, ci_hi, marker_sizes)):
            ax.plot([lo, hi], [y, y], color='#4a4a4a', linewidth=0.8, zorder=1)
            ax.scatter(g, y, s=ms, color='#2166ac', edgecolors='white', linewidth=0.5, zorder=2)

        # Combined estimate as diamond
        y_combined = len(gdata) + 0.8
        meta_lo = gmeta.ci_lo
        meta_hi = gmeta.ci_hi
        meta_g = gmeta.hedges_g

        diamond_x = [meta_lo, meta_g, meta_hi, meta_g]
        diamond_y = [y_combined, y_combined + 0.35, y_combined, y_combined - 0.35]
        ax.fill(diamond_x, diamond_y, color='#b2182b', alpha=0.85, zorder=3)

        # Reference line at zero
        ax.axvline(0, color='grey', linestyle='--', linewidth=0.6, alpha=0.6, zorder=0)

        # Labels
        labels = list(gdata.dataset)
        yticks = list(y_positions) + [y_combined]
        yticklabels = labels + ['Combined']

        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize=7.5)
        ax.set_xlabel("Hedges' g")
        ax.set_title(f'{gene_name}', fontweight='bold', fontsize=11)

        # Add annotation box
        txt = f"g = {meta_g:.2f}  [{meta_lo:.2f}, {meta_hi:.2f}]\n"
        txt += f"I² = {gmeta.I2:.0%}   τ² = {gmeta.tau2:.2f}\n"
        txt += f"q = {gmeta.q_value:.1e}   k = {int(gmeta.datasets)}"
        ax.text(0.97, 0.02, txt, transform=ax.transAxes, fontsize=7,
                verticalalignment='bottom', horizontalalignment='right',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#f7f7f7', edgecolor='#cccccc', alpha=0.9))

        # Add n_uc / n_ctrl on the right
        for i, (y, row) in enumerate(zip(y_positions, gdata.itertuples())):
            ax.text(1.0, y, f'  {int(row.n_uc)}/{int(row.n_ctrl)}',
                    transform=ax.get_yaxis_transform(), fontsize=6.5,
                    verticalalignment='center', color='#666666')

        ax.set_ylim(-0.8, y_combined + 1.0)

    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/{filename}', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {filename}")


# Metabolic / mitochondrial genes (downregulated)
forest_plot_multi(
    ['PPARGC1A', 'OPA1', 'ACSF2'],
    per_gene, meta_gene,
    'forest_metabolic.pdf'
)

# Top upregulated + one more for contrast
forest_plot_multi(
    ['LPCAT1', 'PEX1'],
    per_gene, meta_gene,
    'forest_other.pdf'
)


# ═════════════════════════════════════════════════════════════════════════════
# 2. I² DIAGNOSTIC PLOTS
# ═════════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# Panel A: Histogram of I²
ax = axes[0]
i2_vals = meta_gene.I2.values
ax.hist(i2_vals, bins=50, color='#4393c3', edgecolor='white', linewidth=0.3, alpha=0.85)
ax.axvline(np.median(i2_vals), color='#b2182b', linestyle='--', linewidth=1.2,
           label=f'Median = {np.median(i2_vals):.2f}')
ax.set_xlabel('I²')
ax.set_ylabel('Number of genes')
ax.set_title('A.  Heterogeneity distribution', fontweight='bold', loc='left')
ax.legend(frameon=False)

# Panel B: |g| vs I², colored by -log10(q)
ax = axes[1]
abs_g = np.abs(meta_gene.hedges_g.values)
i2 = meta_gene.I2.values
neg_log_q = -np.log10(meta_gene.q_value.values + 1e-300)
# Cap for color scale
neg_log_q_cap = np.clip(neg_log_q, 0, 30)

sc = ax.scatter(abs_g, i2, c=neg_log_q_cap, cmap='RdYlBu_r', s=4, alpha=0.4,
                edgecolors='none', rasterized=True)
cbar = plt.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
cbar.set_label('−log₁₀(q)')

ax.set_xlabel("|Hedges' g|")
ax.set_ylabel('I²')
ax.set_title('B.  Effect size vs heterogeneity', fontweight='bold', loc='left')

# Add Spearman correlation
rho, pval = stats.spearmanr(abs_g, i2)
ax.text(0.97, 0.05, f'Spearman ρ = {rho:.2f}',
        transform=ax.transAxes, fontsize=8, ha='right',
        bbox=dict(facecolor='white', edgecolor='#cccccc', alpha=0.9, boxstyle='round,pad=0.3'))

plt.tight_layout()
plt.savefig(f'{OUTDIR}/i2_diagnostics.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved i2_diagnostics.pdf")

# Print summary stats for the text
print(f"\nI² summary statistics:")
print(f"  Median: {np.median(i2_vals):.2f}")
print(f"  Mean: {np.mean(i2_vals):.2f}")
print(f"  % genes with I²=0: {(i2_vals == 0).sum() / len(i2_vals):.1%}")
print(f"  % genes with I²<0.25: {(i2_vals < 0.25).sum() / len(i2_vals):.1%}")
print(f"  % genes with I²>0.75: {(i2_vals > 0.75).sum() / len(i2_vals):.1%}")
print(f"  Spearman |g| vs I²: rho={rho:.3f}, p={pval:.2e}")

# Dual ranking analysis
meta_gene_sorted_q = meta_gene.sort_values('q_value')
meta_gene_sorted_g = meta_gene.sort_values('hedges_g', key=abs, ascending=False)

print(f"\nTop 20 by q-value — median I²: {meta_gene_sorted_q.head(20).I2.median():.2f}, median |g|: {meta_gene_sorted_q.head(20).hedges_g.abs().median():.2f}")
print(f"Top 20 by |g|     — median I²: {meta_gene_sorted_g.head(20).I2.median():.2f}, median |g|: {meta_gene_sorted_g.head(20).hedges_g.abs().median():.2f}")
print(f"Top 100 by q-value — median I²: {meta_gene_sorted_q.head(100).I2.median():.2f}")
print(f"Top 100 by |g|     — median I²: {meta_gene_sorted_g.head(100).I2.median():.2f}")


# ═════════════════════════════════════════════════════════════════════════════
# 3. STOUFFER'S SENSITIVITY CHECK
# ═════════════════════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("STOUFFER'S SENSITIVITY CHECK")
print("="*70)

# For each pathway, combine per-dataset p-values using Stouffer's method
# weighted by sqrt(sample_size)

# Get sample sizes per dataset from per_gene (use first gene per dataset)
dataset_sizes = per_gene.groupby('dataset').agg(
    n_total=('n_uc', 'first')  # just need to identify datasets
).reset_index()
# Actually get total n per dataset
ds_n = per_gene.groupby('dataset').apply(
    lambda x: x.iloc[0].n_uc + x.iloc[0].n_ctrl
).to_dict()

# Get unique pathways from the per-dataset GSEA
pathway_terms = per_gsea.Term.unique()

stouffer_results = []
for term in pathway_terms:
    pdata = per_gsea[per_gsea.Term == term].copy()
    if len(pdata) < 2:
        continue

    # Determine direction from NES signs
    n_pos = (pdata.NES > 0).sum()
    n_neg = (pdata.NES < 0).sum()
    direction = 'up' if n_pos >= n_neg else 'down'

    # Convert to one-sided p-values in the majority direction
    z_scores = []
    weights = []
    for _, row in pdata.iterrows():
        p = row['NOM p-val']
        # Clamp p-values
        p = max(p, 1e-10)
        p = min(p, 1 - 1e-10)

        nes = row.NES
        # Convert to one-sided: if gene set is enriched in the majority direction
        if direction == 'up':
            if nes > 0:
                p_one = p / 2
            else:
                p_one = 1 - p / 2
        else:
            if nes < 0:
                p_one = p / 2
            else:
                p_one = 1 - p / 2

        z = stats.norm.ppf(1 - p_one)
        w = np.sqrt(ds_n.get(row.dataset, 50))  # sqrt(n) weight
        z_scores.append(z)
        weights.append(w)

    z_scores = np.array(z_scores)
    weights = np.array(weights)

    # Weighted Stouffer
    z_combined = np.sum(weights * z_scores) / np.sqrt(np.sum(weights**2))
    p_combined = 2 * (1 - stats.norm.cdf(abs(z_combined)))  # two-sided

    stouffer_results.append({
        'term': term,
        'z_stouffer': z_combined,
        'p_stouffer': p_combined,
        'direction': direction,
        'n_datasets': len(pdata)
    })

stouffer_df = pd.DataFrame(stouffer_results)

# BH correction
from statsmodels.stats.multitest import multipletests
_, stouffer_df['q_stouffer'], _, _ = multipletests(stouffer_df.p_stouffer, method='fdr_bh')

# Merge with REML results
# Need the full inflamed pathway meta for comparison
pw_all = pw_all_infl.copy()
pw_all['term_clean'] = pw_all.term
merged = stouffer_df.merge(pw_all[['term', 'combined_nes', 'q_value', 'direction']],
                            left_on='term', right_on='term', how='inner',
                            suffixes=('_stouffer', '_reml'))

# Concordance analysis
merged['sig_reml'] = merged.q_value < 0.05
merged['sig_stouffer'] = merged.q_stouffer < 0.05

both_sig = (merged.sig_reml & merged.sig_stouffer).sum()
reml_only = (merged.sig_reml & ~merged.sig_stouffer).sum()
stouffer_only = (~merged.sig_reml & merged.sig_stouffer).sum()
neither = (~merged.sig_reml & ~merged.sig_stouffer).sum()
total = len(merged)

print(f"\nPathway-level concordance (q < 0.05):")
print(f"  Both significant:     {both_sig} ({both_sig/total:.1%})")
print(f"  REML only:            {reml_only} ({reml_only/total:.1%})")
print(f"  Stouffer only:        {stouffer_only} ({stouffer_only/total:.1%})")
print(f"  Neither:              {neither} ({neither/total:.1%})")
print(f"  Total pathways:       {total}")

# Among REML significant, what fraction also Stouffer significant?
reml_sig = merged[merged.sig_reml]
concordance = reml_sig.sig_stouffer.mean()
print(f"\n  Of {len(reml_sig)} REML-significant pathways, {reml_sig.sig_stouffer.sum()} ({concordance:.1%}) also Stouffer-significant")

# Direction concordance among both-significant
both = merged[merged.sig_reml & merged.sig_stouffer]
same_dir = (np.sign(both.combined_nes) == np.sign(both.z_stouffer)).mean()
print(f"  Direction concordance among jointly significant: {same_dir:.1%}")

# Correlation of z-scores
r_z, p_z = stats.pearsonr(merged.combined_nes, merged.z_stouffer)
print(f"  Pearson r (NES vs Stouffer z): {r_z:.3f} (p={p_z:.2e})")

# Focus on the knee-filtered pathways from the paper
pw_sig_terms = set(pw_sig.term.values)
merged['in_knee'] = merged.term.isin(pw_sig_terms)
knee_merged = merged[merged.in_knee]
knee_concordance = knee_merged.sig_stouffer.mean()
print(f"\n  Of {len(knee_merged)} knee-filtered pathways, {knee_merged.sig_stouffer.sum()} ({knee_concordance:.1%}) also Stouffer-significant")

# Plot: NES vs Stouffer z
fig, ax = plt.subplots(figsize=(5.5, 5))
# Non-significant in grey
nonsig = merged[~merged.sig_reml & ~merged.sig_stouffer]
ax.scatter(nonsig.combined_nes, nonsig.z_stouffer, s=6, alpha=0.25, color='#cccccc',
           edgecolors='none', rasterized=True, label='Neither significant')

# Significant in both
both_df = merged[merged.sig_reml & merged.sig_stouffer]
ax.scatter(both_df.combined_nes, both_df.z_stouffer, s=10, alpha=0.5, color='#2166ac',
           edgecolors='none', rasterized=True, label='Both significant')

# Discordant
reml_only_df = merged[merged.sig_reml & ~merged.sig_stouffer]
stouffer_only_df = merged[~merged.sig_reml & merged.sig_stouffer]
if len(reml_only_df) > 0:
    ax.scatter(reml_only_df.combined_nes, reml_only_df.z_stouffer, s=15, alpha=0.7,
               color='#d6604d', edgecolors='none', marker='s', label=f'REML only (n={len(reml_only_df)})')
if len(stouffer_only_df) > 0:
    ax.scatter(stouffer_only_df.combined_nes, stouffer_only_df.z_stouffer, s=15, alpha=0.7,
               color='#4dac26', edgecolors='none', marker='^', label=f'Stouffer only (n={len(stouffer_only_df)})')

ax.axhline(0, color='grey', linewidth=0.5, linestyle='--')
ax.axvline(0, color='grey', linewidth=0.5, linestyle='--')
ax.set_xlabel('REML combined NES')
ax.set_ylabel("Stouffer's weighted z")
ax.set_title(f'Pathway meta-analysis: REML vs Stouffer\n(Pearson r = {r_z:.2f})', fontsize=10)
ax.legend(fontsize=7, loc='upper left', frameon=True, fancybox=True)

plt.tight_layout()
plt.savefig(f'{OUTDIR}/stouffer_concordance.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("\nSaved stouffer_concordance.pdf")


# ═════════════════════════════════════════════════════════════════════════════
# 4. CONSTITUTIVE GENE ANALYSIS — summary for paper text
# ═════════════════════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("CONSTITUTIVE GENE ANALYSIS")
print("="*70)

comp = pd.read_csv('/mnt/user-data/uploads/comparison_stats.csv')

const_down = comp[(comp.category == 'constitutive') & (comp.sig_uninfl == True) & (comp.direction_uninfl == 'down')]
const_down = const_down.sort_values('hedges_g_uninfl')

print(f"\n27 constitutively downregulated genes (sorted by uninflamed g):")
print(f"{'Gene':<12} {'g_infl':>8} {'g_uninfl':>10} {'Δg':>8} {'q_Δ':>8}")
print("-" * 50)
for _, r in const_down.iterrows():
    print(f"{r.gene:<12} {r.hedges_g_infl:>8.2f} {r.hedges_g_uninfl:>10.2f} {r.delta_g:>8.2f} {r.delta_q:>8.3f}")

# Mitochondrial/metabolic subset
mito_genes = ['OPA1', 'DNAJA3', 'PEX1', 'ENDOG', 'ETHE1', 'PDP2', 'PANK4', 'TMEM65']
print(f"\nMitochondrial/metabolic subset:")
for g in mito_genes:
    row = const_down[const_down.gene == g]
    if len(row) > 0:
        r = row.iloc[0]
        print(f"  {g}: g_infl={r.hedges_g_infl:.2f}, g_uninfl={r.hedges_g_uninfl:.2f}")

print("\nDone.")