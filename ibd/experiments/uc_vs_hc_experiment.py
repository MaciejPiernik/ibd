"""
UC vs Healthy Control meta-analysis experiment.

Per-dataset: Welch's t-test + Hedges' g for each gene, GSEA prerank for pathways.
Cross-dataset: REML random-effects meta-analysis on effect sizes + Robust Rank Aggregation.
Filtering: minimum dataset representation, direction consistency, FDR control.

Outputs
-------
dataset_summary.csv                 One row per included dataset.
per_dataset_gene_stats.csv          Gene-level stats for every dataset.
per_dataset_gsea.csv                GSEA results for every dataset (with leading edge).
gene_meta_all.csv                   Random-effects meta-analysis for all genes.
gene_meta_significant.csv           Filtered gene table (paper Table / Supplement).
gene_rra.csv                        Robust Rank Aggregation results for all genes.
pathway_meta_all.csv                Random-effects meta-analysis for all pathways.
pathway_meta_significant.csv        Filtered pathway table.
pathway_leading_edge_detail.csv     Per-gene leading-edge frequency for significant pathways.
"""

import logging
import os
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import gseapy as gp
from scipy import stats
from scipy.optimize import minimize_scalar
from statsmodels.stats.multitest import multipletests
from tqdm.auto import tqdm

from ibd.core.data.utils import read_all_datasets, read_all_metadata
from ibd.core.data.processing.normalization import robust_zscore_normalization_per_dataset
from ibd.core.utils.genes import entrez_id_to_gene_symbol


DEFAULT_GENESET_LIBRARIES = [
    'WikiPathway_2023_Human',
    'KEGG_2021_Human',
    'GO_Biological_Process_2023',
    'GO_Molecular_Function_2023',
    'GO_Cellular_Component_2023',
    'Reactome_2022',
]


# ===================================================================
# Configuration
# ===================================================================

@dataclass
class UcVsHcConfig:
    data_dir: str
    output_dir: str

    # Per-dataset testing
    alpha: float = 0.05
    min_samples: int = 5
    min_per_group: int = 3

    # Normalisation
    normalize: str = 'none'  # 'none' | 'robust_zscore'
    skip_gene_mapping: bool = False

    # Meta-analysis filtering
    min_datasets: int = 4             # gene/pathway must appear in >= this many datasets
    direction_threshold: float = 0.8  # fraction of datasets that must agree on direction

    # GSEA
    geneset_libraries: List[str] = field(default_factory=lambda: list(DEFAULT_GENESET_LIBRARIES))
    gsea_permutations: int = 1000
    gsea_min_size: int = 5
    gsea_max_size: int = 1000
    gsea_seed: int = 23


# ===================================================================
# Data loading & sample selection
# ===================================================================

def _standardize_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    result = metadata.copy()
    if 'inflammation' not in result.columns:
        result['inflammation'] = None
    if 'time_of_biopsy' not in result.columns:
        result['time_of_biopsy'] = None
    result['inflammation'] = np.where(
        result['inflammation'].isna() & (result['disease'] == 'UC'),
        'Inflamed', result['inflammation'],
    )
    result['inflammation'] = np.where(
        result['inflammation'].isna() & (result['disease'] == 'Ctrl'),
        'Uninflamed', result['inflammation'],
    )
    return result


def _select_uc_vs_ctrl(metadata: pd.DataFrame) -> pd.Series:
    time_ok = metadata['time_of_biopsy'].isna() | metadata['time_of_biopsy'].isin(['W0', 'Before'])
    mask = (
        ((metadata['disease'] == 'UC') & (metadata['inflammation'] == 'Inflamed'))
        | ((metadata['disease'] == 'Ctrl') & (metadata['inflammation'] == 'Uninflamed'))
    ) & time_ok
    return metadata.loc[mask, 'disease']


def _collapse_duplicate_genes(data: pd.DataFrame) -> pd.DataFrame:
    if data.columns.duplicated().any():
        data = data.T.groupby(level=0).mean().T
    return data


# ===================================================================
# Per-dataset gene-level statistics
# ===================================================================

def _hedges_g(mean_uc, mean_ctrl, var_uc, var_ctrl, n_uc, n_ctrl):
    """Bias-corrected standardised mean difference (Hedges' g)."""
    if n_uc < 2 or n_ctrl < 2:
        return np.nan
    pooled_var = ((n_uc - 1) * var_uc + (n_ctrl - 1) * var_ctrl) / (n_uc + n_ctrl - 2)
    if pooled_var <= 0:
        return np.nan
    d = (mean_uc - mean_ctrl) / np.sqrt(pooled_var)
    J = 1 - 3 / (4 * (n_uc + n_ctrl) - 9)
    return d * J


def _hedges_g_variance(g, n_uc, n_ctrl):
    """Approximate sampling variance of Hedges' g."""
    n = n_uc + n_ctrl
    return n / (n_uc * n_ctrl) + g ** 2 / (2 * n)


def _compute_gene_stats(data: pd.DataFrame, labels: pd.Series) -> pd.DataFrame:
    ctrl_mask = labels == 'Ctrl'
    uc_mask = labels == 'UC'
    results = []

    for gene in tqdm(data.columns, desc='Genes', unit='gene', leave=False):
        ctrl_vals = data.loc[ctrl_mask, gene].dropna().astype(float)
        uc_vals = data.loc[uc_mask, gene].dropna().astype(float)
        n_ctrl, n_uc = len(ctrl_vals), len(uc_vals)
        if n_ctrl < 2 or n_uc < 2:
            continue

        t_stat, p_value = stats.ttest_ind(uc_vals, ctrl_vals, equal_var=False)
        mean_ctrl, mean_uc = ctrl_vals.mean(), uc_vals.mean()
        g = _hedges_g(mean_uc, mean_ctrl, uc_vals.var(ddof=1), ctrl_vals.var(ddof=1), n_uc, n_ctrl)
        g_var = _hedges_g_variance(g, n_uc, n_ctrl) if not np.isnan(g) else np.nan

        results.append({
            'gene': gene,
            'n_uc': n_uc,
            'n_ctrl': n_ctrl,
            'mean_uc': mean_uc,
            'mean_ctrl': mean_ctrl,
            'mean_diff': mean_uc - mean_ctrl,
            't_stat': t_stat,
            'p_value': p_value,
            'hedges_g': g,
            'hedges_g_var': g_var,
        })

    df = pd.DataFrame(results)
    if df.empty:
        return df
    df['q_value'] = multipletests(df['p_value'], method='fdr_bh')[1]
    # Rank within dataset by |t_stat| (1 = strongest)
    df['rank'] = df['t_stat'].abs().rank(ascending=False, method='min').astype(int)
    return df


# ===================================================================
# Random-effects meta-analysis (REML)
# ===================================================================

def _cochran_Q_and_I2(effects: np.ndarray, variances: np.ndarray):
    """Compute Cochran's Q statistic and I² from fixed-effects model."""
    w = 1.0 / variances
    mu_fixed = np.sum(w * effects) / w.sum()
    Q = np.sum(w * (effects - mu_fixed) ** 2)
    df = len(effects) - 1
    I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    return Q, I2


def _reml_tau2(effects: np.ndarray, variances: np.ndarray) -> float:
    """
    Estimate between-study variance τ² via REML.

    Maximises the restricted (residual) log-likelihood directly
    using bounded scalar optimisation.
    """
    k = len(effects)
    if k < 2:
        return 0.0

    def neg_reml_loglik(tau2):
        w = 1.0 / (variances + tau2)
        w_sum = w.sum()
        mu = np.sum(w * effects) / w_sum
        resid = effects - mu
        ll = (-0.5 * np.sum(np.log(variances + tau2))
              - 0.5 * np.log(w_sum)
              - 0.5 * np.sum(w * resid ** 2))
        return -ll

    upper = max(10.0 * np.var(effects), 1.0)
    result = minimize_scalar(neg_reml_loglik, bounds=(0, upper), method='bounded')
    return max(0.0, result.x)


def _random_effects_meta(effects: np.ndarray, variances: np.ndarray):
    """
    REML random-effects meta-analysis.

    Returns: (combined_effect, se, p_value, tau2, I2)
    """
    mask = np.isfinite(effects) & np.isfinite(variances) & (variances > 0)
    effects, variances = effects[mask], variances[mask]
    k = len(effects)

    if k == 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    if k == 1:
        se = np.sqrt(variances[0])
        z = effects[0] / se if se > 0 else np.nan
        p = 2 * stats.norm.sf(abs(z)) if np.isfinite(z) else np.nan
        return effects[0], se, p, 0.0, 0.0

    # Estimate τ² via REML
    tau2 = _reml_tau2(effects, variances)

    # Cochran's Q and I² (computed from fixed-effects weights, standard definition)
    _, I2 = _cochran_Q_and_I2(effects, variances)

    # Combined estimate with random-effects weights
    w_re = 1.0 / (variances + tau2)
    w_re_sum = w_re.sum()
    mu_re = np.sum(w_re * effects) / w_re_sum
    se_re = np.sqrt(1.0 / w_re_sum)

    z = mu_re / se_re
    p = 2 * stats.norm.sf(abs(z))

    return mu_re, se_re, p, tau2, I2


# ===================================================================
# Gene-level meta-analysis
# ===================================================================

def _meta_analyse_genes(per_dataset: pd.DataFrame, config: UcVsHcConfig) -> pd.DataFrame:
    records = []
    for gene, grp in per_dataset.groupby('gene'):
        effects = grp['hedges_g'].values
        variances = grp['hedges_g_var'].values
        mu, se, p, tau2, I2 = _random_effects_meta(effects, variances)

        valid = effects[np.isfinite(effects)]
        n_up = int(np.sum(valid > 0))
        n_down = int(np.sum(valid < 0))
        n_valid = n_up + n_down
        direction_ratio = max(n_up, n_down) / n_valid if n_valid > 0 else 0
        direction = 'up' if n_up >= n_down else 'down'
        sign_p = stats.binomtest(max(n_up, n_down), n_valid, 0.5).pvalue if n_valid >= 2 else np.nan

        records.append({
            'gene': gene,
            'datasets': len(grp),
            'hedges_g': mu,
            'se': se,
            'ci_lo': mu - 1.96 * se if np.isfinite(se) else np.nan,
            'ci_hi': mu + 1.96 * se if np.isfinite(se) else np.nan,
            'p_value': p,
            'tau2': tau2,
            'I2': I2,
            'direction': direction,
            'n_up': n_up,
            'n_down': n_down,
            'direction_ratio': direction_ratio,
            'sign_test_p': sign_p,
        })

    df = pd.DataFrame(records)
    if df.empty:
        return df
    df['q_value'] = multipletests(df['p_value'].fillna(1), method='fdr_bh')[1]
    return df.sort_values(['q_value', 'p_value']).reset_index(drop=True)


def _filter_genes(meta: pd.DataFrame, config: UcVsHcConfig) -> pd.DataFrame:
    return meta[
        (meta['datasets'] >= config.min_datasets)
        & (meta['direction_ratio'] >= config.direction_threshold)
        & (meta['q_value'] < config.alpha)
    ].sort_values(['q_value', 'p_value']).reset_index(drop=True)


# ===================================================================
# Robust Rank Aggregation (RRA)
# ===================================================================

def _rra_rho(normalized_ranks: np.ndarray) -> float:
    """
    Compute the RRA rho statistic (Kolde et al., 2012).

    For sorted normalized ranks u_(1) <= ... <= u_(k), compute:
        rho = min_j  Beta_CDF(u_(j); j, k-j+1)

    Small rho means the gene appears near the top of ranked lists
    more often than expected under uniform null.
    """
    k = len(normalized_ranks)
    if k == 0:
        return np.nan
    u = np.sort(normalized_ranks)
    betas = np.array([
        stats.beta.cdf(u[j], j + 1, k - j)
        for j in range(k)
    ])
    return float(np.min(betas))


def _rra_pvalue(rho: float, k: int) -> float:
    """
    Approximate p-value for the rho statistic.
    P(rho <= x) ≈ x * k  for small x (Bonferroni correction on k beta tests).
    Exact for the minimum of independent uniform RVs; conservative otherwise.
    """
    if np.isnan(rho):
        return np.nan
    return min(1.0, rho * k)


def _compute_rra(per_dataset_genes: pd.DataFrame, config: UcVsHcConfig) -> pd.DataFrame:
    """
    Robust Rank Aggregation on gene rankings across datasets.

    Each gene is ranked by |t_stat| within its dataset (1 = strongest).
    Ranks are normalised to [0, 1] by dividing by total genes in that dataset.
    RRA is computed on the normalised ranks.

    Additionally, directional RRA is computed:
      - rra_up: ranks by t_stat descending (most upregulated = rank 1)
      - rra_down: ranks by -t_stat descending (most downregulated = rank 1)
    """
    # Precompute total genes per dataset for normalisation
    genes_per_dataset = per_dataset_genes.groupby('dataset')['gene'].count().to_dict()

    # Compute directional ranks per dataset
    all_ranks = []
    for ds_name, grp in per_dataset_genes.groupby('dataset'):
        n_genes = genes_per_dataset[ds_name]
        sub = grp[['gene', 't_stat', 'rank']].copy()
        # |t_stat| rank already exists as 'rank'
        sub['norm_rank'] = sub['rank'] / n_genes
        # Directional ranks
        sub['rank_up'] = sub['t_stat'].rank(ascending=False, method='min').astype(int)
        sub['rank_down'] = sub['t_stat'].rank(ascending=True, method='min').astype(int)
        sub['norm_rank_up'] = sub['rank_up'] / n_genes
        sub['norm_rank_down'] = sub['rank_down'] / n_genes
        sub['dataset'] = ds_name
        all_ranks.append(sub[['gene', 'dataset', 'norm_rank', 'norm_rank_up', 'norm_rank_down']])

    rank_df = pd.concat(all_ranks, ignore_index=True)

    records = []
    for gene, grp in rank_df.groupby('gene'):
        k = len(grp)
        nr = grp['norm_rank'].values
        nr_up = grp['norm_rank_up'].values
        nr_down = grp['norm_rank_down'].values

        rho = _rra_rho(nr)
        rho_up = _rra_rho(nr_up)
        rho_down = _rra_rho(nr_down)

        records.append({
            'gene': gene,
            'datasets': k,
            'rra_rho': rho,
            'rra_p': _rra_pvalue(rho, k),
            'rra_rho_up': rho_up,
            'rra_p_up': _rra_pvalue(rho_up, k),
            'rra_rho_down': rho_down,
            'rra_p_down': _rra_pvalue(rho_down, k),
        })

    df = pd.DataFrame(records)
    if df.empty:
        return df

    df['rra_q'] = multipletests(df['rra_p'].fillna(1), method='fdr_bh')[1]
    df['rra_q_up'] = multipletests(df['rra_p_up'].fillna(1), method='fdr_bh')[1]
    df['rra_q_down'] = multipletests(df['rra_p_down'].fillna(1), method='fdr_bh')[1]

    return df.sort_values('rra_p').reset_index(drop=True)


# ===================================================================
# Per-dataset GSEA
# ===================================================================

def _parse_tag_percent(tag_str: str) -> Tuple[int, int]:
    """Parse '35/156' -> (n_leading_edge, gene_set_size)."""
    try:
        parts = str(tag_str).split('/')
        return int(parts[0]), int(parts[1])
    except (ValueError, IndexError):
        return 0, 0


def _run_gsea_per_dataset(
    stats_df: pd.DataFrame,
    dataset_name: str,
    libraries: List[str],
    config: UcVsHcConfig,
) -> pd.DataFrame:
    ranking = stats_df[['gene', 't_stat']].dropna().copy()
    if ranking.empty:
        return pd.DataFrame()

    if ranking['t_stat'].duplicated().any():
        rng = np.random.default_rng(config.gsea_seed)
        ranking['t_stat'] += rng.normal(0, 1e-9, len(ranking))

    ranking = ranking.sort_values('t_stat', ascending=False).set_index('gene')

    all_res = []
    for i, library in enumerate(libraries, 1):
        logging.info('  GSEA library %d/%d: %s', i, len(libraries), library)
        try:
            pre_res = gp.prerank(
                rnk=ranking,
                gene_sets=library,
                threads=1,
                min_size=config.gsea_min_size,
                max_size=config.gsea_max_size,
                permutation_num=config.gsea_permutations,
                outdir=None,
                seed=config.gsea_seed,
                verbose=False,
            )
        except Exception as exc:
            logging.warning('GSEA failed for %s / %s: %s', dataset_name, library, exc)
            continue

        res = pre_res.res2d.copy()
        res['NES'] = pd.to_numeric(res['NES'], errors='coerce')
        res['NOM p-val'] = pd.to_numeric(res['NOM p-val'], errors='coerce')
        res['FDR q-val'] = pd.to_numeric(res['FDR q-val'], errors='coerce')

        # Parse Tag % -> n_leading_edge / gene_set_size
        parsed = res['Tag %'].apply(_parse_tag_percent)
        res['n_leading_edge'] = parsed.apply(lambda x: x[0])
        res['gene_set_size'] = parsed.apply(lambda x: x[1])

        if 'Lead_genes' not in res.columns:
            res['Lead_genes'] = ''

        res['dataset'] = dataset_name
        res['library'] = library
        res['rank'] = res['NES'].abs().rank(ascending=False, method='min').astype(int)

        all_res.append(res)

    if not all_res:
        return pd.DataFrame()

    combined = pd.concat(all_res, ignore_index=True)

    keep_cols = [
        'dataset', 'library', 'Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val',
        'gene_set_size', 'n_leading_edge', 'Lead_genes', 'rank',
    ]
    return combined[[c for c in keep_cols if c in combined.columns]]


# ===================================================================
# NES variance estimation
# ===================================================================

def _nes_variance_from_pvalue(nes: float, p: float) -> float:
    """
    Approximate NES sampling variance from nominal p-value.
    Under permutation null NES ~ N(0,1), so z = Phi^-1(1-p/2), var ~ (NES/z)^2.
    """
    p = np.clip(p, 1e-10, 1 - 1e-10)
    z = stats.norm.isf(p / 2)
    if abs(z) < 1e-6:
        return np.nan
    return (nes / z) ** 2


# ===================================================================
# Pathway-level meta-analysis
# ===================================================================

def _aggregate_leading_edge(grp: pd.DataFrame) -> Tuple[Dict[str, int], int]:
    """
    Aggregate leading-edge genes across datasets for one pathway.
    Returns (gene -> count, n_datasets_with_leading_edge).
    """
    gene_counts: Dict[str, int] = Counter()
    n_datasets = 0

    for _, row in grp.iterrows():
        lead_str = row.get('Lead_genes', '')
        if pd.isna(lead_str) or lead_str == '':
            continue
        genes = [g.strip() for g in str(lead_str).split(';') if g.strip()]
        if genes:
            gene_counts.update(genes)
            n_datasets += 1

    return dict(gene_counts), n_datasets


def _meta_analyse_pathways(per_dataset: pd.DataFrame, config: UcVsHcConfig) -> pd.DataFrame:
    records = []

    for (library, term), grp in per_dataset.groupby(['library', 'Term']):
        k = len(grp)
        nes_vals = grp['NES'].values.astype(float)
        p_vals = grp['NOM p-val'].values.astype(float)
        variances = np.array([_nes_variance_from_pvalue(n, p) for n, p in zip(nes_vals, p_vals)])

        mu, se, p, tau2, I2 = _random_effects_meta(nes_vals, variances)

        valid = nes_vals[np.isfinite(nes_vals)]
        n_up = int(np.sum(valid > 0))
        n_down = int(np.sum(valid < 0))
        n_valid = n_up + n_down
        direction_ratio = max(n_up, n_down) / n_valid if n_valid > 0 else 0
        direction = 'up' if n_up >= n_down else 'down'
        sign_p = stats.binomtest(max(n_up, n_down), n_valid, 0.5).pvalue if n_valid >= 2 else np.nan

        # Gene set size and leading edge stats
        gs_sizes = grp['gene_set_size'].values if 'gene_set_size' in grp.columns else np.array([])
        n_leads = grp['n_leading_edge'].values if 'n_leading_edge' in grp.columns else np.array([])
        valid_size_mask = gs_sizes > 0
        median_gs_size = int(np.median(gs_sizes[valid_size_mask])) if np.any(valid_size_mask) else 0
        lead_fracs = n_leads[valid_size_mask] / gs_sizes[valid_size_mask] if np.any(valid_size_mask) else np.array([])
        mean_lead_fraction = float(np.mean(lead_fracs)) if len(lead_fracs) > 0 else 0.0

        # Aggregate leading edge
        gene_counts, n_ds_with_lead = _aggregate_leading_edge(grp)
        core_threshold = max(1, n_ds_with_lead / 2)
        core_genes = sorted([g for g, c in gene_counts.items() if c >= core_threshold])
        all_lead_genes = sorted(gene_counts.keys())

        records.append({
            'library': library,
            'term': term,
            'datasets': k,
            'combined_nes': mu,
            'se': se,
            'ci_lo': mu - 1.96 * se if np.isfinite(se) else np.nan,
            'ci_hi': mu + 1.96 * se if np.isfinite(se) else np.nan,
            'p_value': p,
            'tau2': tau2,
            'I2': I2,
            'direction': direction,
            'n_up': n_up,
            'n_down': n_down,
            'direction_ratio': direction_ratio,
            'sign_test_p': sign_p,
            'median_gene_set_size': median_gs_size,
            'mean_lead_fraction': mean_lead_fraction,
            'lead_genes_all': ';'.join(all_lead_genes),
            'lead_genes_all_size': len(all_lead_genes),
            'lead_genes_core': ';'.join(core_genes),
            'lead_genes_core_size': len(core_genes),
        })

    df = pd.DataFrame(records)
    if df.empty:
        return df
    df['q_value'] = multipletests(df['p_value'].fillna(1), method='fdr_bh')[1]
    return df.sort_values(['q_value', 'p_value']).reset_index(drop=True)


def _filter_pathways(meta: pd.DataFrame, config: UcVsHcConfig) -> pd.DataFrame:
    return meta[
        (meta['datasets'] >= config.min_datasets)
        & (meta['direction_ratio'] >= config.direction_threshold)
        & (meta['q_value'] < config.alpha)
    ].sort_values(['q_value', 'p_value']).reset_index(drop=True)


# ===================================================================
# Pathway leading-edge detail table
# ===================================================================

def _build_leading_edge_detail(
    sig_pathways: pd.DataFrame,
    per_dataset_gsea: pd.DataFrame,
    gene_meta: pd.DataFrame,
) -> pd.DataFrame:
    """
    For each significant pathway, list every leading-edge gene
    with its cross-dataset frequency and gene-level meta-analysis results.
    """
    gene_lookup = {}
    if gene_meta is not None and not gene_meta.empty:
        gene_lookup = gene_meta.set_index('gene')[
            ['hedges_g', 'q_value', 'direction', 'datasets']
        ].to_dict('index')

    records = []
    for _, pw_row in sig_pathways.iterrows():
        lib, term = pw_row['library'], pw_row['term']
        mask = (per_dataset_gsea['library'] == lib) & (per_dataset_gsea['Term'] == term)
        pw_ds = per_dataset_gsea.loc[mask]
        n_tested = len(pw_ds)

        gene_counts = Counter()
        for _, ds_row in pw_ds.iterrows():
            lead_str = ds_row.get('Lead_genes', '')
            if pd.isna(lead_str) or lead_str == '':
                continue
            for g in str(lead_str).split(';'):
                g = g.strip()
                if g:
                    gene_counts[g] += 1

        for gene, count in gene_counts.items():
            gm = gene_lookup.get(gene, {})
            records.append({
                'library': lib,
                'term': term,
                'gene': gene,
                'n_datasets_in_lead': count,
                'n_datasets_tested': n_tested,
                'lead_frequency': count / n_tested if n_tested > 0 else 0,
                'gene_hedges_g': gm.get('hedges_g', np.nan),
                'gene_q_value': gm.get('q_value', np.nan),
                'gene_direction': gm.get('direction', ''),
                'gene_datasets': gm.get('datasets', 0),
            })

    df = pd.DataFrame(records)
    if not df.empty:
        df = df.sort_values(['library', 'term', 'lead_frequency'], ascending=[True, True, False])
    return df.reset_index(drop=True)


# ===================================================================
# Main experiment
# ===================================================================

def run_uc_vs_hc_experiment(config: UcVsHcConfig) -> None:
    if config.normalize not in ('none', 'robust_zscore'):
        raise ValueError("normalize must be 'none' or 'robust_zscore'")

    # ---- Load & prepare ----
    logging.info('Loading datasets from %s', config.data_dir)
    dataset, dataset_labels = read_all_datasets(config.data_dir, dropna=False)
    metadata = read_all_metadata(config.data_dir)

    common_index = dataset.index.intersection(metadata.index)
    dataset = dataset.loc[common_index]
    dataset_labels = dataset_labels.loc[common_index]
    metadata = metadata.loc[common_index]
    metadata = _standardize_metadata(metadata)

    if not config.skip_gene_mapping:
        try:
            dataset = entrez_id_to_gene_symbol(dataset)
        except Exception as exc:
            logging.warning('Gene ID mapping failed: %s', exc)

    dataset = _collapse_duplicate_genes(dataset)

    labels = _select_uc_vs_ctrl(metadata)
    if labels.empty:
        logging.warning('No samples matched the UC vs Ctrl selection criteria')
        return

    dataset = dataset.loc[labels.index]
    dataset_labels = dataset_labels.loc[labels.index]

    if config.normalize == 'robust_zscore':
        dataset = robust_zscore_normalization_per_dataset(dataset, dataset_labels)

    os.makedirs(config.output_dir, exist_ok=True)

    # ---- Per-dataset analysis ----
    all_gene_stats = []
    all_gsea = []
    dataset_summaries = []

    dataset_names = sorted(dataset_labels.unique())
    for idx, ds_name in enumerate(dataset_names, 1):
        ds_mask = dataset_labels == ds_name
        ds_data = dataset.loc[ds_mask]
        ds_labels = labels.loc[ds_mask]

        n_total = len(ds_data)
        n_uc = int((ds_labels == 'UC').sum())
        n_ctrl = int((ds_labels == 'Ctrl').sum())

        if n_total < config.min_samples or n_uc < config.min_per_group or n_ctrl < config.min_per_group:
            logging.info(
                'Skipping %s (%d/%d): n=%d UC=%d Ctrl=%d',
                ds_name, idx, len(dataset_names), n_total, n_uc, n_ctrl,
            )
            continue

        logging.info(
            'Processing %s (%d/%d): n=%d UC=%d Ctrl=%d',
            ds_name, idx, len(dataset_names), n_total, n_uc, n_ctrl,
        )

        gene_stats = _compute_gene_stats(ds_data, ds_labels)
        if gene_stats.empty:
            logging.warning('No gene stats for %s', ds_name)
            continue
        gene_stats['dataset'] = ds_name
        all_gene_stats.append(gene_stats)

        gsea_df = _run_gsea_per_dataset(gene_stats, ds_name, config.geneset_libraries, config)
        if not gsea_df.empty:
            all_gsea.append(gsea_df)

        dataset_summaries.append({
            'dataset': ds_name,
            'n_total': n_total,
            'n_uc': n_uc,
            'n_ctrl': n_ctrl,
            'n_genes_tested': len(gene_stats),
            'n_sig_genes_fdr05': int((gene_stats['q_value'] < 0.05).sum()),
        })

    if not all_gene_stats:
        logging.warning('No datasets produced results')
        return

    # ---- Save per-dataset results ----
    per_dataset_genes = pd.concat(all_gene_stats, ignore_index=True)
    per_dataset_genes.to_csv(os.path.join(config.output_dir, 'per_dataset_gene_stats.csv'), index=False)

    summary_df = pd.DataFrame(dataset_summaries)
    summary_df.to_csv(os.path.join(config.output_dir, 'dataset_summary.csv'), index=False)

    per_dataset_gsea = pd.DataFrame()
    if all_gsea:
        per_dataset_gsea = pd.concat(all_gsea, ignore_index=True)
        per_dataset_gsea.to_csv(os.path.join(config.output_dir, 'per_dataset_gsea.csv'), index=False)

    # ---- Gene meta-analysis ----
    logging.info('Running gene-level meta-analysis')
    gene_meta = _meta_analyse_genes(per_dataset_genes, config)
    gene_meta.to_csv(os.path.join(config.output_dir, 'gene_meta_all.csv'), index=False)

    gene_sig = _filter_genes(gene_meta, config)
    gene_sig.to_csv(os.path.join(config.output_dir, 'gene_meta_significant.csv'), index=False)
    logging.info('Genes: %d total, %d significant', len(gene_meta), len(gene_sig))

    # ---- Gene RRA ----
    logging.info('Running Robust Rank Aggregation')
    gene_rra = _compute_rra(per_dataset_genes, config)
    gene_rra.to_csv(os.path.join(config.output_dir, 'gene_rra.csv'), index=False)
    logging.info('RRA: %d genes, %d with FDR<0.05',
                 len(gene_rra), int((gene_rra['rra_q'] < 0.05).sum()))

    # ---- Pathway meta-analysis ----
    if not per_dataset_gsea.empty:
        logging.info('Running pathway-level meta-analysis')
        pathway_meta = _meta_analyse_pathways(per_dataset_gsea, config)
        pathway_meta.to_csv(os.path.join(config.output_dir, 'pathway_meta_all.csv'), index=False)

        pathway_sig = _filter_pathways(pathway_meta, config)
        pathway_sig.to_csv(os.path.join(config.output_dir, 'pathway_meta_significant.csv'), index=False)
        logging.info('Pathways: %d total, %d significant', len(pathway_meta), len(pathway_sig))

        lead_detail = _build_leading_edge_detail(pathway_sig, per_dataset_gsea, gene_meta)
        lead_detail.to_csv(os.path.join(config.output_dir, 'pathway_leading_edge_detail.csv'), index=False)

    logging.info('All results saved to %s', config.output_dir)


# ===================================================================
# CLI helper
# ===================================================================

def build_config(
    data_dir: str,
    output_dir: str,
    alpha: float = 0.05,
    min_samples: int = 5,
    min_per_group: int = 3,
    normalize: str = 'none',
    skip_gene_mapping: bool = False,
    min_datasets: int = 4,
    direction_threshold: float = 0.8,
    geneset_libraries: Optional[str] = None,
    gsea_permutations: int = 1000,
    gsea_min_size: int = 5,
    gsea_max_size: int = 1000,
    gsea_seed: int = 23,
) -> UcVsHcConfig:
    libraries = DEFAULT_GENESET_LIBRARIES
    if geneset_libraries:
        libraries = [lib.strip() for lib in geneset_libraries.split(',') if lib.strip()]

    return UcVsHcConfig(
        data_dir=data_dir,
        output_dir=output_dir,
        alpha=alpha,
        min_samples=min_samples,
        min_per_group=min_per_group,
        normalize=normalize,
        skip_gene_mapping=skip_gene_mapping,
        min_datasets=min_datasets,
        direction_threshold=direction_threshold,
        geneset_libraries=libraries,
        gsea_permutations=gsea_permutations,
        gsea_min_size=gsea_min_size,
        gsea_max_size=gsea_max_size,
        gsea_seed=gsea_seed,
    )