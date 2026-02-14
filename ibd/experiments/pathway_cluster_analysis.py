"""
Pathway cluster analysis — post-processing for UC vs HC meta-analysis.

Groups significant pathways by overlap of their significant leading-edge genes,
producing deduplicated cluster-level summaries for reporting.

Methodology
-----------
1. For each significant pathway, collect leading-edge genes that are also
   individually significant in the gene-level meta-analysis (FDR < gene_alpha).
2. Compute pairwise overlap coefficient between pathway gene sets:
       overlap(A, B) = |A ∩ B| / min(|A|, |B|)
   This captures subset relationships (e.g. a GO child term fully contained
   in its parent) better than Jaccard, which penalises size asymmetry.
3. Convert to distance (1 - overlap) and apply agglomerative hierarchical
   clustering with average linkage.
4. Automatically select the number of clusters by scanning distance thresholds
   and maximising the mean silhouette score on the distance matrix.
5. For each cluster, elect a representative pathway (highest |combined_nes|)
   and aggregate statistics.

Outputs
-------
pathway_clusters_summary.csv
    One row per cluster: representative pathway, aggregate NES, direction,
    gene union, and top genes ranked by |hedges_g|.
pathway_clusters_detail.csv
    All significant pathways with their cluster assignment, representative
    label, and similarity to the cluster representative.

Dependencies
------------
numpy, pandas, scipy, scikit-learn
"""

import logging
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import average, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score


# ===================================================================
# Configuration
# ===================================================================

@dataclass
class PathwayClusterConfig:
    """Configuration for pathway clustering."""

    # Input paths
    pathway_significant_path: str
    leading_edge_detail_path: str
    gene_significant_path: str
    output_dir: str

    # Gene significance threshold for building pathway gene sets
    gene_alpha: float = 0.05

    # Clustering
    min_genes_per_pathway: int = 3    # pathways with fewer significant genes become singletons
    distance_threshold: Optional[float] = None  # if None, auto-select via silhouette
    silhouette_min_threshold: float = 0.01  # lowest distance threshold to try
    silhouette_max_threshold: float = 0.99  # highest distance threshold to try
    silhouette_steps: int = 100             # number of thresholds to scan

    # Output
    top_genes_per_cluster: int = 20  # how many top genes (by |hedges_g|) to list in summary


# ===================================================================
# Gene set construction
# ===================================================================

def _build_pathway_gene_sets(
    pathways: pd.DataFrame,
    lead_detail: pd.DataFrame,
    sig_genes: pd.DataFrame,
    gene_alpha: float,
    min_genes: int,
) -> Dict[Tuple[str, str], Set[str]]:
    """
    For each pathway, collect significant leading-edge genes.

    A gene is included if it appears in the pathway's leading edge AND
    is significant in the gene-level meta-analysis (q_value < gene_alpha).

    Returns dict mapping (library, term) -> set of gene symbols.
    """
    # Set of significant genes
    sig_gene_set = set(sig_genes.loc[sig_genes['q_value'] < gene_alpha, 'gene'])
    logging.info('Significant gene pool: %d genes (FDR < %g)', len(sig_gene_set), gene_alpha)

    # Build per-pathway gene sets from leading edge detail
    pathway_genes: Dict[Tuple[str, str], Set[str]] = {}

    for (lib, term), grp in lead_detail.groupby(['library', 'term']):
        key = (lib, term)
        # Only keep leading-edge genes that are themselves significant
        genes = set(grp['gene']) & sig_gene_set
        if len(genes) >= min_genes:
            pathway_genes[key] = genes

    # Also handle pathways that might not be in lead_detail (edge case)
    pathway_keys = set(zip(pathways['library'], pathways['term']))
    missing = pathway_keys - set(pathway_genes.keys())
    if missing:
        logging.info('%d pathways had fewer than %d significant leading-edge genes (will be singletons)',
                     len(missing), min_genes)

    return pathway_genes


# ===================================================================
# Overlap coefficient and distance matrix
# ===================================================================

def _overlap_coefficient(set_a: Set[str], set_b: Set[str]) -> float:
    """
    Overlap coefficient: |A ∩ B| / min(|A|, |B|).

    Returns 1.0 when the smaller set is entirely contained in the larger.
    Returns 0.0 when the sets are disjoint or either is empty.
    """
    if not set_a or not set_b:
        return 0.0
    return len(set_a & set_b) / min(len(set_a), len(set_b))


def _compute_distance_matrix(
    pathway_keys: List[Tuple[str, str]],
    pathway_genes: Dict[Tuple[str, str], Set[str]],
) -> np.ndarray:
    """
    Compute pairwise distance matrix (1 - overlap_coefficient).

    Pathways not in pathway_genes (too few significant genes) get
    distance 1.0 to everything (they will become singletons).
    """
    n = len(pathway_keys)
    dist = np.ones((n, n))
    np.fill_diagonal(dist, 0.0)

    for i in range(n):
        genes_i = pathway_genes.get(pathway_keys[i])
        if genes_i is None:
            continue
        for j in range(i + 1, n):
            genes_j = pathway_genes.get(pathway_keys[j])
            if genes_j is None:
                continue
            overlap = _overlap_coefficient(genes_i, genes_j)
            dist[i, j] = 1.0 - overlap
            dist[j, i] = 1.0 - overlap

    return dist


# ===================================================================
# Clustering
# ===================================================================

def _auto_select_threshold(
    dist_matrix: np.ndarray,
    min_t: float,
    max_t: float,
    steps: int,
    output_dir: str,
) -> Tuple[float, float]:
    """
    Scan distance thresholds and return the one with highest silhouette score.

    Returns (best_threshold, best_silhouette_score).
    """
    condensed = squareform(dist_matrix, checks=False)
    Z = average(condensed)
    n = dist_matrix.shape[0]

    best_t, best_score = 0.5, -1.0  # default fallback
    thresholds = np.linspace(min_t, max_t, steps)
    scores = []

    for t in thresholds:
        labels = fcluster(Z, t=t, criterion='distance')
        n_clusters = len(set(labels))

        # Silhouette needs at least 2 clusters and at most n-1
        if n_clusters < 2 or n_clusters >= n:
            continue

        # Skip if almost all singletons (silhouette is meaningless)
        sizes = np.bincount(labels)
        non_singleton = np.sum(sizes[sizes > 0] > 1)
        if non_singleton < 2:
            continue

        try:
            score = silhouette_score(dist_matrix, labels, metric='precomputed')
            scores.append(score)
            if score > best_score:
                best_score = score
                best_t = t
        except ValueError:
            continue

    # plot silhouette scores vs thresholds for diagnostics
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.plot(thresholds, scores)
    plt.xlabel('Distance Threshold')
    plt.ylabel('Silhouette Score')
    plt.title('Silhouette Scores vs Distance Thresholds')
    plt.savefig(os.path.join(output_dir, 'silhouette_threshold_scan.png'))
    plt.close()

    return best_t, best_score


def _cluster_pathways(
    dist_matrix: np.ndarray,
    threshold: float,
) -> np.ndarray:
    """
    Hierarchical clustering with average linkage, cut at given distance threshold.

    Returns cluster labels (1-indexed).
    """
    condensed = squareform(dist_matrix, checks=False)
    Z = average(condensed)
    return fcluster(Z, t=threshold, criterion='distance')


# ===================================================================
# Cluster aggregation
# ===================================================================

def _aggregate_clusters(
    pathways: pd.DataFrame,
    cluster_labels: np.ndarray,
    pathway_keys: List[Tuple[str, str]],
    pathway_genes: Dict[Tuple[str, str], Set[str]],
    gene_meta: pd.DataFrame,
    top_n_genes: int,
) -> pd.DataFrame:
    """
    Produce one-row-per-cluster summary with aggregate statistics.
    """
    # Gene lookup for ranking
    gene_lookup = {}
    if gene_meta is not None and not gene_meta.empty:
        gene_lookup = gene_meta.set_index('gene')[['hedges_g', 'q_value', 'direction']].to_dict('index')

    # Add cluster labels to pathways
    pw = pathways.copy()
    pw['_key'] = list(zip(pw['library'], pw['term']))
    key_to_cluster = dict(zip(pathway_keys, cluster_labels))
    pw['cluster'] = pw['_key'].map(key_to_cluster)

    records = []
    for cluster_id, grp in pw.groupby('cluster'):
        n_pathways = len(grp)

        # Representative: highest |combined_nes|
        rep_idx = grp['combined_nes'].abs().idxmax()
        rep = grp.loc[rep_idx]

        # Direction: majority vote
        n_up = int((grp['direction'] == 'up').sum())
        n_down = int((grp['direction'] == 'down').sum())
        direction = 'up' if n_up >= n_down else 'down'

        # NES statistics
        nes_vals = grp['combined_nes'].values
        abs_nes = np.abs(nes_vals)

        # Gene union across all pathways in cluster
        gene_union: Set[str] = set()
        for _, row in grp.iterrows():
            key = (row['library'], row['term'])
            genes = pathway_genes.get(key, set())
            gene_union |= genes

        # Rank genes by |hedges_g|
        gene_stats = []
        for g in gene_union:
            info = gene_lookup.get(g)
            if info:
                gene_stats.append((g, info['hedges_g'], info['q_value'], info['direction']))

        gene_stats.sort(key=lambda x: abs(x[1]), reverse=True)
        top_genes = gene_stats[:top_n_genes]

        records.append({
            'cluster': int(cluster_id),
            'n_pathways': n_pathways,
            'direction': direction,
            'n_up': n_up,
            'n_down': n_down,
            'representative_library': rep['library'],
            'representative_term': rep['term'],
            'representative_nes': float(rep['combined_nes']),
            'representative_p_value': float(rep['p_value']),
            'representative_q_value': float(rep['q_value']),
            'representative_I2': float(rep['I2']),
            'mean_abs_nes': float(np.mean(abs_nes)),
            'median_abs_nes': float(np.median(abs_nes)),
            'min_abs_nes': float(np.min(abs_nes)),
            'max_abs_nes': float(np.max(abs_nes)),
            'n_significant_genes': len(gene_union),
            'top_genes': '; '.join(
                f"{g} ({hg:+.2f})" for g, hg, _, _ in top_genes
            ),
            'all_significant_genes': ';'.join(sorted(gene_union)),
            'all_terms': ' | '.join(
                f"[{row['library']}] {row['term']}"
                for _, row in grp.sort_values('combined_nes', key=abs, ascending=False).iterrows()
            ),
        })

    df = pd.DataFrame(records)
    df = df.sort_values('mean_abs_nes', ascending=False).reset_index(drop=True)

    # Re-number clusters by mean |NES| rank (1 = strongest)
    df['cluster'] = range(1, len(df) + 1)

    return df


def _build_detail_table(
    pathways: pd.DataFrame,
    cluster_labels: np.ndarray,
    pathway_keys: List[Tuple[str, str]],
    cluster_summary: pd.DataFrame,
    dist_matrix: np.ndarray,
) -> pd.DataFrame:
    """
    Original pathway table augmented with cluster assignment and metadata.
    """
    pw = pathways.copy()

    # Map old cluster labels to new sorted cluster IDs
    key_to_old_cluster = dict(zip(pathway_keys, cluster_labels))

    # Build old_cluster -> new_cluster mapping from summary
    # Summary was re-numbered, so we need to trace back
    # Re-do: assign clusters based on representative membership
    old_to_new = {}
    for new_id, row in cluster_summary.iterrows():
        new_cluster = row['cluster']
        rep_key = (row['representative_library'], row['representative_term'])
        old_cluster = key_to_old_cluster[rep_key]
        old_to_new[old_cluster] = new_cluster

    pw['_key'] = list(zip(pw['library'], pw['term']))
    pw['cluster'] = pw['_key'].map(lambda k: old_to_new.get(key_to_old_cluster.get(k, -1), -1))

    # Add representative info
    cluster_rep = cluster_summary.set_index('cluster')[['representative_library', 'representative_term']].to_dict('index')
    pw['cluster_representative'] = pw['cluster'].map(
        lambda c: cluster_rep.get(c, {}).get('representative_term', '')
    )
    pw['is_representative'] = pw.apply(
        lambda r: (r['library'] == cluster_rep.get(r['cluster'], {}).get('representative_library', ''))
                  and (r['term'] == cluster_rep.get(r['cluster'], {}).get('representative_term', '')),
        axis=1,
    )

    # Similarity to cluster representative
    key_to_idx = {k: i for i, k in enumerate(pathway_keys)}
    rep_idx_map = {}
    for c, info in cluster_rep.items():
        rep_key = (info['representative_library'], info['representative_term'])
        if rep_key in key_to_idx:
            rep_idx_map[c] = key_to_idx[rep_key]

    def _sim_to_rep(row):
        cluster = row['cluster']
        key = (row['library'], row['term'])
        if key not in key_to_idx or cluster not in rep_idx_map:
            return np.nan
        i = key_to_idx[key]
        j = rep_idx_map[cluster]
        if i == j:
            return 1.0
        return 1.0 - dist_matrix[i, j]

    pw['similarity_to_representative'] = pw.apply(_sim_to_rep, axis=1)

    pw = pw.drop(columns=['_key'])
    pw = pw.sort_values(['cluster', 'combined_nes'], key=lambda s: s.abs() if s.name == 'combined_nes' else s,
                        ascending=[True, False]).reset_index(drop=True)

    return pw


# ===================================================================
# Main entry point
# ===================================================================

def run_pathway_clustering(config: PathwayClusterConfig) -> None:
    """Run pathway clustering analysis."""

    logging.info('Loading inputs')
    pathways = pd.read_csv(config.pathway_significant_path)
    lead_detail = pd.read_csv(config.leading_edge_detail_path)
    gene_meta = pd.read_csv(config.gene_significant_path)

    logging.info('Input: %d significant pathways, %d leading-edge records, %d significant genes',
                 len(pathways), len(lead_detail), len(gene_meta))

    # ---- Build gene sets ----
    pathway_genes = _build_pathway_gene_sets(
        pathways, lead_detail, gene_meta,
        gene_alpha=config.gene_alpha,
        min_genes=config.min_genes_per_pathway,
    )
    logging.info('Pathways with sufficient significant genes: %d / %d',
                 len(pathway_genes), len(pathways))

    # ---- Compute distance matrix ----
    pathway_keys = list(zip(pathways['library'], pathways['term']))

    logging.info('Computing pairwise overlap distances')
    dist_matrix = _compute_distance_matrix(pathway_keys, pathway_genes)

    # ---- Cluster ----
    if config.distance_threshold is not None:
        threshold = config.distance_threshold
        logging.info('Using user-specified distance threshold: %.2f', threshold)
    else:
        logging.info('Auto-selecting distance threshold via silhouette analysis')
        threshold, score = _auto_select_threshold(
            dist_matrix,
            config.silhouette_min_threshold,
            config.silhouette_max_threshold,
            config.silhouette_steps,
            config.output_dir,
        )
        logging.info('Selected threshold: %.3f (silhouette score: %.3f)', threshold, score)

    cluster_labels = _cluster_pathways(dist_matrix, threshold)
    n_clusters = len(set(cluster_labels))
    sizes = np.bincount(cluster_labels)
    n_singletons = int(np.sum(sizes[sizes > 0] == 1))

    logging.info('Clusters: %d total (%d singletons, %d multi-pathway)',
                 n_clusters, n_singletons, n_clusters - n_singletons)

    # ---- Aggregate ----
    logging.info('Aggregating cluster statistics')
    summary = _aggregate_clusters(
        pathways, cluster_labels, pathway_keys, pathway_genes,
        gene_meta, config.top_genes_per_cluster,
    )

    detail = _build_detail_table(
        pathways, cluster_labels, pathway_keys, summary, dist_matrix,
    )

    # ---- Save ----
    os.makedirs(config.output_dir, exist_ok=True)

    summary_path = os.path.join(config.output_dir, 'pathway_clusters_summary.csv')
    detail_path = os.path.join(config.output_dir, 'pathway_clusters_detail.csv')

    summary.to_csv(summary_path, index=False)
    detail.to_csv(detail_path, index=False)

    logging.info('Saved %d clusters to %s', len(summary), summary_path)
    logging.info('Saved %d pathway assignments to %s', len(detail), detail_path)

    # ---- Log top clusters ----
    logging.info('Top 10 clusters by mean |NES|:')
    for _, row in summary.head(10).iterrows():
        logging.info(
            '  Cluster %d: %d pathways, direction=%s, mean|NES|=%.2f, '
            'representative=[%s] %s, %d genes',
            row['cluster'], row['n_pathways'], row['direction'],
            row['mean_abs_nes'], row['representative_library'],
            row['representative_term'], row['n_significant_genes'],
        )


# ===================================================================
# CLI helper
# ===================================================================

def build_cluster_config(
    pathway_significant_path: str,
    leading_edge_detail_path: str,
    gene_significant_path: str,
    output_dir: str,
    gene_alpha: float = 0.05,
    min_genes_per_pathway: int = 3,
    distance_threshold: Optional[float] = None,
    top_genes_per_cluster: int = 20,
) -> PathwayClusterConfig:
    return PathwayClusterConfig(
        pathway_significant_path=pathway_significant_path,
        leading_edge_detail_path=leading_edge_detail_path,
        gene_significant_path=gene_significant_path,
        output_dir=output_dir,
        gene_alpha=gene_alpha,
        min_genes_per_pathway=min_genes_per_pathway,
        distance_threshold=distance_threshold,
        top_genes_per_cluster=top_genes_per_cluster,
    )