"""
Generate pathway overlap network graph data (JSON) from meta-analysis CSVs.

Usage:
    python ibd/scripts/build_pathway_graph.py \
        --pathways results/new/uc_vs_hc/pathway_meta_significant.csv \
        --genes results/new/uc_vs_hc/gene_meta_significant.csv \
        --output results/new/uc_vs_hc/pathway_graph.json \
        --gene-alpha 0.05

The output JSON can be loaded into the pathway_network_viewer.html.
"""

import argparse
import json
import logging
import sys
from typing import Dict, Optional, Set, Tuple

import numpy as np
import pandas as pd


def build_pathway_graph(
    pathways: pd.DataFrame,
    gene_meta: pd.DataFrame,
    gene_alpha: float = 0.05,
    min_genes: int = 2,
) -> dict:
    """
    Build a pathway overlap network graph.

    Parameters
    ----------
    pathways : DataFrame
        Significant pathways with columns: library, term, combined_nes, q_value, etc.
    gene_meta : DataFrame
        Gene-level meta-analysis with columns: gene, q_value.
    gene_alpha : float
        FDR threshold for gene significance.
    min_genes : int
        Minimum significant leading-edge genes for a pathway to be included.

    Returns
    -------
    dict with 'nodes' and 'edges' lists, ready for JSON serialization.
    """
    sig_genes = set(gene_meta.loc[gene_meta['q_value'] < gene_alpha, 'gene'])
    logging.info('Significant gene pool: %d genes (FDR < %g)', len(sig_genes), gene_alpha)

    pathway_genes: Dict[Tuple[str, str], Set[str]] = {}
    gene_col = None
    for col in ('lead_genes_all', 'lead_genes_core'):
        if col in pathways.columns:
            gene_col = col
            break

    if gene_col is not None:
        logging.info('Building gene sets from pathways[%s]', gene_col)
        for _, r in pathways.iterrows():
            key = (r['library'], r['term'])
            raw = r[gene_col]
            if pd.isna(raw) or not str(raw).strip():
                continue
            all_genes = set(str(raw).split(';'))
            genes = all_genes & sig_genes
            if len(genes) >= min_genes:
                pathway_genes[key] = genes
    else:
        logging.warning('No gene set source available — all nodes will have 0 genes')

    n_with_genes = sum(1 for _, r in pathways.iterrows()
                       if (r['library'], r['term']) in pathway_genes)
    logging.info('Pathways with gene sets: %d / %d', n_with_genes, len(pathways))

    # Build nodes
    nodes = []
    for i, (_, r) in enumerate(pathways.iterrows()):
        key = (r['library'], r['term'])
        genes = pathway_genes.get(key, set())

        # Shorten label
        term_short = r['term']
        # Strip Reactome IDs for display
        if ' R-HSA-' in term_short:
            term_short = term_short.split(' R-HSA-')[0]
        if len(term_short) > 50:
            term_short = term_short[:47] + '...'

        node = {
            'id': i,
            'label': term_short,
            'fullTerm': r['term'],
            'library': r['library'],
            'nes': round(float(r['combined_nes']), 2),
            'absNes': round(float(abs(r['combined_nes'])), 2),
            'direction': r['direction'],
            'nGenes': len(genes),
            'qValue': float(r['q_value']),
        }

        # Optional columns
        if 'datasets' in r.index:
            node['datasets'] = int(r['datasets'])
        if 'I2' in r.index:
            node['I2'] = round(float(r['I2']), 2)

        nodes.append(node)

    # Build edges (overlap coefficient)
    keys = [(r['library'], r['term']) for _, r in pathways.iterrows()]
    edges = []
    for i in range(len(keys)):
        gi = pathway_genes.get(keys[i], set())
        if not gi:
            continue
        for j in range(i + 1, len(keys)):
            gj = pathway_genes.get(keys[j], set())
            if not gj:
                continue
            overlap = len(gi & gj) / min(len(gi), len(gj))
            edges.append({
                'source': i,
                'target': j,
                'weight': round(overlap, 3),
                'shared': len(gi & gj),
            })

    logging.info('Graph: %d nodes, %d edges', len(nodes), len(edges))

    return {'nodes': nodes, 'edges': edges}


def main():
    parser = argparse.ArgumentParser(description='Build pathway overlap network JSON')
    parser.add_argument('--pathways', required=True, help='Pathway meta-analysis CSV')
    parser.add_argument('--genes', required=True, help='Gene meta-analysis CSV')
    parser.add_argument('--output', default='pathway_graph.json', help='Output JSON path')
    parser.add_argument('--gene-alpha', type=float, default=0.05, help='Gene FDR threshold')
    parser.add_argument('--min-genes', type=int, default=2, help='Min significant genes per pathway')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    pathways = pd.read_csv(args.pathways)
    gene_meta = pd.read_csv(args.genes)

    graph = build_pathway_graph(pathways, gene_meta, gene_alpha=args.gene_alpha, min_genes=args.min_genes)

    with open(args.output, 'w') as f:
        json.dump(graph, f, separators=(',', ':'))

    logging.info('Saved to %s (%d bytes)', args.output, len(json.dumps(graph)))


if __name__ == '__main__':
    main()