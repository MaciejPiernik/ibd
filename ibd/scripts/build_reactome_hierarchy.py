"""
Build a Reactome pathway hierarchy tree with enrichment overlay.

Downloads the Reactome hierarchy files, maps enriched pathways onto the tree,
and outputs a JSON suitable for the collapsible tree viewer.

Usage:
    python ibd/scripts/build_reactome_hierarchy.py \
        --pathways results/paper/uc_vs_hc/pathway_meta_significant.csv \
        --output results/paper/uc_vs_hc/reactome_tree.json

Dependencies: pandas, requests (or urllib)
"""

import argparse
import json
import logging
import os
import re
import urllib.request
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

# ── Reactome download URLs ──────────────────────────────────────────
# "current" redirects to the latest release
REACTOME_PATHWAYS_URL = 'https://reactome.org/download/current/ReactomePathways.txt'
REACTOME_RELATIONS_URL = 'https://reactome.org/download/current/ReactomePathwaysRelation.txt'
CACHE_DIR = os.path.join(os.path.expanduser('~'), '.cache', 'reactome')


def _download_or_cache(url: str, filename: str) -> str:
    """Download file if not cached, return local path."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, filename)
    if os.path.exists(path):
        logging.info('Using cached %s', path)
        return path
    logging.info('Downloading %s ...', url)
    urllib.request.urlretrieve(url, path)
    logging.info('Saved to %s', path)
    return path


# ── Parse Reactome files ────────────────────────────────────────────

def _load_pathway_names(path: str) -> Dict[str, Tuple[str, str]]:
    """
    ReactomePathways.txt: tab-separated (stId, name, species).
    Returns dict: stId -> (name, species).
    """
    names = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                names[parts[0]] = (parts[1], parts[2])
    return names


def _load_relations(path: str) -> List[Tuple[str, str]]:
    """
    ReactomePathwaysRelation.txt: tab-separated (parent_stId, child_stId).
    Returns list of (parent, child) tuples.
    """
    relations = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                relations.append((parts[0], parts[1]))
    return relations


# ── Build tree ──────────────────────────────────────────────────────

def _extract_reactome_id(term: str) -> Optional[str]:
    """Extract R-HSA-XXXXXXX from a term string like 'Pathway Name R-HSA-1234567'."""
    match = re.search(r'R-HSA-\d+', term)
    return match.group(0) if match else None


def build_hierarchy_tree(
    pathways_df: pd.DataFrame,
    species: str = 'Homo sapiens',
) -> dict:
    """
    Build a Reactome hierarchy tree with enrichment data overlaid.

    Parameters
    ----------
    pathways_df : DataFrame
        Enriched pathways with columns: term, combined_nes, direction, q_value, etc.
    species : str
        Filter Reactome pathways to this species.

    Returns
    -------
    dict : Tree structure ready for JSON serialization.
    """
    # Download / cache Reactome files
    names_path = _download_or_cache(REACTOME_PATHWAYS_URL, 'ReactomePathways.txt')
    rels_path = _download_or_cache(REACTOME_RELATIONS_URL, 'ReactomePathwaysRelation.txt')

    all_pathways = _load_pathway_names(names_path)
    all_relations = _load_relations(rels_path)

    # Filter to species
    species_ids = {sid for sid, (name, sp) in all_pathways.items() if sp == species}
    logging.info('Reactome %s pathways: %d', species, len(species_ids))

    # Build parent -> children and child -> parents maps
    children_map: Dict[str, List[str]] = defaultdict(list)
    parent_map: Dict[str, List[str]] = defaultdict(list)

    for parent, child in all_relations:
        if parent in species_ids and child in species_ids:
            children_map[parent].append(child)
            parent_map[child].append(parent)

    # Find root nodes (no parents)
    all_children = set()
    for _, child in all_relations:
        if child in species_ids:
            all_children.add(child)
    roots = species_ids - all_children
    logging.info('Root pathways (top-level): %d', len(roots))

    # Build enrichment lookup from input data
    enriched: Dict[str, dict] = {}
    for _, row in pathways_df.iterrows():
        rid = _extract_reactome_id(str(row['term']))
        if rid and rid in species_ids:
            enriched[rid] = {
                'nes': round(float(row['combined_nes']), 2),
                'direction': row['direction'],
                'q_value': float(row['q_value']) if 'q_value' in row.index else None,
                'I2': round(float(row['I2']), 2) if 'I2' in row.index else None,
                'datasets': int(row['datasets']) if 'datasets' in row.index else None,
            }

    logging.info('Enriched pathways matched to Reactome IDs: %d / %d',
                 len(enriched), len(pathways_df))

    # Determine which nodes are needed: enriched pathways + all ancestors to roots
    needed: Set[str] = set(enriched.keys())
    queue = list(needed)
    while queue:
        node = queue.pop()
        for parent in parent_map.get(node, []):
            if parent not in needed:
                needed.add(parent)
                queue.append(parent)

    # Also include direct children of needed nodes so we see the full context
    # (shows siblings of enriched pathways at each level)
    expanded_needed = set(needed)
    for node in needed:
        for child in children_map.get(node, []):
            expanded_needed.add(child)

    logging.info('Tree nodes (enriched + ancestors + siblings): %d', len(expanded_needed))

    # Recursive tree builder.
    # Reactome is a DAG (pathways can have multiple parents), so nodes may
    # appear under more than one branch. We allow this — the only guard is
    # a per-path `ancestors` set that prevents infinite loops from cycles.
    # Counts are computed bottom-up from actually-rendered children so they
    # always match what's visible in the viewer.

    def _build_node(node_id: str, depth: int = 0, ancestors: Set[str] = None) -> Optional[dict]:
        if ancestors is None:
            ancestors = set()
        if node_id in ancestors:
            return None  # cycle guard only
        path_ancestors = ancestors | {node_id}

        name = all_pathways.get(node_id, (node_id, ''))[0]

        # Build children recursively
        child_ids = [c for c in children_map.get(node_id, []) if c in expanded_needed]
        child_nodes = []
        for cid in sorted(child_ids, key=lambda x: all_pathways.get(x, (x, ''))[0]):
            cn = _build_node(cid, depth + 1, path_ancestors)
            if cn is not None:
                child_nodes.append(cn)

        # Compute counts bottom-up from what was actually rendered
        is_self_enriched = 1 if node_id in enriched else 0
        if child_nodes:
            n_enriched = is_self_enriched + sum(c['enriched'] for c in child_nodes)
            n_total = 1 + sum(c['total'] for c in child_nodes)
        else:
            n_enriched = is_self_enriched
            n_total = 1

        node = {
            'id': node_id,
            'name': name,
            'enriched': n_enriched,
            'total': n_total,
            'depth': depth,
        }

        if node_id in enriched:
            node['data'] = enriched[node_id]

        if child_nodes:
            node['children'] = child_nodes

        return node

    # Build from roots — no global visited, each root gets its own traversal
    tree_roots = []
    for rid in sorted(roots, key=lambda x: all_pathways.get(x, (x, ''))[0]):
        if rid not in expanded_needed:
            continue
        node = _build_node(rid, depth=0)
        if node is not None:
            tree_roots.append(node)

    tree = {
        'name': 'Reactome',
        'children': tree_roots,
        'enriched': sum(r.get('enriched', 0) for r in tree_roots),
        'total': sum(r.get('total', 0) for r in tree_roots),
    }

    logging.info('Tree built: %d top-level categories', len(tree_roots))

    # Report DAG duplication
    def _count_tree_nodes(n):
        c = 1
        for ch in n.get('children', []):
            c += _count_tree_nodes(ch)
        return c
    total_tree_nodes = sum(_count_tree_nodes(r) for r in tree_roots)
    logging.info('Tree size: %d rendered nodes (%d unique Reactome IDs — '
                 'difference is DAG duplication)',
                 total_tree_nodes, len(expanded_needed))

    return tree


def main():
    parser = argparse.ArgumentParser(description='Build Reactome hierarchy with enrichment overlay')
    parser.add_argument('--pathways', required=True,
                        help='Pathway meta-analysis CSV (must have "term" column with R-HSA IDs)')
    parser.add_argument('--output', default='reactome_tree.json', help='Output JSON path')
    parser.add_argument('--species', default='Homo sapiens', help='Species filter')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    pathways_df = pd.read_csv(args.pathways)
    tree = build_hierarchy_tree(pathways_df, species=args.species)

    with open(args.output, 'w') as f:
        json.dump(tree, f, indent=1)

    logging.info('Saved to %s', args.output)


if __name__ == '__main__':
    main()