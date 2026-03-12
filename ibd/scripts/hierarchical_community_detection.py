from operator import index
import pandas as pd
import numpy as np
import networkx as nx

MAX_DEPTH = 1  # Maximum depth of the hierarchy
MIN_SIZE = 6   # Minimum cluster size to consider for further splitting
THRESHOLD = 0.33  # Minimum overlap coefficient to consider an edge
DIR = 'results/new/cd_vs_hc/'  # Output directory for results

# 1. Load Data
df = pd.read_csv(f'{DIR}pathway_meta_significant_knee.csv')

# 2. Build the Graph
gene_sets = [set(x.split(';')) for x in df['lead_genes_core']]
n = len(gene_sets)
similarity_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(i + 1, n):
        set_a, set_b = gene_sets[i], gene_sets[j]
        if not set_a or not set_b: continue
        
        # Overlap Coefficient
        sim = len(set_a.intersection(set_b)) / min(len(set_a), len(set_b))
        
        if sim > THRESHOLD:
            similarity_matrix[i, j] = sim
            similarity_matrix[j, i] = sim

G = nx.Graph()
G.add_nodes_from(range(n))
rows, cols = np.where(np.triu(similarity_matrix, k=1) > 0)
G.add_weighted_edges_from(zip(rows, cols, similarity_matrix[rows, cols]))

# 3. Recursive Clustering Function
def get_representative(node_indices, df):
    """Selects the term with the highest Impact Score: |NES| * -log10(p)"""
    if not node_indices: return "Empty"
    subset = df.iloc[node_indices].copy()
    subset['impact_score'] = subset['combined_nes'].abs() * -np.log10(subset['p_value'] + 1e-100)
    return subset.loc[subset['impact_score'].idxmax(), 'term']

def recursive_louvain(graph, node_indices, current_id="1", min_size=6, max_depth=None, depth=1):
    """
    Recursively splits clusters until they are smaller than min_size 
    or cannot be split further by modularity optimization.
    """
    # Base case: Cluster too small or max depth reached
    if len(node_indices) < min_size or (max_depth is not None and depth >= max_depth):
        return [{
            'Cluster_ID': current_id,
            'Node_Indices': node_indices,
            'Size': len(node_indices),
            'Leaf': True
        }]
    
    # Extract subgraph for this cluster
    subG = graph.subgraph(node_indices).copy()
    
    # Run Louvain on this isolated subgraph
    try:
        communities = nx.community.louvain_communities(subG, weight='weight', resolution=1.0, seed=23)
    except:
        communities = [set(node_indices)] # Fail-safe

    # Base case: No split found (returns itself)
    if len(communities) <= 1:
        return [{
            'Cluster_ID': current_id,
            'Node_Indices': node_indices,
            'Size': len(node_indices),
            'Leaf': True
        }]
    
    # Sort communities by size (largest first)
    communities = sorted(communities, key=len, reverse=True)
    
    results = []
    for i, comm in enumerate(communities):
        # Generate ID (e.g., "1.1", "1.2")
        new_id = f"{current_id}.{i+1}"
        # Recurse
        results.extend(recursive_louvain(graph, list(comm), new_id, min_size, max_depth, depth + 1))
        
    return results

# 4. Run Hierarchy
# Start with global clusters (Level 1)
global_communities = nx.community.louvain_communities(G, weight='weight', resolution=1.0, seed=23)
global_communities = sorted(global_communities, key=len, reverse=True)

final_clusters = []
for i, comm in enumerate(global_communities):
    cid = str(i + 1)
    # Run recursion on each global cluster
    final_clusters.extend(recursive_louvain(G, list(comm), current_id=cid, min_size=MIN_SIZE, max_depth=MAX_DEPTH, depth=1))

# 5. Export Results
results_data = []
for c in final_clusters:
    rep_term = get_representative(c['Node_Indices'], df)
    mean_nes = df.iloc[c['Node_Indices']]['combined_nes'].mean()
    
    results_data.append({
        'Hierarchy_ID': c['Cluster_ID'],
        'Representative': rep_term,
        'Size': c['Size'],
        'Mean_NES': mean_nes,
        'Pathways': "; ".join(df.iloc[c['Node_Indices']]['term'].tolist())
    })

results_df = pd.DataFrame(results_data)
results_df.to_csv(f'{DIR}hierarchical_pathway_clusters.csv', index=False)

# Add cluster labels back to original df for reference and save
cluster_labels = {}
for c in final_clusters:
    for idx in c['Node_Indices']:
        cluster_labels[idx] = c['Cluster_ID']
        
df['Hierarchy_ID'] = df.index.map(cluster_labels)

df.to_csv(f'{DIR}pathway_meta_significant_knee_with_clusters.csv', index=False)

print(f"Recursion Complete. Found {len(results_df)} granular clusters.")
print(results_df[['Hierarchy_ID', 'Size', 'Representative', 'Mean_NES']].head(10))