import pandas as pd
import numpy as np
import igraph as ig
import leidenalg


THRESHOLD = 0.3

# 1. Load Data
df = pd.read_csv('results/uc_vs_hc/pathway_meta_significant.csv')

# 2. Compute Overlap Coefficient Matrix
gene_sets = [set(x.split(';')) for x in df['lead_genes_core']]
n = len(gene_sets)
similarity_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(i + 1, n):
        set_a, set_b = gene_sets[i], gene_sets[j]
        if not set_a or not set_b: continue
        sim = len(set_a.intersection(set_b)) / min(len(set_a), len(set_b))
        if sim > THRESHOLD:
            similarity_matrix[i, j] = sim
            similarity_matrix[j, i] = sim

# 3. Build Graph (igraph)
rows, cols = np.where(np.triu(similarity_matrix, k=1) > 0)
edges = list(zip(rows, cols))
weights = similarity_matrix[rows, cols]

G = ig.Graph(n, edges=edges)
G.es['weight'] = weights

# 4. Run Leiden Algorithm
partition = leidenalg.find_partition(G, leidenalg.ModularityVertexPartition, weights=weights)

print(f"Found {len(partition)} clusters.")

# 5. Assign Labels to DataFrame
df['cluster_id'] = -1
for cid, members in enumerate(partition):
    df.loc[members, 'cluster_id'] = cid

# 6. Find Representative (Highest absolute NES) for Each Cluster
df['abs_combined_nes'] = df['combined_nes'].abs()
reps = df.loc[df.groupby('cluster_id')['abs_combined_nes'].idxmax()][['cluster_id', 'term']]
reps.rename(columns={'term': 'rep_term'}, inplace=True)
df = df.merge(reps, on='cluster_id', how='left')

df.to_csv('results/uc_vs_hc/pathway_clusters_leiden.csv', index=False)

print(f"Total Clusters: {len(partition)}")