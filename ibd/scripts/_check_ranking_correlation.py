import pandas as pd
from scipy.stats import spearmanr

meta = pd.read_csv('results/paper/uc_vs_hc/gene_meta_significant.csv', encoding='latin1')
rra = pd.read_csv('results/paper/uc_vs_hc/gene_rra.csv', encoding='latin1')

merged = meta.merge(rra, on='gene', how='inner')
merged = merged[merged['datasets_x'] >=7]
rank_g = merged['hedges_g'].abs().rank(ascending=False)
rank_rho = merged['rra_rho'].rank(ascending=True)
r, p = spearmanr(rank_g, rank_rho)
print(f"Spearman r={r:.3f}, p={p:.2e}")

# diff = abs(rank_g - rank_rho)
# merged['rank_g'] = rank_g
# merged['rank_rho'] = rank_rho
# merged['diff'] = diff

# print(f"Top 23 genes with largest rank difference that rank in top 50:")
# top_diff = merged[(merged['rank_g'] <= 50) | (merged['rank_rho'] <= 50)].nlargest(23, 'diff')
# print(top_diff[['gene', 'rank_g', 'rank_rho', 'diff']])