import pandas as pd

from kneefinder import KneeFinder


DIR = 'results/new/cd_vs_hc/'

df = pd.read_csv(f'{DIR}pathway_meta_all.csv')

df['abs_combined_nes'] = df['combined_nes'].abs()
df = df[(df['abs_combined_nes'] > 1) & (df['q_value'] < 0.05)]
sorted_nes = df['abs_combined_nes'].sort_values(ascending=False)

kf = KneeFinder(range(len(sorted_nes)), sorted_nes)

knee_x, knee_y = kf.find_knee()

df_filtered = df[df['abs_combined_nes'] >= knee_y]
df_filtered.to_csv(f'{DIR}pathway_meta_significant_knee.csv', index=False)