import pandas as pd

from kneefinder import KneeFinder


DIRS = ['results/paper/uc_vs_cd/', 'results/paper/uc_vs_hc/', 'results/paper/uninf_vs_hc/', 'results/paper/cd_vs_hc/']

for dir in DIRS:
    df = pd.read_csv(f'{dir}pathway_meta_all.csv')

    min_datasets = df['datasets'].max() / 2
    df['abs_combined_nes'] = df['combined_nes'].abs()
    df = df[(df['q_value'] < 0.05)
            & (df['direction_ratio'] >= 0.8)
            & (df['datasets'] >= min_datasets)]
    sorted_nes = df['abs_combined_nes'].sort_values(ascending=False)

    kf = KneeFinder(range(len(sorted_nes)), sorted_nes)

    knee_x, knee_y = kf.find_knee()

    df_filtered = df[df['abs_combined_nes'] >= knee_y]
    df_filtered.to_csv(f'{dir}pathway_meta_significant_knee.csv', index=False)