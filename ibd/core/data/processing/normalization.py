import pandas as pd
import numpy as np

from sklearn.preprocessing import RobustScaler


# Function to normalize datasets individually
def robust_zscore_normalization_per_dataset(df, labels=None):
    if labels is None:
        labels = pd.Series([0] * len(df), index=df.index)

    # Initialize a DataFrame to hold the normalized data
    df_normalized = pd.DataFrame(index=df.index, columns=df.columns)

    # Loop over each unique dataset label and apply RobustScaler individually
    for dataset in labels.unique():
        # Filter the data for the current dataset
        dataset_mask = labels == dataset
        data_subset = df.loc[dataset_mask]

        if data_subset.empty:
            continue

        # Initialize RobustScaler object and fit to the current dataset
        scaler = RobustScaler()
        scaled_subset = scaler.fit_transform(data_subset)

        # Assign the normalized data back to the respective positions in the normalized DataFrame
        df_normalized.loc[dataset_mask] = scaled_subset

    return df_normalized


def quantile_normalize(df):
    # Step 1: Rank each column
    ranks = df.rank(method='min').stack().astype(int)

    # Step 2: Sort the entire DataFrame and average the values with the same rank
    sorted_df = pd.DataFrame(np.sort(df.values, axis=0), index=df.index, columns=df.columns)
    rank_means = sorted_df.stack().groupby(ranks).mean()

    # Step 3: Assign the average value to each original rank position
    normalized_df = ranks.map(rank_means).unstack()
    
    return normalized_df