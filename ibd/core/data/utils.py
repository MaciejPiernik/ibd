import logging
import os
import pandas as pd
import numpy as np


def read_all_datasets(data_dir: str, dropna=True) -> pd.DataFrame:
    db = pd.DataFrame()
    dataset_label = pd.Series(name='dataset')
    for filename in os.listdir(data_dir):
        if not filename.endswith('.parquet') or filename.endswith('_metadata.parquet'):
            continue

        df = pd.read_parquet(os.path.join(data_dir, filename))
        dataset_name = filename[:-8]

        # Auto-detect and fix non-log-transformed expression data
        median_expr = df.median(axis=1).median()
        if median_expr > 100:
            logging.warning('%s appears non-log-transformed (median=%.1f). Applying log2(x+1).', dataset_name, median_expr)
            df = np.log2(df + 1)

        db = pd.concat([db, df])
        tmp = pd.Series([dataset_name] * len(df), index=df.index, name='dataset')
        dataset_label = pd.concat([dataset_label, tmp])

    if dropna:
        db = db.dropna(axis=1)

    return db, dataset_label


def read_all_metadata(data_dir: str) -> pd.DataFrame:
    metadata = pd.DataFrame()
    for filename in os.listdir(data_dir):
        if not filename.endswith('_metadata.csv'):
            continue

        # read csv
        df = pd.read_csv(os.path.join(data_dir, filename), index_col=0)

        # append to dfs
        metadata = pd.concat([metadata, df])

    return metadata