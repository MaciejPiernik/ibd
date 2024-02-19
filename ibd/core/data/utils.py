import os
import pandas as pd


def read_all_datasets(data_dir: str, dropna=True) -> pd.DataFrame:
    db = pd.DataFrame()
    dataset_label = pd.Series(name='dataset')
    for filename in os.listdir(data_dir):
        if not filename.endswith('.parquet') or filename.endswith('_metadata.parquet'):
            continue

        df = pd.read_parquet(os.path.join(data_dir, filename))

        db = pd.concat([db, df])
        tmp = pd.Series([filename[:-8]] * len(df), index=df.index, name='dataset')
        dataset_label = pd.concat([dataset_label, tmp])

    if dropna:
        db = db.dropna(axis=1)

    return db, dataset_label


def read_all_metadata(data_dir: str) -> pd.DataFrame:
    metadata = pd.DataFrame()
    for filename in os.listdir(data_dir):
        if not filename.endswith('_metadata.parquet'):
            continue

        # read csv
        df = pd.read_parquet(os.path.join(data_dir, filename))

        # append to dfs
        metadata = pd.concat([metadata, df])

    return metadata