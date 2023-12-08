import GEOparse
import pandas as pd
import numpy as np

# gse = GEOparse.get_GEO(filepath='path_to_your_file')

# Load the dataset from GEO
gse = GEOparse.get_GEO(geo='GSE11223', destdir="./data/geo_cache")

# Get the platform ID (if needed)
platform_id = gse.metadata['platform_id'][0]
print(f"Platform ID: {platform_id}")

# Get expression matrix
expr_matrix = gse.pivot_samples('VALUE')

# Transpose the matrix and convert it to a DataFrame
expr_matrix_t = expr_matrix.T
expr_df = pd.DataFrame(expr_matrix_t)

# Get gene symbols
gene_symbols = gse.gpls[next(iter(gse.gpls))].table['Gene Symbol']

# Assign gene symbols as columns names
expr_df.columns = gene_symbols

# Remove columns with empty gene symbols
expr_df = expr_df.loc[:, expr_df.columns != ""]

# Show the DataFrame
print(expr_df.head())