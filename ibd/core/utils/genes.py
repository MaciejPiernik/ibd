import mygene


def entrez_id_to_gene_symbol(X, drop_na=True):
    X.columns = X.columns.map(str)  # Make sure ENTREZ IDs are strings

    # Initialize MyGeneInfo
    mg = mygene.MyGeneInfo()

    # Create a list of ENTREZ IDs
    entrez_ids = list(X.columns)

    # Query MyGene.Info for the gene symbols
    gene_info = mg.querymany(entrez_ids, scopes='entrezgene', fields='symbol', species='human', as_dataframe=True)

    # Check for any missing symbols and handle them appropriately
    gene_info['symbol'] = gene_info['symbol'].dropna()

    # Create a mapping from ENTREZ IDs to symbols where 'notfound' is not True and drop the duplicates
    entrez_to_symbol = gene_info.loc[gene_info['notfound'] != True, 'symbol'].drop_duplicates()

    # Replace columns in the dataframe with their gene symbols
    X.columns = X.columns.map(entrez_to_symbol)

    X.columns = X.columns.astype(str)

    # drop columns with nan as name
    if drop_na and ('nan' in X.columns):
        X = X.drop(['nan'], axis=1)

    return X