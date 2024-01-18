import os

import numpy as np
import pandas as pd

from sklearn.base import BaseEstimator, TransformerMixin


class CorrelationFilter(BaseEstimator, TransformerMixin):
    """ A transformer that removes highly correlated features.
        Can be used in a sklearn Pipeline, but doesn't make sense if the dataset contains too many features, as will take forever.
        In that case use the cache_file parameter to save the correlation matrix to a file and reuse it later and simply run the filter before the pipeline.
        It's not 'pure', but at least it won't take a year.
    """
    def __init__(self, threshold: int = 1, cache_file: str = None):
        self.threshold = threshold
        self.cache_file = cache_file
        self.correlation_matrix = None
        self.correlated_pairs = None
        self.columns_to_drop = None

    def fit(self, X, y=None):
        if not isinstance(X, pd.DataFrame):
            if self.cache_file is not None:
                raise TypeError('X must be a pandas DataFrame')
            else:
                X = pd.DataFrame(X)

        if self.cache_file is not None and os.path.exists(self.cache_file):
            self.correlated_pairs = pd.read_csv(self.cache_file)
        else:
            self.correlation_matrix = X.corr().abs()
            upper_triangle = np.triu(np.ones(self.correlation_matrix.shape), k=1).astype(bool)
            tmp = self.correlation_matrix.where(upper_triangle).stack().rename('Correlation')
            tmp.index = tmp.index.rename(['Feature 1', 'Feature 2'])
            self.correlated_pairs = tmp.reset_index()

            if self.cache_file is not None:
                self.correlated_pairs.to_csv(self.cache_file, index=False)
        
        self.columns_to_drop = self.correlated_pairs[self.correlated_pairs['Correlation'] > self.threshold]['Feature 2'].unique()

        return self

    def transform(self, X):
        X = X.drop(columns=self.columns_to_drop, axis=1)

        return X
    
    def fit_transform(self, X, y=None):
        self.fit(X)

        return self.transform(X)