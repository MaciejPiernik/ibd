import logging
import pandas as pd

from ibd.experiments.Config import Config

class Experiment:
    name: str = None
    config: Config = None
    data: pd.DataFrame = None

    def __init__(self, name: str):
        self.name = name

    def load_data(self):
        logging.info(f'Loading data for experiment: {self.name}')

        data = pd.DataFrame()
        
        self.data = data

    def run(self):
        logging.info('Running experiment: %s', self.__class__.__name__)
        results = {}

        return results