import logging
import click
import pandas as pd

from ibd.experiments.Experiment import Experiment
from ibd.core.data_processing.processing import process_all

logging.basicConfig(level=logging.INFO)


@click.group()
def cli():
    pass


@cli.command()
@click.option('--datasets', '-d', help='The datasets to process', default='all')
def run_pipeline(datasets):
    logging.info('Running data processing pipeline')

    if datasets == 'all':
        daq = pd.read_csv('./data/DAQ.csv', sep=';')
        datasets = daq['Accession number']
    else:
        datasets = datasets.split(',')

    process_all(datasets)


@cli.command()
@click.option('--experiment_name', '-e', help='The experiment to run')
def run_experiment(experiment_name):
    experiment = Experiment(experiment_name)
    
    experiment.load_data()

    results = experiment.run()

    print(results)


if __name__ == '__main__':
    cli()