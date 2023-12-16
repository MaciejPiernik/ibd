import logging
import click
import pandas as pd
from ibd.core.data_structures.Dataset import Dataset

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

    daq = pd.read_csv('./data/DAQ.csv', sep=';')

    if datasets == 'all':
        dataset_ids = daq['Accession number']
    else:
        dataset_ids = datasets.split(',')

    datasets = []
    for dataset_id in dataset_ids:
        if dataset_id not in daq['Accession number'].values.tolist():
            raise Exception(f'Dataset {dataset_id} not found in database')
        
        ensembl_release = daq[daq['Accession number'] == dataset_id]['Ensembl release'].iloc[0]
        new_dataset = Dataset(dataset_id, ensembl_release)

        datasets.append(new_dataset)

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