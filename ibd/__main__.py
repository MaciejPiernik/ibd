import logging
import click
import pandas as pd

from dotenv import load_dotenv

from ibd.core.data.structures.Dataset import Dataset
from ibd.experiments.Experiment import Experiment
from pathlib import Path

logging.basicConfig(level=logging.INFO)

load_dotenv()


@click.group()
def cli():
    pass


@cli.command()
@click.option('--datasets', '-d', help='The datasets to process', default='all')
@click.option('--recompute_data', '-r', help='Recompute data', is_flag=True, default=False)
@click.option('--recompute_metadata', '-m', help='Recompute metadata', is_flag=True, default=False)
def run_pipeline(datasets, recompute_data, recompute_metadata):
    db_path = Path("./db")
    if not db_path.exists():
        db_path.mkdir(parents=True)
    logging.info('Running data processing pipeline')

    daq = pd.read_csv('./data/DAQ.csv', sep=';')

    if datasets == 'all':
        dataset_ids = daq['Accession number']
    else:
        dataset_ids = datasets.split(',')

    for dataset_id in dataset_ids:
        if dataset_id not in daq['Accession number'].values.tolist():
            raise Exception(f'Dataset {dataset_id} not found in database')

        dataset = Dataset(dataset_id)

        logging.info('Processing dataset {}'.format(dataset.id))
        try:
            if dataset.data_in_db() and not recompute_data and dataset.metadata_in_db() and not recompute_metadata:
                logging.info('Dataset {} data already in database, skipping'.format(dataset.id))
                continue
            
            dataset.process(recompute_data, recompute_metadata)

        except Exception as e:
            logging.error('Error processing dataset {}: {}'.format(dataset.id, e))


@cli.command()
@click.option('--experiment_name', '-e', help='The experiment to run')
def run_experiment(experiment_name):
    experiment = Experiment(experiment_name)
    
    experiment.load_data()

    results = experiment.run()

    print(results)


if __name__ == '__main__':
    cli()