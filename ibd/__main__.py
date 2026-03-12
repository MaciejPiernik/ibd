import logging
import os
import click
import pandas as pd

from dotenv import load_dotenv

from ibd.core.data.structures.Dataset import Dataset
from ibd.experiments.experiment import build_config, experiment

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


@cli.command(name='run-experiment')
@click.option('--data-dir', help='Directory with dataset parquet and metadata CSV files', default='./db')
@click.option('--output-dir', help='Directory to write CSV results')
@click.option('--alpha', help='FDR threshold', default=0.05, type=float)
@click.option('--min-samples', help='Minimum samples in a dataset', default=6, type=int)
@click.option('--min-per-group', help='Minimum samples per group', default=3, type=int)
@click.option('--skip-gene-mapping', help='Skip Entrez to gene symbol mapping', is_flag=True, default=False)
@click.option('--geneset-libraries', help='Comma-separated GSEA libraries (defaults to standard set)', default=None)
@click.option('--gsea-permutations', help='Number of GSEA permutations', default=1000, type=int)
@click.option('--gsea-min-size', help='Minimum gene set size', default=5, type=int)
@click.option('--gsea-max-size', help='Maximum gene set size', default=1000, type=int)
@click.option('--gsea-seed', help='Random seed for GSEA', default=23, type=int)
def run_experiment(
    data_dir,
    output_dir,
    alpha,
    min_samples,
    min_per_group,
    skip_gene_mapping,
    geneset_libraries,
    gsea_permutations,
    gsea_min_size,
    gsea_max_size,
    gsea_seed,
):
    config = build_config(
        data_dir=data_dir,
        output_dir=output_dir,
        alpha=alpha,
        min_samples=min_samples,
        min_per_group=min_per_group,
        skip_gene_mapping=skip_gene_mapping,
        geneset_libraries=geneset_libraries,
        gsea_permutations=gsea_permutations,
        gsea_min_size=gsea_min_size,
        gsea_max_size=gsea_max_size,
        gsea_seed=gsea_seed,
    )

    experiment(config)


if __name__ == '__main__':
    cli()