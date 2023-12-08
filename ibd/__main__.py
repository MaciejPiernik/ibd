import logging
import click

from ibd.experiments.Experiment import Experiment


logging.basicConfig(level=logging.INFO)


@click.group()
def cli():
    pass


@cli.command()
def run_pipeline():
    logging.info('Running data pipeline')


@cli.command()
@click.option('--experiment_name', '-e', help='The experiment to run')
def run_experiment(experiment_name):
    experiment = Experiment(experiment_name)
    
    experiment.load_data()

    results = experiment.run()

    print(results)


if __name__ == '__main__':
    cli()