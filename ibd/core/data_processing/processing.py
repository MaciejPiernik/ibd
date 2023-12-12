import logging

from ibd.core.data_processing.loaders import load_raw_dataset
from ibd.core.platforms.platforms import get_platform


def process_all(datasets):
    for dataset in datasets:
        process(dataset)


def process(dataset_id):
    logging.info('Processing dataset {}'.format(dataset_id))

    dataset = load_raw_dataset(dataset_id)

    dataset.process()

    dataset.persist()
