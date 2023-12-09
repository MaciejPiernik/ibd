import logging

from ibd.core.data_processing.loaders import load_raw_dataset
from ibd.core.platforms.Platform import Platform


def process_all(datasets):
    for dataset in datasets:
        process(dataset)


def process(dataset_id):
    logging.info('Processing dataset {}'.format(dataset_id))

    raw_dataset = load_raw_dataset(dataset_id)

    platform_id = raw_dataset.metadata['platform_id'][0]
    platform = Platform.get_platform(platform_id)

    dataset = platform.process(raw_dataset)

    dataset.save()
