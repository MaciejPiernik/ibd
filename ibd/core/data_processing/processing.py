import logging


def process_all(datasets):
    for dataset in datasets:
        process(dataset)


def process(dataset):
    logging.info('Processing dataset {}'.format(dataset.id))

    dataset.load_raw_data()

    dataset.process()

    dataset.persist()
