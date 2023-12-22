import logging


def process_all(datasets):
    for dataset in datasets:
        try:
            process(dataset)
        except Exception as e:
            logging.error('Error processing dataset {}: {}'.format(dataset.id, e))


def process(dataset):
    logging.info('Processing dataset {}'.format(dataset.id))

    dataset.load_raw_data()

    dataset.process()

    dataset.persist()
