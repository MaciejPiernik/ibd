
import logging

class GPL1708(Platform):
    def process(self, dataset):
        logging.info('Processing dataset {}'.format(dataset))

        return dataset