import GEOparse
import logging


def load_raw_dataset(dataset_id, cache_dir='./data/geo_cache'):
    logging.info('Loading dataset {}'.format(dataset_id))
    
    gse = GEOparse.get_GEO(geo=dataset_id, destdir=cache_dir)

    return gse