import logging

from ibd.core.platforms.GPL1708 import GPL1708
from ibd.core.platforms.GPL6244 import GPL6244


def get_platform(platform_id):
    if platform_id == 'GPL1708':
        return GPL1708()
    elif platform_id == 'GPL6244':
        return GPL6244()
    else:
        raise Exception('Unknown platform {}'.format(platform_id))
    
