import logging

from ibd.core.platforms.GPL1708 import GPL1708

class Platform:
    @staticmethod
    def get_platform(platform_id):
        if platform_id == 'GPL1708':
            return GPL1708()
        else:
            raise Exception('Unknown platform {}'.format(platform_id))