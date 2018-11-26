from datetime import date
import logging
import logging.config
import yaml
import os

logging.getLogger(__name__).addHandler(logging.NullHandler())

__version__ = '0.1.2'
__copyright__ = 'read2tree (C) 2017-{:d} David Dylus' \
                .format(date.today().year)

path = './log.yaml'
if os.path.exists(path):
    with open(path, 'rt') as f:
        config = yaml.load(f.read())
    logging.config.dictConfig(config)
