from datetime import date
import logging
import logging.config
import yaml
import os
from pkg_resources import resource_string
logging.getLogger(__name__).addHandler(logging.NullHandler())

__version__ = '2.0.1'
__copyright__ = 'read2tree (C) 2017-{:d} David Dylus,  Adrian M. Altenhoff, Sina Majidian  ' \
                .format(date.today().year)


# path = './log.yaml'
# if os.path.exists(path):
#     with open(path, 'rt') as f:
#         config = yaml.load(f.read())
#     logging.config.dictConfig(config)

#conf = resource_string(__name__, 'logging/log.yaml')

# D = yaml.load(conf, Loader=yaml.FullLoader)
# D.setdefault('version', 1)
# logging.config.dictConfig(D)
# del D
