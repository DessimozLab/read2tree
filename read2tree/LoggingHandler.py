import logging
from multiprocessing_logging import install_mp_handler
import os


class LoggingHandler(object):

    def __init__(self, args, name, logger=None):
        if not logger:
            self.logger = logging.getLogger(name)
            self.formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
            self.file_handler = logging.FileHandler(os.path.join(args.output_path,
                                                                 'info.log'))
            self.file_handler.setFormatter(self.formatter)
            self.stream_handler = logging.StreamHandler()
            self.stream_handler.setFormatter(self.formatter)

            if args.debug:
                self.logger.setLevel(logging.DEBUG)
                self.file_handler.setLevel(logging.DEBUG)
                # stream_handler.setLevel(logging.DEBUG)
            else:
                self.logger.setLevel(logging.INFO)
                self.file_handler.setLevel(logging.INFO)
                # stream_handler.setLevel(logging.INFO)

            self.logger.addHandler(self.file_handler)
            # logger.addHandler(stream_handler)
        else:
            self.logger = logger
            self.formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
            self.file_handler = logging.FileHandler(os.path.join(args.output_path,
                                                                 'info.log'))
            self.file_handler.setFormatter(self.formatter)
            self.stream_handler = logging.StreamHandler()
            self.stream_handler.setFormatter(self.formatter)

            if args.debug:
                self.logger.setLevel(logging.DEBUG)
                self.file_handler.setLevel(logging.DEBUG)
                # stream_handler.setLevel(logging.DEBUG)
            else:
                self.logger.setLevel(logging.INFO)
                self.file_handler.setLevel(logging.INFO)
                # stream_handler.setLevel(logging.INFO)

            self.logger.addHandler(self.file_handler)
            # logger.addHandler(stream_handler)
