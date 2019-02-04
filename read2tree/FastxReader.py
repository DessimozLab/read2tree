from __future__ import division
import logging
import gzip
import mimetypes
# from memory_profiler import memory_usage

class FastxReader(object):

    def __init__(self, file):

        self._file = file
        guessed_type = mimetypes.guess_type(file)[1]
        if guessed_type:
            if 'gzip' in guessed_type:
                self._file_handle = 'gzip'
        else:
            self._file_handle = 'txt'

    def open_fastx(self):
        if self._file_handle in 'gzip':
            return gzip.open(self._file, 'rt')
        else:
            return open(self._file, 'rt')

    def readfq_id(self, file_handle):
        for l in file_handle:
            name = l.rstrip()
            seq = next(file_handle).rstrip()
            tmp = next(file_handle).rstrip()
            qual = next(file_handle).rstrip()
            yield name.split(' ')[0]

    def readfq(self, file_handle):
        for l in file_handle:
            name = l.rstrip()
            seq = next(file_handle).rstrip()
            tmp = next(file_handle).rstrip()
            qual = next(file_handle).rstrip()
            yield name, seq, qual

    def readfa(self, file_handle):
        for l in file_handle:
            name = l.rstrip()
            seq = next(file_handle).rstrip()
            yield name, seq
