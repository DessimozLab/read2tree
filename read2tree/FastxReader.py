from __future__ import division
import logging
import gzip
import mimetypes
# from memory_profiler import memory_usage

logger = logging.getLogger(__name__)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('info.log')
file_handler.setFormatter(formatter)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)


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

    def readfx(self, file_handle):  # this is a generator function
        '''
        This function was copy and pasted from https://github.com/lh3/readfq
        Readfq is a fast implementation of a read iterator and provides a
        massive spead up compared to regular
        implementations
        :param file_handle: is a filehandle
        :return: name, seq, quality
        '''
        last = None  # this is a buffer keeping the last unprocessed line
        while True:  # mimic closure; is it a bad idea?
            if not last:  # the first record or a record following a fastq
                for l in file_handle:  # search for the start of the next record
                    if l[0] in '>@':  # fasta/q header line
                        last = l[:-1]  # save this line
                        break
            if not last:
                break
            name, seqs, last = last, [], None
            for l in file_handle:  # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+':  # this is a fasta record
                yield name, ''.join(seqs), None  # yield a fasta record
                if not last:
                    break
            else:  # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in file_handle:  # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq):  # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs)  # yield a fastq record
                        break
                if last:  # reach EOF before reading enough quality
                    yield name, seq, None  # yield a fasta record instead
                    break
