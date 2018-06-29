from __future__ import division
import logging
import gzip
import mimetypes
import tqdm
import time
import tempfile
import random
import os

from math import ceil
import shutil
# from memory_profiler import memory_usage

logger = logging.getLogger(__name__)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('info.log')
file_handler.setFormatter(formatter)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)


class Reads(object):

    def __init__(self, args, load=True):

        self.args = args
        self.coverage = self.args.coverage
        self.genome_len = self.args.genome_len
        self.elapsed_time = 0
        self.total_reads = 0
        if args.debug:
            logger.setLevel(logging.DEBUG)
            file_handler.setLevel(logging.DEBUG)
            # stream_handler.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)
            file_handler.setLevel(logging.INFO)
            # stream_handler.setLevel(logging.INFO)

        logger.addHandler(file_handler)
        # logger.addHandler(stream_handler)

        self.split_len = args.split_len
        self.split_overlap = args.split_overlap
        self.split_min_read_len = args.split_min_read_len

        if self.args.reads:
            if mimetypes.guess_type(self.args.reads[0])[1] in 'gzip':
                self._file_handle = 'gzip'
            else:
                self._file_handle = 'txt'

            if len(self.args.reads) == 1:
                self._reads = self.args.reads[0]
                self._species_name = self._reads.split("/")[-1].split(".")[0]
            elif len(self.args.reads) == 2:
                self._reads = self.args.reads
                self._species_name = self._reads[0].split("/")[-1] \
                    .split(".")[0]

        if self.args.species_name:
            self._species_name = self.args.species_name

        if load:
            if self.args.split_reads:
                print('--- Splitting reads from {} ---'.format(self._reads))
                # print(memory_usage(self.process_reads))
                # self.split_reads = self._write_to_tmp_file(self\
                #       .process_reads())
                self.reads = self.process_reads()
            else:
                self.reads = self._reads

            if self.args.sample_reads:
                print('--- Sampling reads from {} ---'.format(self.reads))
                self.reads = self.sample_from_reads(self.reads)
            else:
                self.reads = self.reads
        else:
            self.reads = self._reads

    def process_reads(self):
        '''
        Function taking in the reads of the object and processing it
        given the provided parameters
        :return: string that contains all the read sequences separated by '\n'
        '''
        # out = ''
        total_new_reads = 0
        total_reads = 0
        start = time.time()
        out_file = tempfile.NamedTemporaryFile(mode='at', suffix='.fq',
                                               delete=False)
        with self._open_reads(self._reads) as f:
            for name, seq, qual in tqdm.tqdm(self._readfq(f),
                                             desc='Splitting reads',
                                             unit=' reads'):
                # total_reads += 1
                read_id = name[1:].split(" ")[0]
                # logger.debug("Process read {}".format(read_id))
                if len(seq) > self.split_min_read_len:
                    x = 0
                    try:
                        new_seq, new_qual = \
                            self._split_len_overlap(seq, self.split_len,
                                                    self.split_overlap), \
                            self._split_len_overlap(qual, self.split_len,
                                                    self.split_overlap)
                    except ValueError:
                        logger.debug('Reads were not split properly!')
                    for i in zip(new_seq, new_qual):
                        # out += self._get_4_line_fastq_string(read_id, x,
                        #                                      i[0],
                        #                                      i[1])
                        out_file.write(self._get_4_line_fastq_string(read_id,
                                                                     x, i[0],
                                                                     i[1]))
                        x += 1
                    total_new_reads += x
                else:
                    # out += self._get_4_line_fastq_string(read_id, None, seq,
                    #                                      qual)
                    out_file.write(self._get_4_line_fastq_string(read_id,
                                                                 None, seq,
                                                                 qual))
                    total_new_reads += 1

        end = time.time()
        self.elapsed_time = end - start

        logger.info('Reads larger than {} were split into {} bp long '
                    'fragments with an overlap of {} bp.'.format(
                        self.split_min_read_len, self.split_len,
                        self.split_overlap))

        logger.info('{} reads were split into {} reads.'.format(
                    total_reads, total_new_reads))

        logger.info('Splitting of reads took {}.'.format(
                    self.elapsed_time))
        out_file.close()
        self.total_reads = total_new_reads
        self._file_handle = 'txt'
        # shutil.move(out_file.name,
        #             '/Users/daviddylus/Research/read2tree/tests\
        #             /data/reads/split.fq')
        return out_file.name

    def sample_from_reads(self, reads):
        '''
        Main function taking in the reads of the object and processing it
        given the provided parameters
        :return: string that contains all the read sequences separated by '\n'
        '''
        sampled_reads = []
        start = time.time()
        if len(self.args.reads) == 1:
            idx_random = self._get_vector_random_reads(reads)
            sampled_reads = self._sample_read_file(reads,
                                                   idx_random)
        elif len(self.args.reads) == 2:
            idx_random = self._get_vector_random_reads(reads[0])
            sampled_reads.append(self._sample_read_file(reads[0],
                                                        idx_random))
            sampled_reads.append(self._sample_read_file(reads[1],
                                                        idx_random))

        end = time.time()
        elapsed_time = end - start
        logger.info('Sampling of reads took {}.'.format(
                    elapsed_time))
        if self.args.debug:
            shutil.copy(sampled_reads,
                        '/Volumes/Untitled/reserach/r2t/test/split.fq')
        return sampled_reads

    def _open_reads(self, file):
        # type = 'txt'
        # try:
        #     type = mimetypes.guess_type(file)[1]
        #     if type in 'gzip':
        #         return gzip.open(file, 'rt')
        #     else:
        #         return open(file, 'rt')
        # except TypeError:
        #     logger.debug("Type of input could not be determined!")
        # else:
        #     type = 'txt'
        #     return open(file, 'rt')
        #     print('second attempt {}'.format(type))
        if self._file_handle in 'gzip':
            return gzip.open(file, 'rt')
        else:
            return open(file, 'rt')

    def _get_num_reads(self, file):
        if self.total_reads > 0:
            return self.total_reads
        else:
            with self._open_reads(file) as f:
                num_lines = sum([1 for l in f])
            return int(num_lines / 4)

    def _get_read_len(self, file):
        if self.total_reads > 0:
            return self.args.split_len
        else:
            with self._open_reads(file) as f:
                head = [next(f) for x in range(2)]
            return len(head[-1])

    def _get_num_reads_by_coverage(self, file):
        read_len = self._get_read_len(file)
        logger.info('Read length estimated to {}.'.format(
                    read_len))
        return int(ceil(self.args.genome_len * self.args.coverage /
                        (len(self.args.reads) * read_len)))

    def _get_vector_random_reads(self, file):
        num_reads_by_coverage = self._get_num_reads_by_coverage(file)
        logger.info('Number of reads {} for {}X coverage.'.format(
                    num_reads_by_coverage, self.coverage))
        total_records = self._get_num_reads(file)
        logger.info('Total number of reads {}.'.format(
                    total_records))
        return set(random.sample(range(total_records + 1),
                                 num_reads_by_coverage))

    def _sample_read_file(self, file, output_sequence_sets):
        initial_length = 0
        sampling_length = 0

        record_number = 0
        # print(os.path.getsize(file))
        out_file = tempfile.NamedTemporaryFile(mode='at', suffix='.fq',
                                               delete=False)
        with self._open_reads(file) as read_input:
            for line1 in tqdm.tqdm(read_input, desc='Sampling reads from {}'
                                   .format(os.path.basename(file)),
                                   unit=' reads'):
                line2 = read_input.readline()
                initial_length += len(line2)
                line3 = read_input.readline()
                line4 = read_input.readline()
                if record_number in output_sequence_sets:
                    out_file.write(line1)
                    out_file.write(line2)
                    out_file.write(line3)
                    out_file.write(line4)
                    sampling_length += len(line2)
                record_number += 1
        logger.info('Cummulative length of all reads {}bp. Cummulative '
                    'length of sampled reads {}bp'.format(initial_length,
                                                          (sampling_length*2
                                                           if len(self._reads)
                                                           == 2 else
                                                           sampling_length)))
        out_file.close()
        return out_file.name

    # def write_split_reads(self, read_string):
    #     outfile = self._reads.replace('.fq', '-split.fq')
    #     with gzip.open(outfile, "wt") as f:
    #         f.write(read_string)

    def _get_4_line_fastq_string(self, read_id, x, seq, qual):
        '''
        Transform 4 lines of read string to new read string providing the
        split information
        :param read_id: Read ID in the form of SRR00001
        :param read_num: Number of read usually after the read ID
        :param x: Numerical iterator
        :param seq: Sequence string
        :param qual: Quality string
        :return: 4 lines that correspond to one read with adapted ID
        '''
        out = ''
        if x:
            new_name = "@" + read_id + "_" + str(x) + ' length=' + \
                str(len(seq))
        else:
            new_name = "@" + read_id + ' length=' + str(len(seq))
        out += new_name + "\n"
        out += seq + "\n"
        out += new_name.replace("@", "+") + '\n'
        out += qual + "\n"
        return out

    def _split_len(self, seq, length):
        split_seqs = [seq[i:i + length] for i in range(0, len(seq), length)]
        if len(split_seqs[-1]) < length:
            split_seqs[-1] = seq[-length:]
        return split_seqs

    def _split_len_overlap(self, seq, length, overlap):
        split_seqs = [seq[i:i + length] for i in range(0, len(seq),
                                                       length - overlap)]
        last_short_value = next((index for index, value in enumerate(
            split_seqs) if len(value) < length), None)
        if last_short_value:
            split_seqs[last_short_value] = seq[-length:]
            return split_seqs[:last_short_value + 1]
        else:
            return split_seqs

    def _readfq(self, fp):  # this is a generator function
        '''
        This function was copy and pasted from https://github.com/lh3/readfq
        Readfq is a fast implementation of a read iterator and provides a
        massive spead up compared to regular
        implementations
        :param fp: is a filehandle
        :return: name, seq, quality
        '''
        last = None  # this is a buffer keeping the last unprocessed line
        while True:  # mimic closure; is it a bad idea?
            if not last:  # the first record or a record following a fastq
                for l in fp:  # search for the start of the next record
                    if l[0] in '>@':  # fasta/q header line
                        last = l[:-1]  # save this line
                        break
            if not last:
                break
            name, seqs, last = last, [], None
            for l in fp:  # read the sequence
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
                for l in fp:  # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq):  # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs)  # yield a fastq record
                        break
                if last:  # reach EOF before reading enough quality
                    yield name, seq, None  # yield a fasta record instead
                    break

    def _write_to_tmp_file(self, split_reads):
        with tempfile.NamedTemporaryFile(mode='wt') as filehandle:
            filehandle.write(split_reads)
            filehandle.seek(0)
        return filehandle.name
