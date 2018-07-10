from __future__ import division
import logging
import gzip
import mimetypes
import time
import tempfile
import random
import os
import numpy as np

from math import ceil
from tqdm import tqdm

from read2tree.FastxReader import FastxReader

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
            guessed_type = mimetypes.guess_type(self.args.reads[0])[1]
            if guessed_type:
                if 'gzip' in guessed_type:
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
                logger.info('{}: --- Splitting reads from {} ---'
                            .format(self._species_name, self._reads))
                # print(memory_usage(self.process_reads))
                # self.split_reads = self._write_to_tmp_file(self\
                #       .process_reads())
                self.reads = self.process_reads()
            else:
                self.reads = self._reads

            if self.args.sample_reads:
                print('--- Sampling reads from {} ---'.format(self.reads))
                logger.info('{}: --- Sampling reads from {} ---'
                            .format(self._species_name, self.reads))
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
        fastq_reader = FastxReader(self._reads)
        with fastq_reader.open_fastx() as f:
            for name, seq, qual in tqdm(fastq_reader.readfx(f),
                                        desc='Splitting reads',
                                        unit=' reads'):
                total_reads += 1
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
                        # write output directly to file to reduce memory
                        # footprint
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

        logger.info('{}: Reads larger than {} were split into {} bp long '
                    'fragments with an overlap of {} bp.'.format(
                        self._species_name, self.split_min_read_len,
                        self.split_len,
                        self.split_overlap))

        logger.info('{}: {} reads were split into {} reads.'
                    .format(self._species_name, total_reads, total_new_reads))

        logger.info('{}: Splitting of reads took {}.'
                    .format(self._species_name, self.elapsed_time))
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
            if idx_random:
                sampled_reads = self._sample_read_file(reads,
                                                       idx_random)
            else:
                sampled_reads = reads
        elif len(self.args.reads) == 2:
            idx_random = self._get_vector_random_reads(reads[0])
            if idx_random:
                sampled_reads.append(self._sample_read_file(reads[0],
                                                            idx_random))
                sampled_reads.append(self._sample_read_file(reads[1],
                                                            idx_random))
            else:
                sampled_reads = reads

        end = time.time()
        elapsed_time = end - start
        logger.info('{}: Sampling of reads took {}.'
                    .format(self._species_name, elapsed_time))
        # if self.args.debug:
        #     shutil.copy(sampled_reads,
        #                 '/Volumes/Untitled/reserach/r2t/test/split.fq')
        return sampled_reads

    def _get_num_reads(self, file):
        if self.total_reads > 0:
            return self.total_reads
        else:
            fastq_reader = FastxReader(file)
            with fastq_reader.open_fastx() as f:
                num_lines = sum([1 for l in f])
            return int(num_lines / 4)

    def _get_read_len(self, file):
        if self.total_reads > 0:
            return self.args.split_len
        else:
            fastq_reader = FastxReader(file)
            with fastq_reader.open_fastx() as f:
                collect_line_len = [len(l) for i, l in enumerate(f)
                                    if (i) % 4 == 1]
                mean_len = np.mean(collect_line_len)
                median_len = np.median(collect_line_len)
            logger.info('{}: The reads have a mean length of {} '
                        'and a median length of {}.'
                        .format(self._species_name, mean_len, median_len))
            return mean_len

    def _get_num_reads_by_coverage(self, file):
        read_len = self._get_read_len(file)
        logger.info('{}: Average read length estimated to {}.'
                    .format(self._species_name, read_len))
        return int(ceil(self.args.genome_len * self.args.coverage /
                        (len(self.args.reads) * read_len)))

    def _get_vector_random_reads(self, file):
        num_reads_by_coverage = self._get_num_reads_by_coverage(file)
        total_records = self._get_num_reads(file)
        logger.info('{}: Sampling {} / {} reads for {}X coverage.'
                    .format(self._species_name, num_reads_by_coverage,
                            total_records, self.coverage))
        if num_reads_by_coverage > total_records:
            logger.info("{}: Not enough reads available for sampling, using "
                        "them all.".format(self._species_name))
            return None
        else:
            return set(random.sample(range(total_records + 1),
                                     num_reads_by_coverage))

    def _sample_read_file(self, file, output_sequence_sets):
        initial_length = 0
        sampling_length = 0

        record_number = 0
        # print(os.path.getsize(file))
        out_file = tempfile.NamedTemporaryFile(mode='at', suffix='.fq',
                                               delete=False)
        fastq_reader = FastxReader(file)
        with fastq_reader.open_fastx() as read_input:
            for line1 in tqdm(read_input, desc='Sampling reads from {}'
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
        logger.info('{}: Cummulative length of all reads {}bp. Cummulative '
                    'length of sampled reads {}bp'
                    .format(self._species_name, initial_length,
                            sampling_length))
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

    def _write_to_tmp_file(self, split_reads):
        with tempfile.NamedTemporaryFile(mode='wt') as filehandle:
            filehandle.write(split_reads)
            filehandle.seek(0)
        return filehandle.name
