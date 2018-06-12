import logging
import gzip
import mimetypes
import tqdm
import time
import tempfile
# from memory_profiler import memory_usage
from Bio import SeqIO

logger = logging.getLogger(__name__)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('info.log')
file_handler.setFormatter(formatter)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)


class Reads(object):

    def __init__(self, args, load=True):

        self.args = args
        self.elapsed_time = 0
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
            if len(self.args.reads) == 1:
                self._reads = self.args.reads[0]
                self._species_name = self._reads.split("/")[-1].split(".")[0]
                if mimetypes.guess_type(self._reads)[1] in 'gzip':
                    self._file_handle = gzip.open(self._reads, 'rt')
                else:
                    self._file_handle = open(self._reads, 'rt')
            elif len(self.args.reads) == 2:
                self._reads = self.args.reads
                self._species_name = self._reads[0].split("/")[-1].split(".")[0]

        if self.args.species_name:
            self._species_name = self.args.species_name

        if load and self.args.split_reads:
            print('--- Splitting reads from {} ---'.format(self._reads))
            # print(memory_usage(self.process_reads))
            # self.split_reads = self._write_to_tmp_file(self.process_reads())
            self.split_reads = self.process_reads()
        else:
            self.split_reads = self._reads

    def process_reads(self):
        '''
        Main function taking in the reads of the object and processing it
        given the provided parameters
        :return: string that contains all the read sequences separated by '\n'
        '''
        # out = ''
        total_new_reads = 0
        total_reads = 0
        start = time.time()
        out_file = tempfile.NamedTemporaryFile(mode='at', delete=False)
        with self._file_handle as f:
            for name, seq, qual in tqdm.tqdm(self._readfq(f),
                                             desc='Splitting reads',
                                             unit=' read'):
                total_reads += 1
                read_id = name[1:].split(" ")[0]
                # logger.debug("Process read {}".format(read_id))
                if len(seq) > self.split_min_read_len:
                    x = 1
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
