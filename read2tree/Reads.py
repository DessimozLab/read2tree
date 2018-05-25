import logging
from Bio import SeqIO, SeqRecord, Seq
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('input.log')
file_handler.setLevel(logging.ERROR)
file_handler.setFormatter(formatter)

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

logger.addHandler(file_handler)
logger.addHandler(stream_handler)

class Reads(object):

    def __init__(self, args, load=False):

        self.args = args

        if self.args.reads:
            if len(self.args.reads) == 1:
                self._reads = self.args.reads[0]
                self._species_name = self._reads.split("/")[-1].split(".")[0]
                self._read_records = list(SeqIO.parse(self._reads, format='fastq'))

        if self.args.species_name:
            self._species_name = self.args.species_name

    # def __call__(self, *args, **kwargs):
    #     split_seq = []
    #     for r in self._reads:
    #         if len(r.seq) > 500:
    #             tmp_seq = self.split_len(r.seq, 400)
    #             tmp_score = self.split_len(r.letter_annotations['phred_quality'], 400)
    #             split_seq += [SeqRecord(s[0], id=r.id + "_" + str(i), name=r.name, description=r.description,
    #                                     letter_annotations={'phred_quality': s[1]}) for i, s in
    #                           enumerate(zip(tmp_seq, tmp_score))]
    #         else:
    #             split_seq.append(r)


    def split_len(self, seq, length):
        split_seqs = [seq[i:i + length] for i in range(0, len(seq), length)]
        if len(split_seqs[-1]) < length:
            split_seqs[-1] = seq[-length:]
        return split_seqs

    def split_len_overlap(self, seq, length, overlap):
        split_seqs = [seq[i:i + length] for i in range(0, len(seq), length - overlap)]
        last_short_value = next(index for index, value in enumerate(split_seqs) if len(value) < length)
        if last_short_value:
            split_seqs[last_short_value] = seq[-length:]
        return split_seqs[:last_short_value + 1]

    def readfq(self, fp):  # this is a generator function
        '''
        This function was copy and pasted from https://github.com/lh3/readfq
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
            if not last: break
            name, seqs, last = last[1:].partition(" ")[0], [], None
            for l in fp:  # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+':  # this is a fasta record
                yield name, ''.join(seqs), None  # yield a fasta record
                if not last: break
            else:  # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp:  # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq):  # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs);  # yield a fastq record
                        break
                if last:  # reach EOF before reading enough quality
                    yield name, seq, None  # yield a fasta record instead
                    break