import unittest
import os
import gzip
import argparse
from Bio import SeqIO
from read2tree.Reads import Reads
from read2tree.FastxReader import FastxReader
from read2tree.main import parse_args
from read2tree._utils import exe_name
dirname = os.path.dirname(__file__)


class ReadTest(unittest.TestCase):

    def setup_long_reads(self, split=False):
        if split:
            argv = ['--output_path', 'data/output', '--reads', 'data/reads/test.fq.gz', '--split_reads',
                    '--split_overlap', '50', '--split_len', '400', '--sample_reads', '--coverage', '10',
                    '--genome_len', '1000']
        else:
            argv = ['--output_path', 'data/output', '--reads', 'data/reads/test.fq.gz']

        args = parse_args(argv, exe_name(), '')
        # args = arg_parser.parse_args(argv)
        return Reads(args)

    def setup_reads_paired(self, sampling=False):

        if sampling:
            argv = ['--output_path', 'data/output', '--reads', 'data/reads/test_1a.fq.gz',
                    'data/reads/test_2a.fq.gz', '--sample_reads', '--coverage', '10', '--genome_len', '1000']
        else:
            argv = ['--output_path', 'data/output', '--reads', 'data/reads/test_1a.fq.gz',
                    'data/reads/test_2a.fq.gz']
        args = parse_args(argv, exe_name(), '')
        return Reads(args)

    def test_split(self):
        test_seq = 'ACGTTTTTTGGAAGAGTTAGAGATTTTTAGAGAGGAGGGGT'
        expected = ['ACGTTTTTTG', 'GAAGAGTTAG', 'AGATTTTTAG', 'AGAGGAGGGG',
                    'GAGGAGGGGT']
        reads = self.setup_long_reads()
        # obtained = reads._split_len(test_seq, 10)
        obtained = reads._split_len_overlap(test_seq, 10, 0)
        self.assertEqual(expected, obtained)

    def test_splitOverlap(self):
        test_seq = 'ACGTTTTTTGGAAGAGTTAGAGATTTTTAGAGAGGAGGGGTTT'
        expected = ['ACGTTTTTTG', 'TTTTGGAAGA', 'GAAGAGTTAG', 'GTTAGAGATT',
                    'AGATTTTTAG', 'TTTAGAGAGG', 'AGAGGAGGGG', 'GGAGGGGTTT']
        reads = self.setup_long_reads()
        obtained = reads._split_len_overlap(test_seq, 10, 5)
        # print(reads._split_len_overlap('TTTTTAGAGAGGAGGGGTTT', 10, 5))
        self.assertEqual(expected, obtained)

    def test_get_4_line_fastq_string(self):
        reads = self.setup_long_reads()
        expected = '@SRR00001 length=16\nACGTTTGGGAAGGTTT\n+SRR00001 ' \
                   'length=16\n????????????????\n'
        read_id = 'SRR00001'
        seq = 'ACGTTTGGGAAGGTTT'
        qual = '????????????????'
        name = reads._get_4_line_fastq_string(read_id, seq, qual, x=0)
        self.assertEqual(name, expected)

    def test_read_num_split(self):
        reads = self.setup_long_reads(split=True)
        num_reads = reads._get_num_reads('data/reads/test.fq.gz')
        self.assertEqual(num_reads, 18)

    def test_read_len_split(self):
        reads = self.setup_long_reads(split=True)
        len_reads = reads._get_read_len('data/reads/test.fq.gz',1000)
        self.assertEqual(len_reads, 400)

    def test_read_num_paired(self):
        reads = self.setup_reads_paired()
        num_reads = reads._get_num_reads('data/reads/test_1a.fq.gz')
        self.assertEqual(num_reads, 1000)

    def test_read_len_paired(self):
        reads = self.setup_reads_paired()
        num_reads = reads._get_read_len('data/reads/test_1a.fq.gz', 1000)
        self.assertEqual(num_reads, 151.0)

    def test_read_num_by_coverage_paired(self):
        reads = self.setup_reads_paired(sampling=True)
        num_reads = reads._get_num_reads_by_coverage(
            'data/reads/test_1a.fq.gz', 1000)
        self.assertEqual(num_reads, 34)

    def test_read_num_by_coverage_split(self):
        reads = self.setup_long_reads(split=True)
        num_reads = reads._get_num_reads_by_coverage(['data/reads/test.fq.gz'],1000)
        self.assertEqual(num_reads, 25)

    def test_read_vec_paired(self):
        reads = self.setup_reads_paired(sampling=True)
        num_reads = reads._get_vector_random_reads(
            'data/reads/test_1a.fq.gz')
        self.assertEqual(len(num_reads), 34)


if __name__ == "__main__":
    unittest.main()
