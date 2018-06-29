import unittest
import os
import gzip
import argparse
from Bio import SeqIO
from read2tree.Reads import Reads

dirname = os.path.dirname(__file__)


class ReadTest(unittest.TestCase):

    def setup_long_reads(self, split=False):
        arg_parser = argparse.ArgumentParser(prog='read2tree')

        arg_parser.add_argument('--reads', nargs='+', default=None,
                                help='Reads to be mapped to reference.'
                                'If paired end add separated by space.')
        arg_parser.add_argument('--split_len', type=int, default=400,
                                help='Set read split length.')
        arg_parser.add_argument('--split_overlap', type=int, default=50,
                                help='Set read split length overlap.')
        arg_parser.add_argument('-s', '--species_name', default=None,
                                help='[Default is name of read] Name of '
                                'species for mapped sequence.')
        arg_parser.add_argument('--debug', action='store_true',
                                help='Changes to debug mode: '
                                     '* bam files are saved!'
                                     '* reads are saved by mapping to OG')
        arg_parser.add_argument('--split_reads', action='store_true',
                                help='Splits reads as defined by split_len '
                                '(400) and split_overlap (0) parameters. ')
        arg_parser.add_argument('--split_min_read_len', type=int, default=500,
                                help='[Default is 500] Reads longer than this '
                                'value are cut into smaller values as defined '
                                'by --split_len.')
        arg_parser.add_argument('--sample_reads', action='store_true',
                                help='Splits reads as defined by split_len (400) '
                                'and split_overlap (0) parameters. ')
        arg_parser.add_argument('--coverage', type=float, default=10,
                                help='[Default is 10] coverage in X.')
        arg_parser.add_argument('--genome_len', type=int, default=1000,
                                help='[Default is 1000] Genome size in bp.')

        if split:
            argv = ['--reads', 'tests/data/reads/test.fq.gz', '--split_reads']
        else:
            argv = ['--reads', 'tests/data/reads/test.fq.gz', '']

        args = arg_parser.parse_args(argv)
        return Reads(args)

    def setup_reads_paired(self, sampling=False):
        arg_parser = argparse.ArgumentParser(prog='read2tree')

        arg_parser.add_argument('--reads', nargs='+', default=None,
                                help='Reads to be mapped to reference.'
                                'If paired end add separated by space.')
        arg_parser.add_argument('--split_len', type=int, default=400,
                                help='Set read split length.')
        arg_parser.add_argument('--split_overlap', type=int, default=50,
                                help='Set read split length overlap.')
        arg_parser.add_argument('-s', '--species_name', default=None,
                                help='[Default is name of read] Name of '
                                'species for mapped sequence.')
        arg_parser.add_argument('--debug', action='store_true',
                                help='Changes to debug mode: '
                                     '* bam files are saved!'
                                     '* reads are saved by mapping to OG')
        arg_parser.add_argument('--split_reads', action='store_true',
                                help='Splits reads as defined by split_len '
                                '(400) and split_overlap (0) parameters. ')
        arg_parser.add_argument('--split_min_read_len', type=int, default=500,
                                help='[Default is 500] Reads longer than this '
                                'value are cut into smaller values as defined '
                                'by --split_len.')
        arg_parser.add_argument('--sample_reads', action='store_true',
                                help='Splits reads as defined by split_len (400) '
                                'and split_overlap (0) parameters. ')
        arg_parser.add_argument('--coverage', type=float, default=10,
                                help='[Default is 10] coverage in X.')
        arg_parser.add_argument('--genome_len', type=int, default=1000,
                                help='[Default is 1000] Genome size in bp.')

        if sampling:
            argv = ['--reads', 'tests/data/reads/test_1a.fq.gz',
                    'tests/data/reads/test_2a.fq.gz', '--sample_reads']
        else:
            argv = ['--reads', 'tests/data/reads/test_1a.fq.gz',
                    'tests/data/reads/test_2a.fq.gz']

        args = arg_parser.parse_args(argv)
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

    def test_readfq(self):
        reads = self.setup_long_reads()
        entry = []
        with gzip.open(os.path.join(dirname, "data/reads/test.fq.gz"), 'rt') \
                as f:
            for name, seq, qual in reads._readfq(f):
                entry.append(name)
        self.assertEqual(entry[-1], '@SRR5314792.33 length=1041')

    def test_get_4_line_fastq_string(self):
        reads = self.setup_long_reads()
        expected = '@SRR00001 length=16\nACGTTTGGGAAGGTTT\n+SRR00001 ' \
                   'length=16\n????????????????\n'
        read_id = 'SRR00001'
        x = 0
        seq = 'ACGTTTGGGAAGGTTT'
        qual = '????????????????'
        name = reads._get_4_line_fastq_string(read_id, x, seq, qual)
        self.assertEqual(name, expected)

    def test_read_num_split(self):
        reads = self.setup_long_reads(split=True)
        num_reads = reads._get_num_reads('tests/data/reads/test.fq.gz')
        self.assertEqual(num_reads, 18)

    def test_read_len_split(self):
        reads = self.setup_long_reads(split=True)
        len_reads = reads._get_read_len('tests/data/reads/test.fq.gz')
        self.assertEqual(len_reads, 400)

    def test_read_num_paired(self):
        reads = self.setup_reads_paired()
        num_reads = reads._get_num_reads('tests/data/reads/test_1a.fq.gz')
        self.assertEqual(num_reads, 1000)

    def test_read_len_paired(self):
        reads = self.setup_reads_paired()
        num_reads = reads._get_read_len('tests/data/reads/test_1a.fq.gz')
        self.assertEqual(num_reads, 152)

    def test_read_num_by_coverage_paired(self):
        reads = self.setup_reads_paired()
        num_reads = reads._get_num_reads_by_coverage(
            'tests/data/reads/test_1a.fq.gz')
        self.assertEqual(num_reads, 33)

    def test_read_num_by_coverage_split(self):
        reads = self.setup_long_reads(split=True)
        num_reads = reads._get_num_reads_by_coverage(
            'tests/data/reads/test.fq.gz')
        self.assertEqual(num_reads, 25)

    def test_read_vec_paired(self):
        reads = self.setup_reads_paired()
        num_reads = reads._get_vector_random_reads(
            'tests/data/reads/test_1a.fq.gz')
        self.assertEqual(len(num_reads), 33)

    # def test_sample_reads_file(self):
    #     reads = self.setup_reads_paired(sampling=True)
    #     # output_sequence_sets = reads._get_vector_random_reads(
    #     #     'tests/data/reads/test_1a.fq.gz')
    #     # read_file = reads._sample_read_file('tests/data/reads/test_1a.fq.gz',
    #     #                                     output_sequence_sets)
    #     read_file = reads.reads
    #     # print(read_file)
    #     read_rec = list(SeqIO.parse(read_file, format="fastq"))
    #     print(read_rec)
    #     num_reads = reads._get_num_reads(read_file)
    #     self.assertEqual(num_reads, 25)


if __name__ == "__main__":
    unittest.main()
