import unittest
import os
import gzip
import argparse
from read2tree.Reads import Reads

dirname = os.path.dirname(__file__)


class ReadTest(unittest.TestCase):

    def setUp(self):
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

        argv = ['--reads', 'tests/data/reads/test.fq.gz', '']

        args = arg_parser.parse_args(argv)
        return Reads(args)

    def test_split(self):
        test_seq = 'ACGTTTTTTGGAAGAGTTAGAGATTTTTAGAGAGGAGGGGT'
        expected = ['ACGTTTTTTG', 'GAAGAGTTAG', 'AGATTTTTAG', 'AGAGGAGGGG',
                    'GAGGAGGGGT']
        reads = self.setUp()
        # obtained = reads._split_len(test_seq, 10)
        obtained = reads._split_len_overlap(test_seq, 10, 0)
        self.assertEqual(expected, obtained)

    def test_splitOverlap(self):
        test_seq = 'ACGTTTTTTGGAAGAGTTAGAGATTTTTAGAGAGGAGGGGTTT'
        expected = ['ACGTTTTTTG', 'TTTTGGAAGA', 'GAAGAGTTAG', 'GTTAGAGATT',
                    'AGATTTTTAG', 'TTTAGAGAGG', 'AGAGGAGGGG', 'GGAGGGGTTT']
        reads = self.setUp()
        obtained = reads._split_len_overlap(test_seq, 10, 5)
        # print(reads._split_len_overlap('TTTTTAGAGAGGAGGGGTTT', 10, 5))
        self.assertEqual(expected, obtained)

    def test_readfq(self):
        reads = self.setUp()
        entry = []
        with gzip.open(os.path.join(dirname, "data/reads/test.fq.gz"), 'rt') \
                as f:
            for name, seq, qual in reads._readfq(f):
                entry.append(name)
        self.assertEqual(entry[-1], '@SRR5314792.33 length=1041')

    def test_get_4_line_fastq_string(self):
        reads = self.setUp()
        expected = '@SRR00001 length=16\nACGTTTGGGAAGGTTT\n+SRR00001 ' \
                   'length=16\n????????????????\n'
        read_id = 'SRR00001'
        x = 0
        seq = 'ACGTTTGGGAAGGTTT'
        qual = '????????????????'
        name = reads._get_4_line_fastq_string(read_id, x, seq, qual)
        self.assertEqual(name, expected)


if __name__ == "__main__":
    unittest.main()
