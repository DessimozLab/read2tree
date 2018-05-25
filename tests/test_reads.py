import unittest
import os
import argparse
from read2tree.Reads import Reads

dirname = os.path.dirname(__file__)


class ReadTest(unittest.TestCase):

    def setUp(self):
        arg_parser = argparse.ArgumentParser(prog='read2tree')

        arg_parser.add_argument('--reads', nargs='+', default=None,
                                help='Reads to be mapped to reference. If paired end '
                                     'add separated by space.')
        arg_parser.add_argument('--read_split_length', type=int, default=400,
                                help='Set read split length.')
        arg_parser.add_argument('--read_split_overlap', type=int, default=50,
                                help='Set read split length overlap.')
        arg_parser.add_argument('-s', '--species_name', default=None,
                                help='[Default is name of read] Name of species '
                                     'for mapped sequence.')

        argv = ['--reads', 'tests/data/reads/test.fq']

        args = arg_parser.parse_args(argv)
        return Reads(args)


    def test_Split(self):
        test_seq = 'ACGTTTTTTGGAAGAGTTAGAGATTTTTAGAGAGGAGGGGTTT'
        expected = ['ACGTTTTTTG', 'GAAGAGTTAG', 'AGATTTTTAG', 'AGAGGAGGGG', 'GGAGGGGTTT']
        reads = self.setUp()
        obtained = reads.split_len(test_seq, 10)
        self.assertEqual(expected, obtained)


    def test_SplitOverlap(self):
        test_seq = 'ACGTTTTTTGGAAGAGTTAGAGATTTTTAGAGAGGAGGGGTTT'
        expected = ['ACGTTTTTTG', 'TTTTGGAAGA', 'GAAGAGTTAG', 'GTTAGAGATT', 'AGATTTTTAG', 'TTTAGAGAGG', 'AGAGGAGGGG', 'GGAGGGGTTT']
        reads = self.setUp()
        obtained = reads.split_len_overlap(test_seq, 10, 5)
        self.assertEqual(expected, obtained)

    def test_readfq(self):
        reads = self.setUp()
        entry = []
        with open(os.path.join(dirname, "data/reads/test.fq")) as f:
            for name, seq, qual in reads.readfq(f):
                entry.append(name)
        self.assertEqual(entry[-1], 'SRR5892449.267')

if __name__ == "__main__":
    unittest.main()
