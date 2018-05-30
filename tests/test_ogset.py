import os
import unittest
from read2tree import OGSet

API_URL = 'http://omabrowser.org/api'

class OGSetTest(unittest.TestCase):
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
        return OGSet(args)

    def test_OGSet(self):
        raise NotImplementedError

    def test_marker_genes_input(self):
        raise NotImplementedError

    def test_omastandalone_input(self):
        raise NotImplementedError

    def test_output_folder_structure(self):
        raise NotImplementedError

    def test_species_removal(self):
        raise NotImplementedError

    def test_species_removal_after_mapping(self):
        raise NotImplementedError

    def test_rest_api_connection(self):
        OGSet._read

    def test_rest_api_dna_downlaod(self):
        raise NotImplementedError


if __name__ == "__main__":
    unittest.main()
