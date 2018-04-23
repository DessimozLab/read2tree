import unittest

from read2tree.Progress import Progress
from read2tree.stats.Coverage import Coverage
from read2tree.stats.SeqCompleteness import SeqCompleteness
import os

class OGSetTest(unittest.TestCase):

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
        raise NotImplementedError

    def test_rest_api_dna_downlaod(self):
        raise NotImplementedError


if __name__ == "__main__":
    unittest.main()
