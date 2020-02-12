import unittest
import os
from Bio import SeqIO
from read2tree.OGSet import OG

dirname = os.path.dirname(__file__)


class OGTest(unittest.TestCase):

    def setup(self):
        aa = list(SeqIO.parse('data/OG4.aa', format='fasta'))
        dna = list(SeqIO.parse('data/OG4.dna', format='fasta'))
        og = OG()
        og.aa = aa
        og.dna = dna
        return og

    def test_init(self):
        og = self.setup()
        self.assertEqual(og.dna[0].id, 'MOUSE21964_OG4')

    def test_get_og_dict(self):
        og = self.setup()
        dna_dict = og._get_og_dict(og)
        self.assertEqual(dna_dict['MOUSE21964'].name, 'MOUSE21964_OG4')

    def test_remove_species_records(self):
        og = self.setup()
        og_wo_mouse = og.remove_species_records('MOUSE')
        self.assertEqual(len(og_wo_mouse[0]), 4)
        self.assertEqual(len(og_wo_mouse[1]), 4)

    def test_get_species_id(self):
        og = self.setup()
        dna = og.dna[0]
        aa = og.aa[0]
        self.assertEqual(og._get_species_id(dna), 'MOUSE')
        self.assertEqual(og._get_species_id(aa), 'MOUSE')


if __name__ == "__main__":
    unittest.main()
