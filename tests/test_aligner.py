import unittest
import os
import gzip
import argparse
from Bio import SeqIO
from read2tree.Reads import Reads
from read2tree.FastxReader import FastxReader
dirname = os.path.dirname(__file__)


class ReadTest(unittest.TestCase):

    def setup_reads_paired(self, sampling=False):
        arg_parser = argparse.ArgumentParser(prog='read2tree')

        arg_parser.add_argument('--standalone_path', default='.',
                                help='[Default is current directory] Path to '
                                     'oma standalone directory.')

        arg_parser.add_argument('--reads', nargs='+', default=None,
                                help='Reads to be mapped to reference. If paired '
                                'end add separated by space.')

        arg_parser.add_argument('--read_type', default='short',
                                help='[Default is short reads] Type of reads to '
                                'use for mapping. Either ngm for short reads or '
                                'ngmlr for long will be used.')

        arg_parser.add_argument('--dna_reference', default='',
                                help='Reference file that contains nucleotide '
                                'sequences (fasta, hdf5). If not given it will use'
                                'the RESTapi and retrieve sequences '
                                'from http://omabrowser.org directly. '
                                'NOTE: internet connection required!')
                                
        arg_parser.add_argument('--keep_all_ogs', action='store_true',
                                help='Keep all orthologs after addition of '
                                'mapped seq, which means also the groups that '
                                'have no mapped sequence. Otherwise only groups '
                                'are used that have the mapped sequence for '
                                'alignment and tree inference.')

        arg_parser.add_argument('-r', '--reference', action='store_true',
                                help='Just generate the reference dataset for '
                                'mapping.')

        arg_parser.add_argument('--remove_species_ogs', default=None,
                                help='[Default is none] Remove species present '
                                'in data set after mapping step completed to '
                                'build OGs. Input is comma separated list '
                                'without spaces, e.g. XXX,YYY,AAA.')

        arg_parser.add_argument('-s', '--species_name', default=None,
                                help='[Default is name of read] Name of species '
                                     'for mapped sequence.')

        arg_parser.add_argument('--output_path', default='.', required=True,
                                help='[Default is current directory] Path to '
                                'output directory.')

        argv = ['--standalone_path', 'tests/data/marker_genes/',
                '--dna_reference', 'tests/data/dna.fa', '--reads',
                'tests/data/mapper/test3/test_1b.fq',
                'tests/data/mapper/test3/test_2b.fq',
                '--output_path', 'tests/data/output', '--read_type',
                'short', '--keep_all_ogs', '--reference',
                '--remove_species_ogs', 'CIOIN', '--species_name', 'ass']

        args = arg_parser.parse_args(argv)
        return alignments = Aligner(args, ogset.ogs, load=True)
