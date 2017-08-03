#!/usr/bin/env python
'''
    pore2tree: from reads to trees using the OMA standalone Output and a set of reads
    1) DNA sequences are retrieved from OMA standalone Output based on the predicted
        ogs and a threshold of minimum species to be included in an OG
    2) Reads are mapped to the DNA sequences (using NGM)
    3) Consensus sequences are build from read mappings
    4) OGs are reconstructed and the most probable aa translation is found
    5) OGs are aligned
    6) Trees are build using the alignments

    This is the main entry point to the program.

    -- David V Dylus, July--XX 2017
'''
import os
from datetime import date
from timeit import default_timer as timer
import pore2tree
from pore2tree.OGSet import OGSet
from pore2tree.ReferenceSet import ReferenceSet
from pore2tree.Mapper import Mapper
from pore2tree.Aligner import Aligner
from pore2tree.Progress import Progress
from pore2tree.TreeInference import TreeInference
import argparse


COPYRIGHT = '(C) 2017-{:d} David V Dylus'.format(date.today().year)


def parse_args(argv, exe_name, desc):
    '''
        Parses the arguments from the terminal.
    '''
    is_standalone = (exe_name == 'pore2tree')
    arg_parser = argparse.ArgumentParser(prog=exe_name,
                                         description=desc,
                                         epilog=pore2tree.__copyright__)

    # # Add standard arguments
    # if not is_standalone:
    #     # If standalone, set in parser.
    arg_parser.add_argument('--output_path', default='.',
                                help='[Default is current directory] Path to '
                                     'output directory.')
    arg_parser.add_argument('--version', action='version',
                            help='Show programme\'s version number and exit.',
                            version=pore2tree.__version__)

    arg_parser.add_argument('--remove_species', default=None,
                            help='Remove species present in dataset and only do analysis on '
                                 'subset of species')

    arg_parser.add_argument('--standalone_path', default='.',
                            help='[Default is current directory] Path to '
                                 'oma standalone directory.', required=True)

    arg_parser.add_argument('--ref_og_aa_folder', default='.',
                            help='Path to preselected og_aa folder')

    # Arguments to generate the reference
    arg_parser.add_argument('-r', '--reference', action='store_true',
                            help='Just generate the reference dataset for mapping')
    arg_parser.add_argument('--min_species', type=int, default=None,
                            help='Min number of species in selected orthologous groups. \
                            If not selected it will be estimated such that around 1000 OGs are available.')
    arg_parser.add_argument('--dna_reference', default=None,
                            help='Reference fasta file that contains nucleotide sequences.')

    # Arguments to map the reads
    arg_parser.add_argument('--ref_folder', default=None,
                            help='Folder containing reference files with sequences sorted by species.')
    arg_parser.add_argument('--reads', default=None,
                            help='Reads to be mapped to reference.')


    # Parse the arguments.
    args = arg_parser.parse_args(argv)

    return args

def check_execution_status():
    raise NotImplementedError

def main(argv, exe_name, desc=''):
    '''
        Main function.
    '''

    t1 = timer()
    # Parse
    args = parse_args(argv, exe_name, desc)

    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    # TODO: Check whether all the necessary binaries are available
    # TODO: Check all given files and throw error if faulty

    # Read in orthologous groups
    progress = Progress(args)

    if progress.status >= 1:
        ogset = OGSet(args, load=False)
    else:
        ogset = OGSet(args)

    if progress.status >= 2:
        reference = ReferenceSet(args, load=False)
    else:
        reference = ReferenceSet(args, ogset=ogset.ogs, load=True)

    mapper = Mapper(args, ref_set=reference.ref, og_set=ogset.ogs)
    ogset.add_mapped_seq(mapper.og_records)
    alignments = Aligner(args, ogset.mapped_ogs)
    concat_alignment = alignments.concat_alignment()
    tree = TreeInference(args, concat_alignment=concat_alignment)
    print(tree)
    #
    # # Map sequences to reference
    # if reference:
    #     pass

    print('Time taken {}'.format(timer() - t1))
