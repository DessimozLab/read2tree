#!/usr/bin/env python
'''
    read2tree: from reads to trees using the OMA standalone Output and a set of reads
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
import read2tree
from read2tree.OGSet import OGSet
from read2tree.ReferenceSet import ReferenceSet
from read2tree.Mapper import Mapper
from read2tree.Aligner import Aligner
from read2tree.Progress import Progress
from read2tree.TreeInference import TreeInference
from read2tree.parser import OMAOutputParser
import argparse


COPYRIGHT = '(C) 2017-{:d} David V Dylus'.format(date.today().year)


def parse_args(argv, exe_name, desc):
    '''
        Parses the arguments from the terminal.
    '''
    is_standalone = (exe_name == 'read2tree')
    arg_parser = argparse.ArgumentParser(prog=exe_name,
                                         description=desc,
                                         epilog=read2tree.__copyright__)

    # # Add standard arguments
    # if not is_standalone:
    #     # If standalone, set in parser.
    arg_parser.add_argument('--output_path', default='.',
                                help='[Default is current directory] Path to '
                                     'output directory.')
    arg_parser.add_argument('--version', action='version',
                            help='Show programme\'s version number and exit.',
                            version=read2tree.__version__)

    arg_parser.add_argument('--remove_species', default=None,
                            help='Remove species present in dataset and only do analysis on '
                                 'subset of species')

    arg_parser.add_argument('--standalone_path', default='.',
                            help='[Default is current directory] Path to '
                                 'oma standalone directory.', required=True)

    # arg_parser.add_argument('--ref_og_aa_folder', default='.',
    #                         help='Path to preselected og_aa folder')

    # Arguments to generate the reference
    arg_parser.add_argument('-r', '--reference', action='store_true',
                            help='Just generate the reference dataset for mapping')
    arg_parser.add_argument('--min_species', type=int, default=None,
                            help='Min number of species in selected orthologous groups. \
                            If not selected it will be estimated such that around 1000 OGs are available.')
    arg_parser.add_argument('--dna_reference', default='',
                            help='Reference fasta file that contains nucleotide sequences.')
    arg_parser.add_argument('--single_mapping', default=None,
                            help='Single species file allowing to map in a job array.')
    arg_parser.add_argument('--threads', type=int, default=None,
                            help='Number of threads for the mapping using ngm / ngmlr!')

    # Arguments to map the reads
    arg_parser.add_argument('--ref_folder', default=None,
                            help='Folder containing reference files with sequences sorted by species.')

    arg_parser.add_argument('--reads', nargs='+', default=None, help='Reads to be mapped to reference.')

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

    oma_output = OMAOutputParser(args)
    args.oma_output_path = oma_output.oma_output_path

    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    # TODO: Check whether all the necessary binaries are available
    # TODO: Check all given files and throw error if faulty

    # Read in orthologous groups
    progress = Progress(args, oma_output)

    if progress.status >= 1:
        ogset = OGSet(args, oma_output, load=False)
    else:
        ogset = OGSet(args, oma_output)

    if progress.status >= 2:
        reference = ReferenceSet(args, load=False)
    else:
        reference = ReferenceSet(args, og_set=ogset.ogs, load=True)

    if not args.reference:
        if progress.status >= 3:
            mapper = Mapper(args, og_set=ogset.ogs, load=False)
        else:
            mapper = Mapper(args, ref_set=reference.ref, og_set=ogset.ogs)

        if args.single_mapping is None:
            ogset.add_mapped_seq(mapper.og_records)
            alignments = Aligner(args, ogset.mapped_ogs)
            concat_alignment = alignments.concat_alignment()
            tree = TreeInference(args, concat_alignment=concat_alignment)
            print(tree.tree)
    else:
        print('--- Finished generating references for mapping! ---')
    #
    # # Map sequences to reference
    # if reference:
    #     pass

    print('Time taken {}'.format(timer() - t1))
