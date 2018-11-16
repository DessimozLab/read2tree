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
import glob
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
    arg_parser.add_argument('--version', action='version',
                            help='Show programme\'s version number and exit.',
                            version=read2tree.__version__)

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

    arg_parser.add_argument('--split_reads', action='store_true',
                            help='Splits reads as defined by split_len (400) '
                            'and split_overlap (0) parameters. ')

    arg_parser.add_argument('--split_len', type=int, default=400,
                            help='[Default is 400] Parameter for selection of '
                            'read split length can only be used in combination'
                            'with with long read option. ')

    arg_parser.add_argument('--split_overlap', type=int, default=0,
                            help='[Default is 0] Reads are split with an '
                            'overlap defined by this argument.')

    arg_parser.add_argument('--split_min_read_len', type=int, default=500,
                            help='[Default is 500] Reads longer than this '
                            'value are cut into smaller values as defined '
                            'by --split_len. ')

    arg_parser.add_argument('--sample_reads', action='store_true',
                            help='Splits reads as defined by split_len (400) '
                            'and split_overlap (0) parameters. ')

    arg_parser.add_argument('--coverage', type=float, default=10,
                            help='[Default is 10] coverage in X.')

    arg_parser.add_argument('--genome_len', type=int, default=2000000,
                            help='[Default is 2000000] Genome size in bp.')

    arg_parser.add_argument('--output_path', default='.', required=True,
                            help='[Default is current directory] Path to '
                            'output directory.')

    arg_parser.add_argument('--dna_reference', default='',
                            help='Reference file that contains nucleotide '
                            'sequences (fasta, hdf5). If not given it will use'
                            'the RESTapi and retrieve sequences '
                            'from http://omabrowser.org directly. '
                            'NOTE: internet connection required!')

    arg_parser.add_argument('--ignore_species', default=None,
                            help='[Default is none] Ignores species part of '
                            'the OMA standalone pipeline. Input is comma '
                            'separated list without spaces, e.g. XXX,YYY,AAA.')

    arg_parser.add_argument('--sc_threshold', type=float, default=0.0,
                            help='[Default is 0.0; Range 0-1] Parameter for '
                            'selection of sequences from mapping by '
                            'completeness compared to its reference sequence '
                            '(number of ACGT basepairs vs length of sequence). '
                            'By default, all sequences are selected.')

    arg_parser.add_argument('--remove_species_mapping', default=None,
                            help='[Default is none] Remove species present in '
                            'data set after mapping step completed and only '
                            'do analysis on subset. Input is comma separated '
                            'list without spaces, e.g. XXX,YYY,AAA.')

    arg_parser.add_argument('--remove_species_ogs', default=None,
                            help='[Default is none] Remove species present '
                            'in data set after mapping step completed to '
                            'build OGs. Input is comma separated list '
                            'without spaces, e.g. XXX,YYY,AAA.')

    arg_parser.add_argument('--ngmlr_parameters', default=None,
                            help='[Default is none] In case this parameters '
                            'need to be changed all 3 values have to be '
                            'changed [x,subread-length,R]. The standard '
                            'is: ont,256,0.25. Possibilities for these '
                            'parameter can be found in the original '
                            'documentation of ngmlr.')

    arg_parser.add_argument('--keep_all_ogs', action='store_true',
                            help='Keep all orthologs after addition of '
                            'mapped seq, which means also the groups that '
                            'have no mapped sequence. Otherwise only groups '
                            'are used that have the mapped sequence for '
                            'alignment and tree inference.')

    arg_parser.add_argument('--check_mate_pairing', action='store_true',
                            help='Check whether in case of paired end '
                            'reads we have consistent mate pairing. Setting '
                            'this option will automatically select the '
                            'overlapping reads and do not consider single '
                            'reads.')

    arg_parser.add_argument('--debug', action='store_true',
                            help='Changes to debug mode: '
                                 '* bam files are saved!'
                                 '* reads are saved by mapping to OG')

    arg_parser.add_argument('-s', '--species_name', default=None,
                            help='[Default is name of read] Name of species '
                                 'for mapped sequence.')

    arg_parser.add_argument('--merge_all_mappings', action='store_true',
                            help='In case multiple species were mapped to '
                            'the same reference this allows to merge this '
                            'mappings and build a tree with all included '
                            'species!')

    # arg_parser.add_argument('--ref_og_aa_folder', default='.',
    #                         help='Path to preselected og_aa folder')

    # Arguments to generate the reference
    arg_parser.add_argument('-r', '--reference', action='store_true',
                            help='Just generate the reference dataset for '
                            'mapping.')

    arg_parser.add_argument('--min_species', type=int, default=None,
                            help='Min number of species in selected '
                            'orthologous groups. If not selected it will be '
                            'estimated such that around 1000 OGs '
                            'are available.')

    arg_parser.add_argument('--single_mapping', default=None,
                            help='Single species file allowing to map in a '
                            'job array.')

    arg_parser.add_argument('--threads', type=int, default=1,
                            help='Number of threads for the mapping '
                            'using ngm / ngmlr!')

    # Arguments to map the reads
    arg_parser.add_argument('--ref_folder', default=None,
                            help='Folder containing reference files with '
                            'sequences sorted by species.')

    # Parse the arguments.
    args = arg_parser.parse_args(argv)

    if not args.split_reads and (args.split_len != 400 or
                                 args.split_overlap != 0 or
                                 args.split_min_read_len != 500):
        arg_parser.error(
            'Arguments --split_len, --split_overlap and --split_min_read_len'
            'can only be set if --split_reads is set.')

    if args.split_reads and len(args.reads) == 2:
        arg_parser.error(
            'Splitting reads does not work for paired end reads.')

    if args.read_type is 'short' and args.ngmlr_parameters:
        arg_parser.error(
            'Arguments for --ngmlr_parameters only work if --read_type is set '
            'to "long".')

    if not args.sample_reads and (args.coverage != 10 or
                                  args.genome_len != 2000000):
        arg_parser.error(
            'Arguments --coverage and --genome_len'
            'can only be set if --sample_reads is set.')

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

    progress = Progress(args)
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    # TODO: Check whether all the necessary binaries are available
    # TODO: Check all given files and throw error if faulty

    if args.species_name:
        species_name = args.species_name
    elif args.reads:
        species_name = args.reads[0].split("/")[-1].split(".")[0]
    else:
        species_name = 'merge'
    progress.get_status(species_name=species_name)

    if (not progress.ref_ogs_01 and not progress.ref_dna_02 and
        not progress.mapping_03 and not progress.append_ogs_04 and
            not progress.align_05):
        oma_output = OMAOutputParser(args)
        args.oma_output_path = oma_output.oma_output_path
        ogset = OGSet(args, oma_output=oma_output)  # Generate the OGs with their DNA sequences
        reference = ReferenceSet(args, og_set=ogset.ogs, load=True)
        alignments = Aligner(args, ogset.ogs, load=True)
        if not args.reference:
            mapper = Mapper(args, og_set=ogset.ogs, ref_set=reference.ref)
            alignments.add_mapped_seq(mapper)
            alignments.write_added_aligns_aa()
            alignments.write_added_aligns_dna()
            # progress.set_status("re_ogs")
            # alignments = Aligner(args, ogset.mapped_ogs, load=True)
            progress.set_status("og_align")
            concat_alignment = alignments.concat_alignment()
            #tree = TreeInference(args, concat_alignment=concat_alignment[0])
            #print(tree.tree)
    # elif (progress.ref_ogs_01 and not progress.ref_dna_02 and
    #       not progress.mapping_03 and not progress.append_ogs_04 and
    #         not progress.align_05):
    #     ogset = OGSet(args, load=False)
    #     reference = ReferenceSet(args, og_set=ogset.ogs, load=True)  # Generate the reference
    #     if not args.reference:  # just generate reference
    #         mapper = Mapper(args, og_set=ogset.ogs, ref_set=reference.ref)
    #         ogset.add_mapped_seq(mapper)
    #         ogset.write_added_ogs_aa()
    #         ogset.write_added_ogs_dna()
    #         progress.set_status("re_ogs")
    #         alignments = Aligner(args, ogset.mapped_ogs, load=True)
    #         progress.set_status("og_align")
    #         concat_alignment = alignments.concat_alignment()
    #         #tree = TreeInference(args, concat_alignment=concat_alignment[0])
    #         #print(tree.tree)
    # elif (progress.ref_ogs_01 and progress.ref_dna_02 and
    #       not progress.mapping_03 and not progress.append_ogs_04 and
    #         not progress.align_05):
    #     if args.single_mapping:
    #         reference = ReferenceSet(args, load=False)
    #         Mapper(args, ref_set=reference.ref)  # Run the mapping
    #     else:
    #         ogset = OGSet(args, load=False)
    #         reference = ReferenceSet(args, load=False)
    #         mapper = Mapper(args, og_set=ogset.ogs, ref_set=reference.ref)  # Run the mapping
    #         ogset.add_mapped_seq(mapper)
    #         ogset.write_added_ogs_aa()
    #         ogset.write_added_ogs_dna()
    #         progress.set_status("re_ogs")
    #         alignments = Aligner(args, ogset.mapped_ogs, load=True)
    #         progress.set_status("og_align")
    #         concat_alignment = alignments.concat_alignment()
    #         #tree = TreeInference(args, concat_alignment=concat_alignment[0])
    #         #print(tree.tree)
    # elif (progress.ref_ogs_01 and progress.ref_dna_02 and
    #       progress.mapping_03 and not progress.append_ogs_04 and
    #         not progress.align_05):
    #     ogset = OGSet(args, load=False)
    #     if not args.merge_all_mappings:
    #         mapper = Mapper(args, og_set=ogset.ogs, load=False)
    #         ogset.add_mapped_seq(mapper)
    #         ogset.write_added_ogs_aa()
    #         ogset.write_added_ogs_dna()
    #         progress.set_status("re_ogs")
    #         alignments = Aligner(args, ogset.mapped_ogs, load=True)
    #         progress.set_status("og_align")
    #         concat_alignment = alignments.concat_alignment()
    #         #tree = TreeInference(args, concat_alignment=concat_alignment[0])
    #         #print(tree.tree)
    #     else:
    #         for folder in glob.glob(os.path.join(args.output_path,
    #                                              "03_mapping_*")):
    #             species_name = folder.split("03_mapping_")[-1]
    #             species_progress = Progress(args)
    #             species_progress.get_status(species_name=species_name)
    #             if species_progress.mapping_03:
    #                 print('--- Addition of {} to all ogs '
    #                       '---'.format(species_name))
    #                 mapper = Mapper(args, og_set=ogset.ogs,
    #                                 species_name=species_name, load=False)
    #                 ogset.add_mapped_seq(mapper, species_name=species_name)
    #         ogset.write_added_ogs_aa(folder_name="04_merge_OGs_aa")
    #         ogset.write_added_ogs_dna(folder_name="04_merge_OGs_dna")
    #         progress.set_status("re_ogs")
    #         alignments = Aligner(args, ogset.mapped_ogs, load=True)
    #         progress.set_status("og_align")
    #         concat_alignment = alignments.concat_alignment()
    #         #tree = TreeInference(args, concat_alignment=concat_alignment[0])
    #         #print(tree.tree)
    # elif (progress.ref_ogs_01 and progress.ref_dna_02 and
    #       progress.mapping_03 and progress.append_ogs_04 and
    #         not progress.align_05):
    #     ogset = OGSet(args, load=False)
    #     alignments = Aligner(args, ogset.mapped_ogs, load=True)
    #     progress.set_status("og_align")
    #     concat_alignment = alignments.concat_alignment()
    #     #tree = TreeInference(args, concat_alignment=concat_alignment[0])
    #     #print(tree.tree)
    # elif (progress.ref_ogs_01 and progress.ref_dna_02 and
    #       progress.mapping_03 and progress.append_ogs_04 and
    #         progress.align_05):
    #     alignments = Aligner(args,  load=False)
    #     progress.set_status("og_align")
    #     concat_alignment = alignments.concat_alignment()
    #     #tree = TreeInference(args, concat_alignment=concat_alignment[0])
    #     #print(tree.tree)

    print('Time taken {}'.format(timer() - t1))
