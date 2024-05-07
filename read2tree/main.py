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
import logging
from datetime import date
from timeit import default_timer as timer
import read2tree
from read2tree.OGSet import OGSet
from read2tree.ReferenceSet import ReferenceSet
from read2tree.Mapper import Mapper
from read2tree.Aligner import Aligner
# from read2tree.Progress import Progress
from read2tree.TreeInference import TreeInference
from read2tree.parser import OMAOutputParser
import argparse
import glob

import sys

COPYRIGHT = '(C) 2017-{:d} David V Dylus'.format(date.today().year)

# logger = logging.getLogger(__name__)
logger_level = "DEBUG"  # DEBUG INFO  # TRACE  DEBUG INFO  WARN  ERROR  FATAL
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)
if logger_level == "INFO":
    logger.setLevel(logging.INFO)


# logger.disabled = True

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

    arg_parser.add_argument('--output_path', default='.',
                            help='[Default is current directory] Path to '
                                 'output directory.')

    arg_parser.add_argument('--standalone_path', default='.', required=True,  # todo change name to marker_gene
                            help='[Default is current directory] Path to '
                                 'oma standalone directory.')

    arg_parser.add_argument('--reads', nargs='+', default=None,
                            help='[Default is none] Reads to be mapped to reference. If paired '
                                 'end add separated by space.')

    arg_parser.add_argument('--read_type', default='short',
                            help='[Default is short] Type of reads to '
                                 'use for mapping: short, long-hifi or long-ont corresponding to  sr,  map-hifi, or map-ont   in minimap2.')

    arg_parser.add_argument('--threads', type=int, default=1,
                            help='[Default is 1] Number of threads for the mapping ')

    # arg_parser.add_argument('--split_reads', action='store_true',
    #                         help='[Default is off] Splits reads as defined by split_len (200) '
    #                         'and split_overlap (0) parameters. ')

    # arg_parser.add_argument('--split_len', type=int, default=200,
    #                         help='[Default is 200] Parameter for selection of '
    #                         'read split length can only be used in combination'
    #                         'with with long read option. ')
    #
    # arg_parser.add_argument('--split_overlap', type=int, default=0,
    #                         help='[Default is 0] Reads are split with an '
    #                         'overlap defined by this argument.')

    # arg_parser.add_argument('--split_min_read_len', type=int, default=200,
    #                         help='[Default is 200] Reads longer than this '
    #                         'value are cut into smaller values as defined '
    #                         'by --split_len. ')

    arg_parser.add_argument('--sample_reads', action='store_true',
                            help='[Default is off] Splits reads as defined by split_len (200) '
                                 'and split_overlap (0) parameters. ')

    arg_parser.add_argument('--genome_len', type=int, default=2000000,
                            help='[Default is 2000000] Genome size in bp.')

    arg_parser.add_argument('--coverage', type=float, default=10,
                            help='[Default is 10] coverage in X. Only considered if --sample reads is selected.')

    arg_parser.add_argument('--min_cons_coverage', type=int, default=1,
                            help='[Default is 1] Minimum number of nucleotides at column.')

    arg_parser.add_argument('--dna_reference', default='',  # todo make it mandatory dna_reference no api
                            help='[Default is None] Reference file that contains nucleotide '
                                 'sequences (fasta, hdf5). If not given it will use'
                                 'the RESTapi and retrieve sequences '
                                 'from http://omabrowser.org directly. '
                                 'NOTE: internet connection required!')

    arg_parser.add_argument('--sc_threshold', type=float, default=0.25,
                            help='[Default is 0.25; Range 0-1] Parameter for '
                                 'selection of sequences from mapping by '
                                 'completeness compared to its reference sequence '
                                 '(number of ACGT basepairs vs length of sequence). '
                                 'By default, all sequences are selected.')

    # arg_parser.add_argument('--ngmlr_parameters', default=None, # todo this could be used for minimap2 options
    #                         help='[Default is none] In case this parameters '
    #                         'need to be changed all 3 values have to be '
    #                         'changed [x,subread-length,R]. The standard '
    #                         'is: ont,256,0.25. Possibilities for these '
    #                         'parameter can be found in the original '
    #                         'documentation of ngmlr.')

    # arg_parser.add_argument('--check_mate_pairing', action='store_true',
    #                         help='Check whether in case of paired end '
    #                         'reads we have consistent mate pairing. Setting '
    #                         'this option will automatically select the '
    #                         'overlapping reads and do not consider single '
    #                         'reads.')

    arg_parser.add_argument('--debug', action='store_true',  # todo make it active always  otherwise change it back
                            help='[Default is false] Changes to debug mode: '
                                 '* bam files are saved!'
                                 '* reads are saved by mapping to OG')

    arg_parser.add_argument('--sequence_selection_mode', default="sc",
                            help='[Default is sc] Possibilities are cov and cov_sc '
                                 'for mapped sequence.')

    arg_parser.add_argument('-s', '--species_name', default="",
                            help='[Default is name of read 1st file] Name of species '
                                 'for mapped sequence.')

    arg_parser.add_argument('--tree', action='store_true',
                            help='[Default is false] Compute tree, otherwise just '
                                 'output concatenated alignment!')

    arg_parser.add_argument('--step', default="all",
                            help='[Default is all  1marker 2map 3combine ')

    # arg_parser.add_argument('--merge_all_mappings', action='store_true',
    #                         help='[Default is off] In case multiple species were mapped to '
    #                         'the same reference this allows to merge this '
    #                         'mappings and build a tree with all included '
    #                         'species!')

    # Arguments to generate the reference
    # arg_parser.add_argument('-r', '--reference', action='store_true',
    #                         help='[Default is off] Just generate the reference dataset for '
    #                         'mapping.')

    arg_parser.add_argument('--min_species', type=int, default=None,
                            help='Min number of species in selected '
                                 'orthologous groups. If not selected it will be '
                                 'estimated such that around 1000 OGs '
                                 'are available.')

    # arg_parser.add_argument('--single_mapping', default=None,
    #                         help='[Default is none] Single species file allowing to map in a '
    #                         'job array.')

    # Arguments to map the reads
    # arg_parser.add_argument('--ref_folder', default=None,
    #                         help='[Default is none] Folder containing reference files with '
    #                         'sequences sorted by species.')

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

    arg_parser.add_argument('--keep_all_ogs', action='store_true', default=True,
                            help='[Default is on] Keep all orthologs after addition of '
                                 'mapped seq, which means also the OGs that '
                                 'have no mapped sequence. Otherwise only OGs '
                                 'are used that have the mapped sequence for '
                                 'alignment and tree inference.')

    arg_parser.add_argument('--ignore_species', default=None,
                            help='[Default is none] Ignores species part of '
                                 'the OMA standalone pipeline. Input is comma '
                                 'separated list without spaces, e.g. XXX,YYY,AAA.')

    # Parse the arguments.
    args = arg_parser.parse_args(argv)

    _reads = ""
    _species_name = ""

    if args.reads:
        # print(args.reads)
        if len(args.reads) == 2:
            _reads = args.reads
            _species_name = _reads[0].split("/")[-1].split(".")[0]
        else:
            _reads = args.reads[0]
            _species_name = _reads.split("/")[-1].split(".")[0]

    if args.species_name:
        _species_name = args.species_name

    if args.step == "3combine":  # todo why is needed?
        _species_name = 'merge'

    args.reads = _reads
    args.species_name = _species_name

    # if not args.split_reads and (args.split_len != 200 or
    #                              args.split_overlap != 0 or
    #                              args.split_min_read_len != 200):
    #     arg_parser.error(
    #         'Arguments --split_len, --split_overlap and --split_min_read_len'
    #         'can only be set if --split_reads is set.')

    # if args.split_reads and len(args.reads) == 2:
    #     arg_parser.error(
    #         'Splitting reads does not work for paired end reads.')

    # if args.read_type == 'short' and args.ngmlr_parameters:
    #     arg_parser.error(
    #         'Arguments for --ngmlr_parameters only work if --read_type is set '
    #         'to "long".')

    if not args.sample_reads and (args.coverage != 10 or
                                  args.genome_len != 2000000):
        arg_parser.error(
            'Arguments --coverage and --genome_len '
            'can only be set if --sample_reads is set.')

    # progress = Progress(args) # todo why calling Progress twice?
    # if progress.num_completed_mappings <= 1 and args.merge_all_mappings:
    #     arg_parser.error('The number of completed mappings ({}) is too '
    #                      'little to perform a merge.'.format(progress.num_completed_mappings))

    return args


def get_finished_mapping_folders2(output_path):
    mapping_folders_finished = []
    # num_expected_mappings = self._get_number_of_references()
    mapping_folders = [x for x in os.listdir(output_path) if '04' in x]
    # logging.debug("Number of mapping expected is " + str(num_expected_mappings) + " we are checking folders in  " + str(path))
    for folder in mapping_folders:
        # NOTE: we are calculating the number of completed mappings as the number of existing cov files,
        # because these are written even if the mapping step did not find any reads to map to a particular reference
        computed_cov = [f for f in
                        glob.glob(os.path.join(output_path,
                                               folder + '/*cov.txt'))]
        # it is finished if the number of generated coverage files is the same as the number of references
        # if self.args.merge_all_mappings:
        #     if num_expected_mappings >= len(computed_cov) and len(computed_cov) >= 1:
        #         mapping_folders_finished.append(folder)
        # elif self._species_name in folder:
        #     if (num_expected_mappings - len(computed_cov)) == 0:
        #         mapping_folders_finished.append(folder)
        if len(computed_cov) >= 1:
            mapping_folders_finished.append(folder)  # todo needs testing

    return mapping_folders_finished


def main(argv, exe_name, desc=''):
    '''
        Main function.
    '''
    from . import __version__ as r2t_version

    logger.info(' ------- Read2Tree version: {} -------'.format(r2t_version))
    # print(' ------- Read2Tree version: '+str(r2t_version)+' -------')

    t1 = timer()
    # Parse
    args = parse_args(argv, exe_name, desc)


    x = ', '.join("{!s}={!r}".format(key, val) for (key, val) in vars(args).items())  # todo why calling Progress twice?
    logger.info('{}: read2tree was run with: {}'.format(args.species_name, x))

    logger.info("Running read2tree in mode " + args.step)

    if os.path.exists(args.output_path):
        if args.step == "all" or args.step == "1marker":
            logger.error(
                "the output folder exist " + args.output_path + ". Since you are running read2tree in all mode, you need to specify output folder which will created by read2tree.")
            sys.exit()
    else:
        os.makedirs(args.output_path)

    if args.step == "all" or args.step == "2map":
        if not args.reads:
            logger.error("reads are not provided in mode " + str(args.step))
            sys.exit()
        if isinstance(args.reads,
                      list):  # in parse_args it is input read is converted to a list, if there are 2 - paired end
            for read_file in args.reads:
                if not os.path.isfile(read_file):
                    logger.error("read file doesn't exist " + read_file)
                    sys.exit()
        elif isinstance(args.reads, str):
            if not os.path.isfile(args.reads):
                logger.error("read file doesn't exist")
                sys.exit()

    if args.step == "all":
        logger.info('{}: ------- NEW RUN -------'.format(args.species_name))
        oma_output = OMAOutputParser(args)
        args.oma_output_path = oma_output.oma_output_path
        ogset = OGSet(args, oma_output=oma_output, step=args.step)  # Generate the OGs with their DNA sequences.   # Write 01_
        # ogset.ogs['OG1188079'].aa[0] = SeqRecord(seq=Seq('FLGMT ...
        # ogset.ogs['OG1188079'].dna[0] = SeqRecord(seq=Seq('TTT

        reference = ReferenceSet(args, og_set=ogset.ogs, step=args.step)                            #  write 02_
        # aa and dna of each species
        # reference.ref['MNELE'].aa[0] =  SeqRecord(seq=Seq('FLGM

        alignments = Aligner(args, ogset.ogs, step=args.step)  # multiple sequence alignment of OGs     #  write 03_
        # alignments.alignments['OG1008242'].aa  <<class 'Bio.Align.MultipleSeqAlignment'> instance (5 records of length 573) at 7f9189695c00>
        # alignments.alignments['OG1008242'].aa[0]  SeqRecord(seq=Seq('---------------MTDFDKL

        mapper = Mapper(args, og_set=ogset.ogs, ref_set=reference.ref, step=args.step)  # map reads onto OG
        alignments.remove_species_from_alignments()
        ogset.remove_species_from_ogs()
        ogset.add_mapped_seq(mapper)
        ogset.write_added_ogs_aa()
        ogset.write_added_ogs_dna()
        # alignments = Aligner(args, ogset.mapped_ogs, load=True)
        alignments.add_mapped_seq(ogset.mapped_ogs)
        alignments.write_added_align_aa()
        alignments.write_added_align_dna()
        concat_alignment = alignments.concat_alignment()
        if args.tree:
            tree = TreeInference(args, concat_alignment=concat_alignment[0])
            logger.info(str(tree.tree))

        logger.info(' ------- Read2Tree finished -*- -------')
        # print("done - all")


    if args.step == "1marker":
        logger.info('{}: ------- NEW RUN -------'.format(args.species_name))
        oma_output = OMAOutputParser(args)
        args.oma_output_path = oma_output.oma_output_path
        ogset = OGSet(args, oma_output=oma_output, step=args.step)  # Generate the OGs with their DNA sequences
        reference = ReferenceSet(args, og_set=ogset.ogs, step=args.step)
        alignments = Aligner(args, ogset.ogs, step=args.step)
        print("done- 1marker")
        logger.info(' ------- Read2Tree step 1marker finished -*- -------')

    if args.step == "2map":
        logger.info('{}: ------- Read2tree RUN step 2map (after running 1marker) -------'.format(args.species_name))
        oma_output = OMAOutputParser(args)
        args.oma_output_path = oma_output.oma_output_path
        ogset = OGSet(args, oma_output=oma_output, step=args.step)  # Generate the OGs with their DNA sequences
        reference = ReferenceSet(args, og_set=ogset.ogs, step=args.step)
        alignments = Aligner(args, ogset.ogs, step=args.step)  # multiple sequence alignment of OGs
        mapper = Mapper(args, og_set=ogset.ogs, ref_set=reference.ref, step=args.step)
        alignments.remove_species_from_alignments()
        ogset.remove_species_from_ogs()
        ogset.add_mapped_seq(mapper)
        ogset.write_added_ogs_aa()
        ogset.write_added_ogs_dna()
        # alignments = Aligner(args, ogset.mapped_ogs, load=True)
        # alignments.add_mapped_seq(ogset.mapped_ogs)
        # alignments.write_added_align_aa()
        # alignments.write_added_align_dna()

        print("done- 2map")
        logger.info(' ------- Read2Tree step 2map finished -*- -------')

    if args.step == "3combine":
        ogset = OGSet(args, step=args.step)
        reference = ReferenceSet(args, step=args.step)
        alignments = Aligner(args, step=args.step)  # og set is not input of this, becuase we only read the aligned seq
        alignments.remove_species_from_alignments()
        ogset.remove_species_from_ogs()
        mapping_folders_finished = get_finished_mapping_folders2(args.output_path)
        for mapping in mapping_folders_finished:
            species_name = mapping.split("04_mapping_")[-1]
            logger.info('--- Addition of {} to all ogs '
                        '---'.format(species_name))
            mapper = Mapper(args, og_set=ogset.ogs, ref_set=reference.ref,
                            species_name=species_name, step=args.step)
            ogset.add_mapped_seq(mapper, species_name=species_name)
            alignments.add_mapped_seq(ogset.mapped_ogs, species_name=species_name)
        ogset.write_added_ogs_aa(folder_name="05_merge_OGs_aa")
        ogset.write_added_ogs_dna(folder_name="05_merge_OGs_dna")
        alignments.write_added_align_aa()
        alignments.write_added_align_dna()
        concat_alignment = alignments.concat_alignment()
        if args.tree:
            tree = TreeInference(args, concat_alignment=concat_alignment[0])
            logger.info(str(tree.tree))

        print("done- 3combine")
        logger.info(' ------- Read2Tree step 3combine finished -*- -------')

        logger.info(' ------- Read2Tree finished -*- -------')

    print("done- main ")

    # TODO: Check whether all the necessary binaries are available
    # TODO: Check all given files and throw error if faulty

    logger.info('{}: Time taken {}'.format(args.species_name, timer() - t1))
    # todo remove sam file