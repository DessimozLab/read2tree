#!/usr/bin/env python
'''
    This file contains definitions of a class that allows to map reads to a
    reference file!
    Importantly the mapper function relies heavy on the following:
    ngm v0.5.4
    ngmlr v0.2.6
    samtools
    bcftools
    vcfutils.pl


    -- David Dylus, July--XXX 2017
'''

import pysam
import tempfile
import os
import shutil
import glob
import subprocess
import logging
import time
from tqdm import tqdm

import functools
from Bio import SeqIO, SeqRecord, Seq
from Bio.Alphabet import generic_dna
from Bio.SeqIO.FastaIO import FastaWriter
from read2tree.OGSet import OG
from read2tree.Reads import Reads
from read2tree.ReferenceSet import Reference
from read2tree.wrappers.read_mappers import NGM
from read2tree.wrappers.read_mappers import NGMLR
from read2tree.Progress import Progress
from read2tree.stats.Coverage import Coverage
from read2tree.stats.SeqCompleteness import SeqCompleteness


logger = logging.getLogger(__name__)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('info.log')
file_handler.setFormatter(formatter)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)


# def sizeof(obj):
#     """APPROXIMATE memory taken by some Python objects in
#     the current 32-bit CPython implementation.
#
#     Excludes the space used by items in containers; does not
#     take into account overhead of memory allocation from the
#     operating system, or over-allocation by lists and dicts.
#     """
#     T = type(obj)
#     if T is int:
#         kind = "fixed"
#         container = False
#         size = 4
#     elif T is list or T is tuple:
#         kind = "variable"
#         container = True
#         size = 4*len(obj)
#     elif T is dict:
#         kind = "variable"
#         container = True
#         size = 144
#         if len(obj) > 8:
#             size += 12*(len(obj)-8)
#     elif T is str:
#         kind = "variable"
#         container = False
#         size = len(obj) + 1
#     else:
#         raise TypeError("don't know about this kind of object")
#     if kind == "fixed":
#         overhead = 8
#     else:  # "variable"
#         overhead = 12
#     if container:
#         garbage_collector = 8
#     else:
#         garbage_collector = 0
#     malloc = 8  # in most cases
#     size = size + overhead + garbage_collector + malloc
#     # Round to nearest multiple of 8 bytes
#     x = size % 8
#     if x != 0:
#         size += 8-x
#         size = (size + 8)
#     return size


class Mapper(object):
    """
    Structure for reference
    """

    def __init__(self, args, ref_set=None, og_set=None, species_name=None,
                 load=True):
        self.args = args
        self.elapsed_time = 0

        if args.debug:
            logger.setLevel(logging.DEBUG)
            file_handler.setLevel(logging.DEBUG)
            # stream_handler.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)
            file_handler.setLevel(logging.INFO)
            # stream_handler.setLevel(logging.INFO)

        logger.addHandler(file_handler)
        # logger.addHandler(stream_handler)

        if self.args.reads:
            if len(self.args.reads) == 2:
                self._reads = self.args.reads
                self._species_name = \
                    self._reads[0].split("/")[-1].split(".")[0]
            else:
                self._reads = self.args.reads[0]
                self._species_name = self._reads.split("/")[-1].split(".")[0]

        if self.args.species_name:
            self._species_name = self.args.species_name

        if species_name:
            self._species_name = species_name

        # load pyopa related stuff
        # self.defaults = pyopa.load_default_environments()
        # self.envs = self.defaults['environments']
        # self.env = self.envs[515]
        self.progress = Progress(args)
        self.all_cov = {}
        self.all_sc = {}

        self.read_og_set = {}

        if load:  # compute mapping
            if ref_set is not None:
                self.mapped_records = \
                    self._map_reads_to_references(ref_set)
                if self.progress.get_mapping_status():
                    self.progress.set_status('map')
            if self.mapped_records and og_set is not None:
                self.og_records = self._sort_by_og(og_set)
        else:  # re-load already computed mapping
            if og_set is not None and not self.args.merge_all_mappings:
                self.mapped_records = self._read_mapping_from_folder()
                self.og_records = self._sort_by_og(og_set)
            elif (og_set is not None and
                  self.args.merge_all_mappings and species_name is not None):
                self.mapped_records = \
                    self._read_mapping_from_folder(species_name=species_name)
                self.og_records = self._sort_by_og(og_set)

    def _call_wrapper(self, ref_file_handle, reads, tmp_output_folder):

        if len(self._reads) is 2:
            ngm_wrapper = NGM(ref_file_handle, reads, tmp_output_folder.name)
            if self.args.threads is not None:
                ngm_wrapper.options.options['-t'].set_value(self.args.threads)
            ngm = ngm_wrapper()
            bam_file = ngm['file']
        elif len(self._reads) is not 2 and 'short' in self.args.read_type:
            ngm_wrapper = NGM(ref_file_handle, reads, tmp_output_folder.name)
            if self.args.threads is not None:
                ngm_wrapper.options.options['-t'].set_value(self.args.threads)
            ngm = ngm_wrapper()
            bam_file = ngm['file']
        elif len(self._reads) is not 2 and 'long' in self.args.read_type:
            ngm_wrapper = NGMLR(ref_file_handle, reads, tmp_output_folder.name)
            if self.args.threads is not None:
                ngm_wrapper.options.options['-t'].set_value(self.args.threads)
            if self.args.ngmlr_parameters is not None:
                par = self.args.ngmlr_parameters.split(',')
                ngm_wrapper.options.options['-x'].set_value(str(par[0]))
                ngm_wrapper.options \
                           .options['--subread-length'].set_value(int(par[1]))
                ngm_wrapper.options.options['-R'].set_value(float(par[2]))
            ngm = ngm_wrapper()
            bam_file = ngm['file']
        logger.info('{}: Mapped {} / {} reads '
                    ''.format(self._species_name, ngm['reads_mapped'],
                              ngm['total_reads']+ngm['reads_mapped']))
        self._rm_file(ref_file_handle + "-enc.2.ngm", ignore_error=True)
        self._rm_file(ref_file_handle + "-ht-13-2.2.ngm", ignore_error=True)
        self._rm_file(ref_file_handle + "-ht-13-2.3.ngm", ignore_error=True)

        return self._post_process_read_mapping(ref_file_handle, bam_file)

    def _read_mapping_from_folder(self, species_name=None):
        """
        Retrieve all the mapped consensus files from folder and add to mapper
        object
        :return: dictionary with key og_name and value sequences mapped to each
                 species
        """
        print('--- Retrieve mapped consensus sequences ---')
        map_reads_species = {}
        if not species_name:
            species_name = self._species_name
        in_folder = os.path.join(self.args.output_path,
                                 "03_mapping_"+species_name)
        for file in tqdm(glob.glob(os.path.join(in_folder, "*_consensus.fa")),
                         desc='Loading consensus read mappings ',
                         unit=' species'):
            species = file.split("/")[-1].split("_")[0]
            map_reads_species[species] = Reference()
            map_reads_species[species].dna = list(SeqIO.parse(file, "fasta"))

            cov = Coverage()
            cov_file_name = os.path.join(in_folder, species + "_OGs_cov.txt")
            for line in open(cov_file_name, "r"):
                if "#" not in line:
                    values = line.split(",")
                    cov.add_coverage(values[2]+"_"+values[1],
                                     [float(values[3]),
                                      float(values[4].replace("\n", ""))])
            self.all_cov.update(cov.coverage)

            seqC = SeqCompleteness()
            seqC_file_name = os.path.join(in_folder, species + "_OGs_sc.txt")
            for line in open(seqC_file_name, "r"):
                if "#" not in line:
                    values = line.split(",")
                    seqC.add_seq_completeness(values[2] + "_" + values[1],
                                              [float(values[3]),
                                               float(values[4]),
                                               int(values[5]),
                                               int(values[6]),
                                               int(values[7].replace("\n",
                                                                     ""))])
            self.all_sc.update(seqC.seq_completeness)

        return map_reads_species

    def _make_tmpdir(self):
        '''
        Make tmpdir for analysis. This is important to run the code on the
        cluster as usually scratch space can be limited and computation is
        accelerated when parsing large files from the node drive.

        '''
        try:
            tmp_output_folder = tempfile.TemporaryDirectory(
                prefix='ngm', dir=os.environ.get("TMPDIR"))
        except NotADirectoryError:
            logger.info('{}: Environmental variable TMPDIR not set, will use \
                        native python tmpdir location.'
                        .format(self._species_name))
        else:
            tmp_output_folder = tempfile.TemporaryDirectory(prefix='ngm_')
            logger.debug('--- Creating tmp directory on local node ---')
        return tmp_output_folder

    def _map_reads_to_references(self, ref):
        """
        Map all the reads to reference
        :param reference:
        :return: dictionary with key og_name and value sequences
                 mapped to each species
        """
        start = time.time()
        print('--- Mapping of reads to reference sequences ---')
        mapped_reads_species = {}
        reference_path = os.path.join(self.args.output_path, "02_ref_dna")

        output_folder = os.path.join(self.args.output_path,
                                     "03_mapping_"+self._species_name)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        tmp_output_folder = self._make_tmpdir()

        # Get reads and perform splitting if requested
        read_container = Reads(self.args)
        reads = read_container.reads
        # print(os.path.getsize(reads[0]))

        if self.args.single_mapping:
            references = [self.args
                          .single_mapping.split("/")[-1].split("_")[0]]
        else:
            references = list(ref.keys())

        # Going through provided references and starting mapping
        for species in tqdm(references,
                            desc='Mapping reads to species',
                            unit=' species'):
            logger.info('{}: --- Mapping of reads to {} reference species '
                        '---'.format(self._species_name, species))

            # write reference into temporary file
            ref_file_handle = os.path.join(reference_path, species+'_OGs.fa')
            ref_tmp_file_handle = os.path.join(tmp_output_folder.name,
                                               species + '_OGs.fa')
            shutil.copy(ref_file_handle, ref_tmp_file_handle)

            # call the WRAPPER here
            processed_reads = self._call_wrapper(ref_tmp_file_handle, reads,
                                                 tmp_output_folder)

            # postprocess mapping and build consensus
            try:
                mapped_reads = list(SeqIO.parse(processed_reads, 'fasta'))
                mapped_reads_species[species] = Reference()
                mapped_reads_species[species].dna = mapped_reads

                # save some general statistics for mapped gene
                if self.args.remove_species_ogs:
                    seqC = SeqCompleteness(mapped_ref=ref[species].dna,
                                           tested_ref=ref[self.args.remove_species_ogs].dna)
                else:
                    seqC = SeqCompleteness(mapped_ref=ref[species].dna)
                seqC.get_seq_completeness(mapped_reads)
                seqC.write_seq_completeness(os
                                            .path.join(output_folder,
                                                       species+"_OGs_sc.txt"))
                self.all_sc.update(seqC.seq_completeness)
            except AttributeError or ValueError:
                logger.debug('Reads not properly processed for further steps.')
            else:
                mapped_reads = []

            # self.progress.set_status('single_map', ref=species)
            self._rm_file(ref_file_handle+".fai", ignore_error=True)

        tmp_output_folder.cleanup()
        end = time.time()
        self.elapsed_time = end - start
        logger.info('{}: Mapping to {} references took {}.'
                    .format(self._species_name, len(references),
                            self.elapsed_time))
        return mapped_reads_species

    def _write_read_query_aling(self, read, og_name_file, write_mode):
        """

        :param read: pysam read object
        :param og_name_file: filename to collect reads
        :param write_mode: either a or wt
        """
        if read.is_paired:
            if read.is_read1:
                read_id = read.query_name + '/1'
            elif read.is_read2:
                read_id = read.query_name + '/2'
        else:
            read_id = read.query_name
        record = "@" + read_id + "\n"
        record += read.query_alignment_sequence + "\n"
        record += "+" + read_id + "\n"
        record += read.qqual + "\n"
        with open(og_name_file, write_mode) as f:
            f.write(record)

    def _write_read_full(self, read, og_name_file, write_mode):
        """

        :param read: pysam read object
        :param og_name_file: filename to collect reads
        :param write_mode: either a or wt
        """
        if read.is_paired:
            if read.is_read1:
                read_id = read.query_name + '/1'
            elif read.is_read2:
                read_id = read.query_name + '/2'
        else:
            read_id = read.query_name
        record = "@" + read_id + "\n"
        record += read.seq + "\n"
        record += "+" + read_id + "\n"
        record += read.qual + "\n"
        with open(og_name_file.replace(".fq", "_full.fq"), write_mode) as f:
            f.write(record)

    def _get_mapping_stats(self, bam_file):
        """
        Function that bins reads into their orthologous groups
        :param ref_file: Current species reference file
        :param bam_file: Mapped bam file
        """
        sfile_idx = pysam.idxstats(bam_file)
        x = sfile_idx.rstrip('\n').split('\n')
        mapped = functools.reduce(lambda x, y: x + y,
                                  [eval('+'.join(l.split('\t')[2]))
                                   for l in x])
        all_reads = functools.reduce(lambda x, y: x + y,
                                     [eval('+'.join(l.split('\t')[2:]))
                                      for l in x])
        return mapped, all_reads

    def _bin_reads(self, ref_file, bam_file):
        """
        Function that bins reads into their  orthologous groups
        :param ref_file: Current species reference file
        :param bam_file: Mapped bam file
        """
        logger.debug("--- Binning reads to {} ---".format(self._species_name))
        output_folder = os.path.join(self.args.output_path, "03_read_ogs_" +
                                     self._species_name)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        tmp_folder = os.path.dirname(bam_file)
        outfile_name = os.path.join(tmp_folder,
                                    ref_file.split('/')[-1].split('.')[0] +
                                    "_post")

        shutil.copy(bam_file, os.path.join(output_folder,
                                           os.path.basename(outfile_name +
                                                            "_sorted.bam")))
        shutil.copy(bam_file + ".bai",
                    os.path.join(output_folder,
                                 os.path.basename(outfile_name +
                                                  "_sorted.bam.bai")))

        if os.path.exists(bam_file):
            bam = pysam.AlignmentFile(bam_file, "rb")
            for read in bam.fetch():
                if not read.is_unmapped:
                    og_name = read.reference_name.split("_")[-1]
                    og_name_file = os.path.join(output_folder, og_name + ".fq")
                    if os.path.exists(og_name_file):
                        og = list(SeqIO.parse(og_name_file, 'fastq'))
                        og_read_ids = [rec.id for rec in og]
                        write_mode = 'a+'
                    else:
                        og_read_ids = []
                        write_mode = 'wt'

                    read_id = read.query_name
                    if read_id not in og_read_ids:
                        self._write_read_query_aling(read, og_name_file,
                                                     write_mode)
                        self._write_read_full(read, og_name_file, write_mode)

            bam.close()

    def _most_common(self, lst):
        return max(set(lst), key=lst.count)

    def _build_consensus_seq(self, ref_file, bam_file):
        """
        Function to build consensus sequence by taking sequence to be mapped
        :param ref_file:
        :param bam_file:
        :return:
        """
        bam = pysam.AlignmentFile(bam_file)
        records = {rec.id: rec for rec in list(SeqIO.parse(ref_file, "fasta"))}
        new_records = {}
        for read in bam.fetch():
            #     logger.info(read.qual)
            read_seq = list(read.seq)
            pairs = read.get_aligned_pairs(matches_only=True, with_seq=True)
            if read.reference_name in new_records.keys():
                seq = list(new_records[read.reference_name])
            else:
                seq = list('N' * len(records[read.reference_name]))
            for tup in pairs:
                seq[tup[1]] = read_seq[tup[0]]

            new_records[read.reference_name] = ("").join(seq)
        return new_records

    def _build_consensus_seq_v2(self, ref_file, bam_file):
        """
        Function to build consensus sequence by taking sequence to be mapped
        :param ref_file:
        :param bam_file:
        :return:
        """
        bam = pysam.AlignmentFile(bam_file)
        references = list(set([read.reference_name for read in bam.fetch()]))
        records = {rec.id: rec for rec in list(SeqIO.parse(ref_file, "fasta"))}
        new_records = {}
        for ref in references:
            #     logger.info(read.qual)
            seq = list('N' * len(records[ref]))
            for pileup_column in bam.pileup(ref, 0, 100000):
                # TODO: improve the selection of a column by its quality
                # qualities = [pileupread.alignment.query_alignment_qualities[pileupread.query_position] for pileupread in
                #        pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip]
                bases = [pileupread.alignment.query_sequence
                         [pileupread.query_position] for pileupread in
                         pileup_column.pileups if not pileupread.is_del and
                         not pileupread.is_refskip]
                if bases:
                    seq[pileup_column.pos] = self._most_common(bases)
            if len(set(seq)) > 1:  # make sure that mapped sequence contains not only N
                new_records[ref] = ("").join(seq)
        return new_records

    def _post_process_read_mapping(self, ref_file, bam_file):
        """
        Function that will perform postprocessing of finished read mapping
        using the pysam functionality
        :param ngm:
        :param reference_file_handle:
        :return:
        """
        # logger.info("--- Postprocessing reads to {} ---".format(self._species_name))
        output_folder = os.path.join(self.args.output_path,
                                     "03_mapping_"+self._species_name)
        tmp_folder = os.path.dirname(bam_file)
        outfile_name = os.path.join(tmp_folder,
                                    ref_file.split('/')[-1].split('.')[0] +
                                    "_post")
        if self.args.single_mapping:
            logger.debug("--- POSTPROCESSING MAPPING FOR {} "
                         "---".format(self._species_name.upper()))

        # ngmlr doesn't have the option to write in bam file directly
        if 'sam' in bam_file.split(".")[-1]:
            sam_file = bam_file
            bam_file = sam_file.replace(".sam", ".bam")
            if os.path.exists(sam_file):
                self._output_shell(
                    'samtools view -F 4 -bh -S -@ ' + str(self.args.threads) +
                    ' -o ' + bam_file + " " + sam_file)
        if self.args.single_mapping:
            logger.debug("---- Samtools view completed for {}"
                         .format(self._species_name))

        if os.path.exists(bam_file):
            self._output_shell(
                'samtools sort -m 2G  -@ ' + str(self.args.threads) +
                ' -o ' + outfile_name + "_sorted.bam " + bam_file)
        if self.args.single_mapping:
            logger.debug("---- Samtools sort completed for {} "
                         .format(self._species_name))

        if os.path.exists(outfile_name + "_sorted.bam"):
            self._output_shell(
                'samtools index -@ ' + str(self.args.threads) + ' ' +
                outfile_name + "_sorted.bam")
            # self._output_shell('bedtools genomecov -bga -ibam ' + outfile_name + '_sorted.bam | grep -w 0$ > ' + outfile_name + "_sorted.bed")
        if self.args.single_mapping:
            logger.debug("---- Samtools index completed for {} "
                         .format(self._species_name))

        # self._rm_file(bam_file, ignore_error=True)
        if self.args.debug:
            self._bin_reads(ref_file, outfile_name + '_sorted.bam')

        consensus = self._build_consensus_seq_v2(ref_file, outfile_name +
                                                 '_sorted.bam')
        # mapped_reads, all_reads = self._get_mapping_stats(outfile_name +
        #                                                   '_sorted.bam')
        #
        # logger.info("---- Mapping of {} against {} with {} "
        #             "of {} mapped reads completed!---"
        #             .format(self._species_name, ref_file.split("_")[0],
        #                     mapped_reads, all_reads))

        if self.args.single_mapping:
            logger.debug("---- Number of mapped reads {} ")

        all_consensus = []
        if consensus:
            try:
                for key, value in consensus.items():
                    seq = Seq.Seq(value, generic_dna)
                    record = SeqRecord.SeqRecord(seq, id=key, description='')
                    all_consensus.append(record)
                handle = open(os.path.join(output_folder, ref_file.split("/")
                                           [-1].split(".")[0] +
                                           '_consensus.fa'), "w")
                writer = FastaWriter(handle, wrap=None)
                writer.write_file(all_consensus)
                handle.close()
            except ValueError:
                pass

        if os.path.exists(os.path.join(output_folder,
                                       ref_file.split("/")[-1].split(".")[0] +
                                       '_consensus.fa')):
            out_file = os.path.join(output_folder, ref_file.split("/")
                                    [-1].split(".")[0] + '_consensus.fa')
        else:
            out_file = None

        # Get effective coverage of each mapped sequence
        cov = Coverage()
        cov.get_coverage_bam(outfile_name + "_sorted.bam")
        cov.write_coverage_bam(os.path.join(
            output_folder, ref_file.split('/')[-1].split('.')[0] + "_cov.txt"))
        self.all_cov.update(cov.coverage)

        return out_file

    def _rm_file(self, *fns, ignore_error=False):
        for fn in fns:
            try:
                os.remove(fn)
            except FileNotFoundError:
                if not ignore_error:
                    raise

    def _clean_up_tmp_files_single(self, species):
        output_folder = os.path.join(self.args.output_path, "03_mapping_" +
                                     self._species_name)
        fn_ends = ('_post.bam', '_post_consensus_call.fq', '_post_sorted.bam',
                   '_post_sorted.bam.bai', '.fa.fai',
                   '.fa.sam', '.fa-ht-13-2.3.ngm', '.fa-ht-13-2.3.ngm',
                   '.fa', '.fa-enc.2.ngm')
        self._rm_file(*[os.path.join(output_folder, species + fn_end)
                        for fn_end in fn_ends], ignore_error=True)

    def _output_shell(self, line):
        """
        Save output of shell line that has pipes
        taken from: https://stackoverflow.com/questions/7389662/link-several-popen-commands-with-pipes
        :param line:
        :return:
        """
        try:
            shell_command = subprocess.Popen(
                line, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                shell=True)
        except OSError:
            return None
        except ValueError:
            return None

        (output, err) = shell_command.communicate()
        shell_command.wait()
        if shell_command.returncode != 0:
            logger.debug("Shell command failed to execute")
            logger.debug(line)
            return None

        return output

    def _sort_by_og(self, og_set):
        """
        Group the mapped sequences according to their OG name
        :return: dictionary with lists of records as values and og name as keys
        """
        og_records = {}
        for records in self.mapped_records.values():
            for record in records.dna:
                name = record.id.split("_")[-1]
                tmp_id = record.id
                species_name = tmp_id[0:5]
                record.description = tmp_id + " [" + species_name + "]"
                if name in og_records.keys():
                    og_records[name].dna.append(record)
                    aa = self._get_protein(record)
                    og_records[name].aa.append(aa)
                else:
                    og_records[name] = OG()
                    og_records[name].dna.append(record)
                    aa = self._get_protein(record)
                    og_records[name].aa.append(aa)
        return og_records

    def _predict_best_protein_position(self, record):
        """
        Given a list of sequences that are derived from mapped reads to
        multiple seq of a OG we find the best corresponding mapped seq by
        comparing it with a representative sequence of the original OG using
        pyopa local alignment and return the sequence with its highest score!
        :return:
        """
        try:
            frame = record.seq[0:].translate(
                table='Standard', stop_symbol='X', to_stop=False, cds=False)
            best_translation = \
                SeqRecord.SeqRecord(frame, id=self._species_name,
                                    description=record.description,
                                    name=record.name)
        except ValueError:
            raise ValueError("Problem with sequence format!")
        return best_translation

    def _get_protein(self, record):
        '''

        :param record: sequence record
        :return: best translation
        '''
        frame = record.seq[0:].translate(
            table='Standard', stop_symbol='X', to_stop=False, cds=False)
        best_translation = SeqRecord.SeqRecord(frame, id=record.id,
                                               description=record.description,
                                               name=record.name)
        return best_translation

    # def _predict_best_protein_pyopa(self, record, og):
    #     """
    #     Given a list of sequences that are derived from mapped reads to multiple seq of a OG
    #     we find the best corresponding mapped seq by comparing it with a representative sequence of the original OG using
    #     pyopa local alignment and return the sequence with its highest score!
    #     :return:
    #     """
    #     ref_og_seq = og.aa[0]
    #     s1 = pyopa.Sequence(str(ref_og_seq.seq))
    #     best_score = 0
    #     try:
    #         frames = [record.seq[i:].translate(table='Standard', stop_symbol='X', to_stop=False, cds=False) for i
    #                   in range(3)]
    #         best_seq_idx = 0
    #         for i, seq in enumerate(frames):
    #             s2 = pyopa.Sequence(str(seq))
    #             # calculating local and global scores for the given sequences
    #             local_double = pyopa.align_double(s1, s2, self.env)
    #             # logger.info('Local score: %f' % local_double[0])
    #             if local_double[0] > best_score:
    #                 best_score = local_double[0]
    #                 best_seq_idx = i
    #         best_translation = SeqRecord.SeqRecord(frames[best_seq_idx], id=self._species_name, description=record.description, name=record.name)
    #     except:
    #         raise ValueError("Problem with sequence format!", ref_og_seq.seq)
    #     return best_translation

    def write_by_og(self, output_folder):
        '''
        Write for each og all the mapped sequences into separate fasta
        files to a specified folder
        :param output_folder: folder where files should be stored
        '''
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        for key, value in tqdm(self.og_records.items(),
                               desc="Writing DNA seq sorted by OG",
                               unit=" OG"):
            handle = open(os.path.join(output_folder,
                                       'mapped_'+key+'.fa'), "w")
            writer = FastaWriter(handle, wrap=None)
            writer.write_file(value)
            handle.close()
