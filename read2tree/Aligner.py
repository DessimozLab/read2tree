#!/usr/bin/env python
'''
    This file contains definitions of a class which surrounds possible alignment methods

    -- David Dylus, July--XXX 2017
'''
import os
import glob
import time
import logging
from multiprocessing import Pool
from collections import ChainMap
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from tqdm import tqdm
from read2tree.wrappers.aligners import Mafft
from read2tree.utils.seq_utils import concatenate
from read2tree.stats.Coverage import Coverage
from read2tree.stats.SeqCompleteness import SeqCompleteness

logger = logging.getLogger('Aligner.py')
formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')


class Aligner(object):

    def __init__(self, args, og_set=None, load=True):

        self.args = args

        self.elapsed_time = 0

        # looging cannot be initalized with object due to multiprocessing and has to be done as it is here
        file_handler = logging.FileHandler(os.path.join(args.output_path,
                                                        'info.log'))
        file_handler.setFormatter(formatter)
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)

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
                self._species_name = self._reads[0].split("/")[-1].split(".")[0]
            else:
                self._reads = self.args.reads[0]
                self._species_name = self._reads.split("/")[-1].split(".")[0]

        if self.args.species_name:
            self._species_name = self.args.species_name

        if not self.args.reads and not self.args.species_name:
            self._species_name = 'merge'

        self.alignments = Alignment()

        if load and og_set is not None:
            print('--- Alignment of {} OGs ---'.format(len(list(og_set.keys()))))
            self._og_set = og_set
            self.alignments = self._align(og_set)
        else:
            self.alignments = self._reload_alignments_from_folder()

        # print(self._get_codon_dict_og(og_set))

    def add_mapped_seq(self, og, species_name=None):
        """
        Add the sequence given from the read mapping to its corresponding
        OG and retain
        all OGs that do not have the mapped sequence, thus all original OGs
        are used for tree inference
        :param cons_og_set: set of ogs with its mapped sequences
        """
        start = time.time()
        cons_og_set = mapper.og_records  # get sequences from mapping
        cov = Coverage()
        seqC = SeqCompleteness()
        if not species_name:
            species_name = self._species_name

        print('--- Add inferred mapped sequence back to OGs ---')

        # iterate through all existing ogs
        for name_og, align in tqdm(self.alignments.items(),
                                desc='Adding mapped seq to OG', unit=' OGs'):
            #og_sc = {rec.id: mapper.all_sc[rec.id]
            #         for rec in og_filt.aa if rec.id in mapper.all_sc.keys()}
            og = self._og_set[name_og]
            if len(align.aa) > 2:
                # continue only if OG is in mapped OGs
                if name_og in cons_og_set.keys():
                    cons_og_filt = cons_og_set[name_og]
                    og_sc = {cons_og_filt._get_id_rec(rec): mapper.all_sc[cons_og_filt._get_id_rec(rec)]
                             for rec in cons_og_filt.aa
                             if cons_og_filt._get_id_rec(rec) in mapper.all_sc.keys()}
                    if len(cons_og_filt.aa) >= 1:  # we had at least one mapped og even after removal
                        best_records = cons_og_filt \
                            .get_best_consensus_by_seq_completeness(
                                og_sc, threshold=self.args.sc_threshold)
                        if best_records:
                            best_record_aa = best_records[0]
                            best_record_dna = best_records[1]
                            best_record_dna._generate_seq_completeness(seqC, mapper,
                                                            og, best_record_dna)
                            cov.add_coverage(self._get_clean_id(best_record_aa),
                                             mapper.all_cov[self._get_clean_id(best_record_aa)])
                            best_record_aa.id = species_name
                            best_record_dna.id = species_name
                            self.aligned_mapped[name_og] = og_filt
                            all_id = [rec.id
                                      for rec in self.aligned_mapped[name_og].aa]
                            if best_record_aa.id not in all_id:  # make sure that repeated run doesn't add the same sequence multiple times at the end of an OG
                                self.aligned_mapped[name_og] \
                                    .aa.append(best_record_aa)
                                self.aligned_mapped[name_og] \
                                    .dna.append(best_record_dna)
                        else:  # case where no best_record_aa reported because it was smaller than the self.args.sc_threshold
                            self.aligned_mapped[name_og] = align
                    else:  # mapping had only one that we removed
                        if self.args.keep_all_ogs:
                            self.aligned_mapped[name_og] = align
                else:  # nothing was mapped to that og
                    if self.args.keep_all_ogs:
                        self.aligned_mapped[name_og] = align
            else:
                logger.debug('{} was left only with a single entry '
                             'and hence not used for further '
                             'processing'.format(name_og))

        cov.write_coverage_bam(os.path.join(self.args.output_path,
                                            species_name+'_all_cov.txt'))
        seqC.write_seq_completeness(os.path.join(self.args.output_path,
                                                 species_name+'_all_sc.txt'))
        end = time.time()
        self.elapsed_time = end-start
        logger.info('{}: Appending {} reconstructed sequences to present OG '
                    'took {}.'
                    .format(self._species_name,
                            len(list(cons_og_set.keys())),
                            self.elapsed_time))

    def _get_codon_dict(self, aa_seq, dna_seq):
        """
        Function that returns codons per sequence for back translation of alignment
        :param aa_seq: amino acid sequence string
        :param dna_seq: dna sequence string
        :return: dictionary with codons, e.g 'M':'AGT'
        """
        codons = {aa_seq[i]: dna_seq[3*i:3*i+3] for i, x in enumerate(aa_seq)}
        codons['-'] = '---'
        return codons

    def _get_codon_dict_og(self, og):
        """
        Function that computes all the codon dictionaries for one OG
        :param og: object of class OG part of OGSet
        :return: dictionary with codon dictionaries, e.g. 'MOUSE':{'M':'AGT',...}
        """
        return {i[0].id: self._get_codon_dict(str(i[0].seq), str(i[1].seq)) for i in zip(og.aa, og.dna)}

    def _get_translated_alignment(self, codons, alignment):
        """
        Function that computes the back translated alignment using the codon dictionary obtained by mapping sequence to
        its corresponding amino acid
        :param codons: dictionary of codons, e.g. 'MOUSE':{'M':'AGT',...}
        :param alignment: MultipleSeqAlignment coming from Biopython
        :return: MultipleSeqAlignment with DNA sequences
        """
        translated_seq = []
        for rec in alignment:
            codon = codons[rec.id]
            translated_seq.append(
                SeqRecord(Seq("".join([codon[s] for s in str(rec.seq)]), generic_dna), id=rec.id))
        return MultipleSeqAlignment(translated_seq)

    def _align_worker(self, og_set):
        align_dict = {}
        output_folder_aa = os.path.join(
            self.args.output_path, "05_align_" + self._species_name + "_aa")
        output_folder_dna = os.path.join(
            self.args.output_path, "05_align_" + self._species_name + "_dna")
        for key, value in og_set.items():
            mafft_wrapper = Mafft(value.aa, datatype="PROTEIN")
            mafft_wrapper.options.options['--localpair'].set_value(True)
            mafft_wrapper.options.options['--maxiterate'].set_value(1000)
            alignment = mafft_wrapper()
            codons = self._get_codon_dict_og(value)
            align = Alignment()
            align.aa = alignment
            align.dna = self._get_translated_alignment(codons, alignment)
            align_dict[key] = align

            og_name = key.split("/")[-1]
            output_handle = open(os.path.join(output_folder_aa, og_name + ".phy"), "w")
            AlignIO.write(align.aa, output_handle, "phylip-relaxed")

            output_handle2 = open(os.path.join(output_folder_dna, og_name + ".phy"), "w")
            AlignIO.write(align.dna, output_handle2, "phylip-relaxed")
        return align_dict

    def _chunkify(self, dic, n):
        key_chunks = [list(dic.keys())[i::n] for i in range(n)]
        return [{key: dic[key] for key in chunk} for chunk in key_chunks]

    def _align(self, og_set):
        """
        Function that computes the alignment of a set of OGs with appended mapping functions
        :param og_set: Object of class OGSet
        :return: alignment dictionary containing Alignment objects with aa and dna MSAs
        """
        # align_dict = {}
        start = time.time()
        self._adapt_id(og_set)

        output_folder_aa = os.path.join(
            self.args.output_path, "05_align_" + self._species_name + "_aa")
        output_folder_dna = os.path.join(
            self.args.output_path, "05_align_" + self._species_name + "_dna")
        if not os.path.exists(output_folder_aa):
            os.makedirs(output_folder_aa)
        if not os.path.exists(output_folder_dna):
            os.makedirs(output_folder_dna)

        og_chunks = self._chunkify(og_set, self.args.threads)
        p = Pool(self.args.threads)
        res_align = p.map(self._align_worker, og_chunks)
        align_dict = dict(ChainMap(*res_align))
        end = time.time()
        self.elapsed_time = end - start
        logger.info('{}: Alignment of {} OGs took {}.'.format(
                    self._species_name,
                    len(list(og_set.keys())),
                    self.elapsed_time))
        return align_dict

    def _reload_alignments_from_folder(self):
        """
        Function that reloads the alignments from the pre-computed folders
        :return: alignment dictionary containing Alignment objects with aa and dna MSAs
        """
        align_dict = {}
        output_folder_aa = os.path.join(
            self.args.output_path, "05_align_" + self._species_name + "_aa")
        output_folder_dna = os.path.join(
            self.args.output_path, "05_align_" + self._species_name + "_dna")
        for f in tqdm(zip(sorted(glob.glob(os.path.join(output_folder_aa, '*.phy'))), sorted(glob.glob(os.path.join(output_folder_dna, '*.phy')))), desc='Loading alignments ', unit=' Alignment'):
            og_name = os.path.basename(f[0]).split(".")[0]
            align_dict[og_name] = Alignment()
            align_dict[og_name].aa = AlignIO.read(f[0], format='phylip-relaxed')
            align_dict[og_name].dna = AlignIO.read(f[1], format='phylip-relaxed')
        return align_dict

    def _adapt_id(selfs, og_set):
        """
        Function that adapts all sequence ids to just use the first 5 letters
        :param og_set: Object of class OGSet
        :return: tuple of concatinated MSAs (aa, dna)
        """
        for key, value in og_set.items():
            for record in value.aa:
                s = record.id[0:5]
                record.id = s

    def concat_alignment(self):
        alignments_aa = []
        alignments_dna = []
        for key, value in self.alignments.items():
            alignments_aa.append(value.aa)
            alignments_dna.append(value.dna)
        concatination_aa = concatenate(alignments_aa)
        concatination_dna = concatenate(alignments_dna)

        if concatination_aa:
            align_output = open(os.path.join(self.args.output_path,
                                             "concat_" + self._species_name + "_aa.phy"), "w")
            AlignIO.write(concatination_aa, align_output, "phylip-relaxed")
            align_output.close()

        if concatination_dna:
            align_output = open(os.path.join(self.args.output_path,
                                             "concat_" + self._species_name + "_dna.phy"), "w")
            AlignIO.write(concatination_dna, align_output, "phylip-relaxed")
            align_output.close()

        return (concatination_aa, concatination_dna)

    def add_to_alignment(self, mapper):
        """
        Add the sequence given from the read mapping to its corresponding
        OG and retain
        all OGs that do not have the mapped sequence, thus all original OGs
        are used for tree inference
        :param cons_og_set: set of ogs with its mapped sequences
        """
        start = time.time()
        cons_og_set = mapper.og_records  # get sequences from mapping


class Alignment(object):

    def __init__(self):
        self.aa = None
        self.dna = None
