#!/usr/bin/env python
'''
    This file contains definitions of a class which surrounds possible alignment methods

    -- David Dylus, July--XXX 2017
'''
import os
import glob
import time
import logging
#from tables import *
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

logger = logging.getLogger(__name__)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('info.log')
file_handler.setFormatter(formatter)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)


class Aligner(object):

    def __init__(self, args, og_set=None, load=True):

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
            self.alignments = self._align(og_set)
        else:
            self.alignments = self._reload_alignments_from_folder()

        # print(self._get_codon_dict_og(og_set))

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


class Alignment(object):

    def __init__(self):
        self.aa = None
        self.dna = None
