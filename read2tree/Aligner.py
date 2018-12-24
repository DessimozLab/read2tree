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

logger = logging.getLogger(__name__)

class Aligner(object):

    def __init__(self, args, og_set=None, load=True):

        self.args = args
        self.mapped_aligns = {}
        self.elapsed_time = 0

        self._reads = self.args.reads
        self._species_name = self.args.species_name

        if self.args.remove_species_ogs:
            self.species_to_remove_ogs = self.args \
                .remove_species_ogs.split(",")
        else:
            self.species_to_remove_ogs = []

        self.alignments = Alignment()

        if load and og_set is not None:
            print('--- Alignment of {} OGs ---'.format(len(list(og_set.keys()))))
            self._og_set = og_set
            self.alignments = self._align(og_set)
        else:
            self.alignments = self._reload_alignments_from_folder()

        # print(self._get_codon_dict_og(og_set))

    def _remove_species_from_alignment(self, current_align):
        """
                Removes sequence records for a species / set of
                species from the reference OGSet
                :param current_og: The current OG object
                :return: OG object with removed species
                """
        if self.args.remove_species_ogs:
            align = Alignment()
            filtered_align = current_align \
                .remove_species_records(self.species_to_remove_ogs)
            if filtered_align:
                align.dna = filtered_align[0]
                align.aa = filtered_align[1]
        else:
            align = current_align
        return align

    def _make_output_path(self, prefix):
        path = os.path.join(self.args.output_path, prefix)
        if not os.path.exists(path):
            os.makedirs(path)
        return path

    def _get_placement_dic(self, alignment, ref_species):
        placement_dic = {}
        ref_rec = None
        for r in alignment:
            if ref_species in r.id:
                ref_rec = r
                k = 0
                for i, aa in enumerate(list(r.seq)):
                    if "-" not in aa:
                        placement_dic[i] = k
                        k = k + 1
                    else:
                        placement_dic[i] = k
        return placement_dic, ref_rec

    def _add_mapseq_align(self, alignment,  map_record, ref_species, species_name):
        placement_dic, ref_rec = self._get_placement_dic(alignment, ref_species)
        new = [map_record[placement_dic[i]] if '-' not in v else '-' for i, v in enumerate(list(ref_rec.seq))]
        alignment.add_sequence(species_name, ''.join(new))
        return alignment

    def _get_species_id(self, record):
        """
        Sequences in OMA are marked by using the first three letters of genus
        and the first 2 letters of species, e.g. Amphiura filiformis = AMPFI.
        This however is not always the case and we therefore prioritize the id
        (e.g. MOUSE over MUSMU that comes from Mus musculus).
        :param description: SeqRecord description
        :return: species_id
        """
        # TODO: add extension for model identifiers
        # model_identifiers_oma = {'MUSMU': 'MOUSE', 'HOMSA': 'HUMAN',
        #                          'SARCE': 'YEAST'}

        sp_description = record.description
        species = sp_description[sp_description.find("[") +
                                     1:sp_description.find("]")]
        if species:  # [Mus musculus]
            if len(species.split(" ")) > 1:
                new_id = species.split(" ")[0][0:3] + \
                    species.split(" ")[1][0:2]
                return new_id.upper()
            else:  # [MUSMU]
                return species

    def add_mapped_seq(self, ogset_add, species_name=None):
        """
        Add the sequence given from the read mapping to its corresponding
        OG and retain
        all OGs that do not have the mapped sequence, thus all original OGs
        are used for tree inference
        :param cons_og_set: set of ogs with its mapped sequences
        """
        start = time.time()
        num_append_seq = 0
        if not species_name:
            species_name = self._species_name
        print('--- Add inferred mapped sequence back to alignment ---')

        # iterate through all existing ogs
        for name_og, align in tqdm(self.alignments.items(),
                                   desc='Adding mapped seq to alignments', unit=' alignments'):
            # print(align.aa)
            align_filt = self._remove_species_from_alignment(align)
            # print(align_filt.aa)
            if len(align_filt.aa) >= 2:
                # get all species that are not mapped from original alignment
                if name_og in ogset_add.keys():
                    # find mapped records from appended records in OGSet
                    map_record_aa = [r for r in ogset_add[name_og].aa if species_name in r.id]
                    # print(map_records_aa)
                    map_record_dna = [r for r in ogset_add[name_og].dna if species_name in r.id]
                    if map_record_aa and map_record_dna:
                        ref_species = self._get_species_id(map_record_aa[0])
                        self.mapped_aligns[name_og] = Alignment()
                        self.mapped_aligns[name_og].aa = self._add_mapseq_align(align_filt.aa, map_record_aa[0], ref_species, species_name)
                        self.mapped_aligns[name_og].dna = self._add_mapseq_align(align_filt.dna, map_record_dna[0], ref_species, species_name)
                        num_append_seq = num_append_seq + 1
                    elif self.args.keep_all_ogs:
                        self.mapped_aligns[name_og] = Alignment()
                        self.mapped_aligns[name_og].aa = align_filt.aa
                        self.mapped_aligns[name_og].dna = align_filt.dna
                elif self.args.keep_all_ogs:
                    self.mapped_aligns[name_og] = Alignment()
                    self.mapped_aligns[name_og].aa = align_filt.aa
                    self.mapped_aligns[name_og].dna = align_filt.dna
        end = time.time()
        self.elapsed_time = end-start
        logger.info('{}: Appending {} reconstructed sequences to present Alignments '
                    'took {}.'
                    .format(self._species_name,
                            num_append_seq,
                            self.elapsed_time))

    def _get_codon_dict(self, aa_seq, dna_seq):
        """
        Function that returns codons per sequence for back translation of alignment
        :param aa_seq: amino acid sequence string
        :param dna_seq: dna sequence string
        :return: dictionary with codons, e.g 'M':'AGT'
        """
        codons = {i: [aa_seq[i], dna_seq[3 * i:3 * i + 3]] for i, x in enumerate(aa_seq)}
        # codons['-'] = '---'
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
                SeqRecord(Seq("".join(self._get_translated_seq_list(codon, str(rec.seq))), generic_dna), id=rec.id))
        return MultipleSeqAlignment(translated_seq)

    def _get_translated_seq_list(self, codon, sequence):
        k = 0
        new_seq = []
        for i, s in enumerate(sequence):
            if '-' in s:
                new_seq.append('---')
            else:
                new_seq.append(codon[k][-1])
                k = k + 1
        return new_seq

    def write_added_align_aa(self, folder_name=None):
        """

        :param self:
        :param folder_name:
        :return:
        """
        if folder_name is None:
            align_with_mapped_seq = self._make_output_path("06_align_" +
                                                         self._species_name +
                                                         "_aa")
        else:
            align_with_mapped_seq = self._make_output_path(folder_name)

        for name, value in self.alignments.items():
            if name in self.mapped_aligns.keys():
                output_file = os.path.join(align_with_mapped_seq, name + ".fa")
                self._write(output_file, self.mapped_aligns[name].aa)
            elif self.args.keep_all_ogs:
                output_file = os.path.join(align_with_mapped_seq, name + ".fa")
                self._write(output_file, value.aa)

    def write_added_align_dna(self, folder_name=None):
        """

        :param self:
        :param folder_name:
        :return:
        """
        if folder_name is None:
            align_with_mapped_seq = self._make_output_path("06_align_" +
                                                         self._species_name +
                                                         "_dna")
        else:
            align_with_mapped_seq = self._make_output_path(folder_name)

        for name, value in self.alignments.items():
            if name in self.mapped_aligns.keys():
                output_file = os.path.join(align_with_mapped_seq, name + ".fa")
                self._write(output_file, self.mapped_aligns[name].dna)
            elif self.args.keep_all_ogs:
                output_file = os.path.join(align_with_mapped_seq, name + ".fa")
                self._write(output_file, value.dna)

    def _align_worker(self, og_set):
        align_dict = {}
        output_folder_aa = os.path.join(
            self.args.output_path, "03_align_aa")
        output_folder_dna = os.path.join(
            self.args.output_path, "03_align_dna")
        for key, value in og_set.items():
            mafft_wrapper = Mafft(value.aa, datatype="PROTEIN")
            mafft_wrapper.options.options['--localpair'].set_value(True)
            mafft_wrapper.options.options['--maxiterate'].set_value(1000)
            alignment = mafft_wrapper()
            codons = self._get_codon_dict_og(value)
            align = Alignment()
            align.aa = alignment
            try:
                align.dna = self._get_translated_alignment(codons, alignment)
            except ValueError as v:
                logger.info('{} with error {}'.format(key, v))
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
            self.args.output_path, "03_align_aa")
        output_folder_dna = os.path.join(
            self.args.output_path, "03_align_dna")
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
            self.args.output_path, "03_align_aa")
        output_folder_dna = os.path.join(
            self.args.output_path, "03_align_dna")
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

        if self.mapped_aligns:
            use_alignments = self.mapped_aligns
        else:
            use_alignments = self.alignments
        alignments_aa = []
        alignments_dna = []
        for key, value in use_alignments.items():
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

    def _write(self, file, value):
        """
        Write output to fasta file
        :param file: file and location of outputfile
        :param value:
        :return:
        """
        output_handle = open(file, "w")
        AlignIO.write(value, output_handle, "phylip-relaxed")
        output_handle.close()


class Alignment(object):

    def __init__(self):
        self.aa = []
        self.dna = []

    def remove_species_records(self, species_to_remove):
        """
        Remove species from reference sequence set
        :param species_to_remove: list of species to be removed
        :param all_species: list of all species present in analysis
        """
        aa = MultipleSeqAlignment([record for i, record in enumerate(
            self.aa) if self._get_species_id(record) not in species_to_remove])
        dna = MultipleSeqAlignment([record for i, record in enumerate(
            self.dna) if self._get_species_id(record) not in species_to_remove])
        if len(aa) > 0 and len(dna) > 0:
            return [dna, aa]
        else:
            return None

    def _get_id_rec(self, record):
        parts = record.id.split('_')
        if len(parts) > 2:
            return record.id.split('_')[0]+'_'+record.id.split('_')[1]
        else:
            return record.id

    def _get_species_id(self, record):
        """
        Sequences in OMA are marked by using the first three letters of genus
        and the first 2 letters of species, e.g. Amphiura filiformis = AMPFI.
        This however is not always the case and we therefore prioritize the id
        (e.g. MOUSE over MUSMU that comes from Mus musculus).
        :param description: SeqRecord description
        :return: species_id
        """
        # TODO: add extension for model identifiers
        # model_identifiers_oma = {'MUSMU': 'MOUSE', 'HOMSA': 'HUMAN',
        #                          'SARCE': 'YEAST'}

        sp_id = record.id
        if sp_id[0:5].isalpha():  # >MUSMU
            return sp_id[0:5]
        else:
            sp_description = record.description
            species = sp_description[sp_description.find("[") +
                                     1:sp_description.find("]")]
            if species:  # [Mus musculus]
                if len(species.split(" ")) > 1:
                    new_id = species.split(" ")[0][0:3] + \
                             species.split(" ")[1][0:2]
                    return new_id.upper()
                else:  # [MUSMU]
                    return species
