#!/usr/bin/env python
"""
    This file contains definitions of a class which holds the orthologous groups.

    -- David Dylus, July--XXX 2017
"""
import glob
import os
import re
import pyham
import requests
import logging
import random
import time
import numpy as np
import gzip

from tqdm import tqdm
from collections import OrderedDict
from Bio import SeqIO, Seq, SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
#from tables import *
# ----------- only to be used internally; requires hdf5 installation -------------------
#from pyoma.browser import db

from . import __version__ as read2tree_version
from read2tree.stats.Coverage import Coverage
from read2tree.stats.SeqCompleteness import SeqCompleteness
from read2tree.FastxReader import FastxReader

API_URL = 'http://omabrowser.org/api'


class OGSet(object):

    def __init__(self, args, oma_output=None, load=True, progress=None):
        self.args = args

        self.logger = logging.getLogger(__name__)
        self.mapped_ogs = {}
        self._remove_species_mapping = False
        self._marker_genes = False

        self.min_species = 0
        self.elapsed_time = 0

        self._reads = self.args.reads
        self._species_name = self.args.species_name

        self._ham_analysis = None

        self.progress = progress
        # self.progress.get_status(species_name=self._species_name)

        if self.args.remove_species_mapping:
            self.species_to_remove_mapping = self.args \
                .remove_species_mapping.split(",")
        else:
            self.species_to_remove_mapping = []

        if self.args.remove_species_ogs:
            self.species_to_remove_ogs = self.args \
                .remove_species_ogs.split(",")
        else:
            self.species_to_remove_ogs = []

        # if not load and self.progress.append_ogs_05:
        #     print("05_ogs_map_" + self._species_name)
        #     self.mapped_ogs = self._reload_ogs_from_folder(
        #         folder_suffix="05_ogs_map_" + self._species_name)
        #     self.ogs = self.mapped_ogs
        if not load and self.progress.ref_ogs_01:
            self.ogs = self._reload_ogs_from_folder()
        elif load and oma_output is not None:
            self.min_species = oma_output.min_species
            self.oma = oma_output
            self.ogs = oma_output.ogs
            self.oma_output_path = self.args.oma_output_path
            self.ogs = self._load_ogs()

    def _reload_ogs_from_folder(self, folder_suffix='01_ref_ogs'):
        """
        Re-load ogs if selection has finished and already exists in output
        folders
        :return: Dictionary with og name as key and list of SeqRecords
        """
        print('--- Re-load ogs and find their corresponding DNA seq '
              'from output folder ---')
        ogs = {}
        ref_ogs_aa = sorted(glob.glob(os.path.join(os.path.join(
            self.args.output_path, folder_suffix+"_aa"), "*.fa")))
        ref_ogs_dna = sorted(glob.glob(os.path.join(os.path.join(
            self.args.output_path, folder_suffix+"_dna"), "*.fa")))
        for file in tqdm(zip(ref_ogs_aa, ref_ogs_dna),
                         desc='Re-loading files', unit=' OGs'):
            name_og = os.path.basename(file[0]).split(".")[0]
            ogs[name_og] = OG()
            ogs[name_og].aa = list(SeqIO.parse(file[0], format='fasta'))
            ogs[name_og].dna = list(SeqIO.parse(file[1], format='fasta'))
            # ensure backward compatibility
            aa_ids = [r.id for r in ogs[name_og].aa if name_og in r.id]
            # if not aa_ids:
            #     for r in ogs[name_og].aa:
                #         tmp = r.id + q + name_og
            #         r.id = tmp
            # if self._remove_species:
            # ogs[name].remove_species_records(self.species_to_remove,
            #                                  species_in_hog)
        return ogs

    def _make_output_path(self, prefix):
        path = os.path.join(self.args.output_path, prefix)
        if not os.path.exists(path):
            os.makedirs(path)
        return path

    def _load_orthoxml(self):
        if self.oma.mode == 'standalone':
            og_orthoxml = os.path.join(self.oma_output_path,
                                       'OrthologousGroups.orthoxml')
            tree_str = os.path.join(self.oma_output_path,
                                    'EstimatedSpeciesTree.nwk')
            self._ham_analysis = pyham.Ham(tree_str, og_orthoxml,
                                     use_internal_name=False)

    def _load_dna_db(self):
        if '.fa' in self.args.dna_reference or \
           '.fasta' in self.args.dna_reference:
            db = {}
            self.logger.info('--- Load ogs and find their corresponding '
                  'DNA seq from {} ---'.format(self.args.dna_reference))
            self.logger.info('Loading {} into memory. This might take a '
                'while . . . '.format(self.args.dna_reference.split("/")[-1]))

            open_ = gzip.open if self.args.dna_reference.endswith('.gz') else open
            with open_(self.args.dna_reference, 'rt') as f:
                for rec in SeqIO.parse(f, 'fasta'):
                    db[rec.id.strip()] = str(rec.seq)
            source = 'fa'
            return db, source
        # ---------------- only to be used internally ----------------------
        # elif '.h5' in self.args.dna_reference:
        #     print('--- Load ogs and find their corresponding DNA \
        #           seq from {} ---'.format(self.args.dna_reference))
        #     self._db = db.Database(self.args.dna_reference)
        #     self._db_id_map = db.OmaIdMapper(self._db)
        #     self._db_source = 'h5'
        #     # self._db_species_list = [row['UniProtSpeciesCode'] \
        #                                .decode("utf-8") for row in
        #                                self._db_id_map.genome_table]
        #     # print(self._db_species_list)
        else:
            self.logger.info('--- Load ogs and find their corresponding DNA seq using '
                  'the REST api ---')
            source = 'REST_api'
            return None, source

    def _load_ogs(self):
        """
        Using the orthoxml file select only the OGs of interest
        that have more species than the min_species threshold
        :return: Dictionary with og name as key and list of SeqRecords
        """
        db, source = self._load_dna_db()
        self._load_orthoxml()
        start = time.time()
        ogs = {}

        orthologous_groups_aa = self._make_output_path("01_ref_ogs_aa")
        orthologous_groups_dna = self._make_output_path("01_ref_ogs_dna")

        names_og = self.ogs

        for name, records in tqdm(names_og.items(), desc='Loading OGs',
                                  unit=' OGs'):
            # name = file.split("/")[-1].split(".")[0]
            ogs[name] = OG()
            ogs[name].aa = self._get_aa_records(name, records)
            output_file_aa = os.path.join(orthologous_groups_aa,
                                          name + ".fa")
            output_file_dna = os.path.join(orthologous_groups_dna,
                                           name + ".fa")

            if source:
                try:
                    ogs[name].dna = self._get_dna_records(ogs[name].aa,
                                                          db, source, name)
                except (ValueError, TypeError):
                    self.logger.warning('This OG {} did not have any DNA'.format(name))
                    pass
                else:
                    all_len_consistent = self._check_dna_aa_length_consistency(name, ogs[name].aa, ogs[name].dna)
                    if "REST_api" in source and not all_len_consistent:
                        msg = "The returned DNA sequences from the REST API do not match the protein sequences. " \
                              "Most likely this is due to an update of the OMA Browser. Please download the DNA " \
                              "sequences from the download page of the OMA Browser that correspond to the release " \
                              "of your reference groups."
                        self.logger.error(msg)
                        raise Exception(msg)
                    self._write(output_file_dna, ogs[name].dna)
                    self._write(output_file_aa, ogs[name].aa)
            else:
                self.logger.debug('DNA reference was not provided. '
                                  'Only amino acid sequences gathered!')
        # self.progress.set_status('ogs')
        end = time.time()
        self.elapsed_time = end-start
        self.logger.info('{}: Gathering of DNA seq for {} OGs took {}.'
                         .format(self._species_name, len(names_og.keys()), self.elapsed_time))
        if db:
            db.clear()
        if self._ham_analysis:
            self._ham_analysis = None
        return ogs

    def _get_aa_records(self, name, records):
        """

        :param og_ham:
        :param records:
        :return:
        """
        if self.oma.mode == 'standalone':
            og_ham = self._ham_analysis.get_hog_by_id(name[2:])
            prot_ids = [gene.prot_id.split(" | ")[0]
                        for gene in og_ham.get_all_descendant_genes()]
            for record in records:
                mystr = record.description
                record.id = [x for x in prot_ids if mystr[mystr.find(
                    "[") + 1:mystr.find("]")] in x[0:5]][0]+"_"+name
                # Remove of stop codon
                if 'X' in record.seq[-1]:
                    tmp_seq = record.seq[0:-1]
                    record.seq = tmp_seq
        elif self.oma.mode == 'marker_genes':
            records = records
            for record in records:
                if 'X' in record.seq[-1]:
                    tmp_seq = record.seq[0:-1]
                    record.seq = tmp_seq
        return records

    def _get_dna_from_h5(self, record):
        """

        :param record:
        :param name:
        :return:
        """
        try:
            oma_db_nr = self._db_id_map.omaid_to_entry_nr(record.id)
        except ValueError:
            self.logger.debug('DNA not found for {}.'.format(record.id))
            pass
        else:
            seq = self._db.get_cdna(oma_db_nr).decode("utf-8")

            # if 'X' in seq:
            cleaned_seq = self._clean_DNA_seq(seq)
            # else:
            #     cleaned_seq = seq

            return SeqRecord.SeqRecord(Seq.Seq(cleaned_seq),
                                       id=record.id,
                                       description="")

    def _get_dna_from_REST(self, record):
        """

        :param record:
        :param name:
        :return:
        """
        tmp_id = re.sub(r'\..*', '', record.id.split("_")[0])
        use_id = re.sub(r'\W+', '', tmp_id)
        try:
            oma_record = requests.get(API_URL + "/protein/" + use_id + "/")
        except requests.exceptions.RequestException:
            self.logger.debug('DNA not found for {}.'.format(use_id))
            pass
        else:
            seq = oma_record.json()['cdna']
            rec_id = oma_record.json()['omaid']

            cleaned_seq = self._clean_DNA_seq(seq)
            # else:
            #     cleaned_seq = dna_record.seq
            return SeqRecord.SeqRecord(cleaned_seq, record.id,
                                       description="", name="")

    def _get_dna_from_REST_bulk(self, records, og_name):
        """

        :param record:
        :param name:
        :return:
        """
        record_ids = [r.id for r in records]
        dna_records = []
        try:
            reply = requests.post('https://omabrowser.org/api/protein/bulk_retrieve/',
                                  json={"ids": record_ids},
                                  headers={'User-Agent': 'read2tree/'+read2tree_version})
        except requests.exceptions.RequestException as exception_type:
            self.logger.warning('DNA not found probably for '+str(record_ids[0])+'. The reason is '+str(exception_type))
            pass
        else:
            group_members = reply.json()
            for memb in group_members:
                # print(">{}\n{}\n\n".format(memb['omaid'], memb['cdna']))
                seq = memb['target']['cdna']
                rec_id = memb['target']['omaid']+"_"+og_name
                cleaned_seq = self._clean_DNA_seq(seq)
                # print(cleaned_seq)
                dna_records.append(SeqRecord.SeqRecord(cleaned_seq, id=rec_id,
                                       description="", name=""))
        return dna_records

    def _get_dna_from_fasta(self, record, db):
        try:
            if record.id.split("_")[0] not in db.keys():
                return self._get_dna_from_REST(record)
            else:
                dna = db[record.id.split("_")[0]]
        except ValueError:
            self.logger.debug('DNA not found for {}.'.format(record.id))
            pass
        else:
            return SeqRecord.SeqRecord(self._clean_DNA_seq(dna),
                                           id=record.id,
                                           description="")
            # else:
            #     return SeqRecord.SeqRecord(Seq.Seq(dna.upper()),  id=record.id,
            #                                description="")

    def _check_dna_aa_length_consistency(self, og_name, aa, dna):
        dna_dic = {r.id.split("_")[0]: r for r in dna}
        aa_dic = {r.id.split("_")[0]: r for r in aa}
        all_consistent = True
        for k, r_dna in dna_dic.items():
            r_aa = aa_dic[k]
            if abs(len(r_dna.seq) - 3*len(r_aa.seq)) > 3:
                self.logger.warning('{}: {} has aa-length {} and dna-length {}'.format(self._species_name, og_name+" "+k, 3*len(r_aa.seq), len(r_dna.seq)))
                all_consistent = False
        return all_consistent

    def _get_dna_records(self, records, db, source, og_name):
        """

        :param records:
        :return:
        """
        og_cdna = []
        if 'REST_api' in source:
            return self._get_dna_from_REST_bulk(records, og_name)
        else:
            for i, record in enumerate(records):
                if 'h5' in source:
                    og_cdna.append(self._get_dna_from_h5(record))
                elif 'fa' in source:
                    og_cdna.append(self._get_dna_from_fasta(record, db))
                # elif 'REST_api' in source:
                #     og_cdna.append(self._get_dna_from_REST(record))

            return og_cdna

    def _clean_DNA_seq(self, record):
        """
        Exchange all X in sequence with N
        :param record: Biopython SeqRecord object
        :return: Biopython seq object with Ns instead of Xs
        """
        # replace all non GATC chars with N
        if isinstance(record, str):
            seq = record.upper()

        else:
            seq = str(record.seq).upper()
        outseq = re.sub('[^GATC]', 'N', seq)
        # remove stopcodon nucleotides
        if 'TGA' in outseq[-3:] or 'TAA' in outseq[-3:] or 'TAG' in outseq[-3:] or 'NNN' in outseq[-3:]:
            outseq = outseq[:-3]

        return Seq.Seq(outseq)

    # def _translate_dna(self, record):
    #     """
    #     Given a list of sequences that are derived from mapped reads to
    #     multiple seq of a OG we find the best corresponding mapped seq by
    #     comparing it with a representative sequence of the original OG using
    #     pyopa local alignment and return the sequence with its highest score!
    #     :return:
    #     """
    #     try:
    #         frame = record.seq[0:].translate(table='Standard', stop_symbol='X', to_stop=False, cds=False)
    #     except ValueError:
    #         raise ValueError("Problem with sequence format!")
    #     return frame

    def remove_species_from_ogs(self):
        for name_og, og in tqdm(self.ogs.items(),
                                desc='Adding mapped seq to OG', unit=' OGs'):
            og_filt = self._remove_species_from_original_set(og)
            self.ogs[name_og] = og_filt

    #TODO: this has to be moved to OG not to OGSet
    def _remove_species_from_original_set(self, current_og):
        """
        Removes sequence records for a species / set of
        species from the reference OGSet
        :param current_og: The current OG object
        :return: OG object with removed species
        """
        if self.args.remove_species_ogs:
            og = OG()
            filtered_og = current_og \
                .remove_species_records(self.species_to_remove_ogs)
            if filtered_og:
                og.dna = filtered_og[0]
                og.aa = filtered_og[1]
        else:
            og = current_og
        return og

    def _remove_species_from_mapping_only(self, mapped_og):
        """
        Removes sequence records for a species / set of species
        from the mapped OGSet
        :param mapped_og: The mapped OG object
        :return: mapped OG object with removed species
        """
        if self.species_to_remove_mapping:  # in case we decided to remove species from the mapping
            consensus_og = OG()
            filtered_mapping = mapped_og \
                .remove_species_records(self.species_to_remove_mapping)
            if filtered_mapping:
                consensus_og.dna = filtered_mapping[0]
                consensus_og.aa = filtered_mapping[1]
        else:  # nothing to remove
            consensus_og = mapped_og
        return consensus_og

    def _get_best_record(self, cons_og, og_sc, species_name):
        best_records = cons_og \
            .get_best_consensus_by_seq_completeness(
                og_sc, threshold=self.args.sc_threshold)
        if best_records:
            best_record_aa = best_records[0]
            best_record_dna = best_records[1]
            best_record_aa.id = species_name
            best_record_dna.id = species_name
            return (best_record_aa, best_record_dna)
        else:
            return none

    def _generate_seq_completeness(self, seqC, mapper, og, best_record_dna):
        if self.args.remove_species_ogs:
            tested_rec = [
                rec for rec in og.dna if
                self.args.remove_species_ogs
                in rec.id]
            if tested_rec:
                seqC2 = SeqCompleteness(
                    mapped_ref=[best_record_dna],
                    tested_ref=tested_rec)
            else:
                seqC2 = SeqCompleteness(
                    mapped_ref=[best_record_dna])
            seqC2.get_seq_completeness([best_record_dna])
            seqC.add_seq_completeness(
                self._get_clean_id(best_record_dna),
                seqC2.seq_completeness[best_record_dna.id])
        else:
            seqC.add_seq_completeness(
                self._get_clean_id(best_record_dna),
                mapper.all_sc[self._get_clean_id(best_record_dna)])

    def _get_id_rec(self, record):
        parts = record.id.split('_')
        if len(parts) > 2:
            return record.id.split('_')[0]+'_'+record.id.split('_')[1]
        else:
            return record.id

    def add_mapped_seq(self, mapper, species_name=None):
        """
        Add the sequence given from the read mapping to its corresponding
        OG and retain
        all OGs that do not have the mapped sequence, thus all original OGs
        are used for tree inference
        :param cons_og_set: set of ogs with its mapped sequences
        """
        start = time.time()
        cons_og_set = mapper.og_records  # get sequences from mapping
        cov = Coverage(self.args)
        seqC = SeqCompleteness()
        if not species_name:
            species_name = self._species_name

        print('--- Add inferred mapped sequence back to OGs ---')

        # iterate through all existing ogs
        for name_og, og in tqdm(self.ogs.items(),
                                desc='Adding mapped seq to OG', unit=' OGs'):
            # og_filt = self._remove_species_from_original_set(og)
            og_filt = og
            if len(og_filt.aa) >= 2:
                # continue only if OG is in mapped OGs
                if name_og in cons_og_set.keys():
                    cons_og_filt = self \
                        ._remove_species_from_mapping_only(cons_og_set[name_og])
                    og_sc = {self._get_id_rec(rec): mapper.all_sc[self._get_id_rec(rec)]
                             for rec in cons_og_filt.aa
                             if self._get_id_rec(rec) in mapper.all_sc.keys()}
                    og_cov = {self._get_id_rec(rec): mapper.all_cov[self._get_id_rec(rec)]
                             for rec in cons_og_filt.aa
                             if self._get_id_rec(rec) in mapper.all_cov.keys()}
                    if len(cons_og_filt.aa) >= 1:  # we had at least one mapped og even after removal
                        best_records = cons_og_filt \
                            .get_best_consensus_by_seq_completeness(self.args.sequence_selection_mode,
                                sc=og_sc, cov=og_cov, threshold=self.args.sc_threshold)
                        if best_records:
                            best_record_aa = best_records[0]
                            best_record_dna = best_records[1]
                            self._generate_seq_completeness(seqC, mapper,
                                                            og, best_record_dna)
                            cov.add_coverage(self._get_clean_id(best_record_aa),
                                             mapper.all_cov[self._get_clean_id(best_record_aa)])
                            best_record_aa.id = species_name
                            best_record_dna.id = species_name
                            self.mapped_ogs[name_og] = og_filt
                            all_id = [rec.id
                                      for rec in self.mapped_ogs[name_og].aa]
                            if best_record_aa.id not in all_id:  # make sure that repeated run doesn't add the same sequence multiple times at the end of an OG
                                self.mapped_ogs[name_og] \
                                    .aa.append(best_record_aa)
                                self.mapped_ogs[name_og] \
                                    .dna.append(best_record_dna)
                        else:  # case where no best_record_aa reported because it was smaller than the self.args.sc_threshold
                            self.mapped_ogs[name_og] = og_filt
                    else:  # mapping had only one that we removed
                        if self.args.keep_all_ogs:
                            self.mapped_ogs[name_og] = og_filt
                else:  # nothing was mapped to that og
                    if self.args.keep_all_ogs:
                        self.mapped_ogs[name_og] = og_filt
            else:
                self.logger.debug('{} was left only with a single entry '
                             'and hence not used for further '
                             'processing'.format(name_og))

        cov.write_coverage_bam(os.path.join(self.args.output_path,
                                            species_name+'_all_cov.txt'))
        seqC.write_seq_completeness(os.path.join(self.args.output_path,
                                                 species_name+'_all_sc.txt'))
        end = time.time()
        self.elapsed_time = end-start
        self.logger.info('{}: Appending {} reconstructed sequences to present OG '
                    'took {}.'
                    .format(self._species_name,
                            len(list(cons_og_set.keys())),
                            self.elapsed_time))

    def write_added_ogs_aa(self, folder_name=None):
        """

        :param self:
        :param folder_name:
        :return:
        """
        if folder_name is None:
            ogs_with_mapped_seq = self._make_output_path("05_ogs_map_" +
                                                         self._species_name +
                                                         "_aa")
        else:
            ogs_with_mapped_seq = self._make_output_path(folder_name)

        for name, value in self.ogs.items():
            if name in self.mapped_ogs.keys():
                output_file = os.path.join(ogs_with_mapped_seq, name + ".fa")
                self._write(output_file, self.mapped_ogs[name].aa)

    def write_added_ogs_dna(self, folder_name=None):
        """

        :param self:
        :param folder_name:
        :return:
        """
        if folder_name is None:
            ogs_with_mapped_seq = self._make_output_path(
                "05_ogs_map_" + self._species_name + "_dna")
        else:
            ogs_with_mapped_seq = self._make_output_path(folder_name)

        for name, value in self.ogs.items():
            if name in self.mapped_ogs.keys():
                output_file = os.path.join(ogs_with_mapped_seq, name + ".fa")
                self._write(output_file, self.mapped_ogs[name].dna)

    def _get_clean_id(self, record):
        """

        :param record:
        :return:
        """
        des = record.description.split(" ")[0]
        des = des.split("_")
        return des[0]+"_"+des[1]

    def _write(self, file, value):
        """
        Write output to fasta file
        :param file: file and location of outputfile
        :param value:
        :return:
        """
        handle = open(file, "w")
        writer = FastaWriter(handle, wrap=None)
        writer.write_file(value)
        handle.close()

    def write_select_og_aa(self):
        """
        Write for each species all the DNA sequences into separate fasta files
        :param output_folder: folder where files should be stored
        """
        output_folder = os.path.join(self.args.output_path, "reference_ogs_aa")
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
            for key, value in tqdm(self.ogs.items(), desc="Writing OGs sorted by species",
                                   unit=" species"):
                handle = open(os.path.join(output_folder, key + '.fa'), "w")
                writer = FastaWriter(handle, wrap=None)
                writer.write_file(value.aa)
                handle.close()
        elif len(self.ogs) == len(glob.glob(os.path.join(output_folder, '*.fa'))):
            print('Folder with files already exists and will not be overwritten.')

    def write_select_og_dna(self):
        """
        Write for each species all the DNA sequences into separate fasta files
        :param output_folder: folder where files should be stored
        """
        output_folder = os.path.join(self.args.output_path, "reference_ogs_dna")
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
            for key, value in tqdm(self.ogs.items(), desc="Writing OGs sorted by species",
                                   unit=" species"):
                handle = open(os.path.join(output_folder, key + '.fa'), "w")
                writer = FastaWriter(handle, wrap=None)
                writer.write_file(value.dna)
                handle.close()
        elif len(self.ogs_dna_by_species) == len(glob.glob(os.path.join(output_folder, '*.fa'))):
            print('Folder with files already exists and will not be overwritten.')

    def append(self, name, record):
        self.ogs[name].append(record)


class OG(object):

    def __init__(self):
        self.aa = []
        self.dna = []

    def _get_og_dict(self, ref_og):
        dna_dict = {}
        for record in ref_og.dna:
            if '_' in record.id:
                tmp = record.id.split("_")[0]
                record.id = tmp

            dna_dict[record.id] = record
        return dna_dict

    def _get_record_by_id(self, records, idx):
        for rec in records:
            if idx in rec.id:
                return rec

    def get_best_consensus_by_seq_completeness(self, sequence_selection_mode, sc=None, cov=None,
                                               threshold=0.0):
        """
        :param ref_og: OG containing reference sequences
        :param threshold: minimum sequence completeness [0.0]
        :return: best amino acid sequence
        """
        if sequence_selection_mode == 'cov_sc':
            sc_cov = {k: v[0]*cov[k][0] for k,v in sc.items()}
            sc_cov_ordered = OrderedDict(sorted(sc_cov.items(), key=lambda t: t[-1]))
            best_record_id = list(sc_cov_ordered.items())[-1][0]
            seq_completenesses = sc[best_record_id][1]
            if seq_completenesses >= threshold:
                return (self._get_record_by_id(self.aa, best_record_id),
                        self._get_record_by_id(self.dna, best_record_id))
            else:
                return None
        elif sequence_selection_mode == 'cov':
            cov_ordered = OrderedDict(sorted(cov.items(), key=lambda t: t[-1]))
            best_record_id = list(cov_ordered.items())[-1][0]
            seq_completenesses = sc[best_record_id][1]
            if seq_completenesses >= threshold:  # check whether best sequence by coverage is above sc threshold
                return (self._get_record_by_id(self.aa, best_record_id),
                        self._get_record_by_id(self.dna, best_record_id))
            else:
                return None
        elif sequence_selection_mode == 'cov_pure':
            cov_ordered = OrderedDict(sorted(cov.items(), key=lambda t: t[-1]))
            best_record_id = list(cov_ordered.items())[-1][0]
            if best_record_id:  # check whether best sequence by coverage is above sc threshold
                return (self._get_record_by_id(self.aa, best_record_id),
                        self._get_record_by_id(self.dna, best_record_id))
            else:
                return None
        elif sequence_selection_mode == 'cov_no_sc':
            cov_ordered = OrderedDict(sorted(cov.items(), key=lambda t: t[-1]))
            best_record_id = list(cov_ordered.items())[-1][0]
            seq_completenesses = sc[best_record_id][1]
            if seq_completenesses >= 0.0:  # check whether best sequence by coverage is above sc threshold
                return (self._get_record_by_id(self.aa, best_record_id),
                        self._get_record_by_id(self.dna, best_record_id))
            else:
                return None
        elif sequence_selection_mode == 'sc':
            sc_ordered = OrderedDict(sorted(sc.items(), key=lambda t: t[-1]))
            best_record_id = list(sc_ordered.items())[-1][0]
            seq_completenesses = sc[best_record_id][1]
            if seq_completenesses >= threshold:
                return (self._get_record_by_id(self.aa, best_record_id),
                        self._get_record_by_id(self.dna, best_record_id))
            else:
                return None
        elif sequence_selection_mode == 'cov_sc_scaled':
            max_cov = np.max([v[0] for k,v in cov.items()])
            sc_cov = {k: v[0] * cov[k][0]/max_cov for k, v in sc.items()}
            sc_cov_ordered = OrderedDict(sorted(sc_cov.items(), key=lambda t: t[-1]))
            best_record_id = list(sc_cov_ordered.items())[-1][0]
            seq_completenesses = sc[best_record_id][1]
            if seq_completenesses >= threshold:
                return (self._get_record_by_id(self.aa, best_record_id),
                        self._get_record_by_id(self.dna, best_record_id))
            else:
                return None
        elif sequence_selection_mode == 'random':
            best_record_id = random.choice(list(sc.keys()))
            # sc_cov = {k: v[0] * cov[k][0]/max_cov for k, v in sc.items()}
            # sc_cov_ordered = OrderedDict(sorted(sc_cov.items(), key=lambda t: t[-1]))
            # best_record_id = list(sc_cov_ordered.items())[-1][0]
            seq_completenesses = sc[best_record_id][1]
            if seq_completenesses >= threshold:
                return (self._get_record_by_id(self.aa, best_record_id),
                        self._get_record_by_id(self.dna, best_record_id))
            else:
                return None
        else:
            return None

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

    def remove_species_records(self, species_to_remove):
        """
        Remove species from reference sequence set
        :param species_to_remove: list of species to be removed
        :param all_species: list of all species present in analysis
        """
        aa = [record for i, record in enumerate(
            self.aa) if self._get_species_id(record) not in species_to_remove]
        dna = [record for i, record in enumerate(
            self.dna) if self._get_species_id(record) not in species_to_remove]
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
