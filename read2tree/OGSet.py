#!/usr/bin/env python
'''
    This file contains definitions of a class which holds the orthologous groups.

    -- David Dylus, July--XXX 2017
'''
import glob
import os
import re
import pyham
import requests
import logging

from ete3 import Tree
from tqdm import tqdm
from Bio import SeqIO, Seq, SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
from Bio.SeqIO.FastaIO import FastaWriter
#from tables import *
# ----------- only to be used internally; requires hdf5 installation -------------------
#from pyoma.browser import db

from read2tree.Progress import Progress
from read2tree.stats.Coverage import Coverage
from read2tree.stats.SeqCompleteness import SeqCompleteness



OMA_STANDALONE_OUTPUT = 'Output'
OMA_MARKER_GENE_EXPORT = 'marker_genes'
API_URL = 'http://omabrowser.org/api'



logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('input.log')
file_handler.setLevel(logging.ERROR)
file_handler.setFormatter(formatter)

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

logger.addHandler(file_handler)
logger.addHandler(stream_handler)



class OGSet(object):

    def __init__(self, args, oma_output=None, load=True):
        self.args = args


        self.progress = Progress(args)
        self.mapped_ogs = {}
        self._db = None
        self._db_id_map = None
        self._db_source = None
        self._db_species_list = None
        self._ham_analysis = None
        self._tree_str = None
        self._og_orthoxml = None
        self._remove_species_mapping = False
        self._marker_genes = False

        self.min_species = self._estimate_best_number_species()

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
            self._species_name = 'merged'

        if self.args.remove_species_mapping:
            self.species_to_remove_mapping = self.args.remove_species_mapping.split(",")
        else:
            self.species_to_remove_mapping = []

        if self.args.remove_species_ogs:
            self.species_to_remove_ogs = self.args.remove_species_ogs.split(",")
        else:
            self.species_to_remove_ogs = []

        if not load and self.progress.append_ogs_04:
            print('here')
            print("04_ogs_map_" + self._species_name)
            self.mapped_ogs = self._reload_ogs_from_folder(folder_suffix="04_ogs_map_" + self._species_name)
        elif not load and self.progress.ref_ogs_01:
            self.ogs = self._reload_ogs_from_folder()
        elif load and oma_output is not None:
            self.oma = oma_output
            self.ogs = oma_output.ogs
            self.oma_output_path = self.args.oma_output_path
            self.ogs = self._load_ogs()


    def _reload_ogs_from_folder(self, folder_suffix='01_ref_ogs'):
        """
        Re-load ogs if selection has finished and already exists in output folders
        :return: Dictionary with og name as key and list of SeqRecords
        """
        print('--- Re-load ogs and find their corresponding DNA seq from output folder ---')
        ogs = {}
        ref_ogs_aa = sorted(glob.glob(os.path.join(os.path.join(self.args.output_path, folder_suffix+"_aa"), "*.fa")))
        ref_ogs_dna = sorted(glob.glob(os.path.join(os.path.join(self.args.output_path, folder_suffix+"_dna"), "*.fa")))
        for file in tqdm(zip(ref_ogs_aa, ref_ogs_dna), desc='Re-loading files', unit=' OGs'):
            name = os.path.basename(file[0]).split(".")[0]
            ogs[name] = OG()
            ogs[name].aa = list(SeqIO.parse(file[0], format='fasta'))
            ogs[name].dna = list(SeqIO.parse(file[1], format='fasta'))
            # if self._remove_species:
            #     ogs[name].remove_species_records(self.species_to_remove, species_in_hog)
        return ogs

    def _make_output_path(self, prefix):
        path = os.path.join(self.args.output_path, prefix)
        if not os.path.exists(path):
            os.makedirs(path)
        return path


    def _load_ogs(self):
        """
        Using the orthoxml file select only the OGs of interest that have more species than the min_species threshold
        :return: Dictionary with og name as key and list of SeqRecords
        """

        if '.fa' in self.args.dna_reference or '.fasta' in self.args.dna_reference:
            print('--- Load ogs and find their corresponding DNA seq from {} ---'.format(self.args.dna_reference))
            print(
                'Loading {} into memory. This might take a while . . . '.format(self.args.dna_reference.split("/")[-1]))
            self._db = SeqIO.index(self.args.dna_reference, "fasta")
            self._db_source = 'fa'
        # ---------------- only to be used internally ----------------------
        # elif '.h5' in self.args.dna_reference:
        #     print('--- Load ogs and find their corresponding DNA seq from {} ---'.format(self.args.dna_reference))
        #     self._db = db.Database(self.args.dna_reference)
        #     self._db_id_map = db.OmaIdMapper(self._db)
        #     self._db_source = 'h5'
        #     # self._db_species_list = [row['UniProtSpeciesCode'].decode("utf-8") for row in self._db_id_map.genome_table]
        #     # print(self._db_species_list)
        else:
            print('--- Load ogs and find their corresponding DNA seq using the REST api ---')
            self._db_source = 'REST_api'

        if self.oma.mode is 'standalone':
            self._og_orthoxml = os.path.join(self.oma_output_path, 'OrthologousGroups.orthoxml')
            self._tree_str = os.path.join(self.oma_output_path, 'EstimatedSpeciesTree.nwk')
            self._ham_analysis = pyham.Ham(self._tree_str, self._og_orthoxml, use_internal_name=False)

        ogs = {}

        orthologous_groups_aa = self._make_output_path("01_ref_ogs_aa")
        orthologous_groups_dna = self._make_output_path("01_ref_ogs_dna")

        names_og = self.ogs

        for name, records in tqdm(names_og.items(), desc='Loading OGs', unit=' OGs'):
            # name = file.split("/")[-1].split(".")[0]
            ogs[name] = OG()
            ogs[name].aa = self._get_aa_records(name, records)
            output_file_aa = os.path.join(orthologous_groups_aa, name + ".fa")
            output_file_dna = os.path.join(orthologous_groups_dna, name + ".fa")

            if self._db_source:
                try:
                    ogs[name].dna = self._get_dna_records(ogs[name].aa, name)
                except (ValueError, TypeError):
                    logger.debug('This OG {} did not have any DNA'.format(name))
                    pass
                else:
                    self._write(output_file_dna, ogs[name].dna)
                    self._write(output_file_aa, ogs[name].aa)
            else:
                print("DNA reference was not provided. Only amino acid sequences gathered!")
        self.progress.set_status('ogs')
        return ogs


    def _get_aa_records(self, name, records):
        """
        
        :param og_ham: 
        :param records: 
        :return: 
        """
        if self.oma.mode is 'standalone':
            og_ham = self._ham_analysis.get_hog_by_id(name[2:])
            prot_ids = [gene.prot_id.split(" | ")[0] for gene in og_ham.get_all_descendant_genes()]
            for record in records:
                mystr = record.description
                record.id = [x for x in prot_ids if mystr[mystr.find("[") + 1:mystr.find("]")] in x[0:5]][0]
        elif self.oma.mode is 'marker_genes':
            records = records
        return records


    def _get_from_h5(self, record, name):
        '''

        :param record:
        :param name:
        :return:
        '''
        try:
            oma_db_nr = self._db_id_map.omaid_to_entry_nr(record.id)
        except ValueError:
            logger.debug('DNA not found for {}.'.format(record.id))
            pass
        else:
            seq = self._db.get_cdna(oma_db_nr).decode("utf-8")

            if 'X' in seq:
                cleaned_seq = self._clean_DNA_seq(seq)
            else:
                cleaned_seq = seq

            return SeqRecord.SeqRecord(Seq.Seq(cleaned_seq), id=record.id + "_" + name, description="")


    def _get_from_REST(self, record, name):
        '''

        :param record:
        :param name:
        :return:
        '''
        try:
            oma_record = requests.get(API_URL + "/protein/" + record.id + "/")
        except requests.exceptions.RequestException:
            logger.debug('DNA not found for {}.'.format(record.id))
            pass
        else:
            dna_record = SeqRecord.SeqRecord(Seq.Seq(oma_record.json()['cdna']), id=record.id + "_" + name, description="")
            if 'X' in str(dna_record.seq):
                cleaned_seq = self._clean_DNA_seq(dna_record)
            else:
                cleaned_seq = dna_record.seq
            return SeqRecord.SeqRecord(cleaned_seq, id=record.id + "_" + name, description="")


    def _get_from_fasta(self, record):
        try:
            dna = self._db[record.id]
        except ValueError:
            logger.debug('DNA not found for {}.'.format(record.id))
            pass
        else:
            if 'X' in str(dna.seq):
                return SeqRecord.SeqRecord(self._clean_DNA_seq(dna), id=record.id + "_" + name, description="")
            else:
                return SeqRecord.SeqRecord(dna.seq,  id=record.id + "_" + name, description="")


    def _get_dna_records(self, records, name):
        """
        
        :param records: 
        :return: 
        """
        og_cdna = []
        for i, record in enumerate(records):
            if 'h5' in self._db_source:
                og_cdna.append(self._get_from_h5(record, name))
            elif 'fa' in self._db_source:
                og_cdna.append(self._get_from_fasta(record))
            elif 'REST_api' in self._db_source:
                og_cdna.append(self._get_from_REST(record, name))

        return og_cdna

    def _estimate_best_number_species(self):
        """
        Estimate min number of species such that around 1000 OGs are selected
        :return:
        """
        min_species = 0
        if self.args.min_species is not None:
            min_species = self.args.min_species

        return min_species

    def _clean_DNA_seq(self, record):
        '''
        Exchange all X in sequence with N
        :param record: Biopython SeqRecord object
        :return: Biopython seq object with Ns instead of Xs
        '''
        return Seq.Seq(re.sub('[^GATC]', 'N', str(record.seq).upper()), SingleLetterAlphabet())

    def _remove_species_from_original_set(self, current_og):
        '''
        Removes sequence records for a species / set of species from the reference OGSet
        :param current_og: The current OG object
        :return: OG object with removed species
        '''
        if self.args.remove_species_ogs:
            og = OG()
            filtered_og = current_og.remove_species_records(self.species_to_remove_ogs)
            if filtered_og:
                og.dna = filtered_og[0]
                og.aa = filtered_og[1]
        else:
            og = current_og
        return og

    def _remove_species_from_mapping_only(self, mapped_og):
        '''
        Removes sequence records for a species / set of species from the mapped OGSet
        :param mapped_og: The mapped OG object
        :return: mapped OG object with removed species
        '''
        if self.species_to_remove_mapping:  # in case we decided to remove species from the mapping
            mapping_og = OG()
            filtered_mapping = mapped_og.remove_species_records(self.species_to_remove_mapping)
            if filtered_mapping:
                mapping_og.dna = filtered_mapping[0]
                mapping_og.aa = filtered_mapping[1]
        else:  # nothing to remove
            mapping_og = mapped_og
        return mapping_og


    def add_mapped_seq(self, mapper, species_name=None):
        """
        Add the sequence given from the read mapping to its corresponding OG and retain
        all OGs that do not have the mapped sequence, thus all original OGs are used for
        tree inference
        :param mapped_og_set: set of ogs with its mapped sequences
        """
        mapped_og_set = mapper.og_records # get sequences from mapping
        cov = Coverage()
        seqC = SeqCompleteness()
        if not species_name:
            species_name = self._species_name

        print('--- Add inferred mapped sequence back to OGs ---')

        # iterate through all existing ogs
        for name, value in tqdm(self.ogs.items(), desc='Adding mapped seq to OG', unit=' OGs'):
            og = self._remove_species_from_original_set(value)
            if len(og.aa) > 2:
                # continue only if OG is in mapped OGs
                if name in mapped_og_set.keys():
                    mapping_og = self._remove_species_from_mapping_only(mapped_og_set[name])

                    if len(mapping_og.aa) >= 1:  # we had at least one mapped og even after removal
                        best_records = mapping_og.get_best_mapping_by_seq_completeness(ref_og=og, threshold=self.args.sc_threshold)
                        if best_records:
                            best_record_aa = best_records[0]
                            best_record_dna = best_records[1]
                            best_record_aa.id = species_name
                            best_record_dna.id = species_name
                            self.mapped_ogs[name] = og
                            all_id = [rec.id for rec in self.mapped_ogs[name].aa]
                            if best_record_aa.id not in all_id:  # make sure that repeated run doesn't add the same sequence multiple times at the end of an OG
                                #print(mapper.all_cov)
                                cov.add_coverage(self._get_clean_id(best_record_aa), mapper.all_cov[self._get_clean_id(best_record_aa)])
                                seqC.add_seq_completeness(self._get_clean_id(best_record_aa), mapper.all_sc[self._get_clean_id(best_record_aa)])
                                self.mapped_ogs[name].aa.append(best_record_aa)
                                self.mapped_ogs[name].dna.append(best_record_dna)
                        else:  # case where no best_record_aa reported because it was smaller than the self.args.sc_threshold
                            self.mapped_ogs[name] = og
                    else:  # mapping had only one that we removed
                        if self.args.keep_all_ogs:
                            self.mapped_ogs[name] = og
                else:  # nothing was mapped to that og
                    if self.args.keep_all_ogs:
                        self.mapped_ogs[name] = og
            else:
                logger.debug('{} was left only with a single entry and hence not used for further processing'.format(name))

        cov.write_coverage_bam(os.path.join(self.args.output_path, species_name+'_all_cov.txt'))
        seqC.write_seq_completeness(os.path.join(self.args.output_path, species_name+'_all_sc.txt'))

    def write_added_ogs_aa(self, folder_name=None):
        '''

        :param self:
        :param folder_name:
        :return:
        '''
        if folder_name is None:
            ogs_with_mapped_seq = self._make_output_path("04_ogs_map_" + self._species_name + "_aa")
        else:
            ogs_with_mapped_seq = self._make_output_path(folder_name)

        for name, value in self.ogs.items():
            if name in self.mapped_ogs.keys():
                output_file = os.path.join(ogs_with_mapped_seq, name + ".fa")
                self._write(output_file, self.mapped_ogs[name].aa)


    def write_added_ogs_dna(self, folder_name=None):
        '''

        :param self:
        :param folder_name:
        :return:
        '''
        if folder_name is None:
            ogs_with_mapped_seq = self._make_output_path("04_ogs_map_" + self._species_name + "_dna")
        else:
            ogs_with_mapped_seq = self._make_output_path(folder_name)

        for name, value in self.ogs.items():
            if name in self.mapped_ogs.keys():
                output_file = os.path.join(ogs_with_mapped_seq, name + ".fa")
                self._write(output_file, self.mapped_ogs[name].dna)


    def _get_clean_id(self, record):
        '''

        :param record:
        :return:
        '''
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
        '''
        Write for each species all the DNA sequences into separate fasta files
        :param output_folder: folder where files should be stored
        '''
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
        '''
        Write for each species all the DNA sequences into separate fasta files
        :param output_folder: folder where files should be stored
        '''
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

    def get_best_mapping_by_seq_completeness(self, ref_og=None, threshold=0.0):
        """

        :param ref_og: OG containing reference sequences
        :param threshold: minimum sequence completeness [0.0]
        :return: best amino acid sequence
        """
        seq_completenesses = self._get_seq_completeness(ref_og=ref_og)
        best_record = seq_completenesses.index(max(seq_completenesses))
        if seq_completenesses[best_record] >= threshold:
            return (self.aa[best_record], self.dna[best_record])
        else:
            return None

    def _get_og_dict(self, ref_og):
        dna_dict = {}
        for record in ref_og.dna:
            if '_' in record.id:
                tmp = record.id.split("_")[0]
                record.id = tmp

            dna_dict[record.id] = record
        return dna_dict

    def _get_seq_completeness(self, ref_og=None):
        """
        TODO: this has to be changed to incorporate the expected sequence length
        :param gene_code:
        :return:
        """
        ref_og_dna = self._get_og_dict(ref_og)
        full_seq_completeness = []
        for record in self.dna:
            map_seq = str(record.seq).upper()
            ref_seq = str(ref_og_dna[record.name.split("_")[0]].seq).upper()
            full_seq_len = len(ref_seq)
            non_n_len = len(map_seq) - map_seq.count('N')
            full_seq_completeness.append(non_n_len / full_seq_len)
        return full_seq_completeness

    def _get_species_id(self, description):
        '''
        Sequences in OMA are marked by using the first three letters of genus and the first 2 letters of species,
        e.g. Amphiura filiformis = AMPFI. This however is not always the case and we therefore prioritize the id
        (e.g. MOUSE over MUSMU that comes from Mus musculus).
        :param description: SeqRecord description
        :return: species_id
        '''
        species = description[description.find("[") + 1:description.find("]")]
        if len(species.split(" ")) > 1:
            new_id = species.split(" ")[0][0:3] + species.split(" ")[1][0:2]
            species = new_id.upper()
        species_id = description[0:5]
        if species_id in species:
            return species
        else:
            return species_id

    def remove_species_records(self, species_to_remove):
        '''
        Remove species from reference sequence set
        :param species_to_remove: list of species to be removed
        :param all_species: list of all species present in analysis
        '''
        aa = [record for i, record in enumerate(self.aa) if self._get_species_id(record.description) not in species_to_remove]
        dna = [record for i, record in enumerate(self.dna) if self._get_species_id(record.description) not in species_to_remove]
        if len(aa) > 0 and len(dna) > 0:
            return [dna, aa]
        else:
            return None
