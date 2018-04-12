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
API_URL = 'http://omadev.cs.ucl.ac.uk/api'


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
            self._species_name = 'merge'

        if self.args.remove_species_mapping:
            self.species_to_remove_mapping = self.args.remove_species_mapping.split(",")
        else:
            self.species_to_remove_mapping = []

        if self.args.remove_species_ogs:
            self.species_to_remove_ogs = self.args.remove_species_ogs.split(",")
        else:
            self.species_to_remove_ogs = []

        if load and oma_output is not None:
            self.oma = oma_output
            self.ogs = oma_output.ogs
            self.oma_output_path = self.args.oma_output_path
            self.ogs = self._load_ogs()
        else:
            self.ogs = self._reload_ogs_from_folder()

    def og(self, name):
        return self.ogs[name]

    def _check_oma_standalone_path(self):
        """
        :return: return true if oma standalone output path provided
        """
        output_path = os.path.join(self.args.standalone_path, OMA_STANDALONE_OUTPUT)
        ortho_group_xml = os.path.getsize(os.path.join(output_path, "OrthologousGroups.orthoxml"))
        if os.path.getsize(os.path.join(output_path, "OrthologousGroups.orthoxml")) > 0:
            return True
        else:
            return False

    # def _remove_species(self):
    #     """
    #     removes species of ogs that are set by user
    #     :return:
    #     """
    #     new_ogs = []
    #     for name, og in self.ogs:
    #         species = og.aa.description[og.aa.description.find("[") + 1:og.aa.description.find("]")]
    #         if species not in self.species_to_remove_mapping:
    #             new_ogs.apped(og)
    #     return new_ogs

    # def _has_species_to_remove(self):
    #     """
    #     :return: true or false depending whether species to remove are part of all the species present
    #     """
    #     not_included_species = []
    #     for spec in self.species_to_remove_mapping:
    #         if spec not in self.species_list:
    #             not_included_species.append(spec)
    #     return not_included_species

    def _get_species_list(self):
        """
        Use nwk string to return list of species
        :return:
        """
        tree_str = os.path.join(self.oma_output_path, 'EstimatedSpeciesTree.nwk')
        t = Tree(tree_str)
        return [leaf.name for leaf in t]

    def _reload_ogs_from_folder(self):
        """
        Re-load ogs if selection has finished and already exists in output folders
        :return: Dictionary with og name as key and list of SeqRecords
        """
        print('--- Re-load ogs and find their corresponding DNA seq from output folder ---')
        ogs = {}
        ref_ogs_aa = os.path.join(self.args.output_path, "01_ref_ogs_aa")
        ref_ogs_dna = os.path.join(self.args.output_path, "01_ref_ogs_dna")
        for file in tqdm(zip(sorted(glob.glob(os.path.join(ref_ogs_aa, "*.fa"))), sorted(glob.glob(os.path.join(ref_ogs_dna, "*.fa")))),desc='Re-loading files',unit=' OGs'):
            name = os.path.basename(file[0]).split(".")[0]
            ogs[name] = OG()
            ogs[name].aa = list(SeqIO.parse(file[0], format='fasta'))
            ogs[name].dna = list(SeqIO.parse(file[1], format='fasta'))
            # if self._remove_species:
            #     ogs[name].remove_species_records(self.species_to_remove, species_in_hog)
        return ogs

    # def _get_num_species_after_removal(self, species_in_og):
    #     in_og = 0
    #     for spec in self.species_to_remove_mapping:
    #         if spec in species_in_og:
    #             in_og += 1
    #
    #     return len(species_in_og)-in_og

    # def _change_record_id(self):
    #     raise NotImplementedError

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

        orthologous_groups_aa = os.path.join(self.args.output_path, "01_ref_ogs_aa")
        if not os.path.exists(orthologous_groups_aa):
            os.makedirs(orthologous_groups_aa)

        orthologous_groups_dna = os.path.join(self.args.output_path, "01_ref_ogs_dna")
        if not os.path.exists(orthologous_groups_dna):
            os.makedirs(orthologous_groups_dna)

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
                except:
                    print('This OG {} did not have any DNA'.format(name))
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

    def _get_dna_records(self, records, name):
        """
        
        :param records: 
        :return: 
        """
        og_cdna = []
        for i, record in enumerate(records):
            # species = record.description[record.description.find("[") + 1:record.description.find("]")]
            # if len(species.split(" ")) > 1:
            #     new_id = species.split(" ")[0][0:3] + species.split(" ")[1][0:2]
            #     species = new_id.upper()
            #-------------- only to be used internally -----------------
            if 'h5' in self._db_source:
                try:
                    oma_db_nr = self._db_id_map.omaid_to_entry_nr(record.id)
                except:
                    print('DNA not found for {}.'.format(record.id))
                    pass
                else:
                    og_cdna.append(SeqRecord.SeqRecord(Seq.Seq(self._db.get_cdna(oma_db_nr).decode("utf-8")),
                                                     id=record.id + "_" + name, description=""))
            elif 'fa' in self._db_source:
                try:
                    og_cdna.append(self._db[record.id])
                except ValueError:
                    print('DNA not found for {}.'.format(record.id))
                    pass
            elif 'REST_api' in self._db_source:
                try:
                    protein = requests.get(API_URL + "/protein/" + record.id + "/")
                except requests.exceptions.RequestException:
                    pass
                else:
                    protein = protein.json()
                    og_cdna.append(SeqRecord.SeqRecord(Seq.Seq(protein['cdna']),
                                                          id=record.id + "_" + name, description=""))

            if 'X' in str(og_cdna[-1].seq):
                cleaned_seq = self._clean_DNA_seq(og_cdna[-1])
                og_cdna[-1].seq = cleaned_seq


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


    def add_mapped_seq_v2(self, mapper, species_name=None):
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
            # remove species from the original set
            if self.args.remove_species_ogs:
                og = OG()
                filtered_og = value.remove_species_records(self.species_to_remove_ogs)
                if filtered_og:
                    og.dna = filtered_og[0]
                    og.aa = filtered_og[1]
            else:
                og = value
            if len(og.aa) > 2:
                # continue only if OG is in mapped OGs
                if name in mapped_og_set.keys():
                    if self.species_to_remove_mapping:  # in case we decided to remove species from the mapping
                        mapping_og = OG()
                        filtered_mapping = mapped_og_set[name].remove_species_records(self.species_to_remove_mapping)
                        if filtered_mapping:
                            mapping_og.dna = filtered_mapping[0]
                            mapping_og.aa = filtered_mapping[1]
                    else:  # nothing to remove
                        mapping_og = mapped_og_set[name]

                    if len(mapping_og.aa) >= 1:  # we had at least one mapped og even after removal
                        best_record_aa = mapping_og.get_best_mapping_by_seq_completeness(ref_og=og, threshold=self.args.sc_threshold)
                        if best_record_aa:
                            best_record_aa.id = species_name
                            self.mapped_ogs[name] = og
                            all_id = [rec.id for rec in self.mapped_ogs[name].aa]
                            if best_record_aa.id not in all_id:  # make sure that repeated run doesn't add the same sequence multiple times at the end of an OG
                                #print(mapper.all_cov)
                                cov.add_coverage(self._get_clean_id(best_record_aa), mapper.all_cov[self._get_clean_id(best_record_aa)])
                                seqC.add_seq_completeness(self._get_clean_id(best_record_aa), mapper.all_sc[self._get_clean_id(best_record_aa)])
                                self.mapped_ogs[name].aa.append(best_record_aa)
                        else:  # case where no best_record_aa reported because it was smaller than the self.args.sc_threshold
                            self.mapped_ogs[name] = og
                    else:  # mapping had only one that we removed
                        if self.args.keep_all_ogs:
                            self.mapped_ogs[name] = og
                else:  # nothing was mapped to that og
                    if self.args.keep_all_ogs:
                        self.mapped_ogs[name] = og

        cov.write_coverage_bam(os.path.join(self.args.output_path, species_name+'_all_cov.txt'))
        seqC.write_seq_completeness(os.path.join(self.args.output_path, species_name+'_all_sc.txt'))

    def write_added_ogs(self, folder_name=None):
        if folder_name is None:
            ogs_with_mapped_seq = os.path.join(self.args.output_path, "04_ogs_map_" + self._species_name)
        else:
            ogs_with_mapped_seq = os.path.join(self.args.output_path, folder_name)

        if not os.path.exists(ogs_with_mapped_seq):
            os.makedirs(ogs_with_mapped_seq)

        for name, value in self.ogs.items():
            if name in self.mapped_ogs.keys():
                output_file = os.path.join(ogs_with_mapped_seq, name + ".fa")
                self._write(output_file, self.mapped_ogs[name].aa)

    def _get_clean_id(self, record):
        des = record.description.split(" ")[0]
        des = des.split("_")
        return des[0]+"_"+des[1]

    def _write(self, file, value):
        """
        Write output to fasta file
        :param folder: file and location of outputfile
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
        seq_completenesses = self._get_seq_completeness_v2(ref_og=ref_og)
        best_record = seq_completenesses.index(max(seq_completenesses))
        if seq_completenesses[best_record] >= threshold:
            return self.aa[best_record]
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

    def _get_seq_completeness_v2(self, ref_og=None):
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
