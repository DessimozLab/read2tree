#!/usr/bin/env python
'''
    This file contains definitions of a class which holds the orthologous groups.

    -- David Dylus, July--XXX 2017
'''
import glob
import os
import re
import pyham

from tqdm import tqdm
from Bio import SeqIO, Seq, SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
from Bio.SeqIO.FastaIO import FastaWriter
from tables import *
from pyoma.browser import db


OMA_STANDALONE_OUTPUT = 'Output'

class OGSet(object):

    def __init__(self, args, load=True):
        self.args = args
        self.ogs = {}
        self.mapped_ogs = {}
        self.min_species = self._estimate_best_number_species()

        if self.args.standalone_path:
            self.oma_output_path = os.path.join(self.args.standalone_path, OMA_STANDALONE_OUTPUT)
        if '.fa' in self.args.dna_reference or '.fasta' in self.args.dna_reference:
            print(
                'Loading {} into memory. This might take a while . . . '.format(self.args.dna_reference.split("/")[-1]))
            self._db = SeqIO.index(self.args.dna_reference, "fasta")
            self._db_source = 'fa'
        elif '.h5' in self.args.dna_reference:
            self._db = db.Database(self.args.dna_reference)
            self._db_id_map = db.OmaIdMapper(self._db)
            self._db_source = 'h5'
        else:
            raise FileNotFoundError

            # self.records = {}
        if load:
            self.ogs = self._load_ogs()
        else:
            self.ogs = self._load_ogs_from_folder()

    def og(self, name):
        return self.ogs[name]

    def _load_ogs_from_folder(self):
        """
        Re-load ogs if selection has finished and already exists in output folders
        :return: Dictionary with og name as key and list of SeqRecords
        """
        print('--- Load ogs and find their corresponding DNA seq from output folder ---')
        ogs = {}
        ref_ogs_aa = os.path.join(self.args.output_path, "01_ref_ogs_aa")
        ref_ogs_dna = os.path.join(self.args.output_path, "01_ref_ogs_dna")
        for file in tqdm(zip(glob.glob(os.path.join(ref_ogs_aa, "*.fa")), glob.glob(os.path.join(ref_ogs_dna, "*.fa"))), desc='Loading files',
                         unit=' OGs'):
            name = file[0].split("/")[-1].split(".")[0]
            ogs[name] = OG()
            ogs[name].aa = list(SeqIO.parse(file[0], format='fasta'))
            ogs[name].dna = list(SeqIO.parse(file[1], format='fasta'))
        return ogs

    def _load_ogs(self):
        """
        Using the orthoxml file select only the OGs of interest that have more species than the min_species threshold
        :return: Dictionary with og name as key and list of SeqRecords
        """
        print('--- Load ogs and find their corresponding DNA seq from a database ---')
        ogs = {}

        orthologous_groups_aa = os.path.join(self.args.output_path, "01_ref_ogs_aa")
        if not os.path.exists(orthologous_groups_aa):
            os.makedirs(orthologous_groups_aa)

        orthologous_groups_dna = os.path.join(self.args.output_path, "01_ref_ogs_dna")
        if not os.path.exists(orthologous_groups_dna):
            os.makedirs(orthologous_groups_dna)

        orthologous_groups_fasta = os.path.join(self.oma_output_path, "OrthologousGroupsFasta")
        og_orthoxml = os.path.join(self.oma_output_path, 'OrthologousGroups.orthoxml')
        tree_str = os.path.join(self.oma_output_path, 'EstimatedSpeciesTree.nwk')

        ham_analysis = pyham.Ham(tree_str, og_orthoxml, use_internal_name=False)

        for file in tqdm((glob.glob(os.path.join(orthologous_groups_fasta, "*.fa")) or glob.glob(os.path.join(orthologous_groups_fasta, "*.fasta"))), desc='Loading files',
                         unit=' OGs'):
            name = file.split("/")[-1].split(".")[0]
            og_ham = ham_analysis.get_hog_by_id(name[2:])
            #     og_new = [for gene in og_ham]
            if len(og_ham.get_all_descendant_genes()) >= self.min_species:
                records = list(SeqIO.parse(file, 'fasta'))
                ogs[name] = OG()
                ogs[name].aa = self._get_aa_records(og_ham, records)
                output_file_aa = os.path.join(orthologous_groups_aa, name+".fa")
                output_file_dna = os.path.join(orthologous_groups_dna, name+".fa")
                self._write(output_file_aa, ogs[name].aa)
                if self.args.dna_reference:
                    ogs[name].dna = self._get_dna_records(ogs[name].aa, name)
                else:
                    print("DNA reference was not provided. Only amino acid sequences gathered!")
                self._write(output_file_dna, ogs[name].dna)
        return ogs

    def _get_aa_records(self, og_ham, records):
        """
        
        :param og_ham: 
        :param records: 
        :return: 
        """
        prot_ids = [gene.prot_id.split(" | ")[0] for gene in og_ham.get_all_descendant_genes()]
        for record in records:
            mystr = record.description
            record.id = [x for x in prot_ids if mystr[mystr.find("[") + 1:mystr.find("]")] in x[0:5]][0]
        return records

    def _get_dna_records(self, records, name):
        """
        
        :param records: 
        :return: 
        """
        og_cdna = [None] * len(records)
        for i, record in enumerate(records):
            if 'h5' in self._db_source:
                oma_db_nr = self._db_id_map.omaid_to_entry_nr(record.id)
                og_cdna[i] = SeqRecord.SeqRecord(Seq.Seq(self._db.get_cdna(oma_db_nr).decode("utf-8")),
                                                 id=record.id + "_" + name, description="")
            elif 'fa' in self._db_source:
                og_cdna[i] = self._db[record.id]

            if 'X' in str(og_cdna[i].seq):
                cleaned_seq = self._clean_DNA_seq(og_cdna[i])
                og_cdna[i].seq = cleaned_seq

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

    def add_mapped_seq(self, mapped_og_set):
        """
        Add the sequence given from the read mapping to its corresponding OG
        :param mapped_og_set: set of ogs with its mapped sequences
        """
        ogs_with_mapped_seq = os.path.join(self.args.output_path, "04_ogs_map")
        if not os.path.exists(ogs_with_mapped_seq):
            os.makedirs(ogs_with_mapped_seq)

        for name, value in mapped_og_set.items():
            best_record_aa = value.get_best_mapping_by_coverage()
            self.mapped_ogs[name] = self.ogs[name]
            self.mapped_ogs[name].aa.append(best_record_aa)
            output_file = os.path.join(ogs_with_mapped_seq, name+".fa")
            self._write(output_file, self.mapped_ogs[name].aa)

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

    def get_best_mapping_by_coverage(self, gene_code='aa'):
        coverages = self._get_coverage(gene_code=gene_code)
        best_record = coverages.index(max(coverages))
        return self.aa[best_record]

    def _get_coverage(self, gene_code='aa'):
        coverage = []
        if gene_code is 'dna':
            for record in self.dna:
                seq_len = len(record)
                non_n_len = len(record) - record.seq.count('n')
                coverage.append(non_n_len / seq_len)
        elif gene_code is 'aa':
            for record in self.dna:
                seq_len = len(record)
                non_n_len = len(record) - record.seq.count('X')
                coverage.append(non_n_len / seq_len)
        return coverage

    def remove_species_records(self, species_list):
        '''
        Remove species from reference sequence set
        :param species_list: list of species to be removed
        '''
        print('--- Removing {} species of list ---'.format(len(species_list )))
        for key, value in tqdm(self.records.items(), desc="Loading records", unit=" record"):
            for i, record in enumerate(value):
                if record.id[0:5] in species_list:
                    del value[i]