#!/usr/bin/env python
'''
    This file contains definitions of a class which allows to create 
    the reference orthologous groups with their DNA sequences.

    -- David Dylus, July--XXX 2017
'''

import re
import os
import glob
from tqdm import tqdm
from itertools import chain
from Bio import SeqIO, Seq, SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Alphabet import SingleLetterAlphabet
from tables import *
from pyoma.browser import db



class ReferenceSet(object):
    '''
    Structure for reference
    '''

    def __init__(self, args, og_set=None, load=True):
        """
        
        :param args: list of arguments from command line
        :param og_set: set of OGs used to obtain reference DNA sequences
        :param load: set to True when reference loaded from folder/file of list of arguments
        """
        self.ref = {}
        self.load = load
        self.args = args

        if load is False:
            self.ref = self._load_records_folder()
        elif og_set is not None and load is True:
                self.ref = self._generate_reference(og_set)
                self.write()

        # if args.remove_species:
        #     self.ref = self._remove_species()

    def _read_fasta(self, ref_file):
        '''
        
        :param ref_file: file that contains all the DNA sequences from the oma database
        :return: 
        '''
        print('--- Reading DNA reference into memory ---')
        return SeqIO.index(ref_file, "fasta")

    def _load_records_folder(self):
        """
        Parse species with their dna sequences from folder
        :return: 
        """
        ref_dict = {}
        print('--- Generating reference for mapping from folder ---')
        ref_dna = os.path.join(self.args.output_path, '02_ref_dna')
        for file in tqdm(glob.glob(os.path.join(ref_dna, "*.fa")), desc="Loading references for mapping from folder", unit=" species"):
            species_name = file.split("/")[-1].split("_")[0]
            ref_dict[species_name] = Reference()
            ref_dict[species_name].dna = list(SeqIO.parse(file, 'fasta'))

        return ref_dict

    def _generate_reference(self, og_set):
        '''
        Split records into dictionary with keys being species and the values the corresponded sequence records
        '''
        print('--- Generating reference for mapping ---')
        ref_set = {}
        for name, og in tqdm(og_set.items(), desc="Loading records", unit=" record"):
            for record in zip(og.aa, og.dna):
                species = record[0].id[0:5]
                record[0].id = record[0].id+"_"+name
                record[1].id = record[1].id+"_"+name
                if species in ref_set.keys():
                    ref_set[species].aa.append(record[0])
                    ref_set[species].dna.append(record[1])
                else:
                    ref_set[species] = Reference()
                    ref_set[species].aa.append(record[0])
                    ref_set[species].dna.append(record[1])

        return ref_set

    def write(self):
        '''
        Write for each species all the DNA sequences into separate fasta files
        :param output_folder: folder where files should be stored
        '''
        out_dna = os.path.join(self.args.output_path, '02_ref_dna')
        if not os.path.exists(out_dna):
            os.makedirs(out_dna)
        for key, value in self.ref.items():
            value.write_dna(key, out_dna)

    def _remove_species(self):
        raise NotImplementedError


class Reference(object):

    def __init__(self, args=None):
        self.args = args
        self.aa = []
        self.dna = []

    def write_aa(self, species, output_folder):
        handle = open(os.path.join(output_folder, species + '_OGs.fa'), "w")
        writer = FastaWriter(handle, wrap=None)
        writer.write_file(self.aa)
        handle.close()

    def write_dna(self, species, output_folder):
        handle = open(os.path.join(output_folder, species + '_OGs.fa'), "w")
        writer = FastaWriter(handle, wrap=None)
        writer.write_file(self.dna)
        handle.close()