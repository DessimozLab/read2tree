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

    def __init__(self, args, ogset=None, load=True):
        """
        
        :param args: list of arguments from command line
        :param ogset: set of OGs used to obtain reference DNA sequences
        :param load: set to True when reference loaded from folder/file of list of arguments
        """
        self.ref = {}
        self.load = load
        self.args = args

        if not load:
            self.ref = self._load_records_folder()
        else:
            if ogset is not None:
                self.ref = self._generate_reference(ogset)

    def _read_fasta(self, ref_file):
        '''
        
        :param ref_file: file that contains all the DNA sequences from the oma database
        :return: 
        '''
        print('--- Reading DNA reference into memory ---')
        return SeqIO.index(ref_file, "fasta")

    def _load_records_folder(self):
        """
        
        :return: 
        """
        ref_dict = {}
        for file in glob.glob(os.path.join(self.args.ref_folder, "*.fa")):
            species_name = file.split("/")[-1].split("_")[0]
            ref_dict[species_name] = Reference()
            ref_dict[species_name].dna = ref_dict[species_name]

        return ref_dict

    def _generate_reference(self, og_set):
        '''
        Split records into dictionary with keys being species and the values the corresponded sequence records
        '''
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

    def write(self, output_folder):
        '''
        Write for each species all the DNA sequences into separate fasta files
        :param output_folder: folder where files should be stored
        '''
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        for key, value in tqdm(self.ref.items(), desc="Writing OGs sorted by species", unit=" species"):
            self.value.write_aa(key, output_folder)
            self.value.write_dna(key, output_folder)

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

    def write_dna(self):
        handle = open(os.path.join(output_folder, species + '_OGs.fa'), "w")
        writer = FastaWriter(handle, wrap=None)
        writer.write_file(self.dna)
        handle.close()