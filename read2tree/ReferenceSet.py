#!/usr/bin/env python
'''
    This file contains definitions of a class which allows to create
    the reference orthologous groups with their DNA sequences.

    -- David Dylus, July--XXX 2017
'''

import os
import glob
import logging
import time
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import sys


class ReferenceSet(object):
    '''
    Structure for reference
    '''

    def __init__(self, args, og_set=None, step=None):
        """

        :param args: list of arguments from command line
        :param og_set: set of OGs used to obtain reference DNA sequences
        :param load: set to True when reference loaded from folder/file of list of arguments
        """
        self.ref = {}
        #self.load = load
        self.args = args
        self.step = step

        self.logger = logging.getLogger(__name__)
        self._species_name = self.args.species_name

        #if load is False:
        if step == "2map":
            self.ref = self._load_records_folder()
        #elif og_set is not None and load is True:
        elif step == "all" or step == "1marker":  #
            self.ref = self._generate_reference(og_set)
            self.write()
            # self.progress.set_status('ref')

        # if args.remove_species:
        #     self.ref = self._remove_species()

    def _read_fasta(self, ref_file):
        '''

        :param ref_file: file that contains all the DNA sequences from the oma database
        :return:
        '''
        self.logger.info('--- Reading DNA reference into memory ---')
        return SeqIO.index(ref_file, "fasta")

    def _load_records_folder(self):
        """
        Parse species with their dna sequences from folder
        :return:
        """
        ref_dict = {}
        self.logger.info('--- Generating reference for mapping from folder ---')
        ref_dna = os.path.join(self.args.output_path, '02_ref_dna')
        if not os.path.exists(ref_dna):
            self.logger.info('The ref_dna folder {} does not exist. Step 1 was incomplete probably.'.format(str(ref_dna)))
            sys.exit(1)


        for file in tqdm(glob.glob(os.path.join(ref_dna, "*.fa")), desc="Re-loading references for mapping from folder", unit=" species"):
            species_name = file.split("/")[-1].split("_")[0]
            ref_dict[species_name] = Reference()
            ref_dict[species_name].dna = list(SeqIO.parse(file, 'fasta'))

        return ref_dict

    def _generate_reference(self, og_set):
        '''
        Split records into dictionary with keys being species and the values the corresponded sequence records
        '''
        self.logger.info('--- Generating reference for mapping ---')
        start = time.time()
        ref_set = {}
        for name, og in tqdm(og_set.items(), desc="Loading records", unit=" record"):
            for record in og.aa:
                species = record.id[0:5]
                record.id = record.id  # +"_"+name
                if species in ref_set.keys():
                    ref_set[species].aa.append(record)
                else:
                    ref_set[species] = Reference()
                    ref_set[species].aa.append(record)

            for record in og.dna:
                species = record.id[0:5]
                record.id = record.id  # + "_" + name
                if species in ref_set.keys():
                    ref_set[species].dna.append(record)
                else:
                    ref_set[species] = Reference()
                    ref_set[species].dna.append(record)
        end = time.time()
        elapsed_time = end - start
        self.logger.info('{}: Extracted {} reference species form {} ogs took {}'
                       .format(self._species_name, len(ref_set.keys()),
                       len(og_set.keys()), elapsed_time))
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
            if value.dna:  # only write if not empty
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
