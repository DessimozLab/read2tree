#!/usr/bin/env python
'''
    This file contains definitions of a class which surrounds possible alignment methods

    -- David Dylus, July--XXX 2017
'''
import os
from tables import *
from Bio import AlignIO

from tqdm import tqdm
from read2tree.wrappers.aligners import Mafft
from read2tree.utils.seq_utils import concatenate


class Aligner(object):

    def __init__(self, args, og_set=None):
        print('--- Alignment of OGs ---')
        self.args = args

        if " " in args.reads:
            self._reads = args.reads.rstrip().split(" ")
        else:
            self._reads = args.reads

        if len(self._reads) == 2:
            self._species_name = self._reads[0].split("/")[-1].split(".")[0]
        else:
            self._species_name = self._reads.split("/")[-1].split(".")[0]

        self.alignments = {}

        if og_set is not None:
            self.alignments = self._align(og_set)

    def _align(self, og_set):
        align_dict = {}
        self._adapt_id(og_set)
        output_folder = os.path.join(self.args.output_path, "05_align_" + self._species_name)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        for key, value in tqdm(og_set.items(), desc='Aligning OGs ', unit=' OG'):
            mafft_wrapper = Mafft(value.aa, datatype="PROTEIN")
            mafft_wrapper.options.options['--localpair'].set_value(True)
            mafft_wrapper.options.options['--maxiterate'].set_value(1000)
            alignment = mafft_wrapper()
            align_dict[key] = alignment

            og_name = key.split("/")[-1]
            output_handle = open(os.path.join(output_folder, og_name + ".phy"), "w")
            AlignIO.write(alignment, output_handle, "phylip-relaxed")

        return align_dict

    def _adapt_id(selfs, og_set):
        for key, value in og_set.items():
            for record in value.aa:
                s = record.id[0:5]
                record.id = s

    def concat_alignment(self):
        alignments = []
        for key, value in self.alignments.items():
            alignments.append(value)

        return concatenate(alignments)
