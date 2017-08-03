#!/usr/bin/env python
'''
    This file contains definitions of a class which surrounds possible alignment methods

    -- David Dylus, July--XXX 2017
'''
import os
from tables import *
from Bio import AlignIO

from tqdm import tqdm
from pore2tree.wrappers.aligners import Mafft
from pore2tree.utils.seq_utils import concatenate


class Aligner(object):

    def __init__(self, args, og_set=None):
        print('--- Alignment of OGs ---')
        self.args = args
        self.alignments = {}

        if og_set is not None:
            self.alignments = self._align(og_set)

    def _align(self, og_set):
        align_dict = {}
        self._adapt_id(og_set)
        output_folder = os.path.join(self.args.output_path, "05_align")
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
