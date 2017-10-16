#!/usr/bin/env python
'''
    This file contains definitions of a class which surrounds possible alignment methods

    -- David Dylus, July--XXX 2017
'''
import os
from tables import *
from Bio import AlignIO
import re

from tqdm import tqdm
from read2tree.wrappers.aligners import Mafft
from read2tree.utils.seq_utils import concatenate


class Analyzer(object):

    def __init__(self, args, og_set=None):
        print('--- Alignment of OGs ---')
        self.args = args
        self._genome_or_transcriptome_length = args.gt_length

        if " " in args.reads:
            self._reads = args.reads.rstrip().split(" ")
        else:
            self._reads = args.reads

        if len(self._reads) == 2:
            self._species_name = self._reads[0].split("/")[-1].split(".")[0]
        else:
            self._species_name = self._reads.split("/")[-1].split(".")[0]

        self.treeStats = {}
        self.alignmentStats = {}

    # def __call__(self, *args, **kwargs):
    #     raise NotImplementedError

    def _get_coverage_reads(self, args):
        """

        :param args:
        :return: coverage
        """
        with open(args.reads[0]) as input:
            read_length = input.readline().split("length=")[-1]
            num_lines = sum([1 for line in input])

        total_records = int(num_lines / 4)
        coverage = (total_records * read_length * len(args.reads))/self._genome_or_transcriptome_length
        return coverage

    def _get_number_results(self):
        raise NotImplementedError

    def _get_rf_dist(self, ref_tree):
        raise NotImplementedError

    def _get_length_align(self):
        raise NotImplementedError

    def _get_num_OGs(self):
        raise NotImplementedError

    def _get_mean_ACGT(self, args):
        for folder in glob.iglob(args.output + '/05_*', recursive=True):
            print(folder)
            all_coverages = []

            for file in glob.iglob(folder + '/*.phy'):
                align = AlignIO.read(file, "phylip-relaxed")
                for record in align:
                    if self._species_name[0:5] in record.id:
                        seq = re.sub('-', '', str(record.seq))
                        xx = seq.count("X")
                        aa = len(seq) - xx
                        all_coverages.append((aa / len(seq)))
            print(sum(all_coverages) / len(all_coverages))

    def _get_branch_length_mapped_seq(self):
        raise NotImplementedError

    def write_to_csv(self):
        raise NotImplementedError