#!/usr/bin/env python
'''
    This class will check the current status of the pipeline computation and determine how to start once several steps have been finished and
    will save a lot of time when several steps are finished. Also it will allow to submit the mapping (part that takes longest time) as job array.

    -- David Dylus, July--XXX 2017
'''

import pyham
import glob
import os

OMA_STANDALONE_OUTPUT = 'Output'

class Progress(object):

    def __init__(self, args):
        self.args = args
        self.status = None
        self.oma_output_path = os.path.join(self.args.standalone_path, OMA_STANDALONE_OUTPUT)
        self._num_ogs = self._determine_num()[0]
        self._num_species = self._determine_num()[1]

        self.status = self._determine_progress()

    def _determine_progress(self):
        status = 0
        ref_ogs_aa = os.path.join(self.args.output_path, "01_ref_ogs_aa")
        ref_ogs_dna = os.path.join(self.args.output_path, "01_ref_ogs_dna")
        ref_dna = os.path.join(self.args.output_path, '02_ref_dna')
        mapping = os.path.join(self.args.output_path, "03_mapping")
        ogs_map = os.path.join(self.args.output_path, "04_ogs_map")

        # check progress of OG selection
        comp_files = 0
        if os.path.exists(ref_ogs_aa) and os.path.exists(ref_ogs_dna):
            for file in zip(glob.glob(os.path.join(ref_ogs_aa, "*.fa")), glob.glob(os.path.join(ref_ogs_dna, "*.fa"))):
                if os.path.getsize(file[0]) > 0 and os.path.getsize(file[1]) > 0:
                    comp_files += 1

        if comp_files == self._num_ogs:
            status = 1

        # check progress of Ref generation
        ref_files = 0
        if os.path.exists(ref_dna):
            for file in glob.glob(os.path.join(ref_dna, "*.fa")):
                if os.path.getsize(file) > 0:
                    ref_files += 1

        if ref_files == self._num_species and status == 1:
            status = 2

        # check progress of mapping

        return status

    def _determine_num(self):
        og_orthoxml = os.path.join(self.oma_output_path, 'OrthologousGroups.orthoxml')
        tree_str = os.path.join(self.oma_output_path, 'EstimatedSpeciesTree.nwk')

        ham_analysis = pyham.Ham(tree_str, og_orthoxml, use_internal_name=False)

        num_select_ogs = 0
        num_species = len(ham_analysis.get_list_extant_genomes())
        hog_dict = ham_analysis.get_dict_top_level_hogs()
        for hog, value in hog_dict.items():
            if len(value.get_all_descendant_genes()) >= self.args.min_species:
                num_select_ogs += 1

        return [num_select_ogs, num_species]
