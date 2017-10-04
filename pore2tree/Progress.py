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
OMA_MARKER_GENE_EXPORT = 'marker_genes'

class Progress(object):

    def __init__(self, args, oma_output):
        self.args = args
        if " " in args.reads:
            self._reads = args.reads.rstrip().split(" ")
        else:
            self._reads = args.reads

        if len(self._reads) == 2:
            self._species_name = self._reads[0].split("/")[-1].split(".")[0]
        else:
            self._species_name = self._reads.split("/")[-1].split(".")[0]

        if self.args.remove_species:
            self.species_to_remove = self.args.remove_species.split(",")
        else:
            self.species_to_remove = []
        self.status = None
        self.oma_output_path = self.args.oma_output_path
        self._num_ogs = oma_output.num_selected_ogs
        self._num_species = oma_output.num_species

        self.status = self._determine_progress()

    def _determine_progress(self):
        status = 0
        ref_ogs_aa = os.path.join(self.args.output_path, "01_ref_ogs_aa")
        ref_ogs_dna = os.path.join(self.args.output_path, "01_ref_ogs_dna")
        ref_dna = os.path.join(self.args.output_path, '02_ref_dna')
        mapping = os.path.join(self.args.output_path, "03_mapping_"+self._species_name)
        ogs_map = os.path.join(self.args.output_path, "04_ogs_map"+self._species_name)

        # check progress of OG selection
        comp_files = 0
        if os.path.exists(ref_ogs_aa) and os.path.exists(ref_ogs_dna):
            for file in zip(glob.glob(os.path.join(ref_ogs_aa, "*.fa")), glob.glob(os.path.join(ref_ogs_dna, "*.fa"))):
                if os.path.getsize(file[0]) > 0 and os.path.getsize(file[1]) > 0:
                    comp_files += 1

        if comp_files >= self._num_ogs:
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
        map_files = 0
        if os.path.exists(mapping):
            for file in glob.glob(os.path.join(mapping, "*consensus.fa")):
                if os.path.getsize(file) > 0:
                    map_files += 1

        if map_files == self._num_species and status == 2:
            status = 3

        return status

    # def _determine_num(self):
    #     og_orthoxml = os.path.join(self.oma_output_path, 'OrthologousGroups.orthoxml')
    #     tree_str = os.path.join(self.oma_output_path, 'EstimatedSpeciesTree.nwk')
    #
    #     ham_analysis = pyham.Ham(tree_str, og_orthoxml, use_internal_name=False)
    #
    #     num_select_ogs = 0
    #     num_species = len(ham_analysis.get_list_extant_genomes()) - len(self.species_to_remove)
    #     hog_dict = ham_analysis.get_dict_top_level_hogs()
    #     for hog, value in hog_dict.items():
    #         genes = self._remove_species(value.get_all_descendant_genes())
    #         if len(genes) >= self.args.min_species:
    #             num_select_ogs += 1
    #
    #     return [num_select_ogs, num_species]
    #
    # def _remove_species(self, genes_in_hog):
    #     return [gene for gene in genes_in_hog if gene.prot_id[0:5] not in self.species_to_remove]
