#!/usr/bin/env python
'''
    This file contains definitions of a class which surrounds the wrappers to build trees given a set of command line arguments.

    -- David Dylus, July--XXX 2017
'''
import os
from read2tree.wrappers.treebuilders import Fasttree


class TreeInference(object):

    def __init__(self, args, concat_alignment=None):
        print('--- Tree inference ---')
        self.args = args

        if " " in args.reads:
            self._reads = args.reads.rstrip().split(" ")
        else:
            self._reads = args.reads

        if len(self._reads) == 2:
            self._species_name = self._reads[0].split("/")[-1].split(".")[0]
        else:
            self._species_name = self._reads.split("/")[-1].split(".")[0]

        self.tree = None
        if concat_alignment is not None:
            self.tree = self._infer_tree(concat_alignment)

    def _infer_tree(self, concat_alignment):
        output_folder = self.args.output_path
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        fasttree_wrapper = Fasttree(concat_alignment, datatype="PROTEIN")
        tree = fasttree_wrapper()
        with open(os.path.join(output_folder, "tree_" + self._species_name + ".nwk"), "w") as text_file:
            text_file.write("{};".format(tree))
        self.tree = "{};".format(tree)
        return tree
