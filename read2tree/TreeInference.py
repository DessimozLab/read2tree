#!/usr/bin/env python
'''
    This file contains definitions of a class which surrounds the wrappers to build trees given a set of command line arguments.

    -- David Dylus, July--XXX 2017
'''
import os
import time
import logging
from read2tree.wrappers.treebuilders import Fasttree
from read2tree.wrappers.treebuilders import Iqtree


logger = logging.getLogger(__name__)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('info.log')
file_handler.setFormatter(formatter)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)


class TreeInference(object):

    def __init__(self, args, concat_alignment=None):
        print('--- Tree inference ---')

        self.args = args

        self.elapsed_time = 0

        if args.debug:
            logger.setLevel(logging.DEBUG)
            file_handler.setLevel(logging.DEBUG)
            # stream_handler.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)
            file_handler.setLevel(logging.INFO)
            # stream_handler.setLevel(logging.INFO)

        logger.addHandler(file_handler)
        # logger.addHandler(stream_handler)

        if self.args.reads:
            if len(self.args.reads) == 2:
                self._reads = self.args.reads
                self._species_name = self._reads[0].split("/")[-1].split(".")[0]
            else:
                self._reads = self.args.reads[0]
                self._species_name = self._reads.split("/")[-1].split(".")[0]

        if self.args.species_name:
            self._species_name = self.args.species_name

        if not self.args.reads and not self.args.species_name:
            self._species_name = 'merge'

        self.tree = None
        if concat_alignment is not None:
            self.tree = self._infer_tree(concat_alignment)

    def _infer_tree(self, concat_alignment):
        start = time.time()
        output_folder = self.args.output_path
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        #fasttree_wrapper = Fasttree(concat_alignment, datatype="PROTEIN")
        #tree = fasttree_wrapper()
        iqtree_wrapper = Iqtree(concat_alignment, datatype="PROTEIN")
        iqtree_wrapper.options.options['-m'].set_value('LG')
        iqtree_wrapper.options.options['-nt'].set_value(self.args.threads)
        tree = iqtree_wrapper()
        with open(os.path.join(output_folder, "tree_" + self._species_name + ".nwk"), "w") as text_file:
            text_file.write("{}".format(tree))
        self.tree = "{}".format(tree)
        end = time.time()
        self.elapsed_time = end - start
        logger.info('{}: Tree inference took {}.'.format(self._species_name,
                                                         self.elapsed_time))

        return tree
