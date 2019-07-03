#!/usr/bin/env python
'''
    This file contains definitions of a class which surrounds possible alignment methods

    -- David Dylus, July--XXX 2017
'''
import logging

logger = logging.getLogger(__name__)

class Aligner(object):

    def __init__(self, args=None, alignments=None):

        self.args = args
        self.alignments = alignments
        self.placement_dic = alignments.placement_dic


    def