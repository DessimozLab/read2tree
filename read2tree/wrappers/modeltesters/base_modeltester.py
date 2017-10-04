import os, types, itertools
from abc import ABCMeta, abstractmethod
from enum import Enum
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from ...seq_utils import is_dna

from zoo.wrappers import WrapperError
from zoo.wrappers.aligners.base_aligner import identify_input

import logging
logger = logging.getLogger(__name__)

AlignmentInput = Enum('AlignmentInput', 'OBJECT FILENAME')
DataType = Enum('DataType', 'DNA PROTEIN UNKNOWN')


class ModelTester(object):
    """
    Base class for wrappers of model testers for phylogeny inference

    The wrapper is written as a callable class.
    This can hold data (state) to do with the operation it performs, so it can keep results,
    execution times and other metadata, as well as perform the task.

    This is a base implementation to be extended. The important parts are
    __init__ (does the setup) and __call__ (does the work). All
    else are helper methods.

    :Example:

    ::

        callable_wrapper = ConcreteModelTester(aln)
        result = callable_wrapper()
        time_taken = callable_wrapper.elapsed_time
        result_again = callable_wrapper.result
    """
    __metaclass__ = ABCMeta

    def __init__(self, alignment=None, datatype=DataType.UNKNOWN, binary=None):
        """
        ..note::  TODO: this documentation is not correct. it needs to be updateted.

        Should work the same whether you're working with a Biopython object or a file
            but the implementation differs, e.g. a Biopython object will need
            to be written temporarily to disk for the Aligner to work on it.

        alignment is one of 4 things:
            a filename
            a Biopython MSA
            a list of Seq objects
            anything else (throw an exception)

        binary is the alignment's executable file, or None
        """

        if alignment is not None:
            self.input_type = identify_input(alignment)  # Figure out what it is - file or object
            if datatype == DataType.UNKNOWN:
                # dup, input_ = itertools.tee(input_)
                self.datatype = guess_datatype(alignment, from_filename=self.input_type == AlignmentInput.FILENAME)
            else:
                self.datatype = datatype

            self.input = alignment  # store it
        else:
            self.input_type = None
            self.input = None


        self.elapsed_time = None
        self.stdout = None
        self.stderr = None
        self.cli = self._init_cli(binary)
        #TODO: the wrapper error is not compatible with calling a function with java!
        #try:
        #    self.cli = self._init_cli(binary)
        #except IOError as err:
        #     raise WrapperError('Error searching for binary: {}'.format(err))
            # End setup

    @abstractmethod
    def __call__(self, *args, **kwargs):
        """
        How to call the underlying aligner
        """
        pass

    @abstractmethod
    def _init_cli(self, binary):
        """
        Set up the command-line interface to the wrapped software
        :param binary: filename of executable binary file
        :return: concrete CLI type inheriting from AbstractCLI
        """
        pass


def guess_datatype(alignment, from_filename=False):
    logger.warning("Guessing is not recommended - specify the sequence type with option datatype={DNA, PROTEIN}, be more confident")
    if from_filename:
        try:
            alignment = list(SeqIO.parse(alignment, 'fasta'))
        except:
            alignment = list(SeqIO.parse(alignment, 'phylip-relaxed'))
    return DataType.DNA if is_dna(alignment) else DataType.PROTEIN


def identify_input(alignment):
    """
    Work out if we're dealing with an alignment (return True), a file
    (return False), or invalid input (raise error)

    :param alignment: either an Biopython MultipleSequenceAlignment or
        a filename pointing to an existing msa file.
    """
    try:
        if isinstance(alignment, (MultipleSeqAlignment, types.GeneratorType, list)):
            # `alignment` is a Biopython MultipleSequenceAlignment
            return AlignmentInput.OBJECT

        elif isinstance(alignment, str) and os.path.exists(alignment):
            # `alignment` is a filepath
            return AlignmentInput.FILENAME

    except:
        # `alignment` is some other thing we can't handle
        raise ValueError('{} is not an alignment object or a valid filename'.format(alignment))


