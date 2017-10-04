import os, types
from abc import ABCMeta, abstractmethod
from enum import Enum
from Bio.SeqRecord import SeqRecord
from read2tree.wrappers import WrapperError

import logging
logger = logging.getLogger(__name__)

ReferenceInput = Enum('ReferenceInput', 'OBJECT FILENAME')
ReadInput = Enum('ReadInput', 'OBJECT FILENAME')

class ReadMapper(object):
    """
    Base class for wrappers of read mapping software

    The wrapper is written as a callable class.
    This can hold data (state) to do with the operation it performs, so it can keep results,
    execution times and other metadata, as well as perform the task.

    This is a base implementation to be extended. The important parts are
    __init__ (does the setup) and __call__ (does the work). All
    else are helper methods.

    :Example:

    ::

        callable_wrapper = ConcreteAligner(aln)
        result = callable_wrapper()
        time_taken = callable_wrapper.elapsed_time
        result_again = callable_wrapper.result
    """
    __metaclass__ = ABCMeta

    def __init__(self, reference=None, reads=None, binary=None):
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
        if reference is not None:
            self.ref_input_type = identify_reference(reference)  # Figure out what it is - file or object
            self.ref_input = reference  # store it
        else:
            self.ref_input_type = None
            self.ref_input = None

        if reads is not None:
            self.read_input_type = identify_reads(reads)  # Figure out what it is - file or object
            self.read_input = reads  # store it
        else:
            self.read_input_type = None
            self.read_input = None

        self.elapsed_time = None
        self.stdout = None
        self.stderr = None
        try:
            self.cli = self._init_cli(binary)
        except IOError as err:
            raise WrapperError('Error searching for binary: {}'.format(err))
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

def identify_reference(sequence):
    """
    Work out if we're dealing with a fasta (return True), a file
    (return False), or invalid input (raise error)

    :param alignment: either an Biopython MultipleSequenceAlignment or
        a filename pointing to an existing msa file.
    """
    try:
        if isinstance(sequence, (SeqRecord, types.GeneratorType, list)):
            # `sequence` is a Biopython MultipleSequenceAlignment
            return ReferenceInput.OBJECT

        elif isinstance(sequence, str) and os.path.exists(sequence):
            # `sequence` is a filepath
            return ReferenceInput.FILENAME

    except:
        # `sequence` is some other thing we can't handle
        raise ValueError('{} is not an sequence object or a valid filename'.format(sequence))


def identify_reads(reads):
    """
    Work out if we're dealing with a fasta (return True), a file
    (return False), or invalid input (raise error)

    :param alignment: either an Biopython MultipleSequenceAlignment or
        a filename pointing to an existing msa file.
    """
    if isinstance(reads, list):
        read = reads[0]
    else:
        read = reads

    try:
        if isinstance(read, (SeqRecord, types.GeneratorType, list)):
            # `sequence` is a Biopython MultipleSequenceAlignment
            return ReadInput.OBJECT

        elif isinstance(read, str) and os.path.exists(read):
            # `sequence` is a filepath
            return ReadInput.FILENAME

    except:
        # `sequence` is some other thing we can't handle
        raise ValueError('{} is not an sequence object or a valid filename'.format(sequence))


