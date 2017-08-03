# Author: Ivana Pilizota
# Date: 1 November 2016

import logging
import os
import time

from Bio import SeqIO
from pyparsing import ParseException
import tempfile

from .base_treebuilder import TreeBuilder, AlignmentInput, DataType
from .parsers import FasttreeParser

from ..abstract_cli import AbstractCLI
from ..options import OptionSet, StringOption, IntegerOption
from ...file_utils import TempFile, TempDir

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


class FasttreeCLI(AbstractCLI):
    @property
    def _default_exe(self):
        return 'fasttree'


def set_default_dna_options(treebuilder):
    """
    Dummy function as sensible default
    """
    treebuilder.options = get_default_options()


def set_default_protein_options(treebuilder):
    """
    Dummy function as sensible default
    """
    treebuilder.options = get_default_options()


class Fasttree(TreeBuilder):

    def __init__(self, alignment, *args, **kwargs):
        self.options = get_default_options()
        super(Fasttree, self).__init__(alignment=alignment, *args, **kwargs)
        if self.input is not None:
            if self.datatype == DataType.DNA:
                set_default_dna_options(self)
            else:
                set_default_protein_options(self)

    def __call__(self, *args, **kwargs):
        """
        Sets up temporary output file location and calls FastTree using _call() function.
        Writes temporary input file if we're working with SeqIO object
        Saves the stdout and stderr and returns
        """
        start = time.time()  # time the execution
        if self.input_type == AlignmentInput.OBJECT:  # different operation depending on what it is
            with tempfile.NamedTemporaryFile(mode='wt') as fh:
                SeqIO.write(self.input, fh, 'phylip-relaxed') # default interleaved
                fh.seek(0)
                output, error = self._call(fh.name, *args, **kwargs)
                self.result = self._read_result(output, error)  # store result
        else:
            filename = os.path.abspath(self.input)
            output, error = self._call(filename, *args, **kwargs)
            self.result = self._read_result(output, error)  # store result

        end = time.time()
        self.elapsed_time = end - start
        return self.result["tree"]
        # End call

    # Any other accessory methods
    def _call(self, filename, *args, **kwargs):
        """
        Call underlying low level FastTree wrapper.
        Options are passed via *args and **kwargs
        [This only covers the simplest automatic
         case]
        """
        #hard code tmp_output as the output name since we don't save it anyway
        #self.cli('{} -log {log_output} {seqfile} > {tmp_path}'.format(self.command(), tmp_path=os.path.join(tmpd,'tmp_output'), log_output=logfile, seqfile=filename), wait=True)
        self.cli('{} {seq_file}'.format(self.command(), seq_file=filename), wait=True)

        return (self.cli.get_stdout(), self.cli.get_stderr())

    def command(self):
        return str(self.options)

    def _read_result(self, stdout, stderr):
        """
        Read back the result.
        """
        parser = FasttreeParser()

        try:
            parser.parse(tree=stdout, other=stderr)
            result = parser.to_dict()
        except IOError as ioerr:
            logger.error('Error reading results')
            result = None
        except ParseException as parseerr:
            logger.error('Other parse error', parseerr)
            result = None

        return result

    def _init_cli(self, binary):
        return FasttreeCLI(executable=binary)


def get_default_options():

    return OptionSet([
        # Algorithm

        # Set datatype to DNA (nt) or AA alignment: AA by default. If set to True will assume DNA format.
        StringOption('-nt', active=False),

        # Set the WAG model for AA alignment. Default Jones-Taylor-Thorton
        StringOption('-wag', active=False),

        # Set the GTR model for nt alignment. Default Jones-Taylor-Thorton
        StringOption('-gtr', active=False),

        # Set the gamma model. Default Jones-Taylor-Thorton
        StringOption('-gamma', active=False),

        # Specify the number of rate categories of sites. Default 20.
        IntegerOption('-cat', 20, active=False),

        # Specify starting tree
        StringOption('-intree', '', active=False),

        # Speed up the neighbor joining phase & reduce memory usage (recommended for >50,000 sequences)
        StringOption('-fastest', active=False),

        # Set the number of rounds of maximum-likelihood NNIs. Deafault 4*log2(N), N = the number of unique sequences
        IntegerOption('-mlnni', 0, active=False),

    ])
