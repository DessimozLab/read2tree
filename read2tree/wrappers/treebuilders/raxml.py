import os
import time
import logging
import random
from pyparsing import ParseException
import shutil
from Bio import AlignIO, SeqIO

from .base_treebuilder import TreeBuilder, AlignmentInput, DataType
from .parsers import RaxmlParser

from ..abstract_cli import AbstractCLI
from ..options import StringOption, FlagOption, IntegerOption, FloatOption, MultiOption, OptionSet

from ...file_utils import TempFile,TempDir

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


class RaxmlCLI(AbstractCLI):
    @property
    def _default_exe(self):
        return ['raxmlHPC','raxmlHPC-PTHREADS']


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


class Raxml(TreeBuilder):

    def __init__(self, alignment, *args, **kwargs):
        self.options = get_default_options()
        super(Raxml, self).__init__(alignment=alignment, *args, **kwargs)
        if self.input is not None:
            if self.datatype == DataType.DNA:
                set_default_dna_options(self)
            else:
                set_default_protein_options(self)



    def __call__(self, *args, **kwargs):
        """
        Sets up temporary output files and calls raxml using _call() function.
        Writes temporary input file if we're working with SeqIO object
        Saves the stdout and stderr and returns
        """
        start = time.time()  # time the execution

        #Need to create temp directory to put raxml output here
        with TempDir() as tmpd:
            if self.input_type is AlignmentInput.OBJECT:  # different operation depending on what it is
                with TempFile() as filename:
                    SeqIO.write(self.input, filename, 'phylip-relaxed') # default interleaved
                    output, error = self._call(filename,tmpd, *args, **kwargs)
            elif self.input_type is AlignmentInput.FILENAME:
                filename = self.input
                output, error = self._call(filename, tmpd, *args, **kwargs)
            else:
                output, error = self._call(None,tmpd, *args, **kwargs)
            self.result = self._read_result(tmpd)  # store result
        self.stdout = output
        self.stderr = error

        end = time.time()
        self.elapsed_time = end - start
        return self.result
        # End call

    # Any other accessory methods
    def _call(self, filename, tmpd, *args, **kwargs):
        """
        Call underlying low level _Raxml wrapper.
        Options are passed via *args and **kwargs
        [This only covers the simplest automatic
         case]
        """
        #hard code tmp_output as the output name since we don't save it anyway
        self.cli('{} -n tmp_output -w {tmp_path} -s {seqfile}'.format(self.command(), tmp_path=tmpd, seqfile=filename),
                wait=True)
        return self.cli.get_stdout(), self.cli.get_stderr()

    def command(self):
        return str(self.options)

    def _read_result(self, tmpd):
        """
        Read back the result.
        """

        expected_outfiles = [os.path.join(tmpd, 'RAxML_info.tmp_output'), os.path.join(tmpd, 'RAxML_bestTree.tmp_output')]


        parser = RaxmlParser()

        try:
            if self.options['-f'].get_value() is not '':
                f_value = os.path.splitext(os.path.basename(self.options['-f'].get_value()))[0]

                result = parser.to_dict(*expected_outfiles, dash_f=f_value)
            else:
                result = parser.to_dict(*expected_outfiles, dash_f=None)

        except IOError as ioerr:
            logger.error('Error reading results')
            result = None
        except ParseException as parseerr:
            logger.error('Other parse error', parseerr)
            result = None

        return result

    def _init_cli(self, binary):
        return RaxmlCLI(executable=binary)


def get_default_options():
    return OptionSet([
        # Algorithm

        # Set the model for either DNA or AA alignment
        StringOption('-m', 'PROTGAMMAGTR', active=True),

        # Number of replicates
        IntegerOption('-p', 12345, active=True),

        # If set to true will assume sequential format
        FlagOption('-q', False, active=False),

        # Turn on bootstrapping - set seed
        IntegerOption('-b', 0, active=False),

        # Number of replicates
        IntegerOption('-#', 0, active=False),

        # Turn on rapid bootstrap - specify seed
        IntegerOption('-x', 0, active=False),

        # Sed number of bootstrap replicates
        IntegerOption('-N', 0, active=False),

        # Set number of threads
        IntegerOption('-T', 0, active=False),

        # Tree topology search operation option
        StringOption('-s', 'NNI', active=False),

        # Select algorithm
        StringOption('-f', '', active=False),

        # Specify starting tree
        StringOption('-t', '', active=False),

        # Specify filename of file containing multiple trees
        StringOption('-z', '', active=False),

    ])
