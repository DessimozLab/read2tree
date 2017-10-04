import os
import time
import logging
import random
from pyparsing import ParseException
import shutil

from .parsers import IqtreeParser
from .base_treebuilder import TreeBuilder, AlignmentInput, DataType

from ..abstract_cli import AbstractCLI
from ..options import StringOption, FlagOption, IntegerOption, FloatOption, MultiOption, OptionSet

from ...file_utils import TempFile,TempDir

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


class IqtreeCLI(AbstractCLI):
    @property
    def _default_exe(self):
        return 'iqtree-omp'


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


class Iqtree(TreeBuilder):

    def __init__(self, input_):
        super(Iqtree, self).__init__(input_)
        self.options = get_default_options()
        if self.datatype == DataType.DNA:
            set_default_dna_options(self)
        else:
            set_default_protein_options(self)

    def __call__(self, *args, **kwargs):
        """
        Sets up temporary output file location and calls iqtree using _call() function.
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
        Call underlying low level _iqtree wrapper.
        Options are passed via *args and **kwargs
        [This only covers the simplest automatic
         case]
        """
        self.cli('{} -pre {tmp_path} -s {seqfile}'.format(self.command(),
                                                          tmp_path=os.path.join(tmpd,'tmp_output'),
                                                          seqfile=filename),
                 wait=True)
        return self.cli.get_stdout(), self.cli.get_stderr()

    def command(self):
        return str(self.options)

    def _read_result(self, tmpd):
        """
        Read back the result.
        """

        expected_outfiles = [os.path.join(tmpd, 'tmp_output.treefile')]

        parser = IqtreeParser()

        try:
            result = parser.to_dict(*expected_outfiles)

        except IOError as ioerr:
            logger.error('Error reading results')
            result = None
        except ParseException as parseerr:
            logger.error('Other parse error', parseerr)
            result = None

        return result["tree"]

    def _init_cli(self, binary):
        return IqtreeCLI(executable=binary)


def get_default_options():
    return OptionSet([
        # Number of threads
        IntegerOption('-nt', 2, active=True),

        # Set the model for either DNA or AA alignment
        StringOption('-m', '', active=False),

        # If set to true will assume sequential format
        #FlagOption('-q', False, active=False),

        # Ultrafast bootstrap (>=1000)
        IntegerOption('-bb', 0, active=False),

        # SH-like approximate likelihood ratio test (SH-aLRT)
        IntegerOption('-alrt', 0, active=False),

        # Bootstrap + ML tree + consensus tree (>=100)
        IntegerOption('-b', 0, active=False)
    ])
