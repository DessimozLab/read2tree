import os
import time
import tempfile
import logging

from pyparsing import ParseException
from Bio import AlignIO, SeqIO

from .parsers import ProtTestParser
from .base_modeltester import ModelTester, AlignmentInput, DataType

from ..abstract_cli import AbstractCLI
from ..options import StringOption, FlagOption, OptionSet

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


class ProtTestCLI(AbstractCLI):
    """
    Especially in this case it is important that the $PROTTEST_HOME environmental variable is set to the installation directory of the prottest tool
    """
    @property
    def _default_exe(self):
        return 'java -jar ' + os.environ['PROTTEST_HOME'] + '/prottest-3.4.2.jar'


def set_default_dna_options(modeltester):
    """
    Dummy function as sensible default
    """
    modeltester.options = get_default_options()


def set_default_protein_options(modeltester):
    """
    Dummy function as sensible default
    """
    modeltester.options = get_default_options()


class ProtTest(ModelTester):
    """ ProtTest to determine the best model for a specific alignment
    This wrapper can be called to test various models for phylogeny inference.
    """

    def __init__(self, alignment, *args, **kwargs):
        """
        :param alignment: input multiple sequence alignment. This can be either
            a filename or an biopython SeqRecord collection.
        """
        self.options = get_default_options()
        super(ProtTest, self).__init__(alignment=alignment, *args, **kwargs)
        if self.datatype == DataType.DNA:
            set_default_dna_options(self)
        else:
            set_default_protein_options(self)

    def __call__(self, *args, **kwargs):
        """
        Anything to do with calling ProtTest should go here.
        If any extra arguments need to be passed they can
        be specified (listed as *args and **kwargs for now).
        """
        start = time.time()  # time the execution
        if self.input_type == AlignmentInput.OBJECT:  # different operation depending on what it is
            with tempfile.NamedTemporaryFile(mode='wt') as filehandle:
                SeqIO.write(self.input, filehandle, 'fasta')
                filehandle.seek(0)
                output, error = self._call(filehandle.name, *args, **kwargs)
        else:
            output, error = self._call(self.input, *args, **kwargs)

        self.result = self._read_result(output)  # store result
        self.stdout = output
        self.stderr = error

        end = time.time()
        self.elapsed_time = end - start
        return self.result
        # End call

    # Any other accessory methods
    def _call(self, filename, *args, **kwargs):
        """
        Call underlying low level _ProtTest wrapper.
        Options are passed via *args and **kwargs
        [This only covers the simplest automatic
         case]
        """
        self.cli('{} -i {}'.format(self.command(), filename),
                 wait=True)
        return self.cli.get_stdout(), self.cli.get_stderr()

    def command(self):
        return str(self.options)

    def _read_result(self, output):

        parser = ProtTestParser()

        try:
            result = parser.to_dict(output)

        except IOError as ioerr:
            logger.error('Error reading results')
            result = None
        except ParseException as parseerr:
            logger.error('Other parse error', parseerr)
            result = None

        return result


    def _init_cli(self, binary):
        return ProtTestCLI(executable=binary)


def get_default_options():
    return OptionSet([
        # Algorithm

        # Display models sorted by Akaike Information Criterion (AIC)
        FlagOption('-AIC', False, active=False),

        # Display models sorted by Decision Theory Criterion
        FlagOption('-DT', False, active=False),

        # Tree file (optional) [default: NJ tree]
        StringOption('-t', '', active=False),

        # Display models sorted by Corrected Akaike Information Criterion (AICc)
        FlagOption('-AICC', False, active=False),

        #Enables / Disables PhyML logging into log directory(see prottest.properties)
        FlagOption('-log', False, active=False)
    ])
