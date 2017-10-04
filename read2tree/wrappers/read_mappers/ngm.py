import time
import tempfile
from Bio import SeqIO

from .parser import NGMParser
from ..abstract_cli import AbstractCLI
from .base_mapper import ReadMapper, ReferenceInput
import logging
from ..options import StringOption, FlagOption, IntegerOption, FloatOption, MultiOption, OptionSet
from pyparsing import Suppress, SkipTo, Word, Regex, Literal, OneOrMore, \
    Group, LineEnd, CharsNotIn, nums, alphanums, ParseException

logger = logging.getLogger(__name__)


class NGMCLI(AbstractCLI):
    @property
    def _default_exe(self):
        return 'ngm'


def set_default_options(readmapper):
    """
    Dummy function as sensible default
    """
    readmapper.options = get_default_options()


class NGM(ReadMapper):

    def __init__(self, reference, reads, *args, **kwargs):
        """
        :param alignment: input multiple sequence alignment. This can be either
            a filename or an biopython SeqRecord collection.
        """
        super(NGM, self).__init__(reference, reads, *args, **kwargs)
        set_default_options(self)

    def __call__(self, *args, **kwargs):
        """
        Anything to do with calling Mafft should go here.
        If any extra arguments need to be passed they can
        be specified (listed as *args and **kwargs for now).
        """
        start = time.time()  # time the execution
        if self.ref_input_type == ReferenceInput.OBJECT:  # different operation depending on what it is
            with tempfile.NamedTemporaryFile(mode='wt') as filehandle:
                SeqIO.write(self.ref_input, filehandle, 'fasta')
                filehandle.seek(0)
                output, error = self._call(filehandle.name, self.ref_input, *args, **kwargs)
                self.result = self._read_result(error, filehandle.name)
        else:
            output, error = self._call(self.ref_input, self.read_input, *args, **kwargs)
            self.result = self._read_result(error, self.ref_input)  # store result

        self.stdout = output
        self.stderr = error
        #
        end = time.time()
        self.elapsed_time = end - start
        return self.result
        # End call

    # Any other accessory methods
    def _call(self, reference, reads, *args, **kwargs):
        """
        Call underlying low level _ngm wrapper. 
        Options are passed via *args and **kwargs
        [This only covers the simplest automatic
         case]
        """
        self.cli('{} -r {} -1 {} -2 {} -o {}.sam'.format(self.command(), reference, reads[0], reads[1], reference), wait=True)

        return self.cli.get_stdout(), self.cli.get_stderr()

    def command(self):
        return str(self.options)

    def _read_result(self, output, filename):
        """
        Read back the result.
        """

        # TODO: change the output dictionary into a better format
        outfile = '{}.sam'.format(filename)
        parser = NGMParser()

        # Phyml outputs two outfiles, a stats file and a tree file.
        # Sometimes it appends .txt, sometimes not. Seems to be platform-specific.
        # Here we assume they are without .txt, but if we can't find them, try
        # looking for the .txt onees instead
        try:
            # parser.parse(output)
            result = parser.to_dict(outfile, output)
        except IOError as ioerr:
            logger.error('Error reading results')
            result = None
        except ParseException as parseerr:
            logger.error('Other parse error', parseerr)
            result = None

        return result

    def _init_cli(self, binary):
        return NGMCLI(executable=binary)

def get_default_options():
    return OptionSet([
        # Algorithm

        # set number of threads
        IntegerOption('-t', 4, active=True),

        # Automatically selects an appropriate strategy from L-INS-i, FFT-NS-i
        # and FFT-NS-2, according to data size. Default: off (always FFT-NS-2)
        FlagOption('--no-unal', False, active=True)


    ])