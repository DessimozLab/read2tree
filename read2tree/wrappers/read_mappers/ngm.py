import time
import os
import tempfile
from Bio import SeqIO

from .parser import NGMParser
from ..abstract_cli import AbstractCLI
from .base_mapper import ReadMapper, ReadInput
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

    def __init__(self, reference, reads, tmp_folder, *args, **kwargs):
        """
        :param alignment: input multiple sequence alignment. This can be either
            a filename or an biopython SeqRecord collection.
        """
        super(NGM, self).__init__(reference, reads, tmp_folder, *args, **kwargs)
        set_default_options(self)

    def __call__(self, *args, **kwargs):
        """
        Anything to do with calling Mafft should go here.
        If any extra arguments need to be passed they can
        be specified (listed as *args and **kwargs for now).
        """
        start = time.time()  # time the execution
        if self.read_input_type == ReadInput.OBJECT:  # different operation depending on what it is
            with tempfile.NamedTemporaryFile(mode='wt') as filehandle:
                SeqIO.write(self.ref_input, filehandle, 'fastq')
                filehandle.seek(0)
                output, error = self._call(self.read_input, filehandle.name,
                                           self.tmp_folder, *args, **kwargs)
                self.result = self._read_result(error, self.ref_input, self.tmp_folder)
        elif self.read_input_type == ReadInput.STRING:
            with tempfile.NamedTemporaryFile(mode='wt') as filehandle:
                filehandle.write(self.read_input)
                filehandle.seek(0)
                output, error = self._call(self.ref_input, filehandle.name,
                                           self.tmp_folder, *args, **kwargs)
                self.result = self._read_result(error, self.ref_input, self.tmp_folder)
        elif self.read_input_type == ReadInput.FILENAME:
            output, error = self._call(self.ref_input, self.read_input,
                                       self.tmp_folder, *args, **kwargs)
            self.result = self._read_result(error, self.ref_input, self.tmp_folder)  # store result

        self.stdout = output
        self.stderr = error
        #
        end = time.time()
        self.elapsed_time = end - start
        return self.result
        # End call

    # Any other accessory methods
    def _call(self, reference, reads, tmp_folder=None, *args, **kwargs):
        """
        Call underlying low level _ngm wrapper.
        Options are passed via *args and **kwargs
        [This only covers the simplest automatic
         case]
        """
        if tmp_folder is None:
            tmp_file = './' + os.path.basename(reference)+".bam"
        if '/' not in tmp_folder[-1]:
            tmp_file = os.path.join(tmp_folder,
                                    os.path.basename(reference))+".bam"
        if len(reads) is 2:
            self.cli('{} -b -r {} -1 {} -2 {} -o {}'.format(self.command(),
                                                            reference,
                                                            reads[0], reads[1],
                                                            tmp_file),
                     wait=True)
        elif len(reads) is not 2:
            self.cli('{} -b -r {} -q {} -o {}'.format(self.command(),
                                                      reference, reads,
                                                      tmp_file), wait=True)

        return self.cli.get_stdout(), self.cli.get_stderr()

    def command(self):
        return str(self.options)

    def _read_result(self, output, filename, tmp_folder):
        """
        Read back the result.
        """

        # TODO: change the output dictionary into a better format
        # outfile = '{}.sam'.format(filename)
        outfile = os.path.join(tmp_folder, os.path.basename(filename))+".bam"
        parser = NGMParser()

        try:
            # parser.parse(output)
            result = parser.to_dict(outfile, output)
        except IOError as ioerr:
            logger.error('Error reading results', ioerr)
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
        IntegerOption('-t', 1, active=True),

        # makes sure that unmapped reads are not saved in bam file
        FlagOption('--no-unal', True, active=True)
    ])
