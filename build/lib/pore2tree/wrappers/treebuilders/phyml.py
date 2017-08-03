import os
import time
import tempfile
import logging
from pyparsing import ParseException
from Bio import AlignIO, SeqIO

from .base_treebuilder import TreeBuilder, AlignmentInput, DataType
from .parsers import PhymlParser

from ..abstract_cli import AbstractCLI
from ..options import StringOption, FlagOption, IntegerOption, FloatOption, MultiOption, OptionSet


logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


class PhymlCLI(AbstractCLI):
    @property
    def _default_exe(self):
        return 'phyml'


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


class Phyml(TreeBuilder):
    """ Phyml tree reconstruction

    This wrapper can be called to reconstruct a phylogenetic tree
    using PhyML.
    """

    def __init__(self, alignment, *args, **kwargs):
        """
        :param alignment: input multiple sequence alignment. This can be either
            a filename or an biopython SeqRecord collection.
        """
        super(Phyml, self).__init__(alignment, *args, **kwargs)
        self.options = get_default_options()
        if self.datatype == DataType.DNA:
            set_default_dna_options(self)
        else:
            set_default_protein_options(self)

    def __call__(self, *args, **kwargs):
        """
        Anything to do with calling Mafft should go here.
        If any extra arguments need to be passed they can
        be specified (listed as *args and **kwargs for now).
        """
        start = time.time()  # time the execution

        if self.input_type == AlignmentInput.OBJECT:  # different operation depending on what it is
            with tempfile.NamedTemporaryFile(mode='wt') as fh:
                SeqIO.write(self.input, fh, 'phylip-relaxed') # default interleaved
                fh.seek(0)
                output, error = self._call(fh.name, *args, **kwargs)
                self.result = self._read_result(fh.name)  # store result
        else:
            path = os.path.dirname(self.input)
            filename = os.path.basename(self.input)
            os.chdir(path)  # some operations done because phyml can not deal with large filenames that are caused due to a large path
            output, error = self._call(filename, *args, **kwargs)
            self.result = self._read_result(filename)  # store result

        self.stdout = output
        self.stderr = error

        end = time.time()
        self.elapsed_time = end - start
        return self.result["tree"]
        # End call

    # Any other accessory methods
    def _call(self, filename, *args, **kwargs):
        """
        Call underlying low level _Phyml wrapper.
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
        """
        Read back the result.
        """

        #TODO: change the output dictionary into a better format
        expected_outfiles = ['{}_phyml_stats'.format(output), '{}_phyml_tree'.format(output)]
        parser = PhymlParser()

        # Phyml outputs two outfiles, a stats file and a tree file.
        # Sometimes it appends .txt, sometimes not. Seems to be platform-specific.
        # Here we assume they are without .txt, but if we can't find them, try
        # looking for the .txt onees instead
        try:
            # Check if these are the .txt style outfiles
            if not os.path.exists(expected_outfiles[0]):
                expected_outfiles = [x + '.txt' for x in expected_outfiles]
            result = parser.to_dict(*expected_outfiles)

        except IOError as ioerr:
            logger.error('Error reading results')
            result = None
        except ParseException as parseerr:
            logger.error('Other parse error', parseerr)
            result = None

        return result

    def _init_cli(self, binary):
        return PhymlCLI(executable=binary)


def get_default_options():
    return OptionSet([
        # Algorithm

        # Set datatype to nt or aa
        StringOption('-d', 'aa', active=True),

        # Set the model for either DNA or AA alignment
        StringOption('-m', '', active=False),

        # If set to true will assume sequential format
        FlagOption('-q', False, active=False),

        # Set bootstrap value
        IntegerOption('-b', 0, active=False),

        # Tree topology search operation option
        StringOption('-s', 'NNI', active=False)
    ])
