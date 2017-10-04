import time
from Bio import AlignIO, SeqIO
from six import StringIO
from ..abstract_cli import AbstractCLI
from .base_aligner import Aligner, AlignmentInput, DataType
from ..options import StringOption, FlagOption, IntegerOption, FloatOption, MultiOption, OptionSet
import tempfile


class ProbConsCLI(AbstractCLI):
    """
    ProbCons low-level command line interface

    :Example:

    ::

        probcons_cli = _ProbConsCLI()
        process = mafft_cli(cmd='mafft args...')
        stdout = mafft_cli.get_stdout()
    """
    @property
    def _default_exe(self):
        return 'probcons'

    # def _set_help(self):
    #     self(help=True, wait=True)
    #     self._help = self.get_stdout()


def set_default_dna_options(aligner):
    """
    Dummy function as sensible default already provided by mafft --auto
    """
    aligner.options = get_default_options()


def set_default_protein_options(aligner):
    """
    Dummy function as sensible default already provided by mafft --auto
    """
    aligner.options = get_default_options()


class ProbCons(Aligner):
    """
    Convenient wrapper for ProbCons multiple sequence aligner

    The wrapper is written as a callable class.
    This can hold data (state) to do with the operation it performs, so it can keep results,
    execution times and other metadata, as well as perform the task.

    This is a basic implementation that can be extended. The important parts are
    __init__ (does the setup) and __call__ (does the work). All
    else are helper methods.

    :Example:

    ::

        callable_wrapper = ProbCons(aln)
        result = callable_wrapper()
        time_taken = callable_wrapper.elapsed_time
        result_again = callable_wrapper.result


    .. note:: There exists an ipython notebook on how to work with wrappers,
         including dealing with non-default parameters.
    """

    def __init__(self, input_, *args, **kwargs):
        super(ProbCons, self).__init__(input_, *args, **kwargs)
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
            with tempfile.NamedTemporaryFile(mode='wt') as filehandle:
                SeqIO.write(self.input, filehandle, 'fasta')
                filehandle.seek(0)
                output, error = self._call(filehandle.name, *args, **kwargs)
                
        else:
            output, error = self._call(self.input, *args, **kwargs)

        self.result = self._read_result(output) # store result
        self.stdout = output
        self.stderr = error

        end = time.time()
        self.elapsed_time = end - start
        return self.result
        # End call

    # Any other accessory methods 
    def _call(self, filename, *args, **kwargs):
        """
        Call underlying low level _Mafft wrapper. 
        Options are passed via *args and **kwargs
        [This only covers the simplest automatic
         case]
        """
        self.cli('{} {}'.format(self.command(), filename),
                 wait=True)
        return self.cli.get_stdout(), self.cli.get_stderr()

    def command(self):
        return str(self.options)

    def _read_result(self, output):
        """
        Read back the result.
        """
        fileobj = StringIO(output)
        return AlignIO.read(fileobj, 'fasta')

    def _init_cli(self, binary):
        return ProbConsCLI(executable=binary)


def get_default_options():
    return OptionSet([
        # Algorithm

        # use CLUSTALW output format instead of MFA
        FlagOption('-clustalw', False, active=False),

        # use 0 <= REPS <= 5 (default: 2) passes of consistency transformation
        IntegerOption('-c', 0, active=False),

        # use 0 <= REPS <= 1000 (default: 100) passes of iterative-refinement
        IntegerOption('-ir', 100, active=False),

        # use 0 <= REPS <= 20 (default: 0) rounds of pretraining
        IntegerOption('-pre', 0, active=False),

        # generate all-pairs pairwise alignments
        FlagOption('-pairs', False, active=False),

        #use Viterbi algorithm to generate all pairs(automatically enables - pairs)
        FlagOption('-viterbi', False, active=False),

        # write annotation for multiple alignment to FILENAME
        StringOption('-annot', '', active=False),

        # print sequences in alignment order rather than input order (default: off)
        FlagOption('-a', False, active=False)

    ])
