import tempfile
import time
from Bio import AlignIO, SeqIO
from six import StringIO
from ..abstract_cli import AbstractCLI
from .base_aligner import Aligner, AlignmentInput, DataType
from ..options import StringOption, FlagOption, IntegerOption, FloatOption, MultiOption, TreeInputOption, OptionSet


class MuscleCLI(AbstractCLI):
    """
    Muscle low-level command line interface

    example:
    muscle_cli = MuscleCLI()
    process = muscle_cli(cmd='muscle args...')
    stdout = muscle_cli.get_stdout()
    """
    @property
    def _default_exe(self):
        return 'muscle'

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

class Muscle(Aligner):
    """
    Convenient wrapper for Muscle multiple sequence aligner

    The wrapper is written as a callable class.
    This can hold data (state) to do with the operation it performs, so it can keep results,
    execution times and other metadata, as well as perform the task.

    This is a basic implementation that can be extended. The important parts are
    __init__ (does the setup) and __call__ (does the work). All
    else are helper methods.

    :Example:

    ::

        callable_wrapper = Muscle(aln)
        result = callable_wrapper()
        time_taken = callable_wrapper.elapsed_time
        result_again = callable_wrapper.result
    """

    def __init__(self, input_, *args, **kwargs):
        super(Muscle, self).__init__(input_, *args, **kwargs)
        self.options = get_default_options()

        if self.datatype == DataType.DNA:
            set_default_dna_options(self)
        else:
            set_default_protein_options(self)

    def __call__(self, *args, **kwargs):
        """
        Anything to do with calling Muscle should go here.
        If any extra arguments need to be passed they can
        be specified (listed as *args and **kwargs for now).
        """
        start = time.time() # time the execution

        if self.input_type == AlignmentInput.OBJECT: # different operation depending on what it is
            with tempfile.NamedTemporaryFile(mode="wt") as filehandle:
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
        Call underlying low level _MuscleCLI wrapper. 
        Options are passed via *args and **kwargs
        [This only covers the simplest automatic
         case]
        """
        self.cli('{} -in {}'.format(self.command(), filename),
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
        return MuscleCLI(executable=binary)


def get_default_options():
    return OptionSet([
        # Algorithm

        # Find diagonals (faster for similar sequences)
        FlagOption('-diags', False, active=False),

        # Maximum number of iterations(integer, default 16)
        IntegerOption('-maxiters', 16, active=False),

        # Maximum time to iterate in hours (default no limit)
        FloatOption('-maxhours', 0.0, active=False)

        #reeInputOption('-usetree', '', active=False)
    ])
