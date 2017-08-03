import time
from Bio import AlignIO, SeqIO
import tempfile
from six import StringIO
from ..abstract_cli import AbstractCLI
from .base_aligner import Aligner, AlignmentInput, DataType
from ..options import StringOption, FlagOption, IntegerOption, FloatOption, MultiOption, OptionSet


class ProGraphMSACLI(AbstractCLI):
    """
    PrographMSA low-level command line interface

    :Example:

    ::

        prograph_cli = ProGraphMSACLI()
        process = prograph_cli(cmd='mafft args...')
        stdout = prograph_cli.get_stdout()
    """

    @property
    def _default_exe(self):
        return 'ProGraphMSA'


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


class ProGraphMSA(Aligner):
    """
    Convenient wrapper for ProGraphMSA multiple sequence aligner

    The wrapper is written as a callable class.
    This can hold data (state) to do with the operation it performs, so it can keep results,
    execution times and other metadata, as well as perform the task.

    This is a basic implementation that can be extended. The important parts are
    __init__ (does the setup) and __call__ (does the work). All
    else are helper methods.

    :Example:

    ::

        callable_wrapper = Mafft(aln)
        result = callable_wrapper()
        time_taken = callable_wrapper.elapsed_time
        result_again = callable_wrapper.result
    """

    def __init__(self, input_, *args, **kwargs):
        super(ProGraphMSA, self).__init__(input_, *args, **kwargs)
        self.options = get_default_options()
        if self.datatype == DataType.DNA:
            set_default_dna_options(self)
        else:
            set_default_protein_options(self)

    def __call__(self, *args, **kwargs):
        """
        Anything to do with calling ProGraphMSA should go here.
        If any extra arguments need to be passed they can
        be specified (listed as *args and **kwargs for now).
        """
        start = time.time()  # time the execution

        if self.input_type == AlignmentInput.OBJECT:  # different operation depending on what it is
            with tempfile.NamedTemporaryFile(mode="wt") as fh:
                SeqIO.write(self.input, fh, 'fasta')
                fh.seek(0)
                output, error = self._call(fh.name, *args, **kwargs)

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
        Call underlying low level ProGraphMSA wrapper.
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
        return ProGraphMSACLI(executable=binary)


def get_default_options():
    return OptionSet([
        # Algorithm

        # output fasta format (instead of stockholm), better because no tree output is produced
        FlagOption('--fasta', True, active=True),

        # output all ancestral sequences
        FlagOption('--ancestral_seqs', False, active=False),

        # output sequences in input order (default: tree order)
        FlagOption('--input_order', False, active=False),

        # output all intermediate guide trees
        FlagOption('--all_trees', False, active=False),

        # use ML distances with gap
        FlagOption('--mldist_gap', False, active=False),

        # use ML distances
        FlagOption('--mldist', False, active=False),

        # use of guide tree
        StringOption('--tree', '', active=False)

    ])
