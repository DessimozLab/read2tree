import time
import os
from Bio import AlignIO, SeqIO
from six import StringIO
from ..abstract_cli import AbstractCLI
from .base_aligner import Aligner, AlignmentInput, DataType
from ...wrappers import WrapperError
from ..options import StringOption, FlagOption, IntegerOption, FloatOption, MultiOption, OptionSet
import tempfile
import logging
logger = logging.getLogger(__name__)


class MafftCLI(AbstractCLI):
    """
    Mafft low-level command line interface

    :Example:

    ::

        mafft_cli = _MafftCLI()
        process = mafft_cli(cmd='mafft args...')
        stdout = mafft_cli.get_stdout()
    """
    @property
    def _default_exe(self):
        return 'mafft'

    # def _set_help(self):
    #     self(help=True, wait=True)
    #     self._help = self.get_stdout()


def set_default_dna_options(aligner):
    """
    Dummy function as sensible default already provided by mafft --auto
    """
    aligner.options = get_default_options()
    aligner.options['--auto'].set_value(True)


def set_default_protein_options(aligner):
    """
    Dummy function as sensible default already provided by mafft --auto
    """
    aligner.options = get_default_options()
    aligner.options['--auto'].set_value(True)


class Mafft(Aligner):
    """
    Convenient wrapper for Mafft multiple sequence aligner

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


    .. note:: There exists an ipython notebook on how to work with wrappers,
         including dealing with non-default parameters.
    """

    def __init__(self, input_, *args, **kwargs):
        super(Mafft, self).__init__(input_, *args, **kwargs)
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
        
        logger.debug('Output of Mafft: stdout={}; stderr={}'.format(output, error))
        if len(output)==0 and len(error)>0:
            logger.warning('is MAFFT_BINARIES set correctly: {}'.format(os.getenv('MAFFT_BINARIES','')))
            raise WrapperError('Mafft did not compute any alignments. StdErr: {}'.format(error))
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
        return MafftCLI(executable=binary)


def get_default_options():
    return OptionSet([
        # Algorithm

        # Automatically selects an appropriate strategy from L-INS-i, FFT-NS-i
        # and FFT-NS-2, according to data size. Default: off (always FFT-NS-2)
        FlagOption('--auto', False, active=False),

        # Distance is calculated based on the number of shared 6mers. Default: on
        FlagOption('--6merpair', True, active=False),

        # All pairwise alignments are computed with the Needleman-Wunsch algorithm.
        # More accurate but slower than --6merpair. Suitable for a set of globally
        # alignable sequences. Applicable to up to ~200 sequences. A combination
        # with --maxiterate 1000 is recommended (G-INS-i). Default: off
        # (6mer distance is used)
        FlagOption('--globalpair', False, active=False),

        # All pairwise alignments are computed with the Smith-Waterman algorithm.
        # More accurate but slower than --6merpair. Suitable for a set of locally
        # alignable sequences. Applicable to up to ~200 sequences. A combination
        # with --maxiterate 1000 is recommended (L-INS-i). Default: off
        # (6mer distance is used)
        FlagOption('--localpair', False, active=False),

        # All pairwise alignments are computed with a local algorithm with the
        # generalized affine gap cost (Altschul 1998). More accurate but slower than
        # --6merpair. Suitable when large internal gaps are expected. Applicable to
        # up to ~200 sequences. A combination with --maxiterate 1000 is recommended
        # (E-INS-i). Default: off (6mer distance is used)
        FlagOption('--genafpair', False, active=False),

        # All pairwise alignments are computed with FASTA (Pearson and Lipman 1988).
        # FASTA is required. Default: off (6mer distance is used)
        FlagOption('--fastapair', False, active=False),

        # Weighting factor for the consistency term calculated from pairwise
        # alignments. Valid when either of --blobalpair, --localpair, --genafpair,
        # --fastapair or --blastpair is selected. Default: 2.7
        FloatOption('--weighti', 2.7, active=False),

        # Guide tree is built number times in the progressive stage. Valid with 6mer
        # distance. Default: 2
        IntegerOption('--retree', 2, active=False),

        # number cycles of iterative refinement are performed. Default: 0
        IntegerOption('--maxiterate', 0, active=False),

        # Use FFT approximation in group-to-group alignment. Default: on
        FlagOption('--fft', True, active=False),

        # Do not use FFT approximation in group-to-group alignment. Default: off
        FlagOption('--nofft', False, active=False),

        #Alignment score is not checked in the iterative refinement stage. Default:
        # off (score is checked)
        FlagOption('--noscore', False, active=False),

        # Use the Myers-Miller (1988) algorithm. Default: automatically turned on
        # when the alignment length exceeds 10,000 (aa/nt).
        FlagOption('--memsave', False, active=False),

        # Use a fast tree-building method (PartTree, Katoh and Toh 2007) with the
        # 6mer distance. Recommended for a large number (> ~10,000) of sequences are
        # input. Default: off
        FlagOption('--parttree', False, active=False),

        # The PartTree algorithm is used with distances based on DP. Slightly more
        # accurate and slower than --parttree. Recommended for a large number
        # (> ~10,000) of sequences are input. Default: off
        FlagOption('--dpparttree', False, active=False),

        # The PartTree algorithm is used with distances based on FASTA. Slightly
        # more accurate and slower than --parttree. Recommended for a large number
        # (> ~10,000) of sequences are input. FASTA is required. Default: off
        FlagOption('--fastaparttree', False, active=False),

        # The number of partitions in the PartTree algorithm. Default: 50
        IntegerOption('--partsize', 50, active=False),

        # Do not make alignment larger than number sequences. Valid only with the
        # --*parttree options. Default: the number of input sequences
        IntegerOption('--groupsize', 1, active=False),

        # Parameter

        # Gap opening penalty at group-to-group alignment. Default: 1.53
        FloatOption('--op', 1.53, active=False),

        # Offset value, which works like gap extension penalty, for group-to-group
        # alignment. Deafult: 0.123
        FloatOption('--ep', 0.123, active=False),

        # Gap opening penalty at local pairwise alignment. Valid when the
        # --localpair or --genafpair option is selected. Default: -2.00
        FloatOption('--lop', -2.0, active=False),

        # Offset value at local pairwise alignment. Valid when the --localpair or
        # --genafpair option is selected. Default: 0.1
        FloatOption('--lep', 0.1, active=False),

        # Gap extension penalty at local pairwise alignment. Valid when the
        # --localpair or --genafpair option is selected. Default: -0.1
        FloatOption('--lexp', -0.1, active=False),

        # Gap opening penalty to skip the alignment. Valid when the --genafpair
        # option is selected. Default: -6.00
        FloatOption('--LOP', -6.00, active=False),

        # Gap extension penalty to skip the alignment. Valid when the --genafpair
        # option is selected. Default: 0.00
        FloatOption('--LEXP', 0.0, active=False),

        # BLOSUM number matrix (Henikoff and Henikoff 1992) is used. number=30, 45,
        # 62 or 80. Default: 62
        IntegerOption('--bl', 62, active=False),

        # JTT PAM number (Jones et al. 1992) matrix is used. number>0.
        # Default: BLOSUM62
        IntegerOption('--jtt', 100, active=False),

        # Transmembrane PAM number (Jones et al. 1994) matrix is used. number>0.
        # Default: BLOSUM62
        IntegerOption('--tm', 100, active=False),

        # Use a user-defined AA scoring matrix. The format of matrixfile is the same
        # to that of BLAST. Ignored when nucleotide sequences are input.
        # Default: BLOSUM62
        StringOption('--aamatrix', '', active=False),

        # Incorporate the AA/nuc composition information into the scoring matrix.
        # Deafult: off
        FlagOption('--fmodel', False, active=False),

        # Output

        # Output format: clustal format. Default: off (fasta format)
        FlagOption('--clustalout', False, active=False),

        # Output order: same as input. Default: on
        FlagOption('--inputorder', True, active=False),

        # Output order: aligned. Default: off (inputorder)
        FlagOption('--reorder', False, active=False),

        # Guide tree is output to the input.tree file. Default: off
        FlagOption('--treeout', False, active=False),

        # Do not report progress. Default: off
        FlagOption('--quiet', False, active=False),

        # Input

        # Assume the sequences are nucleotide. Default: auto
        FlagOption('--nuc', False, active=False),

        # Assume the sequences are amino acid. Default: auto
        FlagOption('--amino', False, active=False),

        # Seed alignments given in alignment_n (fasta format) are aligned with
        # sequences in input. The alignment within every seed is preserved.
        MultiOption('--seed', None, active=False),

        # Choose analysis mode automatically (default=True)
        FlagOption('--auto', True, active=True)
    ])
