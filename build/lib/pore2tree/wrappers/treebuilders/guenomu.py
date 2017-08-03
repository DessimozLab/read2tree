import os
import time
import random
from pyparsing import ParseException
import shutil
from Bio import AlignIO, SeqIO
import dendropy
from zoo.wrappers import WrapperError
import logging

from .parsers import GuenomuParser
from .base_treebuilder import TreeBuilder, AlignmentInput, DataType

from ..abstract_cli import AbstractCLI
from ..options import StringOption, FlagOption, IntegerOption, FloatOption, MultiOption, OptionSet

from ...file_utils import TempFile, TempDir

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


class GuenomuCLI(AbstractCLI):
    @property
    def _default_exe(self):
        return "guenomu"


class Guenomu(TreeBuilder):
    """
    Species tree reconstruction 
    
    This is a wrapper for the guenomu software, that estimates species trees from (distributions of) gene family trees.
    It does rely on previously built trees for the gene families (that don't need to be orthologs-only, by the way).
    """

    def __init__(self, genetrees, species, *args, **kwargs):
        """
        :param genetrees: a list of dendropy.TreeList objects (which may have one or more trees). Must be a list, even
        for a single dendropy.TreeList
        :param species: a list of species names (as strings), that must be found within the gene leaves
        """
        self.options = get_default_options()
        self.genetrees = genetrees  # TODO: check if is a list of TreeList
        self.species = species  # TODO: check if it's a list of strings
        self.elapsed_time = None
        self.stdout = None
        self.stderr = None
        try:
            self.cli = self._init_cli()
        except IOError as err:
            raise WrapperError('Error searching for binary: {}'.format(err))
            # End setup

    def __call__(self, *args, **kwargs):
        """
        Sets up temporary output files and calls guenomu using _call() function.
        Writes temporary input files since we are working with dendropy.TreeList
        Saves the stdout and stderr and returns
        """
        start = time.time()  # time the execution

        with TempDir() as tmpd:  # creates temporary directory and removes it when out of scope (i.e. leave this block)
            os.chdir(tmpd)
            # create one nexus file per gene family tree
            glistnames = []
            for i in range(len(self.genetrees)):
                print(i)  ## DEBUG
                glistnames.append("gene" + str(i))
                self.genetrees[i].write(path="gene" + str(i), schema="nexus")
            # create file with name of all gene family tree files
            with open("gene_tree_list", 'w') as f:
                for s in glistnames:
                    f.write(s + '\n')
            # create file with names of all species
            with open("species_names", 'w') as f:
                for s in self.species:
                    f.write(s + '\n')

            self.stdout, self.stderr = self._call(tmpd, *args, **kwargs)
            self.result = self._read_result(tmpd)  # store result

        end = time.time()
        self.elapsed_time = end - start
        return self.result["best_topology"]
        # End call

    # Any other accessory methods
    def _call(self, tempdir, *args, **kwargs):
        """
        Call underlying low level wrapper.

        Options are passed via *args and **kwargs [This only
        covers the simplest automatic case]
        """
        self.cli('{} -S species_names -G gene_tree_list -z 0'.format(self.command()), wait=True)  # importance sampling
        self.cli('{} -S species_names -G gene_tree_list -z 1'.format(self.command()), wait=True)  # analyse output
        return self.cli.get_stdout(), self.cli.get_stderr()

    def command(self):
        return str(self.options)

    def _read_result(self, tmpd):
        """
        Read back the result species tree files.
        """
        ## job0.<genetreefilename>.trprobs also available, as well as job0.params.txt
        expected_outfiles = [os.path.join(tmpd, 'species.tre'), os.path.join(tmpd, 'job0.unrooted.trprobs')]
        parser = GuenomuParser()
        result = parser.to_dict(expected_outfiles)
        return result

    def _init_cli(self):
        return GuenomuCLI()


def get_default_options():
    return OptionSet([
        # -L, --lambda=<double>        baseline reconciliation hyperprior that control distances (one or several values, should be << 1)
        FloatOption('-L', 0.0001, active=True),
        # -A, --sa_samples=<int>       number of cycles/samples simulated Annealing
        IntegerOption('-A', 10, active=True),
        # -a, --sa_iter=<int>          iterations per cycle for simulated Annealing
        IntegerOption('-a', 1000, active=True),
        # -t, --temp_i=<double>        initial (inverse) temperature for each cycle of simulated annealing
        FloatOption('-t', 0.1, active=True),
        # -T, --temp_f=<double>        final (inverse) temperature at each cycle of simulated annealing
        FloatOption('-T', 5., active=True),
        # -b, --burnin=<int>           iterations for burnin stage of Bayesian posterior sampling
        IntegerOption('-b', 0, active=True),  # default is no burnin since we run only Simul Annealing
        # -B, --b_iterations=<int>     iterations for main Bayesian posterior sampling stage (zero if you only want simulated annealing)
        IntegerOption('-B', 0, active=True),  # default is no burnin since we run only Simul Annealing
        # -n, --samples=<int>          Number of samples to be stored to file and analysed later
        IntegerOption('-n', 1, active=True),  # under SimulAnneal-only, this will default to number of cycles (-A)
        # -N, --printout=<int>         how many times summary info should be output on screen
        IntegerOption('-N', 1, active=True),
        # -m, --mtm=<int>              sample size (number of trials) for each Multiple-Try Metropolis proposal
        IntegerOption('-m', 4, active=True),
        # -p, --partition=<int>        auxiliary MCMC length to avoid the Partition function calculation (zero to neglect normalization)
        IntegerOption('-p', 0, active=True),
        # -D, --distances=<string>     string with seven  1's and 0's describing distances to be used [dup, los, ils, rf, hdist, hdist_final, spr]
        StringOption('-D', '1111111', active=True),
    ])


"""
guenomu  [-h|--help] [-L|--lambda=<double>]... [-A|--sa_samples=<int>] [-a|--sa_iter=<int>] [-t|--temp_i=<double>] [-T|--temp_f=<double>] [-b|--burnin=<int>] [-B|--b_iterations=<int>] [-n|--samples=<int>] [-N|--printout=<int>] [-m|--mtm=<int>] [-p|--partition=<int>] [-z|--action=<int>] [-D|--distances=<string>] [-S|--species=<file>] [-G|--genetrees=<file>] [<file>]
  -z, --action=<int>           execution mode for program: 0='importance sampling', 1='analisis of output' etc.
  -S, --species=<file>         file with species names, one name per line (nexus-style bracketed comments are allowed)
  -G, --genetrees=<file>       file with list of gene tree file names, one gene tree file name per line
  <file>                       Name of control file with remaining parameters (those defined also here overwrite those from control file)

"""
