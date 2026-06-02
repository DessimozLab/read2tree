import os
import time
import logging
from ..abstract_cli import AbstractCLI
from ..options import IntegerOption, OptionSet
from read2tree.wrappers import WrapperError

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


class AsterCLI(AbstractCLI):
    @property
    def _default_exe(self):
        return ['astral3', 'astral-pro3', 'astral-pro2', 'wastral', 'astral4']


class Aster(object):
    """
    Wrapper for the ASTER suite of coalescent species tree estimators.

    Takes a file of gene trees (one Newick tree per line) and writes a species
    tree to output_file using the selected ASTER algorithm.  Returns the
    species tree in Newick format.

    :Example:

    ::

        aster_wrapper = Aster('gene_trees.nwk', 'species_tree.nwk')
        aster_wrapper.options.options['-t'].set_value(8)
        result = aster_wrapper()
        time_taken = aster_wrapper.elapsed_time
    """

    def __init__(self, gene_tree_file, output_file, binary=None):
        self.gene_tree_file = gene_tree_file
        self.output_file = output_file
        self.options = get_default_options()
        self.elapsed_time = None
        self.stdout = None
        self.stderr = None
        self.result = None
        try:
            self.cli = AsterCLI(executable=binary)
        except IOError as err:
            raise WrapperError('Error searching for ASTER binary: {}'.format(err))

    def __call__(self, *args, **kwargs):
        start = time.time()
        output, error = self._call(self.gene_tree_file, self.output_file)
        self.stdout = output
        self.stderr = error
        self.result = self._read_result(self.output_file)
        end = time.time()
        self.elapsed_time = end - start
        return self.result

    def _call(self, gene_tree_file, output_file):
        self.cli('{} -i {} -o {}'.format(self.command(), gene_tree_file, output_file), wait=True)
        return self.cli.get_stdout(), self.cli.get_stderr()

    def command(self):
        return str(self.options)

    def _read_result(self, output_file):
        try:
            with open(output_file, 'r') as fh:
                return fh.read().strip()
        except IOError:
            logger.error('Error reading ASTER output: {}'.format(output_file))
            return None


def get_default_options():
    return OptionSet([
        IntegerOption('-t', 1, active=True),
    ])
