import os
import time
import logging
from ..abstract_cli import AbstractCLI
from ..options import StringOption, FloatOption, OptionSet
from read2tree.wrappers import WrapperError

logger = logging.getLogger(__name__)


class ClipkitCLI(AbstractCLI):
    @property
    def _default_exe(self):
        return 'clipkit'


class Clipkit(object):
    """
    Wrapper for ClipKIT alignment trimmer.

    Takes an input alignment file and writes a trimmed alignment to output_file.
    Returns the output file path on success, None if the output is empty.

    :Example:

    ::

        clipkit_wrapper = Clipkit('alignment.fasta', 'alignment_trimmed.fasta')
        result = clipkit_wrapper()
        time_taken = clipkit_wrapper.elapsed_time
    """

    def __init__(self, input_file, output_file, binary=None):
        self.input_file = input_file
        self.output_file = output_file
        self.options = get_default_options()
        self.elapsed_time = None
        self.stdout = None
        self.stderr = None
        self.result = None
        try:
            self.cli = ClipkitCLI(executable=binary)
        except IOError as err:
            raise WrapperError('Error searching for clipkit binary: {}'.format(err))

    def __call__(self, *args, **kwargs):
        start = time.time()
        output, error = self._call(self.input_file, self.output_file)
        self.stdout = output
        self.stderr = error
        self.result = self.output_file if os.path.exists(self.output_file) and os.path.getsize(self.output_file) > 0 else None
        end = time.time()
        self.elapsed_time = end - start
        return self.result

    def _call(self, input_file, output_file):
        self.cli('{} {} -o {}'.format(input_file, self.command(), output_file), wait=True)
        return self.cli.get_stdout(), self.cli.get_stderr()

    def command(self):
        return str(self.options)

    def _init_cli(self, binary):
        return ClipkitCLI(executable=binary)


def get_default_options():
    return OptionSet([
        StringOption('-m', 'gappy', active=True),
        FloatOption('-g', 0.8, active=True),
    ])
