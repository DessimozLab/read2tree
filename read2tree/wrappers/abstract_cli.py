import os
import shlex
import sys
import threading
from abc import ABCMeta, abstractmethod, abstractproperty
from locale import getpreferredencoding
from subprocess import PIPE, Popen
import logging


IS_PY3 = sys.version_info[0] == 3
POSIX = 'posix' in sys.builtin_module_names
DEFAULT_ENCODING = getpreferredencoding() or "UTF-8"
logger = logging.getLogger(__name__)

if IS_PY3:
    from queue import Queue
else:
    from Queue import Queue


class ExternalProcessError(Exception):
    pass


def _py2_and_3_joiner(sep, joinable):
    """
    Allow '\n'.join(...) statements to work in Py2 and Py3.
    :param sep:
    :param joinable:
    :return:
    """
    if IS_PY3:
        sep = bytes(sep, DEFAULT_ENCODING)
    joined = sep.join(joinable)
    return joined.decode(DEFAULT_ENCODING) if IS_PY3 else joined


def _kwargs_to_args(kwargs, num_hyphens):
    processed = []
    for k, v in kwargs.items():
        # we're passing a short arg as a kwarg, example:
        # cut(d="\t")
        if v is not False:
            if len(k) == 1:
                processed.append('-' + k)
            else:
                processed.append('-' * num_hyphens + k)
            if v is not True:
                processed.append(v)
    return ' '.join(processed)


class AbstractCLI(object):
    """
    Abstract Base Class:
    Run an external program as a subprocess with non-blocking collection of stdout.
    This class uses subprocess.Popen to handle the system call, and threading.Thread
    to allow non-blocking reads of stderr and stdout.

    Specific program wrappers inherit from this. They should implement two methods:
    1: _default_exe and 2: _set_help.

    :Example:

    ::

       class SomeProgram(AbstractInternal):
           @property
           def _default_exe(self):
               return 'some_program'           # The usual name of the program binary
                                               # or alternatively specify a list of names

           def _set_help(self):
               self('--help', wait=True)       # calls 'some_program --help'
               self._help = self.get_stdout()  # puts the help output into self._help

           #### Possible Extra Requirement ####
           @property
           def _hyphen_policy(self):
               return 1       # If some_program arguments take a single leading hyphen
                              # then set this to 1 (i.e. to the number of leading hyphens)
                              # This only affects calling some_program using named arguments

    That's it!

    :seealso:

        See blog post: `<http://www.zultron.com/2012/06/python-subprocess-example-running-a-background-subprocess-with-non-blocking-output-processing/>`_
        Also, StackOverflow: `<https://stackoverflow.com/questions/375427/non-blocking-read-on-a-subprocess-pipe-in-python>`_

    Keyword handling borrowed from sh: `<https://amoffat.github.io/sh>`_

    """
    __metaclass__ = ABCMeta

    def __init__(self, executable=None, verbose=True):
        """
        Sets up the wrapper. A custom executable can be passed in, otherwise
        it will search the PATH.

        :param executable: Path to a custom executable to use
        :return:
        """
        exe = None
        if executable:
            exe = self._search_for_executable(executable)
        if exe is None:
            if isinstance(self._default_exe, list): 
                for default_name in self._default_exe:
                    exe = self._search_for_executable(default_name)
                    if exe is not None:
                        break
            else:
                exe = self._search_for_executable(self._default_exe)
            if not exe:
                logger.error("executable not found. PATH='{}'".format(os.getenv('PATH', '')))
                raise IOError(executable if executable else self._default_exe)

        self.exe = exe  # The wrapped executable
        self.verbose = False  # Controls printing of output
        self.stdout_q = Queue()  # The output-reading threads pipe output into
        self.stderr_q = Queue()  # these queues.
        self.stdout_l = list()  # Queued output gets placed in these lists,
        self.stderr_l = list()  # making it available to the caller
        self.stdout_thread = None  # monitoring thread that captures stdout
        self.stderr_thread = None  # monitoring thread that captures stderr
        self.process = None  # This holds the running process once the wrapped program is called
        self._help = None  # A place to hold a help string

    def __repr__(self):
        return '{}(executable=\'{}\')'.format(self.__class__.__name__, self.exe)

    # Private
    @abstractproperty
    def _default_exe(self):
        pass

    @property
    def _hyphen_policy(self):
        """
        Returns 'n', where 'n' is the number of hyphens prepended to commandline flags.
        Used internally when constructing command line arguments from parameters
        passed to __call__ as keywords.
        :return: 2 (as in 2 hyphens, '--')
        """
        return 2

    def _set_help(self):
        pass

    def _log_thread(self, pipe, queue):
        """
        Start a thread logging output from pipe
        """

        # thread function to log subprocess output (LOG is a queue)
        def enqueue_output(out, q):
            for line in iter(out.readline, b''):
                q.put(line.rstrip())
            out.close()

        # start thread
        t = threading.Thread(target=enqueue_output,
                                  args=(pipe, queue))
        t.daemon = True  # thread dies with the program
        t.start()
        return t

    def _search_for_executable(self, executable):
        """
        Search for file give in "executable". If it is not found, we try the environment PATH.
        Returns either the absolute path to the found executable, or None if the executable
        couldn't be found.
        """
        # Check in case java wrapper is called where executable is called with "java -j XXX"
        if len(executable.split(" ")) > 1:
            return executable

        if os.path.isfile(executable):
            return os.path.abspath(executable)
        else:
            envpath = os.getenv('PATH')
            if envpath is None:
                return
            for path in envpath.split(os.pathsep):
                exe = os.path.join(path, executable)
                if os.path.isfile(exe):
                    return os.path.abspath(exe)

    # Public
    def __call__(self, cmd=None, wait=False, **flags):
        """
        Spawns the subprocess and the threads used to monitor stdout and stderr without blocking.
        :param cmd: Pass the command line arguments as a string
        :param wait: Block until the process returns
        :param flags: Pass the commandline arguments as a dictionary. Will be appended to any content in cmd.

        :return:
        """
        # Check there is not already a process running
        if self.running():
            self.kill()

        # Wipe any stdout/stderr from previous processes
        self.stderr_l = []
        self.stdout_l = []

        # Assemble command line
        if cmd is None:
            cmd = ''
        if flags:
            cmd = ' '.join([cmd.strip(), _kwargs_to_args(flags, self._hyphen_policy)])
        self.cmd = '{} {}'.format(self.exe, cmd)
        logger.debug('Running following command: {}'.format(self.cmd))

        # spawn
        self.process = Popen(shlex.split(self.cmd),
                             shell=False, stdout=PIPE, stderr=PIPE, bufsize=1, close_fds=POSIX)
        if self.verbose:
            print('Launched {} with PID {}'.format(self.exe, self.process.pid))

        # start stdout and stderr logging threads
        self.stdout_thread = self._log_thread(self.process.stdout, self.stdout_q)
        self.stderr_thread = self._log_thread(self.process.stderr, self.stderr_q)

        if wait:
            self.process.wait()
            self._join_threads()

    @property
    def help(self):
        """
        Returns a helpful string, preferably derived from the wrapped program
        :return:
        """
        if self._help is None:
            self._set_help()
        return self._help

    def _join_threads(self):
        for thread in [self.stdout_thread, self.stderr_thread]:
            if thread is not None and thread.is_alive():
                thread.join()
                thread = None

    def get_stderr(self, tail=None):
        """
        Returns current total output written to standard error.
        :param tail: Return this number of most-recent lines.
        :return: copy of stderr stream
        """
        while not self.stderr_q.empty():
            self.stderr_l.append(self.stderr_q.get_nowait())
        if tail is None:
            tail = len(self.stderr_l)
        return _py2_and_3_joiner('\n', self.stderr_l[:tail])

    def get_stdout(self, tail=None):
        """
        Returns current total output written to standard output.
        :param tail: Return this number of most-recent lines.
        :return: copy of stdout stream
        """
        while not self.stdout_q.empty():
            self.stdout_l.append(self.stdout_q.get_nowait())
        if tail is None:
            tail = len(self.stdout_l)
        return _py2_and_3_joiner('\n', self.stdout_l[:tail])

    def finished(self):
        """
        Check if the running process is finished. Raises an exception if no process has ever been launched.
        :return: bool
        """
        if self.process is None:
            raise ExternalProcessError('No process has been launched from this instance')
        return self.process.poll() is not None

    def kill(self):
        """
        Kill the running process (if there is one)
        :return: void
        """
        if self.running():
            if self.verbose:
                print('Killing {} with PID {}'.format(self.exe, self.process.pid))
            self.process.kill()

            # Threads *should* tidy up by themselves, but we do it explicitly
            self._join_threads()

    def running(self):
        """
        True if there is a running process. False if either no process is associated with this instance,
        or if the associated process has finished.
        :return: bool
        """
        if self.process is None:
            return False
        return self.process.poll() is None
