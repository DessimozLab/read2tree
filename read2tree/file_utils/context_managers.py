import os
import shutil
import tempfile

__all__ = ['TempFile', 'TempDir', 'ChDir', 'MkDir', 'NonDeletingTempDir']


class TempFile(object):
    """ 
    Context manager for working with a temporary file
    that automatically cleans up.

    Usage:

    with TempFile() as tmp:
        # In scope, tmp exists on the disk
        # Do some work with tmp, e.g. tmp.write('something')

    # Out of scope, tmp is deleted

    with TempFile('local_temp_space') as tmp:
        # tmp is created in the directory 'local_temp_space'
        # The specified directory must exist, or an error is thrown

    """

    def __init__(self, dir_=None):
        if dir_ is not None and not os.path.exists(dir_):
            raise IOError('Directory "{}"" does not exist'.format(dir_))
        self.dir = dir_

    def __enter__(self):
        self._fd, self._wrapped_tmp = tempfile.mkstemp(dir=self.dir)
        return os.path.abspath(self._wrapped_tmp)

    def __exit__(self, type, value, tb):
        os.close(self._fd)
        os.remove(self._wrapped_tmp)


class TempDir(object):
    """
    Context manager for working with a temporary file
    that automatically cleans up.

    Usage:

    with TempDir() as tmpd:
        # In scope, tmpd exists on the disk
        # Do some work with tmpd ...

    # Out of scope, tmpd is deleted along with all its content

    Can be nested with TempFile, e.g.

    with TempDir() as tmpd, TempFile(tmpd) as tmpf:
        # tempfile tmpf is created inside temporary directory tmpd
    # On exit, everything is deleted

    """

    def __enter__(self):
        self._wrapped_tmpdir = tempfile.mkdtemp()
        return os.path.abspath(self._wrapped_tmpdir)

    def __exit__(self, type, value, tb):
        shutil.rmtree(self._wrapped_tmpdir)


class NonDeletingTempDir(TempDir):
    def __exit__(self, tpye, value, tb):
        pass


class ChDir(object):
    """
    Context manager to switch to a working directory,
    and return to the current directory (like 'Dir.chdir do' block in Ruby)

    Usage:

    with TempDir() as dir, ChDir(dir):
        # Do some work in the working temp directory 'dir'

    # Exit 'dir'
    """

    def __init__(self, working_dir):
        if not os.path.exists(working_dir):
            raise IOError('Directory "{}"" does not exist'.format(working_dir))
        self._cdir = os.getcwd()
        self._wdir = working_dir

    def __enter__(self):
        os.chdir(self._wdir)

    def __exit__(self, type, value, tb):
        os.chdir(self._cdir)


class MkDir(ChDir):
    """
    Context manager to create and switch to a working directory,
    then return to the current directory.

    Usage:

    with TempDir() as dir, MkDir(dir):
        # Do some work in the working temp directory 'dir'

    # Exit 'dir'
    """

    def __init__(self, working_dir):
        if not os.path.exists(working_dir):
            try:
                os.makedirs(working_dir)
            except OSError as e:
                if e.errno != 17:
                    raise
                pass  # path was created by another thread / process
                # this is a race condition, but probably benign

    def __enter__(self):
        pass

    def __exit__(self, type, value, tb):
        pass
