#!/usr/bin/env python
'''
    Utilities for parsing the annotations files.

    -- Alex Warwick Vesztrocy - March--June 2016
'''
import bz2
import gzip
import os
import sys


# File opening. This is based on the example on SO here:
# http://stackoverflow.com/a/26986344
fmagic = {b'\x1f\x8b\x08': gzip.open,
          b'\x42\x5a\x68': bz2.BZ2File}


def auto_open(fn, *args):
    '''
        Opens files based on their "magic bytes". Supports bz2 and gzip. If it
        finds neither of these, presumption is it is a standard, uncompressed
        file.
    '''
    if os.path.isfile(fn) and os.stat(fn).st_size > 0:
        with open(fn, 'rb') as fp:
            fs = fp.read(max([len(x) for x in fmagic]))
        for (magic, _open) in fmagic.items():
            if fs.startswith(magic):
                return _open(fn, *args)
    else:
        if fn.endswith('gz'):
            return gzip.open(fn, *args)
        elif fn.endswith('bz2'):
            return bz2.BZ2File(fn, *args)

    return open(fn, *args)


def exe_name():
    '''
        Return the executable's basename, for inclusion in the help (with the
        help of argparse).
    '''
    return os.path.basename(sys.argv[0])


class LazyProperty(object):
    '''
        Decorator to evaluate a property only on access.

        Compute the attribute value and caches it in the instance.
        Python Cookbook (Denis Otkidach)
        http://stackoverflow.com/users/168352/denis-otkidach
        This decorator allows you to create a property which can be computed
        once and accessed many times.

        (Include from pyoma.browser.models - Adrian Altenhoff)
    '''
    def __init__(self, method, name=None):
        # record the unbound-method and the name
        self.method = method
        self.name = name or method.__name__
        self.__doc__ = method.__doc__

    def __get__(self, inst, cls):
        if inst is None:
            return self
        # compute, cache and return the instance's attribute value
        result = self.method(inst)
        # setattr redefines the instance's attribute so this doesn't get called
        # again
        setattr(inst, self.name, result)
        return result


def get_job_id():
    '''
        Gets job ID.
    '''
    if 'JOB_ID' in os.environ:
        # SGE
        return int(os.environ['JOB_ID'])
    elif 'LSB_JOBID' in os.environ:
        # LSF
        return int(os.environ['LSB_JOBID'])
    elif 'PBS_JOBID' in os.environ:
        # PBS / Torque
        return int(os.environ['PBS_JOBID'])
    elif 'SLURM_ARRAY_JOB_ID' in os.environ:
        # Slurm
        return int(os.environ['SLURM_ARRAY_JOB_ID'])
    else:
        # No parallelism detected.
        return None


def get_worker_id():
    '''
        Gets worker ID from the array ID in the job handler.
        number of workers.
    '''
    try:
        if 'SGE_TASK_ID' in os.environ:
            # SGE
            return int(os.environ['SGE_TASK_ID'])
        elif 'LSB_JOBINDEX' in os.environ:
            return int(os.environ['LSB_JOBINDEX'])
        elif 'PBS_ARRAYID' in os.environ:
            # PBS / Torque
            return int(os.environ['PBS_ARRAYID'])
        elif 'SLURM_ARRAY_TASK_ID' in os.environ:
            # Slurm
            return int(os.environ['SLURM_ARRAY_TASK_ID'])
    except ValueError:
        # int() to base10 error
        pass

    # No parallelism detected.
    return None


def check_array_ids(args):
    '''
        Checks the IDs added to args for array jobs. Raises errors if not setup
        correctly.
    '''
    if args.worker_id > args.array or args.worker_id == 0:
        raise RuntimeError('Recognised: worker ID {} and array size {}. '
                           'Worker IDs should run from 1-N (N is array size'
                           ').'.format(args.worker_id, args.array))
    if args.job_id is None or args.worker_id is None:
        raise RuntimeError('User requested HOGPROP to run as job array.'
                           'Can\'t find job ID ({}) or array ID ({}).'
                           .format(args.job_id, args.worker_id))
