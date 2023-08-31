# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'check_mpi4py',
    'check_mpi_is_needed',
    'get_mpi_rank',
    'get_mpi_size',
    'mpi_barrier',
]

import importlib
import os
import sys
import warnings


def check_mpi4py():
    """
    Detect when the code was called with MPI and mpi4py module is missing

    Only raise an error when needed (more than one processor required),
    otherwise you might be running multiple runs in parallel but not
    talking to each other.
    """
    size = 1
    # Detect MPI in call (might not be exhaustive)
    if 'OMPI_COMM_WORLD_SIZE' in os.environ:
       size = int(os.environ['OMPI_COMM_WORLD_SIZE'])
    elif 'PMI_SIZE' in os.environ:
       size = int(os.environ['PMI_SIZE'])

    # Detect mpi4pi package is installed
    mpi4py_exists = importlib.util.find_spec('mpi4py') is not None

    # Complain only if necessary:
    if size > 1 and not mpi4py_exists:
        raise ModuleNotFoundError(
            "Attempted to run pyratbay with MPI, but module 'mpi4py' is not "
            "installed. Run 'pip install mpi4py' and try again"
        )

def check_mpi_is_needed(inputs):
    """
    Prevent using parallel processes through MPI when not needed
    (only required for MultiNest runs).
    """
    size = get_mpi_size()
    rank = get_mpi_rank()
    mpi_needed = (
        inputs.runmode == 'retrieval' and
        inputs.sampler == 'multinest'
    )

    if size > 1 and not mpi_needed:
        # Keep only rank-zero process to reach completion
        msg = (
            'Attempting to use MPI, but this is only needed for MultiNest '
            'runs. Subprocesses will be terminated'
        )
        if rank == 0:
            warnings.warn(msg, category=Warning)
        else:
            sys.exit(0)


def get_mpi_rank():
    """
    Get the MPI rank of the current process (intended for MPI runs).
    If mpi4py is not installed, return zero.

    Returns
    -------
    rank: Interger
        The MPI process rank.
    """
    rank = 0
    if 'PBAY_NO_MPI' in os.environ:
        return rank
    mpi_exists = importlib.util.find_spec('mpi4py') is not None
    if mpi_exists:
        from mpi4py import MPI
        rank = MPI.COMM_WORLD.Get_rank()
    return rank


def get_mpi_size():
    """
    Get the size of the current group of process (intended for MPI runs).
    If mpi4py is not installed, return one.

    Returns
    -------
    size: Interger
        The size of the MPI group of processes.
    """
    size = 1
    if 'PBAY_NO_MPI' in os.environ:
        return size
    mpi_exists = importlib.util.find_spec('mpi4py') is not None
    if mpi_exists:
        from mpi4py import MPI
        size = MPI.COMM_WORLD.Get_size()
    return size


def mpi_barrier():
    """
    Make an MPI barrier() call. Ignore it if mpi4py is not installed.
    """
    if 'PBAY_NO_MPI' in os.environ:
        return
    mpi_exists = importlib.util.find_spec('mpi4py') is not None
    if mpi_exists:
        from mpi4py import MPI
        MPI.COMM_WORLD.barrier()
