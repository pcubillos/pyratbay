# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import subprocess
import pytest

from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


@pytest.mark.skip(reason="TBD: configure MPI into github actions")
def test_call_check_mpi4py():
    # But first I need to make it think mpi4py is not installed...
    subprocess.call(
        'mpirun -n 2 pbay -c tests/configs/pt_isothermal.cfg'.split()
    )
    # Assert no warning raised


@pytest.mark.skip(reason="TBD: configure MPI into github actions")
def test_call_check_mpi_is_needed__ignored():
    subprocess.call(
        'mpirun -n 1 pbay -c tests/configs/pt_isothermal.cfg'.split()
    )
    # Assert no warning raised


@pytest.mark.skip(reason="TBD: configure MPI into github actions")
def test_call_check_mpi_is_needed__raise_warning():
    # Assert run to completion
    match = (
        'Attempting to use MPI, but this is only needed for MultiNest '
        'runs. Subprocesses will be terminated'
    )
    with pytest.warns(Warning, match=match):
        # Figure out how to trigger a mpirun call (without subprocess?)
        subprocess.call(
            'mpirun -n 2 pbay -c tests/configs/pt_isothermal.cfg'.split()
        )


@pytest.mark.skip(reason="TBD: configure MPI into github actions")
def test_get_mpi_rank():
    # Figure out how to trigger a mpirun call (without subprocess?)
    pass


@pytest.mark.skip(reason="TBD: configure MPI into github actions")
def test_get_mpi_size():
    # Figure out how to trigger a mpirun call (without subprocess?)
    pass


@pytest.mark.skip(reason="TBD: configure MPI into github actions")
def test_mpi_barrier():
    # Figure out how to trigger a mpirun call (without subprocess?)
    pass

