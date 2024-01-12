# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import subprocess
import pytest

import numpy as np

from conftest import make_config
import pyratbay as pb
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


def test_call_from_command_line1():
    subprocess.call('pbay -c configs/spectrum_transmission_test.cfg'.split())


def test_call_from_command_line2():
    subprocess.call(
        'python pbay -c configs/spectrum_transmission_test.cfg'.split())


def test_call_from_interpreter():
    pyrat = pb.run('configs/spectrum_transmission_test.cfg')
    assert pyrat is not None


def test_no_logfile(tmp_path):
    # Run a spectrum:
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    with open(pyrat.log.logname, 'r') as f:
        log = f.read()
    assert 'Computed transit spectrum' in log

    # Now, initialize without overwritting log:
    pyrat = pb.Pyrat(cfg, log=False)
    assert pyrat.log.logname is None


def test_mute(tmp_path, capfd):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.Pyrat(cfg, mute=True)
    captured = capfd.readouterr()
    assert captured.out == ''


def test_get_ec_lbl(tmp_path):
    cfile = f'{ROOT}pyratbay/data/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'csfile': cfile, 'wllow': '0.55'},
    )
    pyrat = pb.run(cfg)
    layer = 31
    ec, labels = pyrat.get_ec(layer)

    expected_labels = ['H2O', 'Na', 'CIA H2-H2', 'lecavelier', 'deck']
    with np.load(f'{ROOT}tests/expected/expected_get_ec_lbl.npz') as d:
        expected_extinction = d['ec']

    assert labels == expected_labels
    np.testing.assert_allclose(ec, expected_extinction)
    # Cloud deck model does not use the ec, rather post processed during RT.
    # This array contains zeros or ones whether one is above or below cloud:


def test_get_ec_line_sample(tmp_path):
    extfile = ROOT+'tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz'
    reset = {
        'extfile': extfile,
        'chemistry': 'tea',
        'species': 'H2 H He Na K H2O CH4 CO CO2 e- H- H+ H2+ Na- Na+ K+ K-',
        'h_ion': 'h_ion_john1988',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset=reset,
        remove=['tlifile', 'clouds'],
    )
    pyrat = pb.run(cfg)
    layer = 31
    ec, labels = pyrat.get_ec(layer)

    expected_labels = [
        'H2O', 'Na', 'CIA H2-H2', 'CIA H2-He', 'lecavelier',
        'H- bound-free/free-free',
    ]
    with np.load(ROOT+'tests/expected/expected_get_ec_ls.npz') as d:
        expected_extinction = d['ec']

    assert labels == expected_labels
    np.testing.assert_allclose(ec, expected_extinction)
