# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os
import subprocess

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
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg')
    pyrat = pb.run(cfg)
    logfile = pyrat.log.logname

    # Now, initialize without overwritting log:
    pyrat = pb.Pyrat(cfg, no_logfile=True)
    assert pyrat.log.logname is None
    with open(logfile, 'r') as f:
        log = f.read()
    assert 'Computed transmission spectrum' in log

    # Continue running, check it doesn't break:
    pyrat.setup_spectrum()
    pyrat.run()
    assert pyrat is not None


def test_mute(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg')
    pyrat = pb.Pyrat(cfg, mute=True)
    captured = capfd.readouterr()
    assert captured.out == ''


def test_get_ec(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'csfile': ROOT
               + 'pyratbay/data/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'})
    pyrat = pb.run(cfg)
    ec, label = pyrat.get_ec(50)
    assert label == ['H2O', 'CIA H2-H2', 'lecavelier', 'deck', 'Na']
    np.testing.assert_allclose(ec[0], pyrat.ex.ec[50])
    np.testing.assert_allclose(ec[1], pyrat.cs.ec[50])
    np.testing.assert_allclose(ec[2], pyrat.rayleigh.ec[50])
    # Cloud deck model does not use the ec, rather post processed during RT.
    # This array contains zeros or ones whether one is above or below cloud:
    np.testing.assert_allclose(ec[3], np.ones(pyrat.spec.nwave))
    np.testing.assert_allclose(ec[4], pyrat.alkali.ec[50])

