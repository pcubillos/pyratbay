import os
import subprocess

import numpy as np

from conftest import make_config
import pyratbay as pb
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


def test_call_from_command_line1():
    subprocess.call('pbay -c spectrum_transmission_test.cfg'.split())


def test_call_from_command_line2():
    subprocess.call('python pbay -c spectrum_transmission_test.cfg'.split())


def test_call_from_interpreter():
    pyrat = pb.run('spectrum_transmission_test.cfg')
    assert pyrat is not None


def test_get_ec(tmp_path):
    cfg = make_config(tmp_path, ROOT + 'tests/spectrum_transmission_test.cfg',
        reset={'csfile': ROOT
               + 'inputs/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'})
    pyrat = pb.run(cfg)
    ec, label = pyrat.get_ec(50)
    assert label == ['H2O', 'CIA H2-H2', 'lecavelier', 'deck', 'Na']
    np.testing.assert_allclose(ec[0], pyrat.ex.ec[50])
    np.testing.assert_allclose(ec[1], pyrat.cs.ec[50])
    np.testing.assert_allclose(ec[2], pyrat.rayleigh.ec[50])
    np.testing.assert_allclose(ec[3], pyrat.cloud.ec[50])
    np.testing.assert_allclose(ec[4], pyrat.alkali.ec[50])

