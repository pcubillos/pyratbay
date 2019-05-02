import os
import sys
import subprocess

import numpy as np

from conftest import make_config

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay as pb

os.chdir(ROOT+'tests')


def test_call_from_command_line1():
    subprocess.call('../pbay.py -c spectrum_transmission_test.cfg'.split())


def test_call_from_command_line2():
    subprocess.call('python ../pbay.py -c spectrum_transmission_test.cfg'
                    .split())


def test_call_from_interpreter():
    pyrat = pb.run('spectrum_transmission_test.cfg')
    assert pyrat is not None
    assert repr(pyrat) == (
        "Pyrat atmospheric model\n"
        "configuration file:  'spectrum_transmission_test.cfg'\n"
        "Pressure profile (bar):  1.00e-06 -- 1.00e+02 (81 layers)\n"
        "Wavelength range (um):  1.10 -- 1.70 (3209 samples, dwn=1.000 cm-1)\n"
        "Composition:  ['H2' 'He' 'Na' 'H2O' 'CH4' 'CO' 'CO2']\n"
        "Opacity sources:  ['H2O', 'CIA H2-H2', 'CIA H2-He', 'lecavelier', 'deck', 'Na']")


def test_get_ec(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'csfile':'../inputs/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'})
    pyrat = pb.run(cfg)
    ec, label = pyrat.get_ec(50)
    assert label == ['H2O', 'CIA H2-H2', 'lecavelier', 'deck', 'Na']
    np.testing.assert_allclose(ec[0], pyrat.ex.ec[50])
    np.testing.assert_allclose(ec[1], pyrat.cs.ec[50])
    np.testing.assert_allclose(ec[2], pyrat.rayleigh.ec[50])
    np.testing.assert_allclose(ec[3], pyrat.haze.ec[50])
    np.testing.assert_allclose(ec[4], pyrat.alkali.ec[50])

