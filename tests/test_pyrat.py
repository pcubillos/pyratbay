import os
import sys
import subprocess

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)

import pyratbay as pb

os.chdir(ROOT+'tests')


def test_call_from_command_line():
    subprocess.call(['python', '../pbay.py', '-c',
                     'spectrum_transmission_test.cfg'])


def test_call_from_command_line2():
    subprocess.call(['../pbay.py', '-c', 'spectrum_transmission_test.cfg'])


def test_call_from_interpreter():
    pyrat = pb.run('spectrum_transmission_test.cfg')
    assert pyrat is not None
    assert repr(pyrat) == (
        "Pyrat atmospheric model\n"
        "configuration file:  'spectrum_transmission_test.cfg'\n"
        "Pressure profile (bar):  1.00e-06 -- 1.00e+02 (81 layers)\n"
        "Wavelength range (um):  1.10 -- 1.70 (3209 samples, dwn=1.000 cm-1)\n"
        "Composition:  ['H2' 'He' 'Na' 'H2O' 'CH4' 'CO' 'CO2']\n"
        "Opacity sources:  ['H2O', 'CIA H2-H2', 'CIA H2-He', 'lecavelier', 'Na']")

