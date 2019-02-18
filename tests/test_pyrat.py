import os
import sys

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)

import pyratbay as pb

os.chdir(ROOT+'tests')


def test_repr():
    pyrat = pb.pyrat.init(ROOT+'tests/spectrum_transmission_test.cfg')
    assert repr(pyrat) == (
        "Pyrat atmospheric model\n"
        "configuration file:  '/Users/pato/Dropbox/IWF/projects/2014_pyratbay/2016-01-08_develop/pyratbay/tests/spectrum_transmission_test.cfg'\n"
        "Pressure profile (bar):  1.00e-06 -- 1.00e+02 (81 layers)\n"
        "Wavelength range (um):  1.70 -- 1.10 (3209 samples, dwn=1.000 cm-1)\n"
        "Composition:  ['H2' 'He' 'Na' 'H2O' 'CH4' 'CO' 'CO2']\n"
        "Opacity sources:  ['H2O', 'H2-H2', 'H2-He', 'lecavelier']")

