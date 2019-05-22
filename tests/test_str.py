import os
import sys
import pytest

import numpy as np

from conftest import make_config

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay as pb

os.chdir(ROOT+'tests')


def test_pyrat_str(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg')
    pyrat = pb.run(cfg)
    print(tmp_path)
    assert pyrat is not None
    assert str(pyrat) == """\
Pyrat atmospheric model
configuration file:  '{:s}/test.cfg'
Pressure profile (bar):  1.00e-06 -- 1.00e+02 (81 layers)
Wavelength range (um):  1.10 -- 1.70 (3209 samples, dwn=1.000 cm-1)
Composition:  ['H2' 'He' 'Na' 'H2O' 'CH4' 'CO' 'CO2']
Opacity sources:  ['H2O', 'CIA H2-H2', 'CIA H2-He', 'lecavelier', 'deck', 'Na']""".format(str(tmp_path))


def test_alkali_str(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg')
    pyrat = pb.run(cfg)
    assert pyrat is not None
    assert str(pyrat.alkali.models[0]) == """\
Model name (name): 'SodiumVdWst'
Model species (mol): Na
Species index in atmosphere (imol): 2
Detuning parameter (detuning): 30.0
Lorentz-width parameter (lpar): 0.071
Partition function (Z): 2.0
Wavenumber  Wavelength          gf   Lower-state energy
      cm-1          um               cm-1
      (wn)                    (gf)   (elow)
  16960.87    0.589592   6.546e-01   0.000e+00
  16978.07    0.588995   1.309e+00   0.000e+00
Extinction-coefficient (ec, cm-1):
[[ 1.915e-26  1.918e-26 ...  2.621e-24  2.625e-24]
 [ 3.036e-26  3.040e-26 ...  4.154e-24  4.161e-24]
 ...
 [ 1.347e-08  1.349e-08 ...  3.597e-07  3.601e-07]
 [ 2.134e-08  2.136e-08 ...  5.687e-07  5.693e-07]]
"""

def test_cloud__deck_str(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'hazes':'deck', 'hpars':'-3'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    assert str(pyrat.haze.models[0]) == """\
Model name (name): 'deck'
Number of model parameters (npars): 1
Parameter name     Value
  (pnames)         (pars)
  log(p_top)       -3.000e+00
Extinction-coefficient (ec, cm-1):
[[ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 ...
 [ 1.401e-02  1.401e-02  1.401e-02 ...  1.401e-02  1.401e-02  1.401e-02]
 [ 1.764e-02  1.764e-02  1.764e-02 ...  1.764e-02  1.764e-02  1.764e-02]
 [ 2.220e-02  2.220e-02  2.220e-02 ...  2.220e-02  2.220e-02  2.220e-02]]
"""


def test_cloud_ccsgray_str(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'hazes':'ccsgray', 'hpars':'0 -4 2'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    assert str(pyrat.haze.models[0]) == """\
Model name (name): 'ccsgray'
Model species (mol): H2
Number of model parameters (npars): 3
Parameter name     Value
  (pnames)         (pars)
  log(f_gray)       0.000e+00
  log(p_top)       -4.000e+00
  log(p_bot)        2.000e+00
Extinction-coefficient (ec, cm2 molec-1):
[[ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 ...
 [ 5.310e-27  5.310e-27  5.310e-27 ...  5.310e-27  5.310e-27  5.310e-27]
 [ 5.310e-27  5.310e-27  5.310e-27 ...  5.310e-27  5.310e-27  5.310e-27]
 [ 5.310e-27  5.310e-27  5.310e-27 ...  5.310e-27  5.310e-27  5.310e-27]]
"""

