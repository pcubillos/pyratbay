# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os

from conftest import make_config

import pyratbay as pb
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


def test_alkali_str(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'wllow':'0.466 um', 'wlhigh':'0.80252 um',
               'alkali_cutoff':'4500.0'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    print(pyrat.alkali.models[0])
    assert str(pyrat.alkali.models[0]) == """\
Model name (name): 'sodium_vdw'
Model species (mol): Na
Species index in atmosphere (imol): 2
Profile hard cutoff from line center (cutoff, cm-1): 4500.0
Detuning parameter (detuning): 30.0
Lorentz-width parameter (lpar): 0.071
Partition function (Z): 2.0
Wavenumber  Wavelength          gf   Lower-state energy
      cm-1          um               cm-1
      (wn)                    (gf)   (elow)
  16960.87    0.589592   6.546e-01   0.000e+00
  16978.07    0.588995   1.309e+00   0.000e+00
Opacity cross section (ec, cm2 molecule-1):
[[ 0.000e+00  1.021e-29 ...  3.136e-29  3.131e-29]
 [ 0.000e+00  1.285e-29 ...  3.948e-29  3.941e-29]
 ...
 [ 0.000e+00  4.999e-21 ...  1.525e-20  1.523e-20]
 [ 0.000e+00  6.269e-21 ...  1.912e-20  1.910e-20]]
"""

def test_cloud_str(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'clouds':'deck ccsgray', 'cpars':'-3.0  0.0 -4.0 2.0'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    assert str(pyrat.cloud.models[0]) == """\
Model name (name): 'deck'
Number of model parameters (npars): 1
Parameter name     Value
  (pnames)         (pars)
  log(p_top)       -3.000e+00
Index of atmospheric layer at or directly below cloud top: 30
Cloud-top pressure: 1.0000e-03 bar
Cloud-top altitude: 72750.14 km
Cloud-top temperature: 1051.39 K
"""

    assert str(pyrat.cloud.models[1]) == """\
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


def test_rayleigh_str(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'rayleigh':'lecavelier dalgarno_H dalgarno_He dalgarno_H2',
               'rpars':'0.0 -4.0'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    assert str(pyrat.rayleigh.models[0]) == """\
Model name (name): 'lecavelier'
Model species (mol): H2
Number of model parameters (npars): 2
Parameter name     Value
  (pnames)         (pars)
  log(f_ray)        0.000e+00
  alpha_ray        -4.000e+00
Opacity cross section (ec, cm2 molec-1):
    [ 9.540e-30  9.547e-30  9.553e-30 ...  5.436e-29  5.439e-29  5.441e-29]
"""

    assert str(pyrat.rayleigh.models[1]) == """\
Model name (name): 'dalgarno_H'
Model species (mol): H
Number of model parameters (npars): 0
Extinction-coefficient (ec, cm2 molec-1):
    [ 7.002e-30  7.007e-30  7.012e-30 ...  4.038e-29  4.040e-29  4.041e-29]
"""

    assert str(pyrat.rayleigh.models[2]) == """\
Model name (name): 'dalgarno_He'
Model species (mol): He
Number of model parameters (npars): 0
Extinction-coefficient (ec, cm2 molec-1):
    [ 6.577e-31  6.582e-31  6.586e-31 ...  3.757e-30  3.758e-30  3.760e-30]
"""

    assert str(pyrat.rayleigh.models[3]) == """\
Model name (name): 'dalgarno_H2'
Model species (mol): H2
Number of model parameters (npars): 0
Extinction-coefficient (ec, cm2 molec-1):
    [ 9.799e-30  9.806e-30  9.813e-30 ...  5.626e-29  5.629e-29  5.631e-29]
"""


def test_pyrat_transmission_str(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg')
    pyrat = pb.run(cfg)
    assert pyrat is not None
    assert str(pyrat) == """\
Pyrat atmospheric model
configuration file:  '{:s}/test.cfg'
Pressure profile (bar):  1.00e-06 -- 1.00e+02 (81 layers)
Wavelength range (um):  1.10 -- 1.70 (3209 samples, dwn=1.000 cm-1)
Composition:  ['H2' 'He' 'Na' 'H2O' 'CH4' 'CO' 'CO2']
Opacity sources:  ['H2O', 'CIA H2-H2', 'CIA H2-He', 'lecavelier', 'deck', 'Na']""".format(str(tmp_path))

    assert str(pyrat.spec) == """\
Spectral information:
Wavenumber internal units: cm-1
Wavelength internal units: cm
Wavelength display units (wlunits): um
Low wavenumber boundary (wnlow):     5882.353 cm-1  (wlhigh =   1.70 um)
High wavenumber boundary (wnhigh):   9090.909 cm-1  (wllow  =   1.10 um)
Number of samples (nwave): 3209
Sampling interval (wnstep): 1.000 cm-1
Wavenumber array (wn, cm-1):
    [ 5882.353  5883.353  5884.353 ...  9088.353  9089.353  9090.353]
Oversampling factor (wnosamp): 2160

Modulation spectrum, (Rp/Rs)**2 (spectrum):
    [ 6.780e-03  6.781e-03  6.780e-03 ...  6.780e-03  6.780e-03  6.780e-03]
"""

    assert str(pyrat.atm) == """\
Atmospheric model information:
Atmospheric file name (atmfile):
    '{:s}/inputs/atmosphere_uniform_test.atm'
Number of layers (nlayers): 81

Pressure display units (punits): bar
Pressure internal units: barye
Pressure at top of atmosphere (ptop):        1.00e-06 bar
Pressure at bottom of atmosphere (pbottom):  1.00e+02 bar
Reference pressure at rplanet (refpressure): 1.00e-01 bar
Pressure profile (press, bar):
    [ 1.000e-06  1.259e-06  1.585e-06 ...  6.310e+01  7.943e+01  1.000e+02]

Radius display units (runits): rjup
Radius internal units: cm
Radius model name (rmodelname): hydro_m
Radius profile (radius, rjup):
    [1.0434 1.0425 1.0416 ... 0.9665 0.9653 0.9641]

Temperature units (tunits): kelvin
Temperature model name (tmodelname): None
Temperature profile (temp, K):
    [ 1047.045  1047.046  1047.048 ...  1663.181  1664.146  1665.360]

Mean molecular mass (mm, amu):
    [  2.3328   2.3328   2.3328 ...   2.3328   2.3328   2.3328]

Abundance units (qunits): None
Abundance internal units: mole mixing fraction
Number of atmospheric species: 7
Abundance profiles (q, mole mixing fraction):
    species [ 0]:   [ 8.500e-01  8.500e-01 ...  8.500e-01  8.500e-01]
    species [ 1]:   [ 1.490e-01  1.490e-01 ...  1.490e-01  1.490e-01]
    species [ 2]:   [ 3.000e-06  3.000e-06 ...  3.000e-06  3.000e-06]
    species [ 3]:   [ 4.000e-04  4.000e-04 ...  4.000e-04  4.000e-04]
    species [ 4]:   [ 1.000e-04  1.000e-04 ...  1.000e-04  1.000e-04]
    species [ 5]:   [ 5.000e-04  5.000e-04 ...  5.000e-04  5.000e-04]
    species [ 6]:   [ 1.000e-07  1.000e-07 ...  1.000e-07  1.000e-07]
Density profiles (d, molecules cm-3):
    species [ 0]:   [ 5.880e+12  7.402e+12 ...  2.939e+20  3.697e+20]
    species [ 1]:   [ 1.031e+12  1.298e+12 ...  5.151e+19  6.480e+19]
    species [ 2]:   [ 2.075e+07  2.613e+07 ...  1.037e+15  1.305e+15]
    species [ 3]:   [ 2.767e+09  3.483e+09 ...  1.383e+17  1.740e+17]
    species [ 4]:   [ 6.918e+08  8.708e+08 ...  3.457e+16  4.349e+16]
    species [ 5]:   [ 3.459e+09  4.354e+09 ...  1.729e+17  2.175e+17]
    species [ 6]:   [ 6.918e+05  8.708e+05 ...  3.457e+13  4.349e+13]
""".format(os.getcwd())

    assert str(pyrat.mol) == f"""\
Atmospheric species information:
Number of species (nmol): 7

Molecule    Mass       Radius
            g/mol      Angstrom
(name)      (mass)     (radius)
  H2          2.0159       1.445
  He          4.0026       1.400
  Na         22.9898       2.270
  H2O        18.0153       1.600
  CH4        16.0425       2.000
  CO         28.0101       1.690
  CO2        44.0095       1.900
Molecular data taken from (molfile):
    '{os.path.realpath('./../pyratbay/data')}/molecules.dat'
"""

    assert str(pyrat.lt) == f"""\
Line-transition information:
Input TLI files (tlifile):
    ['{os.getcwd()}/outputs/HITRAN_H2O_1.1-1.7um_test.tli']
Number of databases (ndb): 1

Database name (name): HITRAN H2O
Species name (molname):  H2O
Number of isotopes (niso): 9
Isotope correlative index (iiso): 0
Number of temperature samples (ntemp): 503
Temperature (temp, K):
    [1.00e+00 5.00e+00 1.00e+01 ... 4.98e+03 4.99e+03 5.00e+03]
Partition function for each isotope (z):
    [ 1.000e+00  1.010e+00  1.328e+00 ...  8.301e+04  8.358e+04  8.416e+04]
    [ 1.000e+00  1.010e+00  1.332e+00 ...  7.706e+04  7.758e+04  7.811e+04]
    [ 6.000e+00  6.059e+00  7.981e+00 ...  4.623e+05  4.654e+05  4.686e+05]
    [ 6.000e+00  6.213e+00  8.396e+00 ...  4.835e+05  4.865e+05  4.895e+05]
    [ 6.000e+00  6.219e+00  8.445e+00 ...  4.702e+05  4.733e+05  4.763e+05]
    [ 3.600e+01  3.729e+01  5.053e+01 ...  2.719e+06  2.737e+06  2.754e+06]
    [ 6.000e+00  6.343e+00  9.129e+00 ...  9.504e+05  9.578e+05  9.652e+05]
    [ 6.000e+00  6.353e+00  9.217e+00 ...  9.775e+05  9.850e+05  9.927e+05]
    [ 3.600e+01  3.809e+01  5.505e+01 ...  5.784e+06  5.829e+06  5.874e+06]

Total number of line transitions (ntransitions): 47,666
Minimum and maximum temperatures (tmin, tmax): [1.0, 5000.0] K
Line-transition isotope IDs (isoid):
    [0 0 0 0 0 0 0 ... 3 3 3 3 3 3 3]
Line-transition wavenumbers (wn, cm-1):
    [5882.494 5883.065 5883.353 ... 7494.719 7504.756 7513.791]
Line-transition lower-state energy (elow, cm-1):
    [ 1.807e+03  2.106e+03  2.630e+03 ...  1.244e+03  5.201e+02  6.531e+02]
Line-transition gf (gf, cm-1):
    [ 1.399e-08  1.188e-09  1.210e-08 ...  5.498e-06  1.558e-07  1.076e-06]
"""

    assert str(pyrat.iso) == """\
Isotopes information:
Number of isotopes (niso): 9

Isotope  Molecule      Mass    Isotopic   Database   Extinc-coeff
            index     g/mol       ratio      index    index
 (name)    (imol)    (mass)     (ratio)   (dbindex)  (iext)
    161         3   18.0106   9.973e-01          0   None
    181         3   20.0148   1.999e-03          0   None
    171         3   19.0148   3.719e-04          0   None
    162         3   19.0168   3.107e-04          0   None
    182         3   21.0211   6.230e-07          0   None
    172         3   20.0211   1.158e-07          0   None
    262         3   20.0210   2.420e-08          0   None
    282         3   22.0000   0.000e+00          0   None
    272         3   21.0000   0.000e+00          0   None
Partition function evaluated at atmosperic layers (z):
    [ 1.325e+03  1.325e+03  1.325e+03 ...  3.407e+03  3.412e+03  3.417e+03]
    [ 1.337e+03  1.337e+03  1.337e+03 ...  3.426e+03  3.430e+03  3.436e+03]
    [ 7.980e+03  7.980e+03  7.980e+03 ...  2.042e+04  2.045e+04  2.048e+04]
    [ 6.954e+03  6.954e+03  6.954e+03 ...  1.909e+04  1.912e+04  1.915e+04]
    [ 7.066e+03  7.066e+03  7.066e+03 ...  1.940e+04  1.943e+04  1.946e+04]
    [ 4.228e+04  4.228e+04  4.228e+04 ...  1.157e+05  1.158e+05  1.160e+05]
    [ 8.967e+03  8.967e+03  8.967e+03 ...  2.652e+04  2.656e+04  2.661e+04]
    [ 9.136e+03  9.136e+03  9.136e+03 ...  2.709e+04  2.713e+04  2.718e+04]
    [ 5.433e+04  5.433e+04  5.433e+04 ...  1.609e+05  1.611e+05  1.614e+05]
"""

    assert str(pyrat.voigt) == """\
Voigt-profile information:

Number of Doppler-width samples (ndop): 40
Number of Lorentz-width samples (nlor): 40
Doppler HWHM (doppler, cm-1):
    [ 4.963e-03  5.243e-03  5.538e-03 ...  3.765e-02  3.977e-02  4.201e-02]
Lorentz HWMH (lorentz, cm-1):
    [ 2.210e-08  3.702e-08  6.202e-08 ...  4.313e+00  7.225e+00  1.210e+01]
Doppler--Lorentz ratio threshold (dlratio): 1.000e-01

Voigt-profiles' extent (extent, in HWHMs): 100.0
Voigt-profiles' cutoff extent (cutoff in cm-1): 25.0
Voigt-profile half-sizes (size) of shape [ndop, nlor]:
[[ 1072  1132 ...  8590  9074]
 [ 1072  1132 ...  8590  9074]
 ...
 [54000 54000 ... 54000 54000]
 [54000 54000 ... 54000 54000]]
Voigt-profile indices (index) of shape [ndop, nlor]:
[[       0     2145 ...   267160   284341]
 [  302490   304635 ...   569650   586831]
 ...
 [14927707 14927707 ... 14927707 14927707]
 [15035708 15035708 ... 15035708 15035708]]

Voigt profiles:
  profile[ 0, 0]: [ 2.85913e-08  2.86448e-08 ...  2.86448e-08  2.85913e-08]
  ...
  profile[39,39]: [ 4.99389e-03  4.99404e-03 ...  4.99404e-03  4.99389e-03]
"""

    assert str(pyrat.ex) == """\
Extinction-coefficient information:
Line-transition strength threshold (ethresh): 1.00e-15

LBL extinction coefficient for the atmospheric model (ec, cm-1) [layer, wave]:
[[2.08e-21 1.80e-14 3.89e-21 ... 3.94e-16 2.08e-17 3.67e-17]
 [2.62e-21 2.27e-14 4.90e-21 ... 4.96e-16 2.62e-17 4.62e-17]
 [5.52e-21 2.85e-14 1.03e-20 ... 6.24e-16 3.30e-17 5.82e-17]
 ...
 [6.43e-07 1.04e-06 1.15e-06 ... 2.10e-06 1.79e-06 1.51e-06]
 [8.11e-07 1.04e-06 1.18e-06 ... 2.67e-06 2.25e-06 1.85e-06]
 [1.02e-06 1.31e-06 1.48e-06 ... 3.36e-06 2.84e-06 2.33e-06]]
Extinction-coefficient table filename(s) (extfile): None
"""

    assert str(pyrat.cs) == """\
Cross-section extinction information:
Number of cross-section files (nfiles): 2

Cross-section file name (files[0]):
    '{:s}/pyratbay/data/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
Cross-section species (molecules): H2-H2
Number of temperature samples: 20
Number of wavenumber samples: 824
Temperature array (temp, K):
  [  60.  100.  150.  200.  250.  300.  350.  400.  500.  600.  700.  800.  900.
    1000. 2000. 3000. 4000. 5000. 6000. 7000.]
Wavenumber array (wavenumber, cm-1):
  [20.0 40.0 60.0 ... 16440.0 16460.0 16480.0]
Input extinction coefficient (absorption, cm-1 amagat-2):
[[ 3.88e-08  9.56e-08 ...  3.39e-13  3.32e-13]
 [ 3.70e-08  1.09e-07 ...  3.79e-13  3.72e-13]
 ...
 [ 3.51e-09  1.41e-08 ...  1.35e-09  1.33e-09]
 [ 3.87e-09  1.55e-08 ...  2.23e-09  2.20e-09]]

Cross-section file name (files[1]):
    '{:s}/pyratbay/data/CIA/CIA_Borysow_H2He_0050-3000K_0.3-030um.dat'
Cross-section species (molecules): H2-He
Number of temperature samples: 20
Number of wavenumber samples: 1652
Temperature array (temp, K):
  [  50.   75.  100.  150.  200.  250.  300.  350.  400.  500.  600.  700.  800.
    900. 1000. 1250. 1500. 2000. 2500. 3000.]
Wavenumber array (wavenumber, cm-1):
  [320.0 340.0 360.0 ... 33300.0 33320.0 33340.0]
Input extinction coefficient (absorption, cm-1 amagat-2):
[[ 5.14e-07  8.42e-07 ...  6.65e-23  6.49e-23]
 [ 6.02e-07  7.64e-07 ...  4.20e-29  4.05e-29]
 ...
 [ 2.93e-06  3.20e-06 ...  1.25e-13  1.25e-13]
 [ 3.09e-06  3.36e-06 ...  4.06e-13  4.05e-13]]

Minimum and maximum temperatures (tmin, tmax) in K: [60.0, 3000.0]
Atmospheric-model extinction coefficient (ec, cm-1):
[[ 2.67e-20  2.66e-20  2.65e-20 ...  2.56e-21  2.56e-21  2.55e-21]
 [ 4.23e-20  4.21e-20  4.20e-20 ...  4.05e-21  4.05e-21  4.05e-21]
 [ 6.70e-20  6.68e-20  6.66e-20 ...  6.42e-21  6.42e-21  6.41e-21]
 ...
 [ 1.19e-04  1.19e-04  1.19e-04 ...  5.55e-06  5.55e-06  5.55e-06]
 [ 1.89e-04  1.89e-04  1.88e-04 ...  8.79e-06  8.79e-06  8.78e-06]
 [ 3.00e-04  2.99e-04  2.99e-04 ...  1.39e-05  1.39e-05  1.39e-05]]
""".format(os.path.realpath('./..'), os.path.realpath('./..'))

    assert str(pyrat.od) == """\
Optical depth information:
Observing geometry (rt_path): transit
Total atmospheric extinction coefficient (ec, cm-1) [layer, wave]:
[[ 5.613e-17  1.806e-14  5.620e-17 ...  7.134e-16  3.406e-16  3.567e-16]
 [ 7.067e-17  2.274e-14  7.076e-17 ...  8.981e-16  4.288e-16  4.490e-16]
 [ 8.898e-17  2.863e-14  8.911e-17 ...  1.131e-15  5.399e-16  5.653e-16]
 ...
 [ 1.200e-04  1.202e-04  1.201e-04 ...  7.658e-06  7.355e-06  7.066e-06]
 [ 1.901e-04  1.899e-04  1.897e-04 ...  1.148e-05  1.106e-05  1.065e-05]
 [ 3.010e-04  3.006e-04  3.002e-04 ...  1.730e-05  1.677e-05  1.626e-05]]

Distance along the ray path across each layer (outside-in) at each impact
    parameter (raypath, km):
    IP[  1]: [3061.10036832]
    IP[  2]: [1268.99486459 3057.53305829]
    IP[  3]: [ 974.36813165 1267.62126565 3053.61109788]
    ...
    IP[ 80]: [ 164.81603538  165.39273132  165.94143453 ... 1097.06940799
    1426.2286258
 3434.22183058]

Maximum optical depth to calculate (maxdepth): 10.00
Layer index where the optical depth reaches maxdepth (ideep):
    [ 30  30  30  30  30  30  30 ...  30  30  30  30  30  30  30]
Maximum ideep (deepest layer reaching maxdepth): 30

Optical depth at each impact parameter, down to max(ideep) (depth):
[[ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 3.881e-08  1.249e-05  3.887e-08 ...  4.933e-07  2.355e-07  2.466e-07]
 [ 6.490e-08  2.088e-05  6.499e-08 ...  8.248e-07  3.938e-07  4.124e-07]
 ...
 [ 4.661e-05  1.232e-02  4.708e-05 ...  4.890e-04  2.336e-04  2.463e-04]
 [ 6.141e-05  1.550e-02  6.217e-05 ...  6.163e-04  2.944e-04  3.112e-04]
 [ 8.157e-05  1.950e-02  8.275e-05 ...  7.768e-04  3.710e-04  3.931e-04]]
"""

    assert str(pyrat.cloud) == """\
Cloud-opacity models (models):

Model name (name): 'deck'
Number of model parameters (npars): 1
Parameter name     Value
  (pnames)         (pars)
  log(p_top)       -3.000e+00
Index of atmospheric layer at or directly below cloud top: 30
Cloud-top pressure: 1.0000e-03 bar
Cloud-top altitude: 72750.14 km
Cloud-top temperature: 1051.39 K

Patchiness fraction (fpatchy): None
Total atmospheric cloud extinction-coefficient (ec, cm-1):
[[ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 ...
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]]
"""

    assert str(pyrat.rayleigh) == """\
Rayleigh-opacity models (models):

Model name (name): 'lecavelier'
Model species (mol): H2
Number of model parameters (npars): 2
Parameter name     Value
  (pnames)         (pars)
  log(f_ray)        0.000e+00
  alpha_ray        -4.000e+00
Opacity cross section (ec, cm2 molec-1):
    [ 9.540e-30  9.547e-30  9.553e-30 ...  5.436e-29  5.439e-29  5.441e-29]

Total atmospheric Rayleigh extinction-coefficient (ec, cm-1):
[[ 5.610e-17  5.614e-17  5.617e-17 ...  3.197e-16  3.198e-16  3.199e-16]
 [ 7.062e-17  7.067e-17  7.072e-17 ...  4.024e-16  4.026e-16  4.028e-16]
 [ 8.891e-17  8.897e-17  8.903e-17 ...  5.066e-16  5.068e-16  5.071e-16]
 ...
 [ 2.228e-09  2.230e-09  2.231e-09 ...  1.270e-08  1.270e-08  1.271e-08]
 [ 2.804e-09  2.806e-09  2.807e-09 ...  1.598e-08  1.598e-08  1.599e-08]
 [ 3.527e-09  3.529e-09  3.532e-09 ...  2.010e-08  2.011e-08  2.011e-08]]
"""

    assert str(pyrat.alkali) == """\
Alkali-opacity models (models):

Model name (name): 'sodium_vdw'
Model species (mol): Na
Species index in atmosphere (imol): 2
Profile hard cutoff from line center (cutoff, cm-1): 4500.0
Detuning parameter (detuning): 30.0
Lorentz-width parameter (lpar): 0.071
Partition function (Z): 2.0
Wavenumber  Wavelength          gf   Lower-state energy
      cm-1          um               cm-1
      (wn)                    (gf)   (elow)
  16960.87    0.589592   6.546e-01   0.000e+00
  16978.07    0.588995   1.309e+00   0.000e+00
Opacity cross section (ec, cm2 molecule-1):
[[ 0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00]
 ...
 [ 0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00]]

Total atmospheric alkali extinction-coefficient (ec, cm-1):
[[ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 ...
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]]
"""

    assert str(pyrat.obs) == """\
Observing information:
Number of data points (ndata): 0

Number of filter pass bands (nfilters): 0
"""

    assert str(pyrat.phy) == """\
Physical properties information:

Stellar effective temperature (tstar, K): 5800.0
Stellar radius (rstar, Rsun): 1.270
Stellar mass (mstar, Msun):   None
Stellar surface gravity (gstar, cm s-2): 22908.0

Planetary radius (rplanet, Rjup): 1.000
Planetary mass (mplanet, Mjup):   0.600
Planetary surface gravity (gplanet, cm s-2): 1487.2
Planetary internal temperature (tint, K):  100.0
Orbital semi-major axis (smaxis, AU): 0.0450
Planet-to-star radius ratio (rprs):   0.08092
Planetary Hill radius (rhill, Rjup):  inf
Input stellar spectrum is a blackbody at Teff = 5800.0 K.
Stellar spectrum wavenumber (starwn, cm-1):
    [  5882.353   5883.353   5884.353 ...   9088.353   9089.353   9090.353]
Stellar flux spectrum (starflux, erg s-1 cm-2 cm):
    [ 2.306e+06  2.307e+06  2.307e+06 ...  3.293e+06  3.293e+06  3.293e+06]
"""

    assert str(pyrat.ret) == """\
Retrieval information:
No retrieval parameters set.
"""


def test_pyrat_transmission_resolution_str(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'resolution':'5000.0'},
        remove=['clouds'])
    pyrat = pb.run(cfg)
    assert pyrat is not None
    assert str(pyrat) == """\
Pyrat atmospheric model
configuration file:  '{:s}/test.cfg'
Pressure profile (bar):  1.00e-06 -- 1.00e+02 (81 layers)
Wavelength range (um):  1.10 -- 1.70 (2177 samples, R=5000.0)
Composition:  ['H2' 'He' 'Na' 'H2O' 'CH4' 'CO' 'CO2']
Opacity sources:  ['H2O', 'CIA H2-H2', 'CIA H2-He', 'lecavelier', 'Na']""".format(str(tmp_path))

    assert str(pyrat.spec) == """\
Spectral information:
Wavenumber internal units: cm-1
Wavelength internal units: cm
Wavelength display units (wlunits): um
Low wavenumber boundary (wnlow):     5882.353 cm-1  (wlhigh =   1.70 um)
High wavenumber boundary (wnhigh):   9090.909 cm-1  (wllow  =   1.10 um)
Number of samples (nwave): 2177
Spectral resolving power (resolution): 5000.0
Wavenumber array (wn, cm-1):
    [ 5882.353  5883.530  5884.706 ...  9086.201  9088.018  9089.836]
Oversampling factor (wnosamp): 2160

Modulation spectrum, (Rp/Rs)**2 (spectrum):
    [ 6.522e-03  6.540e-03  6.523e-03 ...  6.670e-03  6.500e-03  6.475e-03]
"""


def test_pyrat_emission_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'rt_path': 'emission'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    assert str(pyrat.spec) == """\
Spectral information:
Wavenumber internal units: cm-1
Wavelength internal units: cm
Wavelength display units (wlunits): um
Low wavenumber boundary (wnlow):     5882.353 cm-1  (wlhigh =   1.70 um)
High wavenumber boundary (wnhigh):   9090.909 cm-1  (wllow  =   1.10 um)
Number of samples (nwave): 3209
Sampling interval (wnstep): 1.000 cm-1
Wavenumber array (wn, cm-1):
    [ 5882.353  5883.353  5884.353 ...  9088.353  9089.353  9090.353]
Oversampling factor (wnosamp): 2160

Intensity zenithal angles (raygrid, degree): [ 0. 20. 40. 60. 80.]
raygrid internal units: radian
Intensity spectra (intensity, erg s-1 cm-2 sr-1 cm):
    [ 7.741e+02  7.734e+02  7.727e+02 ...  3.549e+01  3.545e+01  3.542e+01]
    [ 7.741e+02  7.734e+02  7.727e+02 ...  3.549e+01  3.545e+01  3.542e+01]
    [ 7.741e+02  7.734e+02  7.727e+02 ...  3.549e+01  3.545e+01  3.542e+01]
    [ 7.741e+02  7.734e+02  7.727e+02 ...  3.549e+01  3.545e+01  3.542e+01]
    [ 7.741e+02  7.734e+02  7.727e+02 ...  3.549e+01  3.545e+01  3.542e+01]
Emission spectrum (spectrum, erg s-1 cm-2 cm):
    [ 2.432e+03  2.430e+03  2.428e+03 ...  1.115e+02  1.114e+02  1.113e+02]
"""

    assert str(pyrat.od) == """\
Optical depth information:
Observing geometry (rt_path): emission
Total atmospheric extinction coefficient (ec, cm-1) [layer, wave]:
[[ 5.613e-17  1.806e-14  5.620e-17 ...  7.134e-16  3.406e-16  3.567e-16]
 [ 7.067e-17  2.274e-14  7.076e-17 ...  8.981e-16  4.288e-16  4.490e-16]
 [ 8.898e-17  2.863e-14  8.911e-17 ...  1.131e-15  5.399e-16  5.653e-16]
 ...
 [ 1.200e-04  1.202e-04  1.201e-04 ...  7.658e-06  7.355e-06  7.066e-06]
 [ 1.901e-04  1.899e-04  1.897e-04 ...  1.148e-05  1.106e-05  1.065e-05]
 [ 3.010e-04  3.006e-04  3.002e-04 ...  1.730e-05  1.677e-05  1.626e-05]]

Distance across each layer along a normal ray path (raypath, km):
    [62.8 62.7 62.6 62.5 ... 86.0 85.8 85.7 85.5]

Maximum optical depth to calculate (maxdepth): 10.00
Layer index where the optical depth reaches maxdepth (ideep):
    [ 30  30  30  30  30  30  30 ...  30  30  30  30  30  30  30]
Maximum ideep (deepest layer reaching maxdepth): 30

Planck emission down to max(ideep) (B, erg s-1 cm-2 sr-1 cm):
[[ 7.487e+02  7.480e+02  7.474e+02 ...  3.370e+01  3.367e+01  3.363e+01]
 [ 7.487e+02  7.480e+02  7.474e+02 ...  3.371e+01  3.367e+01  3.363e+01]
 [ 7.487e+02  7.480e+02  7.474e+02 ...  3.371e+01  3.367e+01  3.364e+01]
 ...
 [ 7.645e+02  7.639e+02  7.632e+02 ...  3.481e+01  3.478e+01  3.474e+01]
 [ 7.687e+02  7.681e+02  7.674e+02 ...  3.511e+01  3.508e+01  3.504e+01]
 [ 7.741e+02  7.734e+02  7.727e+02 ...  3.549e+01  3.545e+01  3.542e+01]]

Optical depth at each layer along a normal ray path into the planet, down to
    max(ideep) (depth):
[[ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 3.984e-10  1.282e-07  3.989e-10 ...  5.063e-09  2.418e-09  2.531e-09]
 [ 8.992e-10  2.893e-07  9.005e-10 ...  1.143e-08  5.457e-09  5.713e-09]
 ...
 [ 1.083e-06  3.008e-04  1.092e-06 ...  1.193e-05  5.695e-06  5.995e-06]
 [ 1.412e-06  3.784e-04  1.426e-06 ...  1.502e-05  7.172e-06  7.562e-06]
 [ 1.853e-06  4.760e-04  1.875e-06 ...  1.892e-05  9.036e-06  9.544e-06]]
"""


def test_pyrat_exfile_str(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/mcmc_transmission_test.cfg',
        reset={'runmode':'spectrum',
               'specfile':f'{ROOT}tests/outputs/extfile_spectrum_test.dat'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    pyrat.band_integrate()
    assert str(pyrat.lt) == """\
Line-transition information:
No input TLI files.
"""

    assert str(pyrat.ex) == """\
Extinction-coefficient information:
Line-transition strength threshold (ethresh): 1.00e-15

LBL extinction coefficient for the atmospheric model (ec, cm-1) [layer, wave]:
[[1.99e-21 1.67e-14 3.74e-21 ... 5.24e-16 2.13e-17 3.54e-17]
 [3.28e-21 2.11e-14 6.21e-21 ... 6.59e-16 2.69e-17 4.46e-17]
 [5.20e-21 2.65e-14 9.78e-21 ... 8.30e-16 3.38e-17 5.62e-17]
 ...
 [6.41e-07 1.03e-06 1.13e-06 ... 2.07e-06 1.77e-06 1.49e-06]
 [8.07e-07 1.04e-06 1.17e-06 ... 2.64e-06 2.23e-06 1.83e-06]
 [1.02e-06 1.31e-06 1.47e-06 ... 3.32e-06 2.80e-06 2.30e-06]]
Extinction-coefficient table filename(s) (extfile):
    {:s}/outputs/exttable_test_300-3000K_1.1-1.7um.npz
Minimum temperature (tmin, K):  300.0
Maximum temperature (tmax, K): 3000.0
Temperature sampling interval (tstep, K): None

Number of species (nspec):              1
Number of temperatures (ntemp):        10
Number of layers (nlayers):            81
Number of spectral samples (nwave):  3209

Species array (species): ['H2O']
Temperature array (temp, K):
   [ 300.  600.  900. 1200. 1500. 1800. 2100. 2400. 2700. 3000.]
Partition function (z): None
Pressure array (press, bar):
   [1.00e-06 1.26e-06 1.58e-06 2.00e-06 2.51e-06 3.16e-06 3.98e-06 5.01e-06
    6.31e-06 7.94e-06 1.00e-05 1.26e-05 1.58e-05 2.00e-05 2.51e-05 3.16e-05
    3.98e-05 5.01e-05 6.31e-05 7.94e-05 1.00e-04 1.26e-04 1.58e-04 2.00e-04
    2.51e-04 3.16e-04 3.98e-04 5.01e-04 6.31e-04 7.94e-04 1.00e-03 1.26e-03
    1.58e-03 2.00e-03 2.51e-03 3.16e-03 3.98e-03 5.01e-03 6.31e-03 7.94e-03
    1.00e-02 1.26e-02 1.58e-02 2.00e-02 2.51e-02 3.16e-02 3.98e-02 5.01e-02
    6.31e-02 7.94e-02 1.00e-01 1.26e-01 1.58e-01 2.00e-01 2.51e-01 3.16e-01
    3.98e-01 5.01e-01 6.31e-01 7.94e-01 1.00e+00 1.26e+00 1.58e+00 2.00e+00
    2.51e+00 3.16e+00 3.98e+00 5.01e+00 6.31e+00 7.94e+00 1.00e+01 1.26e+01
    1.58e+01 2.00e+01 2.51e+01 3.16e+01 3.98e+01 5.01e+01 6.31e+01 7.94e+01
    1.00e+02]
Wavenumber array (wn, cm-1):
    [5882.3529 5883.3529 5884.3529 ... 9088.3529 9089.3529 9090.3529]
Tabulated extinction coefficient (etable, cm2 molecule-1) of shape
    [nmol, ntemp, nlayers, nwave]:
[[[[1.16e-32 1.09e-26 ... 1.48e-26 3.73e-28]
   [1.92e-32 1.09e-26 ... 1.48e-26 3.73e-28]
   ...
   [3.61e-25 4.10e-25 ... 1.46e-24 1.18e-24]
   [3.34e-25 3.67e-25 ... 1.46e-24 1.24e-24]]

  [[1.88e-31 1.44e-24 ... 1.35e-26 7.00e-27]
   [3.11e-31 1.44e-24 ... 1.35e-26 7.00e-27]
   ...
   [2.95e-24 3.32e-24 ... 4.11e-24 3.39e-24]
   [2.95e-24 3.32e-24 ... 4.11e-24 3.39e-24]]

  ...

  [[2.06e-31 3.91e-24 ... 6.80e-27 1.11e-26]
   [2.06e-31 3.91e-24 ... 6.80e-27 1.11e-26]
   ...
   [3.07e-24 4.90e-24 ... 1.17e-23 9.96e-24]
   [3.23e-24 4.23e-24 ... 1.14e-23 9.48e-24]]

  [[1.64e-31 3.12e-24 ... 6.85e-27 1.13e-26]
   [1.64e-31 3.12e-24 ... 6.85e-27 1.13e-26]
   ...
   [2.47e-24 3.95e-24 ... 9.96e-24 8.48e-24]
   [2.62e-24 3.45e-24 ... 9.71e-24 8.05e-24]]]]
""".format(os.getcwd())

    assert str(pyrat.obs) == """\
Observing information:
Data/bandflux display units (units): none
Data/bandflux internal units: none
Number of data points (ndata): 20
        Data  Uncertainty   Wavenumber  Wavelength
        none         none         cm-1          um
      (data)     (uncert)     (bandwn)
     0.00661      0.00002      8826.64       1.133
     0.00660      0.00002      8635.75       1.158
     0.00660      0.00002      8450.01       1.183
     0.00651      0.00002      8271.44       1.209
     0.00645      0.00002      8097.29       1.235
     0.00641      0.00002      7936.90       1.260
     0.00647      0.00002      7782.21       1.285
     0.00648      0.00002      7631.05       1.310
     0.00666      0.00002      7485.12       1.336
     0.00673      0.00002      7345.18       1.361
     0.00677      0.00002      7207.29       1.387
     0.00674      0.00002      7077.43       1.413
     0.00676      0.00002      6951.76       1.438
     0.00670      0.00002      6830.87       1.464
     0.00667      0.00002      6715.97       1.489
     0.00658      0.00002      6600.71       1.515
     0.00656      0.00002      6493.74       1.540
     0.00646      0.00002      6387.78       1.565
     0.00650      0.00002      6285.57       1.591
     0.00649      0.00002      6188.15       1.616

Number of filter pass bands (nfilters): 20
Wavenumber  Wavelength    Bandflux  Filter name
      cm-1          um        none
  (bandwn)              (bandflux)  (filters)
   8826.64       1.133     0.00661  filter_test_WFC3_G141_1.133um.dat
   8635.75       1.158     0.00659  filter_test_WFC3_G141_1.158um.dat
   8450.01       1.183     0.00658  filter_test_WFC3_G141_1.183um.dat
   8271.44       1.209     0.00651  filter_test_WFC3_G141_1.209um.dat
   8097.29       1.235     0.00648  filter_test_WFC3_G141_1.235um.dat
   7936.90       1.260     0.00646  filter_test_WFC3_G141_1.260um.dat
   7782.21       1.285     0.00651  filter_test_WFC3_G141_1.285um.dat
   7631.05       1.310     0.00654  filter_test_WFC3_G141_1.310um.dat
   7485.12       1.336     0.00666  filter_test_WFC3_G141_1.336um.dat
   7345.18       1.361     0.00672  filter_test_WFC3_G141_1.361um.dat
   7207.29       1.387     0.00672  filter_test_WFC3_G141_1.387um.dat
   7077.43       1.413     0.00672  filter_test_WFC3_G141_1.413um.dat
   6951.76       1.438     0.00671  filter_test_WFC3_G141_1.438um.dat
   6830.87       1.464     0.00668  filter_test_WFC3_G141_1.464um.dat
   6715.97       1.489     0.00664  filter_test_WFC3_G141_1.489um.dat
   6600.71       1.515     0.00660  filter_test_WFC3_G141_1.515um.dat
   6493.74       1.540     0.00657  filter_test_WFC3_G141_1.540um.dat
   6387.78       1.565     0.00652  filter_test_WFC3_G141_1.565um.dat
   6285.57       1.591     0.00651  filter_test_WFC3_G141_1.591um.dat
   6188.15       1.616     0.00651  filter_test_WFC3_G141_1.616um.dat
"""

    assert str(pyrat.ret) == """\
Retrieval information:
  Parameter name        value        pmin        pmax       pstep  Model type
  (pnames)           (params)      (pmin)      (pmax)     (pstep)  (retflag)
  log(kappa')      -5.000e+00  -9.000e+00   5.000e+00   3.000e-01  temp
  log(gamma1)       0.000e+00  -3.000e+00   3.000e+00   3.000e-01  temp
  log(gamma2)       0.000e+00  -3.000e+00   3.000e+00   0.000e+00  temp
  alpha             0.000e+00   0.000e+00   1.000e+00   0.000e+00  temp
  T_irr (K)         1.486e+03   0.000e+00   7.000e+03   5.000e+01  temp
  T_int (K)         1.000e+02   0.000e+00   5.000e+02   0.000e+00  temp
  Rp (km)           7.150e+04   3.000e+04   1.500e+05   1.000e+02  rad
  log(H2O)         -4.000e+00  -9.000e+00  -1.000e+00   5.000e-01  mol

Retrieval algorithm (sampler): snooker
Number of retrieval samples (nsamples): 300
Number of parallel chains (nchains):   21
Number of burned-in samples (burnin):  10
Thinning factor (thinning): 1

Upper boundary for sum of metal abundances (qcap): 1.0
Temperature upper boundary (tlow, K):   300.0
Temperature lower boundary (thigh, K): 3000.0

Retrieval posterior file (mcmcfile):
    {:s}/outputs/MCMC_transmission_test.npz
""".format(os.getcwd())

