# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import pytest

from conftest import make_config

import pyratbay as pb
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


def test_alkali_str(tmp_path):
    reset = {
        'wllow': '0.466 um',
        'wlhigh': '0.80252 um',
        'alkali_cutoff': '4500.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    assert pyrat is not None
    # print(pyrat.alkali.models[0])
    assert str(pyrat.alkali.models[0]) == """\
Model name (name): 'sodium_vdw'
Model species (mol): Na
Species mass (mass, amu): 22.989769
Profile hard cutoff from line center (cutoff, cm-1): 4500.0
Detuning parameter (detuning): 30.0
Lorentz-width parameter (lpar): 0.071
Partition function (Z): 2.0
Wavenumber  Wavelength          gf   Lower-state energy
      cm-1          um               cm-1
      (wn)                    (gf)   (elow)
  16960.87    0.589592   6.546e-01   0.000e+00
  16978.07    0.588995   1.309e+00   0.000e+00
Wavenumber (wn, cm-1):
   [12460.75 12461.75 12462.75 ... 21456.75 21457.75 21458.75]
Pressure (pressure, barye):
[1.000e+00 1.259e+00 1.585e+00 1.995e+00 2.512e+00 3.162e+00 3.981e+00
 5.012e+00 6.310e+00 7.943e+00 1.000e+01 1.259e+01 1.585e+01 1.995e+01
 2.512e+01 3.162e+01 3.981e+01 5.012e+01 6.310e+01 7.943e+01 1.000e+02
 1.259e+02 1.585e+02 1.995e+02 2.512e+02 3.162e+02 3.981e+02 5.012e+02
 6.310e+02 7.943e+02 1.000e+03 1.259e+03 1.585e+03 1.995e+03 2.512e+03
 3.162e+03 3.981e+03 5.012e+03 6.310e+03 7.943e+03 1.000e+04 1.259e+04
 1.585e+04 1.995e+04 2.512e+04 3.162e+04 3.981e+04 5.012e+04 6.310e+04
 7.943e+04 1.000e+05 1.259e+05 1.585e+05 1.995e+05 2.512e+05 3.162e+05
 3.981e+05 5.012e+05 6.310e+05 7.943e+05 1.000e+06 1.259e+06 1.585e+06
 1.995e+06 2.512e+06 3.162e+06 3.981e+06 5.012e+06 6.310e+06 7.943e+06
 1.000e+07 1.259e+07 1.585e+07 1.995e+07 2.512e+07 3.162e+07 3.981e+07
 5.012e+07 6.310e+07 7.943e+07 1.000e+08]
Cross section (cross_section, cm2 molecule-1):
[[0.000e+00 1.089e-29 ... 3.344e-29 3.338e-29]
 [0.000e+00 1.371e-29 ... 4.210e-29 4.203e-29]
 ...
 [0.000e+00 5.273e-21 ... 1.608e-20 1.606e-20]
 [0.000e+00 6.612e-21 ... 2.017e-20 2.014e-20]]
"""

def test_cloud_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'clouds':'deck ccsgray', 'cpars':'-3.0  0.0 -4.0 2.0'},
    )
    pyrat = pb.run(cfg)
    assert pyrat is not None
    assert str(pyrat.cloud.models[0]) == """\
Model name (name): 'deck'
Number of model parameters (npars): 1
Parameter name     Value
  (pnames)         (pars)
  log_p_cl         -3.000e+00
Index of atmospheric layer at or directly below cloud top: 30
Cloud-top pressure: 1.0000e-03 bar
Cloud-top altitude: 72750.08 km
Cloud-top temperature: 1051.39 K
"""

    assert str(pyrat.cloud.models[1]) == """\
Model name (name): 'ccsgray'
Model species (mol): H2
Number of model parameters (npars): 3
Parameter name     Value
  (pnames)         (pars)
  log_k_gray        0.000e+00
  log_p_top        -4.000e+00
  log_p_bot         2.000e+00
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
  log_k_ray         0.000e+00
  alpha_ray        -4.000e+00
Wavenumber (wn, cm-1):
   [5882.35 5883.35 5884.35 ... 9088.35 9089.35 9090.35]
Cross section (cross_section, cm2 molec-1):
   [ 9.540e-30  9.547e-30  9.553e-30 ...  5.436e-29  5.439e-29  5.441e-29]
"""

    assert str(pyrat.rayleigh.models[1]) == """\
Model name (name): 'dalgarno_H'
Model species (mol): H
Number of model parameters (npars): 0
Wavenumber (wn, cm-1):
   [5882.35 5883.35 5884.35 ... 9088.35 9089.35 9090.35]
Cross section (cross_section, cm2 molec-1):
   [7.002e-30 7.007e-30 7.012e-30 ... 4.038e-29 4.040e-29 4.041e-29]
"""

    assert str(pyrat.rayleigh.models[2]) == """\
Model name (name): 'dalgarno_He'
Model species (mol): He
Number of model parameters (npars): 0
Wavenumber (wn, cm-1):
   [5882.35 5883.35 5884.35 ... 9088.35 9089.35 9090.35]
Cross section (cross_section, cm2 molec-1):
   [6.577e-31 6.582e-31 6.586e-31 ... 3.757e-30 3.758e-30 3.760e-30]
"""

    assert str(pyrat.rayleigh.models[3]) == """\
Model name (name): 'dalgarno_H2'
Model species (mol): H2
Number of model parameters (npars): 0
Wavenumber (wn, cm-1):
   [5882.35 5883.35 5884.35 ... 9088.35 9089.35 9090.35]
Cross section (cross_section, cm2 molec-1):
   [9.799e-30 9.806e-30 9.813e-30 ... 5.626e-29 5.629e-29 5.631e-29]
"""


def test_pyrat_transmission_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    assert pyrat is not None
    assert str(pyrat) == f"""\
Pyrat atmospheric model
configuration file:  '{tmp_path}/test.cfg'
Pressure profile:  1.00e-06 -- 1.00e+02 bar (81 layers)
Wavelength range:  1.10 -- 1.70 um (3209 samples, delta-wn=1.000 cm-1)
Composition:
  ['H2' 'He' 'Na' 'H2O' 'CH4' 'CO' 'CO2']
Opacity sources:
  ['H2O', 'CIA H2-H2', 'CIA H2-He', 'lecavelier', 'deck', 'Na']"""


def test_pyrat_transmission_spec_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
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

Gaussian quadrature cos(theta) angles (quadrature_mu):
    [1.    0.94  0.766 0.5   0.174]
Gaussian quadrature weights (quadrature_weights):
    [0.095 0.691 1.058 0.931 0.367]

Transmission spectrum, (Rp/Rs)**2 (spectrum):
    [ 6.780e-03  6.781e-03  6.780e-03 ...  6.780e-03  6.780e-03  6.780e-03]
"""


def test_pyrat_transmission_atm_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.atm) == """\
Atmospheric model information:
Input atmospheric file name (input_atmfile):
    '{:s}/inputs/atmosphere_uniform_test.atm'
Output atmospheric file name (atmfile): 'None'
Number of layers (nlayers): 81

Pressure display units (punits): bar
Pressure internal units: barye
Pressure at top of atmosphere (ptop):        1.00e-06 bar
Pressure at bottom of atmosphere (pbottom):  1.00e+02 bar
Reference pressure at rplanet (refpressure): 1.00e-01 bar
Pressure profile (press, bar):
    [ 1.000e-06  1.259e-06  1.585e-06 ...  6.310e+01  7.943e+01  1.000e+02]
Atmospheric species information:
Number of species (nmol): 7

Molecule    Mass       Radius
            g/mol      Angstrom
(name)      (mass)     (radius)
  H2          2.0160       1.440
  He          4.0026       1.400
  Na         22.9898       2.200
  H2O        18.0150       1.600
  CH4        16.0430       2.000
  CO         28.0100       1.690
  CO2        44.0090       1.900
Molecular data taken from (molfile):
    '/Users/pato/Dropbox/IWF/projects/2014_pyratbay/pyratbay/pyratbay/data/molecules.dat'

Planetary radius (rplanet, Rjup): 1.000
Planetary mass (mplanet, Mjup): 0.600
Planetary surface gravity (gplanet, cm s-2): 1487.3
Planetary internal temperature (tint, K):  100.0
Planetary Hill radius (rhill, Rjup):  inf
Orbital semi-major axis (smaxis, AU): 0.0450

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
    [  2.3329   2.3329   2.3329 ...   2.3329   2.3329   2.3329]

Abundance units (qunits): vmr
Abundance internal units: mole mixing fraction
Number of atmospheric species: 7
Abundance profiles (vmr, mole mixing fraction):
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


def test_pyrat_transmission_lbl_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.lbl) == f"""\
Line-transition information:
Input TLI files (tlifile):
    ['{os.getcwd()}/outputs/HITRAN_H2O_1.1-1.7um_test.tli']
Number of databases (ndb): 1

Database name (name): HITRAN H2O
Species name (molname):  H2O
Number of isotopes (niso): 7
Number of temperature samples (ntemp): 1001
Temperature (temp, K):
    [1.000e+00 5.000e+00 1.000e+01 ... 4.990e+03 4.995e+03 5.000e+03]
Partition function for each isotope (z):
    [ 1.000e+00  1.010e+00  1.328e+00 ...  8.358e+04  8.387e+04  8.416e+04]
    [ 1.000e+00  1.010e+00  1.332e+00 ...  7.758e+04  7.785e+04  7.811e+04]
    [ 6.000e+00  6.058e+00  7.981e+00 ...  4.655e+05  4.670e+05  4.686e+05]
    [ 6.000e+00  6.213e+00  8.396e+00 ...  4.865e+05  4.880e+05  4.895e+05]
    [ 6.000e+00  6.219e+00  8.445e+00 ...  4.733e+05  4.748e+05  4.763e+05]
    [ 3.600e+01  3.729e+01  5.053e+01 ...  2.737e+06  2.746e+06  2.754e+06]
    [ 6.000e+00  6.343e+00  9.129e+00 ...  9.578e+05  9.615e+05  9.652e+05]

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
Line-transition strength threshold (ethresh): 1.00e-15
Isotopes information:
Number of isotopes (niso): 7

Isotope  Molecule      Mass    Isotopic   Database
            index     g/mol       ratio
 (name)    (imol)    (mass)     (ratio)
    161         3   18.0106   9.973e-01   HITRAN H2O
    181         3   20.0148   2.000e-03   HITRAN H2O
    171         3   19.0148   3.719e-04   HITRAN H2O
    162         3   19.0167   3.107e-04   HITRAN H2O
    182         3   21.0210   6.230e-07   HITRAN H2O
    172         3   20.0210   1.159e-07   HITRAN H2O
    262         3   20.0229   2.420e-08   HITRAN H2O
"""


def test_pyrat_transmission_voigt_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.voigt) == """\
Voigt-profile information:

Number of Doppler-width samples (ndop): 40
Number of Lorentz-width samples (nlor): 40
Doppler HWHM (doppler, cm-1):
    [4.963e-03 5.243e-03 5.538e-03 ... 3.765e-02 3.977e-02 4.201e-02]
Lorentz HWMH (lorentz, cm-1):
    [2.210e-08 3.702e-08 6.202e-08 ... 4.313e+00 7.225e+00 1.210e+01]
Doppler--Lorentz ratio threshold (dlratio): 1.000e-01

Voigt-profiles extent (extent, in HWHMs): 100.0
Voigt-profiles cutoff extent (cutoff in cm-1): 25.0
Voigt-profile half-sizes (size) of shape [ndop, nlor]:
[[ 1072  1132 ...  8590  9074]
 [ 1072  1132 ...  8590  9074]
 ...
 [54000 54000 ... 54000 54000]
 [54000 54000 ... 54000 54000]]
Voigt-profile indices (index) of shape [ndop, nlor]:
[[       0     2145 ...   267162   284343]
 [  302492   304637 ...   569654   586835]
 ...
 [14927817 14927817 ... 14927817 14927817]
 [15035818 15035818 ... 15035818 15035818]]

Voigt profiles:
  profile[ 0, 0]: [2.85914e-08 2.86448e-08 ... 2.86448e-08 2.85914e-08]
  ...
  profile[39,39]: [4.99389e-03 4.99404e-03 ... 4.99404e-03 4.99389e-03]
"""


def test_pyrat_transmission_cs_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.cs) == """\
CIA extinction information:
Number of CIA files (nfiles): 2

CIA file name (files[0]):
    '{:s}/pyratbay/data/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
CIA species: H2-H2
Number of temperature samples: 20
Number of wavenumber samples: 3209
Temperature array (temp, K):
  [  60.  100.  150.  200.  250.  300.  350.  400.  500.  600.  700.  800.  900.
    1000. 2000. 3000. 4000. 5000. 6000. 7000.]
Wavenumber array (wavenumber, cm-1):
  [5882.4 5883.4 5884.4 ... 9088.4 9089.4 9090.4]

CIA file name (files[1]):
    '{:s}/pyratbay/data/CIA/CIA_Borysow_H2He_0050-3000K_0.3-030um.dat'
CIA species: H2-He
Number of temperature samples: 20
Number of wavenumber samples: 3209
Temperature array (temp, K):
  [  50.   75.  100.  150.  200.  250.  300.  350.  400.  500.  600.  700.  800.
    900. 1000. 1250. 1500. 2000. 2500. 3000.]
Wavenumber array (wavenumber, cm-1):
  [5882.4 5883.4 5884.4 ... 9088.4 9089.4 9090.4]

Minimum and maximum temperatures (tmin, tmax) in K: [60.0, 3000.0]
CIA extinction coefficient (ec, cm-1):
[[2.78e-20 2.77e-20 2.77e-20 ... 2.57e-21 2.57e-21 2.57e-21]
 [4.41e-20 4.40e-20 4.38e-20 ... 4.08e-21 4.07e-21 4.07e-21]
 [6.99e-20 6.97e-20 6.95e-20 ... 6.46e-21 6.45e-21 6.45e-21]
 ...
 [1.24e-04 1.24e-04 1.23e-04 ... 5.64e-06 5.64e-06 5.63e-06]
 [1.96e-04 1.96e-04 1.95e-04 ... 8.93e-06 8.92e-06 8.92e-06]
 [3.11e-04 3.10e-04 3.10e-04 ... 1.41e-05 1.41e-05 1.41e-05]]
""".format(os.path.realpath('./..'), os.path.realpath('./..'))


def test_pyrat_transmission_od_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.od) == """\
Optical depth information:
Observing geometry (rt_path): transit
Total atmospheric extinction coefficient (ec, cm-1) [layer, wave]:
[[ 5.613e-17  1.806e-14  5.621e-17 ...  7.134e-16  3.406e-16  3.567e-16]
 [ 7.067e-17  2.274e-14  7.077e-17 ...  8.981e-16  4.288e-16  4.490e-16]
 [ 8.898e-17  2.863e-14  8.911e-17 ...  1.131e-15  5.399e-16  5.653e-16]
 ...
 [ 1.245e-04  1.246e-04  1.245e-04 ...  7.748e-06  7.443e-06  7.152e-06]
 [ 1.970e-04  1.969e-04  1.966e-04 ...  1.162e-05  1.120e-05  1.079e-05]
 [ 3.119e-04  3.116e-04  3.112e-04 ...  1.752e-05  1.699e-05  1.647e-05]]

Distance along the ray path across each layer (outside-in) at each impact
    parameter (raypath, km):
    IP[  1]: [3061.02491505]
    IP[  2]: [1268.96353345 3057.45786888]
    IP[  3]: [ 974.34404583 1267.59004145 3053.53618033]
    ...
    IP[ 80]: [ 164.81155312  165.38824754  165.93695005 ... 1097.04770713
    1426.20058392
 3434.15476801]

Maximum optical depth to calculate (maxdepth): 10.00
Layer index where the optical depth reaches maxdepth (ideep):
    [ 30  30  30  30  30  30  30 ...  30  30  30  30  30  30  30]
Maximum ideep (deepest layer reaching maxdepth): 30

Optical depth at each impact parameter, down to max(ideep) (depth):
[[ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 3.881e-08  1.249e-05  3.887e-08 ...  4.933e-07  2.355e-07  2.466e-07]
 [ 6.490e-08  2.088e-05  6.499e-08 ...  8.248e-07  3.938e-07  4.123e-07]
 ...
 [ 4.696e-05  1.232e-02  4.744e-05 ...  4.890e-04  2.336e-04  2.463e-04]
 [ 6.197e-05  1.550e-02  6.274e-05 ...  6.163e-04  2.944e-04  3.112e-04]
 [ 8.247e-05  1.950e-02  8.365e-05 ...  7.768e-04  3.710e-04  3.931e-04]]
"""


def test_pyrat_transmission_cloud_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.cloud) == """\
Cloud-opacity models (models):

Model name (name): 'deck'
Number of model parameters (npars): 1
Parameter name     Value
  (pnames)         (pars)
  log_p_cl         -3.000e+00
Index of atmospheric layer at or directly below cloud top: 30
Cloud-top pressure: 1.0000e-03 bar
Cloud-top altitude: 72750.08 km
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


def test_pyrat_transmission_rayleigh_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.rayleigh) == """\
Rayleigh-opacity models (models):

Model name (name): 'lecavelier'
Model species (mol): H2
Number of model parameters (npars): 2
Parameter name     Value
  (pnames)         (pars)
  log_k_ray         0.000e+00
  alpha_ray        -4.000e+00
Wavenumber (wn, cm-1):
   [5882.35 5883.35 5884.35 ... 9088.35 9089.35 9090.35]
Cross section (cross_section, cm2 molec-1):
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


def test_pyrat_transmission_alkali_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.alkali) == """\
Alkali-opacity models (models):

Model name (name): 'sodium_vdw'
Model species (mol): Na
Species mass (mass, amu): 22.989769
Profile hard cutoff from line center (cutoff, cm-1): 4500.0
Detuning parameter (detuning): 30.0
Lorentz-width parameter (lpar): 0.071
Partition function (Z): 2.0
Wavenumber  Wavelength          gf   Lower-state energy
      cm-1          um               cm-1
      (wn)                    (gf)   (elow)
  16960.87    0.589592   6.546e-01   0.000e+00
  16978.07    0.588995   1.309e+00   0.000e+00
Wavenumber (wn, cm-1):
   [5882.35 5883.35 5884.35 ... 9088.35 9089.35 9090.35]
Pressure (pressure, barye):
[1.000e+00 1.259e+00 1.585e+00 1.995e+00 2.512e+00 3.162e+00 3.981e+00
 5.012e+00 6.310e+00 7.943e+00 1.000e+01 1.259e+01 1.585e+01 1.995e+01
 2.512e+01 3.162e+01 3.981e+01 5.012e+01 6.310e+01 7.943e+01 1.000e+02
 1.259e+02 1.585e+02 1.995e+02 2.512e+02 3.162e+02 3.981e+02 5.012e+02
 6.310e+02 7.943e+02 1.000e+03 1.259e+03 1.585e+03 1.995e+03 2.512e+03
 3.162e+03 3.981e+03 5.012e+03 6.310e+03 7.943e+03 1.000e+04 1.259e+04
 1.585e+04 1.995e+04 2.512e+04 3.162e+04 3.981e+04 5.012e+04 6.310e+04
 7.943e+04 1.000e+05 1.259e+05 1.585e+05 1.995e+05 2.512e+05 3.162e+05
 3.981e+05 5.012e+05 6.310e+05 7.943e+05 1.000e+06 1.259e+06 1.585e+06
 1.995e+06 2.512e+06 3.162e+06 3.981e+06 5.012e+06 6.310e+06 7.943e+06
 1.000e+07 1.259e+07 1.585e+07 1.995e+07 2.512e+07 3.162e+07 3.981e+07
 5.012e+07 6.310e+07 7.943e+07 1.000e+08]
Cross section (cross_section, cm2 molecule-1):
[[0.000e+00 0.000e+00 ... 0.000e+00 0.000e+00]
 [0.000e+00 0.000e+00 ... 0.000e+00 0.000e+00]
 ...
 [0.000e+00 0.000e+00 ... 0.000e+00 0.000e+00]
 [0.000e+00 0.000e+00 ... 0.000e+00 0.000e+00]]

Total atmospheric alkali extinction-coefficient (ec, cm-1):
[[0.000e+00 0.000e+00 0.000e+00 ... 0.000e+00 0.000e+00 0.000e+00]
 [0.000e+00 0.000e+00 0.000e+00 ... 0.000e+00 0.000e+00 0.000e+00]
 [0.000e+00 0.000e+00 0.000e+00 ... 0.000e+00 0.000e+00 0.000e+00]
 ...
 [0.000e+00 0.000e+00 0.000e+00 ... 0.000e+00 0.000e+00 0.000e+00]
 [0.000e+00 0.000e+00 0.000e+00 ... 0.000e+00 0.000e+00 0.000e+00]
 [0.000e+00 0.000e+00 0.000e+00 ... 0.000e+00 0.000e+00 0.000e+00]]
"""

def test_pyrat_transmission_obs_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.obs) == """\
Observing information:
Data/bandflux display units (units): none
Data/bandflux internal units: none
Number of data points (ndata): 0

Number of filter pass bands (nfilters): 0
"""


def test_pyrat_transmission_phy_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.phy) == """\
Physical properties information:

Stellar effective temperature (tstar, K): 5800.0
Stellar radius (rstar, Rsun): 1.270
Stellar mass (mstar, Msun):   None
Stellar surface gravity (log_gstar, cm s-2): 4.36
Input stellar spectrum is a blackbody at Teff = 5800.0 K.
Stellar spectrum wavenumber (starwn, cm-1):
    [  5882.353   5883.353   5884.353 ...   9088.353   9089.353   9090.353]
Stellar flux spectrum (starflux, erg s-1 cm-2 cm):
    [ 2.306e+06  2.307e+06  2.307e+06 ...  3.293e+06  3.293e+06  3.293e+06]
"""

def test_pyrat_transmission_ret_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.ret) == """\
Retrieval information:
No retrieval parameters set.
"""


def test_pyrat_transmission_resolution_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'resolution':'5000.0'},
        remove=['clouds'],
    )
    pyrat = pb.run(cfg)
    assert pyrat is not None
    assert str(pyrat) == f"""\
Pyrat atmospheric model
configuration file:  '{tmp_path}/test.cfg'
Pressure profile:  1.00e-06 -- 1.00e+02 bar (81 layers)
Wavelength range:  1.10 -- 1.70 um (2177 samples, R=5000.0)
Composition:
  ['H2' 'He' 'Na' 'H2O' 'CH4' 'CO' 'CO2']
Opacity sources:
  ['H2O', 'CIA H2-H2', 'CIA H2-He', 'lecavelier', 'Na']"""

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

Gaussian quadrature cos(theta) angles (quadrature_mu):
    [1.    0.94  0.766 0.5   0.174]
Gaussian quadrature weights (quadrature_weights):
    [0.095 0.691 1.058 0.931 0.367]

Transmission spectrum, (Rp/Rs)**2 (spectrum):
    [ 6.525e-03  6.541e-03  6.525e-03 ...  6.670e-03  6.501e-03  6.473e-03]
"""

@pytest.mark.skip(reason="TBD")
def test_pyrat_transmission_wl_step_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'wlstep':'1e-5 um'},
        remove=['clouds'],
    )
    pyrat = pb.run(cfg)
    assert pyrat is not None
    assert str(pyrat) == """\
Pyrat atmospheric model
configuration file:  '{:s}/test.cfg'
Pressure profile (bar):  1.00e-06 -- 1.00e+02 (81 layers)
Wavelength range (um):  1.10 -- 1.70 (2177 samples, dwl=1e-05 um)
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
    [ 6.522e-03  6.540e-03  6.523e-03 ...  6.670e-03  6.500e-03  6.473e-03]
"""


def test_pyrat_emission_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'rt_path': 'emission'},
    )
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

Gaussian quadrature cos(theta) angles (quadrature_mu):
    [1.    0.94  0.766 0.5   0.174]
Gaussian quadrature weights (quadrature_weights):
    [0.095 0.691 1.058 0.931 0.367]
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
[[ 5.613e-17  1.806e-14  5.621e-17 ...  7.134e-16  3.406e-16  3.567e-16]
 [ 7.067e-17  2.274e-14  7.077e-17 ...  8.981e-16  4.288e-16  4.490e-16]
 [ 8.898e-17  2.863e-14  8.911e-17 ...  1.131e-15  5.399e-16  5.653e-16]
 ...
 [ 1.245e-04  1.246e-04  1.245e-04 ...  7.748e-06  7.443e-06  7.152e-06]
 [ 1.970e-04  1.969e-04  1.966e-04 ...  1.162e-05  1.120e-05  1.079e-05]
 [ 3.119e-04  3.116e-04  3.112e-04 ...  1.752e-05  1.699e-05  1.647e-05]]

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
 [ 3.984e-10  1.282e-07  3.989e-10 ...  5.063e-09  2.417e-09  2.531e-09]
 [ 8.992e-10  2.893e-07  9.004e-10 ...  1.143e-08  5.456e-09  5.713e-09]
 ...
 [ 1.090e-06  3.008e-04  1.099e-06 ...  1.193e-05  5.695e-06  5.995e-06]
 [ 1.422e-06  3.784e-04  1.436e-06 ...  1.502e-05  7.172e-06  7.562e-06]
 [ 1.869e-06  4.760e-04  1.891e-06 ...  1.892e-05  9.036e-06  9.544e-06]]
"""


def test_pyrat_exfile_str(tmp_path):
    reset = {
        'runmode': 'spectrum',
        'specfile': f'{ROOT}tests/outputs/extfile_spectrum_test.dat',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/mcmc_transmission_test.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    assert pyrat is not None
    pyrat.band_integrate()
    assert str(pyrat.lbl) == """\
Line-transition information:
No input TLI files.
"""

    assert str(pyrat.opacity) == """\
Opacity extinction information:
Model           type           T_min   T_max
H2O             line_sample    300.0  3000.0
"""

    assert str(pyrat.obs) == """\
Observing information:
Data/bandflux display units (units): none
Data/bandflux internal units: none
Number of data points (ndata): 20
        Data  Uncertainty   Wavenumber  Wavelength
        none         none         cm-1          um
      (data)     (uncert)     (bandwn)
     0.00661      0.00002      8826.31       1.133
     0.00660      0.00002      8635.38       1.158
     0.00660      0.00002      8449.69       1.183
     0.00651      0.00002      8271.12       1.209
     0.00645      0.00002      8097.00       1.235
     0.00641      0.00002      7936.66       1.260
     0.00647      0.00002      7781.94       1.285
     0.00648      0.00002      7630.82       1.310
     0.00666      0.00002      7484.88       1.336
     0.00673      0.00002      7344.98       1.361
     0.00677      0.00002      7207.07       1.388
     0.00674      0.00002      7077.26       1.413
     0.00676      0.00002      6951.56       1.439
     0.00670      0.00002      6830.71       1.464
     0.00667      0.00002      6715.80       1.489
     0.00658      0.00002      6600.55       1.515
     0.00656      0.00002      6493.61       1.540
     0.00646      0.00002      6387.63       1.566
     0.00650      0.00002      6285.45       1.591
     0.00649      0.00002      6188.02       1.616

Number of filter pass bands (nfilters): 20
Wavenumber  Wavelength    Bandflux  Filter name
      cm-1          um        none
  (bandwn)              (bandflux)  (filters)
   8826.31       1.133     0.00657  filter_test_WFC3_G141_1.133um
   8635.38       1.158     0.00655  filter_test_WFC3_G141_1.158um
   8449.69       1.183     0.00654  filter_test_WFC3_G141_1.183um
   8271.12       1.209     0.00648  filter_test_WFC3_G141_1.209um
   8097.00       1.235     0.00646  filter_test_WFC3_G141_1.235um
   7936.66       1.260     0.00645  filter_test_WFC3_G141_1.260um
   7781.94       1.285     0.00649  filter_test_WFC3_G141_1.285um
   7630.82       1.310     0.00651  filter_test_WFC3_G141_1.310um
   7484.88       1.336     0.00664  filter_test_WFC3_G141_1.336um
   7344.98       1.361     0.00670  filter_test_WFC3_G141_1.361um
   7207.07       1.388     0.00670  filter_test_WFC3_G141_1.387um
   7077.26       1.413     0.00671  filter_test_WFC3_G141_1.413um
   6951.56       1.439     0.00671  filter_test_WFC3_G141_1.438um
   6830.71       1.464     0.00667  filter_test_WFC3_G141_1.464um
   6715.80       1.489     0.00662  filter_test_WFC3_G141_1.489um
   6600.55       1.515     0.00657  filter_test_WFC3_G141_1.515um
   6493.61       1.540     0.00655  filter_test_WFC3_G141_1.540um
   6387.63       1.566     0.00651  filter_test_WFC3_G141_1.565um
   6285.45       1.591     0.00651  filter_test_WFC3_G141_1.591um
   6188.02       1.616     0.00651  filter_test_WFC3_G141_1.616um
"""

    assert str(pyrat.ret) == """\
Retrieval information:
  Parameter name        value        pmin        pmax       pstep
  (pnames)           (params)      (pmin)      (pmax)     (pstep)
  log_kappa'       -5.000e+00  -9.000e+00   5.000e+00   3.000e-01
  log_gamma1        0.000e+00  -3.000e+00   3.000e+00   3.000e-01
  log_gamma2        0.000e+00  -3.000e+00   3.000e+00   0.000e+00
  alpha             0.000e+00   0.000e+00   1.000e+00   0.000e+00
  T_irr             1.486e+03   0.000e+00   7.000e+03   5.000e+01
  T_int             1.000e+02   0.000e+00   5.000e+02   0.000e+00
  R_planet          1.020e+00   5.000e-01   4.500e+00   3.000e-02
  log_H2O          -4.000e+00  -9.000e+00  -1.000e+00   5.000e-01

Parameter name     Prior
  log_kappa'       Uniform between     [-9.000e+00,  5.000e+00]
  log_gamma1       Uniform between     [-3.000e+00,  3.000e+00]
  log_gamma2       Fixed at   0.000e+00
  alpha            Fixed at   0.000e+00
  T_irr            Uniform between     [ 0.000e+00,  7.000e+03]
  T_int            Fixed at   1.000e+02
  R_planet         Uniform between     [ 5.000e-01,  4.500e+00]
  log_H2O          Uniform between     [-9.000e+00, -1.000e+00]

Retrieval algorithm (sampler): snooker
Number of retrieval samples (nsamples): 300
Number of parallel chains (nchains):   21
Number of burned-in samples (burnin):  10
Thinning factor (thinning): 1

Upper boundary for sum of metal abundances (qcap): None
Temperature upper boundary (tlow, K):   300.0
Temperature lower boundary (thigh, K): 3000.0

Retrieval posterior file (mcmcfile):
    {:s}/outputs/MCMC_transmission_test.npz
""".format(os.getcwd())

