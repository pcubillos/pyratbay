# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import pytest

import numpy as np
from conftest import make_config

import pyratbay as pb
import pyratbay.atmosphere as pa
import pyratbay.constants as pc
import pyratbay.opacity as op
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')
extfile = '{ROOT}/tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz'


def test_pyrat_transmission_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'extfile': f'{extfile}'},
        remove=['tlifile'],
    )
    pyrat = pb.run(cfg)
    assert str(pyrat) == f"""\
Pyrat atmospheric model
configuration file:  '{tmp_path}/test.cfg'
Pressure profile:  1.00e-06 -- 1.00e+02 bar (51 layers)
Wavelength range:  1.10 -- 1.70 um (3209 samples, delta-wn=1.000 cm-1)
Composition:
  ['H2' 'He' 'H' 'Na' 'K' 'H2O' 'CH4' 'CO' 'CO2']
Opacity sources:
  ['H2O', 'Na', 'CIA H2-H2', 'CIA H2-He', 'lecavelier', 'deck']"""


def test_pyrat_opacity_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'extfile': f'{extfile}'},
        remove=['tlifile'],
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.opacity) == f"""\
Opacity extinction information:
Model           type           T_min   T_max
H2O             line_sample    300.0  3000.0
sodium_vdw      alkali
CIA H2-H2       cia             60.0  3000.0
CIA H2-He       cia             60.0  3000.0
lecavelier      cloud
deck            cloud
"""


def test_opacity_alkali_str():
    pressure = pa.pressure('1e-6 bar', '1e2 bar', nlayers=51)
    temperature = pa.tmodels.Isothermal(pressure)(1500.0)
    wn_min = 1.0 / (0.80252 * pc.um)
    wn_max = 1.0 / (0.466 * pc.um)
    wn = np.arange(wn_min, wn_max, 1.0)

    model = op.alkali.SodiumVdW(pressure, wn=wn, cutoff=4500.0)
    model.calc_cross_section(temperature)
    assert str(model) == """\
Model name (name): 'sodium_vdw'
Model species (species): Na
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
Pressure (pressure, bar):
[1.000e-06 1.445e-06 2.089e-06 3.020e-06 4.365e-06 6.310e-06 9.120e-06
 1.318e-05 1.905e-05 2.754e-05 3.981e-05 5.754e-05 8.318e-05 1.202e-04
 1.738e-04 2.512e-04 3.631e-04 5.248e-04 7.586e-04 1.096e-03 1.585e-03
 2.291e-03 3.311e-03 4.786e-03 6.918e-03 1.000e-02 1.445e-02 2.089e-02
 3.020e-02 4.365e-02 6.310e-02 9.120e-02 1.318e-01 1.905e-01 2.754e-01
 3.981e-01 5.754e-01 8.318e-01 1.202e+00 1.738e+00 2.512e+00 3.631e+00
 5.248e+00 7.586e+00 1.096e+01 1.585e+01 2.291e+01 3.311e+01 4.786e+01
 6.918e+01 1.000e+02]
Cross section (cross_section, cm2 molecule-1):
[[0.000e+00 4.874e-29 ... 1.488e-28 1.486e-28]
 [0.000e+00 7.045e-29 ... 2.151e-28 2.149e-28]
 ...
 [0.000e+00 3.337e-21 ... 1.019e-20 1.018e-20]
 [0.000e+00 4.769e-21 ... 1.456e-20 1.454e-20]]
"""


def test_opacity_cloud_deck_str():
    pressure = pa.pressure('1e-6 bar', '1e2 bar', nlayers=51)
    wn_min = 1.0 / (1.7 * pc.um)
    wn_max = 1.0 / (1.1 * pc.um)
    wn = np.arange(wn_min, wn_max, 1.0)
    model = op.clouds.Deck(pressure, wn)

    pars = [-3.0]
    tmodel = pa.tmodels.Guillot(pressure)
    temperature = tmodel([-4.83, -0.8, 0, 0, 1200, 100])
    p0 = 0.1  # bar
    mu = 2.3
    radius = pa.hydro_m(
        pressure, temperature, mu, 1.0*pc.mjup, p0, 1.0*pc.rjup,
    )
    model.calc_extinction_coefficient(radius, temperature, pars)
    assert str(model) == """\
Model name (name): 'deck'
Number of model parameters (npars): 1
Parameter name     Value
  (pnames)         (pars)
  log_p_cl         -3.000e+00
Index of atmospheric layer at or directly below cloud top: 19
Cloud-top pressure: 1.0000e-03 bar
Cloud-top altitude: 72253.11 km
Cloud-top temperature: 1051.32 K
"""


def test_opacity_cloud_ccsgray_str():
    pressure = pa.pressure('1e-6 bar', '1e2 bar', nlayers=51)
    wn_min = 1.0 / (1.7 * pc.um)
    wn_max = 1.0 / (1.1 * pc.um)
    wn = np.arange(wn_min, wn_max, 1.0)

    model = op.clouds.CCSgray(pressure, wn)
    pars = [0.0, -3.0, 2.0]
    tmodel = pa.tmodels.Guillot(pressure)
    temperature = tmodel([-4.83, -0.8, 0, 0, 1200, 100])
    model.calc_extinction_coefficient(temperature, pars)
    assert str(model) == """\
Model name (name): 'ccsgray'
Number of model parameters (npars): 3
Parameter name     Value
  (pnames)         (pars)
  log_k_gray        0.000e+00
  log_p_top        -3.000e+00
  log_p_bot         2.000e+00
Extinction coefficient (ec, cm-1):
[[0.000e+00 0.000e+00 0.000e+00 ... 0.000e+00 0.000e+00 0.000e+00]
 [0.000e+00 0.000e+00 0.000e+00 ... 0.000e+00 0.000e+00 0.000e+00]
 [0.000e+00 0.000e+00 0.000e+00 ... 0.000e+00 0.000e+00 0.000e+00]
 ...
 [1.108e-06 1.108e-06 1.108e-06 ... 1.108e-06 1.108e-06 1.108e-06]
 [1.600e-06 1.600e-06 1.600e-06 ... 1.600e-06 1.600e-06 1.600e-06]
 [2.310e-06 2.310e-06 2.310e-06 ... 2.310e-06 2.310e-06 2.310e-06]]
"""


def test_opacity_rayleigh_lecavelier_str():
    wn_min = 1.0 / (1.7 * pc.um)
    wn_max = 1.0 / (1.1 * pc.um)
    wn = np.arange(wn_min, wn_max, 1.0)
    model = op.rayleigh.Lecavelier(wn)
    assert str(model) == """\
Model name (name): 'lecavelier'
Model species (species): H2
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


def test_opacity_rayleigh_dalgarno_H_str():
    wn_min = 1.0 / (1.7 * pc.um)
    wn_max = 1.0 / (1.1 * pc.um)
    wn = np.arange(wn_min, wn_max, 1.0)
    model = op.rayleigh.Dalgarno(wn, species='H')
    assert str(model) == """\
Model name (name): 'dalgarno_H'
Model species (species): H
Number of model parameters (npars): 0
Wavenumber (wn, cm-1):
   [5882.35 5883.35 5884.35 ... 9088.35 9089.35 9090.35]
Cross section (cross_section, cm2 molec-1):
   [7.002e-30 7.007e-30 7.012e-30 ... 4.038e-29 4.040e-29 4.041e-29]
"""

def test_opacity_rayleigh_dalgarno_He_str():
    wn_min = 1.0 / (1.7 * pc.um)
    wn_max = 1.0 / (1.1 * pc.um)
    wn = np.arange(wn_min, wn_max, 1.0)
    model = op.rayleigh.Dalgarno(wn, species='He')
    assert str(model) == """\
Model name (name): 'dalgarno_He'
Model species (species): He
Number of model parameters (npars): 0
Wavenumber (wn, cm-1):
   [5882.35 5883.35 5884.35 ... 9088.35 9089.35 9090.35]
Cross section (cross_section, cm2 molec-1):
   [6.577e-31 6.582e-31 6.586e-31 ... 3.757e-30 3.758e-30 3.760e-30]
"""

def test_opacity_rayleigh_dalgarno_H2_str():
    wn_min = 1.0 / (1.7 * pc.um)
    wn_max = 1.0 / (1.1 * pc.um)
    wn = np.arange(wn_min, wn_max, 1.0)
    model = op.rayleigh.Dalgarno(wn, species='H2')
    assert str(model) == """\
Model name (name): 'dalgarno_H2'
Model species (species): H2
Number of model parameters (npars): 0
Wavenumber (wn, cm-1):
   [5882.35 5883.35 5884.35 ... 9088.35 9089.35 9090.35]
Cross section (cross_section, cm2 molec-1):
   [9.799e-30 9.806e-30 9.813e-30 ... 5.626e-29 5.629e-29 5.631e-29]
"""


def test_opacity_cia_borysow_H2H2_str():
    wn_min = 1.0 / (1.7 * pc.um)
    wn_max = 1.0 / (1.1 * pc.um)
    wn = np.arange(wn_min, wn_max, 1.0)
    cia_file = f'{pc.ROOT}pyratbay/data/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
    model = op.Collision_Induced(cia_file, wn=wn)
    assert str(model) == f"""\
CIA file name (cia_file):
    '{ROOT}pyratbay/data/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
CIA species (species): ['H2', 'H2']
Number of temperature samples (ntemp): 20
Number of wavenumber samples (nwave): 3209
Temperature array (temps, K):
[  60.  100.  150.  200.  250.  300.  350.  400.  500.  600.  700.  800.
  900. 1000. 2000. 3000. 4000. 5000. 6000. 7000.]
Wavenumber array (wn, cm-1):
[5882.353 5883.353 5884.353 ... 9088.353 9089.353 9090.353]
Wavelength array (um):
[1.70000 1.69971 1.69942 ... 1.10031 1.10019 1.10007]
Tabulated cross section (tab_cross_section, cm5 molec-2):
[[6.79e-49 6.71e-49 6.64e-49 ... 5.35e-48 5.29e-48 5.22e-48]
 [2.59e-48 2.58e-48 2.57e-48 ... 8.53e-48 8.45e-48 8.37e-48]
 [5.35e-48 5.31e-48 5.27e-48 ... 1.35e-47 1.34e-47 1.33e-47]
 ...
 [1.22e-44 1.22e-44 1.22e-44 ... 3.76e-46 3.76e-46 3.76e-46]
 [1.62e-44 1.61e-44 1.61e-44 ... 5.48e-46 5.47e-46 5.47e-46]
 [1.89e-44 1.89e-44 1.89e-44 ... 7.42e-46 7.41e-46 7.40e-46]]
Minimum and maximum temperatures (tmin, tmax) in K: [60.0, 7000.0]
"""

def test_opacity_cia_borysow_H2He_str():
    wn_min = 1.0 / (1.7 * pc.um)
    wn_max = 1.0 / (1.1 * pc.um)
    wn = np.arange(wn_min, wn_max, 1.0)
    cia_file = f'{pc.ROOT}pyratbay/data/CIA/CIA_Borysow_H2He_0050-3000K_0.3-030um.dat'
    model = op.Collision_Induced(cia_file, wn=wn)
    assert str(model) == f"""\
CIA file name (cia_file):
    '{ROOT}pyratbay/data/CIA/CIA_Borysow_H2He_0050-3000K_0.3-030um.dat'
CIA species (species): ['H2', 'He']
Number of temperature samples (ntemp): 20
Number of wavenumber samples (nwave): 3209
Temperature array (temps, K):
[  50.   75.  100.  150.  200.  250.  300.  350.  400.  500.  600.  700.
  800.  900. 1000. 1250. 1500. 2000. 2500. 3000.]
Wavenumber array (wn, cm-1):
[5882.353 5883.353 5884.353 ... 9088.353 9089.353 9090.353]
Wavelength array (um):
[1.70000 1.69971 1.69942 ... 1.10031 1.10019 1.10007]
Tabulated cross section (tab_cross_section, cm5 molec-2):
[[1.23e-49 1.22e-49 1.22e-49 ... 8.34e-50 8.29e-50 8.24e-50]
 [2.83e-49 2.81e-49 2.80e-49 ... 1.35e-49 1.35e-49 1.34e-49]
 [5.78e-49 5.75e-49 5.72e-49 ... 2.06e-49 2.05e-49 2.04e-49]
 ...
 [1.28e-45 1.28e-45 1.27e-45 ... 3.93e-47 3.92e-47 3.91e-47]
 [1.99e-45 1.99e-45 1.98e-45 ... 5.51e-47 5.50e-47 5.49e-47]
 [2.85e-45 2.85e-45 2.84e-45 ... 7.41e-47 7.40e-47 7.38e-47]]
Minimum and maximum temperatures (tmin, tmax) in K: [50.0, 3000.0]
"""


def test_opacity_lbl_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)
    lbl_index = pyrat.opacity.models_type.index('lbl')
    lbl = pyrat.opacity.models[lbl_index]

    assert str(lbl) == f"""\
Line-transition information:
Input TLI files (tlifile):
    ['{os.getcwd()}/outputs/HITRAN_H2O_1.1-1.7um_test.tli']
Number of databases (ndb): 1

Database name (name): HITRAN H2O
Species name (molname):  H2O
Number of isotopes (niso): 9
Number of temperature samples (ntemp): 1201
Temperature (temp, K):
    [1.000e+00 5.000e+00 1.000e+01 ... 5.990e+03 5.995e+03 6.000e+03]
Partition function for each isotope (z):
    [ 1.000e+00  1.010e+00  1.328e+00 ...  1.564e+05  1.568e+05  1.573e+05]
    [ 1.000e+00  1.010e+00  1.332e+00 ...  1.450e+05  1.454e+05  1.459e+05]
    [ 6.000e+00  6.058e+00  7.981e+00 ...  8.709e+05  8.735e+05  8.760e+05]
    [ 6.000e+00  6.213e+00  8.396e+00 ...  8.406e+05  8.426e+05  8.446e+05]
    [ 6.000e+00  6.219e+00  8.445e+00 ...  8.538e+05  8.561e+05  8.584e+05]
    [ 3.600e+01  3.729e+01  5.053e+01 ...  4.911e+06  4.924e+06  4.938e+06]
    [ 6.000e+00  6.343e+00  9.129e+00 ...  1.949e+06  1.955e+06  1.962e+06]
    [ 6.000e+00  6.353e+00  9.217e+00 ...  2.006e+06  2.013e+06  2.019e+06]
    [ 3.600e+01  3.809e+01  5.505e+01 ...  1.187e+07  1.191e+07  1.195e+07]

Total number of line transitions (ntransitions): 47,658
Minimum and maximum temperatures (tmin, tmax): [1.0, 6000.0] K
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
Number of isotopes (niso): 9

Isotope  Molecule      Mass    Isotopic   Database
            index     g/mol       ratio
 (name)    (imol)    (mass)     (ratio)
    161         5   18.0106   9.973e-01   HITRAN H2O
    181         5   20.0148   2.000e-03   HITRAN H2O
    171         5   19.0148   3.719e-04   HITRAN H2O
    162         5   19.0167   3.107e-04   HITRAN H2O
    182         5   21.0210   6.230e-07   HITRAN H2O
    172         5   20.0210   1.159e-07   HITRAN H2O
    262         5   20.0229   2.420e-08   HITRAN H2O
    282         5   22.0274   4.500e-09   HITRAN H2O
    272         5   21.0273   8.600e-10   HITRAN H2O
"""

    assert str(pyrat.voigt) == """\
Voigt-profile information:

Number of Doppler-width samples (ndop): 50
Number of Lorentz-width samples (nlor): 100
Doppler HWHM (doppler, cm-1):
    [4.963e-03 5.184e-03 5.415e-03 ... 3.850e-02 4.022e-02 4.201e-02]
Lorentz HWMH (lorentz, cm-1):
    [2.210e-08 2.708e-08 3.318e-08 ... 8.061e+00 9.878e+00 1.210e+01]
Doppler--Lorentz ratio threshold (dlratio): 1.000e-01

Voigt-profiles extent (extent, in HWHMs): 100.0
Voigt-profiles cutoff extent (cutoff in cm-1): 25.0
Voigt-profile half-sizes (size) of shape [nlor, ndop]:
[[ 1072  1120 ...  8687  9074]
 [ 1072  1120 ...  8687  9074]
 ...
 [54000 54000 ... 54000 54000]
 [54000 54000 ... 54000 54000]]
Voigt-profile indices (index) of shape [nlor, ndop]:
[[       0     2145 ...   341896   359271]
 [  377420   379565 ...   719318   736693]
 ...
 [46933096 46933096 ... 46933096 46933096]
 [47041097 47041097 ... 47041097 47041097]]

Voigt profiles:
  profile[ 0, 0]: [2.85914e-08 2.86448e-08 ... 2.86448e-08 2.85914e-08]
  ...
  profile[99,49]: [4.99389e-03 4.99404e-03 ... 4.99404e-03 4.99389e-03]
"""


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


def test_atm_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'extfile': f'{extfile}'},
        remove=['tlifile'],
    )
    pyrat = pb.run(cfg)
    print(pyrat.atm)
    assert str(pyrat.atm) == f"""\
Atmospheric model information:
Input atmospheric file name (atmfile):
    '{os.getcwd()}/inputs/atmosphere_uniform_test.atm'
Output atmospheric file name (output_atmfile): 'None'
Number of layers (nlayers): 51

Planetary radius (rplanet, Rjup): 1.000
Planetary mass (mplanet, Mjup): 0.600
Planetary surface gravity (gplanet, cm s-2): 1487.3
Planetary internal temperature (tint, K):  100.0
Planetary Hill radius (rhill, Rjup):  inf
Orbital semi-major axis (smaxis, AU): 0.0450

Pressure display units (punits): bar
Pressure internal units: bar
Pressure at top of atmosphere (ptop):        1.00e-06 bar
Pressure at bottom of atmosphere (pbottom):  1.00e+02 bar
Reference pressure at rplanet (refpressure): 1.00e-01 bar
Pressure profile (press, bar):
    [1.000e-06 1.445e-06 2.089e-06 ... 4.786e+01 6.918e+01 1.000e+02]

Temperature units (tunits): kelvin
Temperature model name (tmodelname): None
Temperature profile (temp, K):
    [ 1046.894  1046.896  1046.899 ...  1662.094  1663.380  1665.234]

Abundance units (qunits): vmr
Abundance internal units: VMR
Abundance model (chemistry): None
Number of species (nmol): 9

Index   Molecule  Mass (g/mol)    Radius (A)
       (species)    (mol_mass)  (mol_radius)
    0         H2       2.01600         1.440
    1         He       4.00260         1.400
    2          H       1.00800         1.100
    3         Na      22.98977         2.200
    4          K      39.09830         2.800
    5        H2O      18.01500         1.600
    6        CH4      16.04300         2.000
    7         CO      28.01000         1.690
    8        CO2      44.00900         1.900
Molecular data taken from (molfile):
    '{ROOT}pyratbay/data/molecules.dat'
Abundance profiles (vmr):
              H2:   [8.500e-01 8.500e-01 ... 8.500e-01 8.500e-01]
              He:   [1.490e-01 1.490e-01 ... 1.490e-01 1.490e-01]
               H:   [1.000e-06 1.000e-06 ... 1.000e-06 1.000e-06]
              Na:   [3.000e-06 3.000e-06 ... 3.000e-06 3.000e-06]
               K:   [5.000e-08 5.000e-08 ... 5.000e-08 5.000e-08]
             H2O:   [4.000e-04 4.000e-04 ... 4.000e-04 4.000e-04]
             CH4:   [1.000e-04 1.000e-04 ... 1.000e-04 1.000e-04]
              CO:   [5.000e-04 5.000e-04 ... 5.000e-04 5.000e-04]
             CO2:   [1.000e-07 1.000e-07 ... 1.000e-07 1.000e-07]
Density profiles (d, molecules cm-3):
              H2:   [5.881e+12 8.500e+12 ... 2.561e+20 3.697e+20]
              He:   [1.031e+12 1.490e+12 ... 4.489e+19 6.481e+19]
               H:   [6.919e+06 1.000e+07 ... 3.012e+14 4.350e+14]
              Na:   [2.076e+07 3.000e+07 ... 9.037e+14 1.305e+15]
               K:   [3.459e+05 5.000e+05 ... 1.506e+13 2.175e+13]
             H2O:   [2.767e+09 4.000e+09 ... 1.205e+17 1.740e+17]
             CH4:   [6.919e+08 1.000e+09 ... 3.012e+16 4.350e+16]
              CO:   [3.459e+09 5.000e+09 ... 1.506e+17 2.175e+17]
             CO2:   [6.919e+05 1.000e+06 ... 3.012e+13 4.350e+13]

Radius display units (runits): rjup
Radius internal units: cm
Radius model name (rmodelname): hydro_m
Radius profile (radius, rjup):
    [1.0433 1.0419 1.0405 ... 0.9679 0.966  0.9641]

Mean molecular mass (mm, amu):
    [  2.3329   2.3329   2.3329 ...   2.3329   2.3329   2.3329]
"""


def test_pyrat_transmission_od_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'extfile': f'{extfile}'},
        remove=['tlifile'],
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.od) == """\
Optical depth information:
Observing geometry (rt_path): transit
Total atmospheric extinction coefficient (ec, cm-1) [layer, wave]:
[[ 5.614e-17  1.705e-14  5.621e-17 ...  8.241e-16  3.412e-16  3.554e-16]
 [ 8.116e-17  2.465e-14  8.127e-17 ...  1.191e-15  4.932e-16  5.138e-16]
 [ 1.174e-16  3.563e-14  1.175e-16 ...  1.722e-15  7.129e-16  7.427e-16]
 ...
 [ 7.173e-05  7.197e-05  7.190e-05 ...  4.776e-06  4.599e-06  4.419e-06]
 [ 1.496e-04  1.496e-04  1.494e-04 ...  9.100e-06  8.734e-06  8.387e-06]
 [ 3.119e-04  3.116e-04  3.112e-04 ...  1.746e-05  1.695e-05  1.645e-05]]

Distance along the ray path across each layer (outside-in) at each impact
    parameter (raypath, km):
    IP[  1]: [3870.14]
    IP[  2]: [1605.36 3862.34]
    IP[  3]: [1233.12 1602.11 3854.56]
    ...
    IP[ 50]: [ 263.89  265.3   266.76 ... 1392.69 1807.63 4345.38]

Maximum optical depth to calculate (maxdepth): 10.00
Layer index where the optical depth reaches maxdepth (ideep):
    [ 19  19  19  19  19  19  19 ...  19  19  19  19  19  19  19]
Maximum ideep (deepest layer reaching maxdepth): 19

Optical depth at each impact parameter, down to max(ideep) (depth):
[[ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 5.313e-08  1.614e-05  5.321e-08 ...  7.800e-07  3.229e-07  3.364e-07]
 [ 9.871e-08  2.998e-05  9.885e-08 ...  1.449e-06  5.998e-07  6.248e-07]
 ...
 [ 3.741e-05  9.556e-03  3.777e-05 ...  4.654e-04  1.922e-04  2.016e-04]
 [ 5.784e-05  1.380e-02  5.855e-05 ...  6.740e-04  2.781e-04  2.925e-04]
 [ 9.147e-05  1.992e-02  9.288e-05 ...  9.775e-04  4.027e-04  4.252e-04]]
"""


def test_pyrat_transmission_obs_str(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'extfile': f'{extfile}'},
        remove=['tlifile'],
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
        reset={'extfile': f'{extfile}'},
        remove=['tlifile'],
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
        reset={'extfile': f'{extfile}'},
        remove=['tlifile'],
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
    assert str(pyrat) == f"""\
Pyrat atmospheric model
configuration file:  '{tmp_path}/test.cfg'
Pressure profile:  1.00e-06 -- 1.00e+02 bar (51 layers)
Wavelength range:  1.10 -- 1.70 um (2177 samples, R=5000.0)
Composition:
  ['H2' 'He' 'H' 'Na' 'K' 'H2O' 'CH4' 'CO' 'CO2']
Opacity sources:
  ['H2O', 'Na', 'CIA H2-H2', 'CIA H2-He', 'lecavelier']"""

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
    [ 6.523e-03  6.540e-03  6.524e-03 ...  6.669e-03  6.500e-03  6.473e-03]
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
Composition:  ['H2' 'He' 'H' 'Na' 'H2O' 'CH4' 'CO' 'CO2']
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
        reset={
            'rt_path': 'emission',
            'extfile': f'{extfile}',
        },
        remove=['tlifile'],
    )
    pyrat = pb.run(cfg)
    assert str(pyrat.spec) == """\
Spectral information:
Wavenumber internal units: cm-1
Wavelength internal units: cm
Wavelength display units (wlunits): um
Low wavenumber boundary (wnlow):     5882.353 cm-1  (wlhigh =   1.70 um)
High wavenumber boundary (wnhigh):   9090.353 cm-1  (wllow  =   1.10 um)
Number of samples (nwave): 3209
Sampling interval (wnstep): 1.000 cm-1
Wavenumber array (wn, cm-1):
    [ 5882.353  5883.353  5884.353 ...  9088.353  9089.353  9090.353]
Oversampling factor (wnosamp): None

Gaussian quadrature cos(theta) angles (quadrature_mu):
    [1.    0.94  0.766 0.5   0.174]
Gaussian quadrature weights (quadrature_weights):
    [0.095 0.691 1.058 0.931 0.367]
Intensity spectra (intensity, erg s-1 cm-2 sr-1 cm):
    [ 7.737e+02  7.730e+02  7.723e+02 ...  3.546e+01  3.542e+01  3.539e+01]
    [ 7.737e+02  7.730e+02  7.723e+02 ...  3.546e+01  3.542e+01  3.539e+01]
    [ 7.737e+02  7.730e+02  7.723e+02 ...  3.546e+01  3.542e+01  3.539e+01]
    [ 7.737e+02  7.730e+02  7.723e+02 ...  3.546e+01  3.542e+01  3.539e+01]
    [ 7.737e+02  7.730e+02  7.723e+02 ...  3.546e+01  3.542e+01  3.539e+01]
Emission spectrum (spectrum, erg s-1 cm-2 cm):
    [ 2.431e+03  2.428e+03  2.426e+03 ...  1.114e+02  1.113e+02  1.112e+02]
"""

    assert str(pyrat.od) == """\
Optical depth information:
Observing geometry (rt_path): emission
Total atmospheric extinction coefficient (ec, cm-1) [layer, wave]:
[[ 5.614e-17  1.705e-14  5.621e-17 ...  8.241e-16  3.412e-16  3.554e-16]
 [ 8.116e-17  2.465e-14  8.127e-17 ...  1.191e-15  4.932e-16  5.138e-16]
 [ 1.174e-16  3.563e-14  1.175e-16 ...  1.722e-15  7.129e-16  7.427e-16]
 ...
 [ 7.173e-05  7.197e-05  7.190e-05 ...  4.776e-06  4.599e-06  4.419e-06]
 [ 1.496e-04  1.496e-04  1.494e-04 ...  9.100e-06  8.734e-06  8.387e-06]
 [ 3.119e-04  3.116e-04  3.112e-04 ...  1.746e-05  1.695e-05  1.645e-05]]

Distance across each layer along a normal ray path (raypath, km):
    [100.5 100.2 99.9 99.7 ... 138.2 137.7 137.3 136.8]

Maximum optical depth to calculate (maxdepth): 10.00
Layer index where the optical depth reaches maxdepth (ideep):
    [ 19  19  19  19  19  19  19 ...  19  19  19  19  19  19  19]
Maximum ideep (deepest layer reaching maxdepth): 19

Planck emission down to max(ideep) (B, erg s-1 cm-2 sr-1 cm):
[[ 7.478e+02  7.471e+02  7.465e+02 ...  3.364e+01  3.361e+01  3.357e+01]
 [ 7.478e+02  7.471e+02  7.465e+02 ...  3.364e+01  3.361e+01  3.357e+01]
 [ 7.478e+02  7.472e+02  7.465e+02 ...  3.365e+01  3.361e+01  3.358e+01]
 ...
 [ 7.612e+02  7.605e+02  7.599e+02 ...  3.458e+01  3.454e+01  3.451e+01]
 [ 7.673e+02  7.666e+02  7.660e+02 ...  3.501e+01  3.497e+01  3.494e+01]
 [ 7.737e+02  7.730e+02  7.723e+02 ...  3.546e+01  3.542e+01  3.539e+01]]

Optical depth at each layer along a normal ray path into the planet, down to
    max(ideep) (depth):
[[ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 6.897e-10  2.095e-07  6.907e-10 ...  1.012e-08  4.192e-09  4.367e-09]
 [ 1.684e-09  5.115e-07  1.687e-09 ...  2.472e-08  1.023e-08  1.066e-08]
 ...
 [ 8.937e-07  2.378e-04  9.008e-07 ...  1.156e-05  4.777e-06  5.002e-06]
 [ 1.362e-06  3.433e-04  1.375e-06 ...  1.673e-05  6.907e-06  7.248e-06]
 [ 2.114e-06  4.956e-04  2.142e-06 ...  2.423e-05  9.995e-06  1.052e-05]]
"""


def test_pyrat_exfile_str(tmp_path):
    reset = {
        'runmode': 'spectrum',
        'specfile': f'{ROOT}tests/outputs/extfile_spectrum_test.dat',
        'extfile': f'{extfile}',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/mcmc_transmission_test.cfg',
        reset=reset,
        remove=['tlifile'],
    )
    pyrat = pb.run(cfg)
    assert pyrat is not None
    pyrat.band_integrate()

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
   6951.56       1.439     0.00670  filter_test_WFC3_G141_1.438um
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

Retrieval posterior file (mcmcfile): None
""".format(os.getcwd())

