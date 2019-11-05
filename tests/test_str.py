import os

from conftest import make_config

import pyratbay as pb
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


def test_alkali_str(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg')
    pyrat = pb.run(cfg)
    assert pyrat is not None
    print(pyrat.alkali.models[0])
    assert str(pyrat.alkali.models[0]) == """\
Model name (name): 'sodium_vdw'
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
Extinction-coefficient (ec, cm2 molecule-1):
[[ 9.230e-34  9.244e-34 ...  1.263e-31  1.265e-31]
 [ 1.162e-33  1.164e-33 ...  1.590e-31  1.593e-31]
 ...
 [ 1.299e-23  1.300e-23 ...  3.468e-22  3.472e-22]
 [ 1.636e-23  1.637e-23 ...  4.358e-22  4.363e-22]]
"""

def test_cloud_str(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'clouds':'deck ccsgray', 'cpars':'-3.0  0.0 -4.0 2.0'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    assert str(pyrat.cloud.models[0]) == """\
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
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
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
  log(f_Ray)        0.000e+00
  alpha_Ray        -4.000e+00
Extinction-coefficient (ec, cm2 molec-1):
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
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg')
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
    [ 6.772e-03  6.773e-03  6.772e-03 ...  6.772e-03  6.772e-03  6.772e-03]
"""

    assert str(pyrat.atm) == """\
Atmospheric model information:
Atmospheric file name (atmfile):
    '{:s}/atmosphere_uniform_test.atm'
Number of layers (nlayers): 81

Pressure display units (punits): bar
Pressure internal units: barye
Pressure at top of atmosphere (ptop):        1.00e-06 bar
Pressure at bottom of atmosphere (pbottom):  1.00e+02 bar
Reference pressure at rplanet (refpressure): 1.00e-01 bar
Pressure profile (press, bar):
    [ 1.000e-06  1.259e-06  1.585e-06 ...  6.310e+01  7.943e+01  1.000e+02]

Radius display units (runits): km
Radius internal units: cm
Radius model name (rmodelname): hydro_m
Radius profile (radius, km):
    [74591.6311 74528.7914 74466.0453 ... 69099.5141 69013.8566 68928.3555]

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

    assert str(pyrat.mol) == """\
Atmospheric species information:
Number of species (nmol): 7

Molecule    ID   Mass      Radius
                 g/mol     Angstrom
(name)     (ID)  (mass)    (radius)
  H2       105    2.0159   1.445
  He         2    4.0026   1.400
  Na        11   22.9898   2.270
  H2O      101   18.0153   1.600
  CH4      102   16.0425   2.000
  CO       103   28.0101   1.690
  CO2      104   44.0095   1.900
Molecular data taken from (molfile):
    '{:s}/molecules.dat'
""".format(os.path.realpath('./../inputs'))

    assert str(pyrat.lt) == """\
Line-transition information:
Input TLI files (tlifile):
    ['{:s}/HITRAN_H2O_1.1-1.7um_test.tli']
Number of databases (ndb): 1

Database name (name): HITRAN H2O
Species name (molname):  H2O
Number of isotopes (niso): 6
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

Total number of line transitions (ntransitions): 47,666
Minimum and maximum temperatures (tmin, tmax): [1.0, 5000.0] K
Line-transition isotope IDs (isoid):
    [0 0 0 0 0 0 0 ... 3 3 3 3 3 3 3]
Line-transition wavenumbers (wn, cm-1):
    [5882.494 5883.065 5883.353 ... 7494.719 7504.756 7513.791]
Line-transition lower-state energy (elow, cm-1):
    [ 1.807e+03  2.106e+03  2.630e+03 ...  1.244e+03  5.201e+02  6.531e+02]
Line-transition gf (gf, cm-1):
    [ 1.511e-08  1.188e-09  1.083e-08 ...  4.975e-06  1.375e-07  9.624e-07]
""".format(os.getcwd())

    assert str(pyrat.iso) == """\
Isotopes information:
Number of isotopes (niso): 6

Isotope  Molecule      Mass    Isotopic   Database   Extinc-coeff
            index     g/mol       ratio      index    index
 (name)    (imol)    (mass)     (ratio)   (dbindex)  (iext)
    161         3   18.0106   9.973e-01          0   None
    181         3   20.0148   1.999e-03          0   None
    171         3   19.0148   3.719e-04          0   None
    162         3   19.0168   3.107e-04          0   None
    182         3   21.0211   6.230e-07          0   None
    172         3   20.0211   1.158e-07          0   None
Partition function evaluated at atmosperic layers (z):
    [ 1.325e+03  1.325e+03  1.325e+03 ...  3.407e+03  3.412e+03  3.417e+03]
    [ 1.337e+03  1.337e+03  1.337e+03 ...  3.426e+03  3.430e+03  3.436e+03]
    [ 7.980e+03  7.980e+03  7.980e+03 ...  2.042e+04  2.045e+04  2.048e+04]
    [ 6.954e+03  6.954e+03  6.954e+03 ...  1.909e+04  1.912e+04  1.915e+04]
    [ 7.066e+03  7.066e+03  7.066e+03 ...  1.940e+04  1.943e+04  1.946e+04]
    [ 4.228e+04  4.228e+04  4.228e+04 ...  1.157e+05  1.158e+05  1.160e+05]
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

Voigt-profiles' extent (extent, in HWHMs): 20.0
Voigt-profile half-sizes (size) of shape [ndop, nlor]:
[[   214    226 ...   1718   1815]
 [   214    226 ...   1718   1815]
 ...
 [312138 312138 ... 312138 312138]
 [522892 522892 ... 522892 522892]]
Voigt-profile indices (index) of shape [ndop, nlor]:
[[      0     429 ...   53470   56907]
 [  60538   60967 ...  114008  117445]
 ...
 [3616511 3616511 ... 3616511 3616511]
 [4240788 4240788 ... 4240788 4240788]]

Voigt profiles:
  profile[ 0, 0]: [ 7.23955e-07  7.30807e-07 ...  7.30807e-07  7.23955e-07]
  ...
  profile[39,39]: [ 6.55809e-05  6.55811e-05 ...  6.55811e-05  6.55809e-05]
"""

    assert str(pyrat.ex) == """\
Extinction-coefficient information:
Line-transition strength threshold (ethresh): 1.00e-15

LBL extinction coefficient for the atmospheric model (ec, cm-1) [layer, wave]:
[[1.42e-21 1.61e-14 7.52e-23 ... 3.33e-16 2.08e-17 2.86e-17]
 [1.78e-21 2.03e-14 9.46e-23 ... 4.19e-16 2.62e-17 3.60e-17]
 [3.76e-21 2.55e-14 2.00e-22 ... 5.28e-16 3.30e-17 4.53e-17]
 ...
 [5.89e-07 9.44e-07 1.05e-06 ... 2.14e-06 1.95e-06 1.64e-06]
 [7.45e-07 9.61e-07 1.08e-06 ... 2.69e-06 2.33e-06 1.94e-06]
 [9.36e-07 1.21e-06 1.36e-06 ... 3.39e-06 2.94e-06 2.44e-06]]
Extinction-coefficient table filename(s) (extfile): None
"""

    assert str(pyrat.cs) == """\
Cross-section extinction information:
Number of cross-section files (nfiles): 2

Cross-section file name (files[0]):
    '{:s}/inputs/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
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
    '{:s}/inputs/CIA/CIA_Borysow_H2He_0050-3000K_0.3-030um.dat'
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
Observing geometry (path): transit
Total atmospheric extinction coefficient (ec, cm-1) [layer, wave]:
[[ 5.613e-17  1.617e-14  5.620e-17 ...  6.528e-16  3.406e-16  3.485e-16]
 [ 7.066e-17  2.035e-14  7.076e-17 ...  8.218e-16  4.288e-16  4.387e-16]
 [ 8.898e-17  2.562e-14  8.910e-17 ...  1.035e-15  5.399e-16  5.524e-16]
 ...
 [ 1.413e-02  1.413e-02  1.413e-02 ...  1.402e-02  1.402e-02  1.402e-02]
 [ 1.783e-02  1.782e-02  1.782e-02 ...  1.765e-02  1.765e-02  1.765e-02]
 [ 2.250e-02  2.250e-02  2.250e-02 ...  2.222e-02  2.222e-02  2.222e-02]]

Distance along the ray path across each layer (outside-in) at each impact
    parameter (raypath, km):
    IP[  1]: [3061.15674278]
    IP[  2]: [1269.01827346 3057.58923561]
    IP[  3]: [ 974.38612723 1267.64459463 3053.66707207]
    ...
    IP[ 80]: [ 164.81938429  165.39608138  165.94478511 ... 1097.08562146
    1426.24957685
 3434.27193531]

Maximum optical depth to calculate (maxdepth): 10.00
Layer index where the optical depth reaches maxdepth (ideep):
    [ 31  31  31  31  31  31  31 ...  31  31  31  31  31  31  31]
Maximum ideep (deepest layer reaching maxdepth): 31

Optical depth at each impact parameter, down to max(ideep) (depth):
[[ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 3.881e-08  1.118e-05  3.886e-08 ...  4.514e-07  2.355e-07  2.410e-07]
 [ 6.490e-08  1.869e-05  6.499e-08 ...  7.547e-07  3.938e-07  4.029e-07]
 ...
 [ 6.113e-05  1.387e-02  6.054e-05 ...  5.631e-04  2.947e-04  3.038e-04]
 [ 8.113e-05  1.745e-02  8.019e-05 ...  7.094e-04  3.716e-04  3.838e-04]
 [ 8.260e+01  8.262e+01  8.260e+01 ...  8.260e+01  8.260e+01  8.260e+01]]
"""

    assert str(pyrat.cloud) == """\
Cloud-opacity models (models):

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

Patchiness fraction (fpatchy): None
Total atmospheric cloud extinction-coefficient (ec, cm-1):
[[ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 ...
 [ 1.401e-02  1.401e-02  1.401e-02 ...  1.401e-02  1.401e-02  1.401e-02]
 [ 1.764e-02  1.764e-02  1.764e-02 ...  1.764e-02  1.764e-02  1.764e-02]
 [ 2.220e-02  2.220e-02  2.220e-02 ...  2.220e-02  2.220e-02  2.220e-02]]
"""

    assert str(pyrat.rayleigh) == """\
Rayleigh-opacity models (models):

Model name (name): 'lecavelier'
Model species (mol): H2
Number of model parameters (npars): 2
Parameter name     Value
  (pnames)         (pars)
  log(f_Ray)        0.000e+00
  alpha_Ray        -4.000e+00
Extinction-coefficient (ec, cm2 molec-1):
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
Detuning parameter (detuning): 30.0
Lorentz-width parameter (lpar): 0.071
Partition function (Z): 2.0
Wavenumber  Wavelength          gf   Lower-state energy
      cm-1          um               cm-1
      (wn)                    (gf)   (elow)
  16960.87    0.589592   6.546e-01   0.000e+00
  16978.07    0.588995   1.309e+00   0.000e+00
Extinction-coefficient (ec, cm2 molecule-1):
[[ 9.230e-34  9.244e-34 ...  1.263e-31  1.265e-31]
 [ 1.162e-33  1.164e-33 ...  1.590e-31  1.593e-31]
 ...
 [ 1.299e-23  1.300e-23 ...  3.468e-22  3.472e-22]
 [ 1.636e-23  1.637e-23 ...  4.358e-22  4.363e-22]]

Total atmospheric alkali extinction-coefficient (ec, cm-1):
[[ 1.915e-26  1.918e-26  1.921e-26 ...  2.617e-24  2.621e-24  2.625e-24]
 [ 3.036e-26  3.040e-26  3.045e-26 ...  4.148e-24  4.154e-24  4.161e-24]
 [ 4.812e-26  4.819e-26  4.826e-26 ...  6.574e-24  6.585e-24  6.595e-24]
 ...
 [ 8.496e-09  8.504e-09  8.513e-09 ...  2.270e-07  2.272e-07  2.275e-07]
 [ 1.347e-08  1.349e-08  1.350e-08 ...  3.593e-07  3.597e-07  3.601e-07]
 [ 2.134e-08  2.136e-08  2.138e-08 ...  5.681e-07  5.687e-07  5.693e-07]]
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
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
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
    [ 5882.459  5883.636  5884.813 ...  9086.365  9088.182  9090.000]
Oversampling factor (wnosamp): 2160

Modulation spectrum, (Rp/Rs)**2 (spectrum):
    [ 6.522e-03  6.538e-03  6.521e-03 ...  6.506e-03  6.462e-03  6.532e-03]
"""


def test_pyrat_emission_str(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
    reset={'path':'eclipse'})
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
    [ 7.815e+02  7.808e+02  7.801e+02 ...  3.602e+01  3.598e+01  3.594e+01]
    [ 7.812e+02  7.805e+02  7.798e+02 ...  3.599e+01  3.596e+01  3.592e+01]
    [ 7.803e+02  7.796e+02  7.790e+02 ...  3.593e+01  3.590e+01  3.586e+01]
    [ 7.789e+02  7.782e+02  7.776e+02 ...  3.583e+01  3.580e+01  3.576e+01]
    [ 7.775e+02  7.768e+02  7.762e+02 ...  3.573e+01  3.570e+01  3.566e+01]
Emission spectrum (spectrum, erg s-1 cm-2 cm):
    [ 2.450e+03  2.448e+03  2.446e+03 ...  1.128e+02  1.127e+02  1.125e+02]
"""

    assert str(pyrat.od) == """\
Optical depth information:
Observing geometry (path): eclipse
Total atmospheric extinction coefficient (ec, cm-1) [layer, wave]:
[[ 5.613e-17  1.617e-14  5.620e-17 ...  6.528e-16  3.406e-16  3.485e-16]
 [ 7.066e-17  2.035e-14  7.076e-17 ...  8.218e-16  4.288e-16  4.387e-16]
 [ 8.898e-17  2.562e-14  8.910e-17 ...  1.035e-15  5.399e-16  5.524e-16]
 ...
 [ 1.413e-02  1.413e-02  1.413e-02 ...  1.402e-02  1.402e-02  1.402e-02]
 [ 1.783e-02  1.782e-02  1.782e-02 ...  1.765e-02  1.765e-02  1.765e-02]
 [ 2.250e-02  2.250e-02  2.250e-02 ...  2.222e-02  2.222e-02  2.222e-02]]

Distance across each layer along a normal ray path (raypath, km):
    [62.8 62.7 62.6 62.5 ... 86.0 85.8 85.7 85.5]

Maximum optical depth to calculate (maxdepth): 10.00
Layer index where the optical depth reaches maxdepth (ideep):
    [ 35  35  35  35  35  35  35 ...  35  35  35  35  35  35  35]
Maximum ideep (deepest layer reaching maxdepth): 35

Planck emission down to max(ideep) (B, erg s-1 cm-2 sr-1 cm):
[[ 7.487e+02  7.480e+02  7.474e+02 ...  3.370e+01  3.367e+01  3.363e+01]
 [ 7.487e+02  7.480e+02  7.474e+02 ...  3.371e+01  3.367e+01  3.363e+01]
 [ 7.487e+02  7.480e+02  7.474e+02 ...  3.371e+01  3.367e+01  3.364e+01]
 ...
 [ 8.003e+02  7.996e+02  7.989e+02 ...  3.736e+01  3.732e+01  3.728e+01]
 [ 8.140e+02  8.134e+02  8.127e+02 ...  3.836e+01  3.832e+01  3.828e+01]
 [ 8.315e+02  8.308e+02  8.301e+02 ...  3.964e+01  3.960e+01  3.956e+01]]

Optical depth at each layer along a normal ray path into the planet, down to
    max(ideep) (depth):
[[ 0.000e+00  0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00  0.000e+00]
 [ 3.984e-10  1.147e-07  3.989e-10 ...  4.633e-09  2.418e-09  2.473e-09]
 [ 8.992e-10  2.590e-07  9.004e-10 ...  1.046e-08  5.457e-09  5.583e-09]
 ...
 [ 5.121e+00  5.121e+00  5.121e+00 ...  5.121e+00  5.121e+00  5.121e+00]
 [ 8.125e+00  8.126e+00  8.125e+00 ...  8.125e+00  8.125e+00  8.125e+00]
 [ 1.191e+01  1.191e+01  1.191e+01 ...  1.191e+01  1.191e+01  1.191e+01]]
"""


def test_pyrat_exfile_str(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/mcmc_transmission_test.cfg',
        reset={'runmode':'spectrum',
               'outspec':'extfiled_spectrum_test.dat'})
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
[[1.35e-21 1.50e-14 7.16e-23 ... 4.43e-16 2.13e-17 2.76e-17]
 [2.21e-21 1.89e-14 1.17e-22 ... 5.57e-16 2.68e-17 3.47e-17]
 [3.54e-21 2.38e-14 1.87e-22 ... 7.02e-16 3.38e-17 4.37e-17]
 ...
 [5.87e-07 9.36e-07 1.04e-06 ... 2.11e-06 1.93e-06 1.62e-06]
 [7.41e-07 9.57e-07 1.08e-06 ... 2.66e-06 2.31e-06 1.91e-06]
 [9.32e-07 1.20e-06 1.36e-06 ... 3.35e-06 2.90e-06 2.41e-06]]
Extinction-coefficient table filename(s) (extfile):
    {:s}/exttable_test_300-3000K_1.1-1.7um.dat
Minimum temperature (tmin, K):  300.0
Maximum temperature (tmax, K): 3000.0
Temperature sampling interval (tstep, K): None

Number of species (nmol):               1
Number of temperatures (ntemp):        10
Number of layers (nlayers):            81
Number of spectral samples (nwave):  3209

Species ID array (molID): [101]
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
[[[[1.26e-32 9.75e-27 ... 1.48e-26 2.90e-28]
   [2.08e-32 9.75e-27 ... 1.48e-26 2.91e-28]
   ...
   [3.26e-25 3.70e-25 ... 1.34e-24 1.09e-24]
   [3.02e-25 3.32e-25 ... 1.33e-24 1.14e-24]]

  [[2.02e-31 1.28e-24 ... 1.35e-26 5.44e-27]
   [3.34e-31 1.28e-24 ... 1.35e-26 5.44e-27]
   ...
   [2.68e-24 3.03e-24 ... 3.87e-24 3.20e-24]
   [2.68e-24 3.03e-24 ... 3.87e-24 3.20e-24]]

  ...

  [[1.59e-31 3.50e-24 ... 6.07e-27 1.01e-26]
   [1.59e-31 3.50e-24 ... 6.07e-27 1.01e-26]
   ...
   [2.82e-24 4.48e-24 ... 1.33e-23 1.12e-23]
   [2.98e-24 3.91e-24 ... 1.22e-23 1.03e-23]]

  [[1.28e-31 2.79e-24 ... 6.08e-27 1.06e-26]
   [1.28e-31 2.79e-24 ... 6.08e-27 1.06e-26]
   ...
   [2.27e-24 3.61e-24 ... 1.13e-23 9.63e-24]
   [2.42e-24 3.18e-24 ... 1.04e-23 8.76e-24]]]]
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
   7485.12       1.336     0.00665  filter_test_WFC3_G141_1.336um.dat
   7345.18       1.361     0.00671  filter_test_WFC3_G141_1.361um.dat
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
  log(kappa)       -1.000e+00  -3.000e+00   3.000e+00   3.000e-01  temp
  log(gamma1)       0.000e+00  -3.000e+00   3.000e+00   3.000e-01  temp
  log(gamma2)       0.000e+00  -3.000e+00   3.000e+00   0.000e+00  temp
  alpha             0.000e+00   0.000e+00   1.000e+00   0.000e+00  temp
  beta              1.000e+00   1.000e-01   2.000e+00   1.000e-01  temp
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
    {:s}/MCMC_transmission_test.npz
""".format(os.getcwd())

