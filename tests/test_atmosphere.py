import os
import sys
import pytest

import numpy as np

sys.path.append('../')
import pyratbay.atmosphere as pa
import pyratbay.constants  as pc


expected_pressure = np.array([1.e-02, 1.e-01, 1.e+00, 1.e+01, 1.e+02, 1.e+03,
                              1.e+04, 1.e+05, 1.e+06, 1.e+07, 1.e+08])
expected_temp = np.array(
      [1175.88704562, 1175.88696601, 1175.88632398, 1175.88144499,
       1175.84806406, 1175.66810064, 1175.38603129, 1185.64394212,
       1286.15511142, 1357.93586563, 1358.82442802])


def test_pressure_floats():
    ptop    = 1e-8
    pbottom = 1e2
    nlayers = 11
    units   = "bar"
    pressure = pa.pressure(ptop, pbottom, nlayers, units)
    np.testing.assert_almost_equal(pressure, expected_pressure, decimal=7)


def test_pressure_units():
    ptop    = "1e-8 bar"
    pbottom = "1e2 bar"
    nlayers = 11
    pressure = pa.pressure(ptop, pbottom, nlayers)
    np.testing.assert_almost_equal(pressure, expected_pressure, decimal=7)


@pytest.mark.parametrize("tinputs",
    [1500.0, [1500.0], (1500,), np.array([1500.0])])
def test_temperature_isothermal(tinputs):
    nlayers = 100
    temp = pa.temp_isothermal(tinputs, nlayers)
    np.testing.assert_equal(temp, np.tile(1500.0, 100))


def test_temperature_TCEA_floats():
    tparams = [-3.0, -0.25, 0.0, 0.0, 1.0]
    pressure = expected_pressure
    rstar = 1.0 * pc.rsun
    tstar = 5800.0
    tint = 100.0
    gplanet = 800.0
    smaxis = 0.05 * pc.au
    temp = pa.temp_TCEA(tparams, pressure, rstar, tstar, tint, gplanet, smaxis)
    np.testing.assert_almost_equal(temp, expected_temp, decimal=7)


def test_temperature_TCEA_units():
    tparams = [-3.0, -0.25, 0.0, 0.0, 1.0]
    pressure = expected_pressure
    rstar = "1.0 rsun"
    tstar = 5800.0
    tint = 100.0
    gplanet = 800.0
    smaxis = "0.05 au"
    temp = pa.temp_TCEA(tparams, pressure, rstar, tstar, tint, gplanet, smaxis)
    np.testing.assert_almost_equal(temp, expected_temp, decimal=7)


def test_uniform():
    atmfile = "atm_test.dat"
    nlayers = 11
    pressure    = np.logspace(-8, 2 , nlayers)
    temperature = np.tile(1500.0, nlayers)
    species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    pa.uniform(atmfile, pressure, temperature, species, abundances)
    # Now what?


