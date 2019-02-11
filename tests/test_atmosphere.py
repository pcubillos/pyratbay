import os
import sys
import pytest

import numpy as np

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay.atmosphere as pa
import pyratbay.constants  as pc


expected_pressure = np.array([1.e-02, 1.e-01, 1.e+00, 1.e+01, 1.e+02, 1.e+03,
                              1.e+04, 1.e+05, 1.e+06, 1.e+07, 1.e+08])
expected_temp = np.array(
      [1175.88704562, 1175.88696601, 1175.88632398, 1175.88144499,
       1175.84806406, 1175.66810064, 1175.38603129, 1185.64394212,
       1286.15511142, 1357.93586563, 1358.82442802])

# Jupiter radius profile at constant g:
radius_g = np.array(
      [1.0563673 , 1.04932138, 1.04227547, 1.03522956, 1.02818365, 1.02113774,
       1.01409182, 1.00704591, 1.        , 0.99295409, 0.98590818])
# Jupiter radius profile with g(r)=GM/r**2:
radius_m = np.array(
      [1.05973436, 1.05188019, 1.04414158, 1.036516  , 1.029001  , 1.02159419,
       1.01429324, 1.00709591, 1.        , 0.99300339, 0.986104  ])

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
    atmfile = "uniform_test.atm"
    nlayers = 11
    pressure    = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    temperature = np.tile(1500.0, nlayers)
    species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    qprofiles = pa.uniform(pressure, temperature, species, abundances)
    assert np.shape(qprofiles) == (nlayers, len(species))
    for q in qprofiles:
        np.testing.assert_equal(q, np.array(abundances))


def test_abundances_uniform():
    atmfile = "atm_test.dat"
    nlayers = 11
    punits  = 'bar'
    pressure    = pa.pressure(1e-8, 1e2, nlayers, punits)
    temperature = pa.temp_isothermal(1500.0, nlayers)
    species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    qprofiles = pa.abundances(atmfile, pressure, temperature, species,
                              quniform=abundances, punits=punits)
    assert np.shape(qprofiles) == (nlayers, len(species))
    for q in qprofiles:
        np.testing.assert_equal(q, np.array(abundances))
    # TBD: Check file is there


def test_hydro_g():
    nlayers = 11
    pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    temperature = pa.temp_isothermal(1500.0, nlayers)
    mu = np.tile(2.3, nlayers)
    g = pc.G * pc.mjup / pc.rjup**2
    r0 = 1.0 * pc.rjup
    p0 = 1.0 * pc.bar
    # Radius profile in Jupiter radii:
    radius = pa.hydro_g(pressure, temperature, mu, g, p0, r0) / pc.rjup
    np.testing.assert_almost_equal(radius, radius_g, decimal=7)


def test_hydro_m():
    nlayers = 11
    pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    temperature = pa.temp_isothermal(1500.0, nlayers)
    mu = np.tile(2.3, nlayers)
    Mp = 1.0 * pc.mjup
    r0 = 1.0 * pc.rjup
    p0 = 1.0 * pc.bar
    radius = pa.hydro_m(pressure, temperature, mu, Mp, p0, r0) / pc.rjup
    # Radius profile in Jupiter radii:
    np.testing.assert_almost_equal(radius, radius_m, decimal=7)


def test_stoich():
    species = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    elements, stoichs = pa.stoich(species)
    assert elements == ['C', 'H', 'He', 'O']
    np.testing.assert_equal(stoichs, np.array([[0, 2, 0, 0],
                                               [0, 0, 1, 0],
                                               [0, 2, 0, 1],
                                               [1, 0, 0, 1],
                                               [1, 0, 0, 2],
                                               [1, 4, 0, 0]]))


@pytest.mark.parametrize("abundances",
    [ [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4],
      [[0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]]])
def test_meanweight(abundances):
    species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    mu = pa.meanweight(abundances, species)
    np.testing.assert_almost_equal(mu, np.array([2.31928918]), decimal=7)
