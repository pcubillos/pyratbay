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
      [1047.04157312, 1047.04189805, 1047.04531644, 1047.08118784,
       1047.45648563, 1051.34469989, 1088.69956369, 1311.86379107,
       1640.12857767, 1660.02396061, 1665.30121021])

# Jupiter radius profile at constant g:
radius_g = np.array(
      [1.0563673 , 1.04932138, 1.04227547, 1.03522956, 1.02818365, 1.02113774,
       1.01409182, 1.00704591, 1.        , 0.99295409, 0.98590818])

# Jupiter radius profile with g(r)=GM/r**2:
radius_m = np.array(
      [1.05973436, 1.05188019, 1.04414158, 1.036516  , 1.029001  , 1.02159419,
       1.01429324, 1.00709591, 1.        , 0.99300339, 0.986104  ])

expected_dens = np.array(
      [4.10241993e+10, 7.24297303e+09, 4.82864869e+06, 4.82864869e+06,
       4.82864869e+02, 4.82864869e+06])

q0 = np.array(
      [[  7.15780000e-01,   1.34580000e-01,   3.48390000e-04,
          4.25490000e-04,   8.00670000e-08,   1.93810000e-22],
       [  8.07420000e-01,   1.41760000e-01,   3.67350000e-04,
          4.48220000e-04,   7.88430000e-08,   2.77920000e-20],
       [  8.38830000e-01,   1.44230000e-01,   3.73770000e-04,
          4.56010000e-04,   7.85590000e-08,   3.11600000e-18],
       [  8.49010000e-01,   1.45030000e-01,   3.75850000e-04,
          4.58530000e-04,   7.84790000e-08,   3.23090000e-16],
       [  8.52260000e-01,   1.45280000e-01,   3.76510000e-04,
          4.59340000e-04,   7.84550000e-08,   3.26810000e-14],
       [  8.53290000e-01,   1.45360000e-01,   3.76720000e-04,
          4.59590000e-04,   7.84480000e-08,   3.27990000e-12],
       [  8.53610000e-01,   1.45390000e-01,   3.76780000e-04,
          4.59670000e-04,   7.84450000e-08,   3.28370000e-10],
       [  8.53720000e-01,   1.45390000e-01,   3.76840000e-04,
          4.59670000e-04,   7.84460000e-08,   3.28440000e-08],
       [  8.53750000e-01,   1.45400000e-01,   3.80050000e-04,
          4.56480000e-04,   7.85620000e-08,   3.23430000e-06],
       [  8.53560000e-01,   1.45440000e-01,   5.31500000e-04,
          3.05290000e-04,   7.34970000e-08,   1.54570000e-04],
       [  8.53190000e-01,   1.45530000e-01,   8.23730000e-04,
          1.36860000e-05,   5.10860000e-09,   4.46510000e-04]])


def test_pressure_default_units():
    ptop    = 1e-8
    pbottom = 1e2
    nlayers = 11
    pressure = pa.pressure(ptop, pbottom, nlayers)
    np.testing.assert_almost_equal(pressure, expected_pressure, decimal=7)


def test_pressure_floats():
    ptop    = 1e-8
    pbottom = 1e2
    nlayers = 11
    units   = "bar"
    pressure = pa.pressure(ptop, pbottom, nlayers, units)
    np.testing.assert_almost_equal(pressure, expected_pressure, decimal=7)


def test_pressure_with_units():
    ptop    = "1e-8 bar"
    pbottom = "1e2 bar"
    nlayers = 11
    pressure = pa.pressure(ptop, pbottom, nlayers)
    np.testing.assert_almost_equal(pressure, expected_pressure, decimal=7)


@pytest.mark.parametrize("tparams",
    [1500.0, [1500.0], (1500,), np.array([1500.0])])
def test_temp_isothermal(tparams):
    nlayers = 100
    temp = pa.temp_isothermal(tparams, nlayers)
    np.testing.assert_equal(temp, np.tile(1500.0, 100))


@pytest.mark.parametrize("tparams",
    [[-1.5, -0.8, -0.8, 0.5, 1.0],
     np.array([-1.5, -0.8, -0.8, 0.5, 1.0])
    ])
def test_temp_TCEA_floats(tparams):
    tparams = [-1.5, -0.8, -0.8, 0.5, 1.0]
    pressure = expected_pressure
    rstar = 0.756 * pc.rsun
    tstar = 5040.0
    tint = 100.0
    gplanet = 2200.0
    smaxis = 0.031 * pc.au
    temp = pa.temp_TCEA(tparams, pressure, rstar, tstar, tint, gplanet, smaxis)
    np.testing.assert_almost_equal(temp, expected_temp, decimal=7)


def test_temp_TCEA_units():
    tparams = [-1.5, -0.8, -0.8, 0.5, 1.0]
    pressure = expected_pressure
    rstar = "0.756 rsun"
    tstar = 5040.0
    tint = 100.0
    gplanet = 2200.0
    smaxis = "0.031 au"
    temp = pa.temp_TCEA(tparams, pressure, rstar, tstar, tint, gplanet, smaxis)
    np.testing.assert_almost_equal(temp, expected_temp, decimal=7)


def test_temperature_isothermal():
    tparams = 1500.0
    nlayers = 11
    temp = pa.temperature("isothermal", tparams=tparams, nlayers=nlayers)
    np.testing.assert_equal(temp, np.tile(tparams, nlayers))


def test_temperature_TCEA():
    nlayers = 11
    pressure = expected_pressure
    tparams = [-1.5, -0.8, -0.8, 0.5, 1.0]
    temp = pa.temperature("TCEA", pressure, rstar="0.756 rsun", tstar=5040,
          tint=100.0, gplanet=2200.0, smaxis="0.031 au", tparams=tparams)
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


@pytest.mark.skip(reason="Skip until implementing the fast TEA")
def test_abundances_TEA():
    atmfile = "atm_test.dat"
    nlayers = 11
    punits  = 'bar'
    pressure    = pa.pressure(1e-8, 1e2, nlayers, punits)
    temperature = pa.temp_isothermal(1500.0, nlayers)
    species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    elements    = ["H", "He", "C", "O"]
    xsolar      = 1.0
    qprofiles = pa.abundances(atmfile, pressure, temperature, species,
                              elements, punits=punits, xsolar=xsolar)
    #solar = pc.ROOT+"inputs/AsplundEtal2009.txt"
    #index, symbol, dex, name, mass = pa.readatomic(solar)


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


def test_IGLdensity():
    atmfile = "uniform_test.atm"
    nlayers = 11
    pressure    = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    temperature = np.tile(1500.0, nlayers)
    species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    qprofiles = pa.uniform(pressure, temperature, species, abundances)
    dens = pa.IGLdensity(qprofiles, pressure, temperature)
    for i,density in enumerate(dens):
        np.testing.assert_allclose(density, expected_dens*10**i, rtol=1e-7)


def test_ratio():
    q = np.tile([0.8, 0.2], (5,1))
    q[4] = 0.5, 0.5
    ibulk = [0, 1]
    bratio, invsrat = pa.ratio(q, ibulk)
    np.testing.assert_equal(bratio[:,0], np.tile(1.0, 5))
    np.testing.assert_equal(bratio[:,1], np.array([0.25,0.25,0.25,0.25,1.0]))
    np.testing.assert_equal(invsrat, np.array([0.8, 0.8, 0.8, 0.8, 0.5]))


def test_balance():
    q = np.tile([0.8, 0.2, 0.5], (5,1))
    q[4] = 0.5, 0.5, 0.5
    ibulk = [0, 1]
    bratio, invsrat = pa.ratio(q, ibulk)
    pa.balance(q, ibulk, bratio, invsrat)
    # Check sum(q) == 1:
    np.testing.assert_equal(np.sum(q,axis=1), np.tile(1.0,5))
    # Check ratio(q) == bratio:
    np.testing.assert_equal(q[:,1]/q[:,0], bratio[:,1]/bratio[:,0])


def test_qscale():
    spec = np.array(["H2", "He", "H2O", "CO", "CO2", "CH4"])
    bulk = np.array(['H2', 'He'])
    molmodel = ['vert', 'scale']
    molfree  = ['H2O', 'CO']
    molpars  = [-4, 1.0]
    q2 = pa.qscale(q0, spec, molmodel, molfree, molpars, bulk)
    nlayers, nspec = np.shape(q0)
    # All H2O abundances set to constant value:
    np.testing.assert_equal(q2[:,2], np.tile(10**molpars[0], nlayers))
    # All CO abundances scaled by value:
    np.testing.assert_allclose(q2[:,3], q0[:,3]*10**molpars[1], rtol=1e-7)


