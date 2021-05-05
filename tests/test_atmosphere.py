# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os
import pytest

import numpy as np

import pyratbay.atmosphere as pa
import pyratbay.constants  as pc

os.chdir(pc.ROOT+'tests')


expected_pressure = np.logspace(-2, 9, 15)

expected_temp_tcea = np.array(
      [1046.89057381, 1046.89075751, 1046.89192532, 1046.89933754,
       1046.94631087, 1047.24341507, 1049.11331707, 1060.60902021,
       1123.15986552, 1339.81840964, 1617.02710403, 1659.45254019,
       1660.78464365, 1668.82931802, 1715.58904031])

expected_temp_madhu_noinv = np.array(
      [ 870.23835532,  875.0069192 ,  888.59773802,  911.24902223,
        942.96082012,  983.73313169, 1033.56595695, 1092.4592959 ,
       1160.41290637, 1236.232834  , 1300.7457372 , 1372.71764431,
       1453.77542722, 1460.59618283, 1460.7314666 ])

expected_temp_madhu_inv = np.array(
      [870.23835532, 875.0069192 , 888.59773802, 911.24902223,
       942.96037048, 981.51425888, 988.43599566, 952.27744069,
       927.46081518, 917.22633753, 921.57466482, 940.50529625,
       971.54884044, 974.31160987, 974.37075759])

# Jupiter radius profile at constant g:
radius_g = np.array([
    1.05636546, 1.04931978, 1.04227409, 1.03522841, 1.02818273,
    1.02113705, 1.01409136, 1.00704568, 1.        , 0.99295432, 0.98590864])

# Jupiter radius profile with g(r)=GM/r**2:
radius_m = np.array([
    1.0597323 , 1.05187841, 1.04414007, 1.03651477, 1.02900003,
    1.02159347, 1.01429277, 1.00709568, 1.        , 0.99300361, 0.98610444])

expected_dens = np.array([
    4.10241850e+10, 7.24297052e+09, 4.82864701e+06, 4.82864701e+06,
    4.82864701e+02, 4.82864701e+06])

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

qtea_expected = np.array([
    [7.3360e-01, 9.2056e-02, 9.0542e-21, 7.2656e-05, 1.7382e-01,
     1.6586e-04, 2.9100e-04, 1.0735e-07, 3.9871e-28],
    [4.1413e-01, 1.1529e-01, 1.7346e-20, 5.5995e-06, 4.6992e-01,
     2.9316e-04, 3.6449e-04, 8.7915e-08, 4.0179e-25],
    [1.7275e-01, 1.3284e-01, 2.5163e-20, 5.2429e-07, 6.9365e-01,
     3.4373e-04, 4.2000e-04, 8.0466e-08, 9.1401e-23],
    [6.3491e-02, 1.4078e-01, 2.8797e-20, 5.7232e-08, 7.9492e-01,
     3.6479e-04, 4.4512e-04, 7.8974e-08, 9.8864e-21],
    [2.2316e-02, 1.4378e-01, 3.0175e-20, 6.5751e-09, 8.3308e-01,
     3.7260e-04, 4.5458e-04, 7.8606e-08, 8.1888e-19],
    [7.7236e-03, 1.4484e-01, 3.0665e-20, 7.6831e-10, 8.4661e-01,
     3.7536e-04, 4.5794e-04, 7.8498e-08, 6.1850e-17],
    [2.6591e-03, 1.4521e-01, 3.0835e-20, 9.0296e-11, 8.5130e-01,
     3.7631e-04, 4.5910e-04, 7.8462e-08, 4.5257e-15],
    [9.1383e-04, 1.4533e-01, 3.0893e-20, 1.0633e-11, 8.5292e-01,
     3.7664e-04, 4.5950e-04, 7.8450e-08, 3.2757e-13],
    [3.1385e-04, 1.4538e-01, 3.0914e-20, 1.2529e-12, 8.5347e-01,
     3.7675e-04, 4.5964e-04, 7.8446e-08, 2.3621e-11],
    [1.0777e-04, 1.4539e-01, 3.0920e-20, 1.4768e-13, 8.5367e-01,
     3.7679e-04, 4.5969e-04, 7.8445e-08, 1.7010e-09],
    [3.7001e-05, 1.4540e-01, 3.0905e-20, 1.7412e-14, 8.5373e-01,
     3.7693e-04, 4.5958e-04, 7.8449e-08, 1.2238e-07],
    [1.2704e-05, 1.4540e-01, 2.9688e-20, 2.0979e-15, 8.5374e-01,
     3.8528e-04, 4.5126e-04, 7.8733e-08, 8.4610e-06],
    [4.3609e-06, 1.4546e-01, 1.0486e-20, 3.8007e-16, 8.5348e-01,
     5.9194e-04, 2.4497e-04, 6.5687e-08, 2.1495e-04],
    [1.4970e-06, 1.4553e-01, 3.0538e-22, 6.2646e-17, 8.5318e-01,
     8.2745e-04, 9.9759e-06, 3.7406e-09, 4.5022e-04],
    [5.1396e-07, 1.4553e-01, 4.3365e-24, 7.4728e-18, 8.5317e-01,
     8.3731e-04, 1.4335e-07, 5.4391e-11, 4.6007e-04]])

q_H2O = [[
     [3.758500e-04, 3.765100e-04, 3.767200e-04, 3.767800e-04, 3.768400e-04,
      3.800500e-04, 5.315000e-04, 8.237300e-04, 8.373100e-04],
     [1.501200e-16, 1.518300e-14, 1.523800e-12, 1.525600e-10, 1.526100e-08,
      1.521900e-06, 1.206600e-04, 7.737300e-04, 8.382700e-04]],
    [[3.736000e-03, 3.742500e-03, 3.744600e-03, 3.745300e-03, 3.745500e-03,
      3.748800e-03, 4.027400e-03, 7.363500e-03, 8.393900e-03],
     [1.389200e-16, 1.405100e-14, 1.410200e-12, 1.411800e-10, 1.412300e-08,
      1.412100e-06, 1.375500e-04, 4.689100e-03, 8.482200e-03]]
]


def test_pressure_default_units():
    ptop    = 1e-8
    pbottom = 1e3
    nlayers = 15
    pressure = pa.pressure(ptop, pbottom, nlayers)
    np.testing.assert_allclose(pressure, expected_pressure)


def test_pressure_floats():
    ptop    = 1e-8
    pbottom = 1e3
    nlayers = 15
    units   = "bar"
    pressure = pa.pressure(ptop, pbottom, nlayers, units)
    np.testing.assert_allclose(pressure, expected_pressure)


def test_pressure_with_units():
    ptop    = "1e-8 bar"
    pbottom = "1e3 bar"
    nlayers = 15
    pressure = pa.pressure(ptop, pbottom, nlayers)
    np.testing.assert_allclose(pressure, expected_pressure)


@pytest.mark.parametrize("params",
    [1500.0, [1500.0], (1500,), np.array([1500.0])])
def test_tmodel_isothermal(params):
    nlayers = 100
    tmodel = pa.tmodels.Isothermal(nlayers)
    np.testing.assert_equal(tmodel(params), np.tile(1500.0, nlayers))


@pytest.mark.parametrize("params",
    [[-4.84, -0.8, -0.8, 0.5, 1200.0, 100.0],
     np.array([-4.84, -0.8, -0.8, 0.5, 1200.0, 100.0])
    ])
def test_tmodel_tcea_floats(params):
    pressure = expected_pressure
    tmodel = pa.tmodels.TCEA(pressure)
    np.testing.assert_allclose(tmodel(params), expected_temp_tcea)


@pytest.mark.parametrize('gravity',
    [None, 2200.0, np.tile(2200.0, len(expected_pressure))])
def test_tmodel_tcea_gravity(gravity):
    params = np.array([-4.84, -0.8, -0.8, 0.5, 1200.0, 100.0])
    pressure = expected_pressure
    if gravity is not None:
        params[0] += np.log10(2200.0)
    tmodel = pa.tmodels.TCEA(pressure, gravity)
    np.testing.assert_allclose(tmodel(params), expected_temp_tcea)


def test_temp_madhu_no_inv():
    pressure = expected_pressure
    tmodel = pa.tmodels.Madhu(pressure)
    params = -0.77, -3.61, 1.45, 0.85, 0.67, 870.0
    np.testing.assert_allclose(tmodel(params), expected_temp_madhu_noinv)


def test_temp_madhu_inv():
    pressure = expected_pressure
    tmodel = pa.tmodels.Madhu(pressure)
    params = -3.61, -0.77, 1.45, 0.85, 0.67, 870.0
    np.testing.assert_allclose(tmodel(params), expected_temp_madhu_inv)


def test_temp_madhu_invalid_params():
    # p1 > p3:
    pressure = expected_pressure
    tmodel = pa.tmodels.Madhu(pressure)
    params = 2.0, -0.77, 1.5, 0.85, 0.67, 870.0
    np.testing.assert_allclose(tmodel(params), np.tile(0.0,len(pressure)))


def test_temperature_isothermal():
    params = 1500.0
    nlayers = 15
    temp = pa.temperature("isothermal", params=params, nlayers=nlayers)
    np.testing.assert_equal(temp, np.tile(params, nlayers))


def test_temperature_tcea():
    pressure = expected_pressure
    params = [-4.84, -0.8, -0.8, 0.5, 1200.0, 100.0]
    temp = pa.temperature('tcea', pressure, params=params)
    np.testing.assert_allclose(temp, expected_temp_tcea)


def test_temperature_madhu():
    pressure = expected_pressure
    params = -0.77, -3.61, 1.45, 0.85, 0.67, 870.0
    temp = pa.temperature('madhu', pressure, params=params)
    np.testing.assert_allclose(temp, expected_temp_madhu_noinv)


def test_uniform():
    nlayers = 11
    pressure    = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    temperature = np.tile(1500.0, nlayers)
    species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    qprofiles = pa.uniform(pressure, temperature, species, abundances)
    assert np.shape(qprofiles) == (nlayers, len(species))
    for q in qprofiles:
        np.testing.assert_equal(q, np.array(abundances))


def test_abundance_uniform():
    atmfile = "outputs/atm_test.dat"
    nlayers = 11
    punits  = 'bar'
    pressure    = pa.pressure(1e-8, 1e2, nlayers, punits)
    tmodel = pa.tmodels.Isothermal(nlayers)
    temperature = tmodel(1500.0)
    species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    qprofiles = pa.abundance(pressure, temperature, species,
        quniform=abundances, atmfile=atmfile, punits=punits)
    assert np.shape(qprofiles) == (nlayers, len(species))
    for q in qprofiles:
        np.testing.assert_equal(q, np.array(abundances))


def test_abundances_tea():
    nlayers = 15
    punits  = 'bar'
    pressure = pa.pressure(1e-10, 1e3, nlayers, punits)
    tmodel = pa.tmodels.Isothermal(nlayers)
    temperature = tmodel(1500.0)
    species     = 'H He C O H2 H2O CO CO2 CH4'.split()
    elements    = 'H He C O'.split()
    qtea = pa.abundance(pressure, temperature, species, elements,
        punits=punits)
    np.testing.assert_allclose(qtea, qtea_expected)


@pytest.mark.parametrize('xsolar', [1.0, 10.0])
@pytest.mark.parametrize('escale', [{}, {'C': 5.0}])
def test_abundances_tea_xsolar_escale(xsolar, escale):
    nlayers = 9
    punits  = 'bar'
    pressure = pa.pressure(1e-5, 1e3, nlayers, punits)
    tmodel = pa.tmodels.Isothermal(nlayers)
    temperature = tmodel(1500.0)
    species     = 'H He C O H2 H2O CO CO2 CH4'.split()
    elements    = 'H He C O'.split()
    q = pa.abundance(pressure, temperature, species, elements,
        punits=punits, xsolar=xsolar, escale=escale)
    np.testing.assert_allclose(q[:,5], q_H2O[xsolar>1]['C' in escale])


def test_hydro_g():
    nlayers = 11
    pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    tmodel = pa.tmodels.Isothermal(nlayers)
    temperature = tmodel(1500.0)
    mu = np.tile(2.3, nlayers)
    g = pc.G * pc.mjup / pc.rjup**2
    r0 = 1.0 * pc.rjup
    p0 = 1.0 * pc.bar
    # Radius profile in Jupiter radii:
    radius = pa.hydro_g(pressure, temperature, mu, g, p0, r0) / pc.rjup
    np.testing.assert_allclose(radius, radius_g)


def test_hydro_m():
    nlayers = 11
    pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    tmodel = pa.tmodels.Isothermal(nlayers)
    temperature = tmodel(1500.0)
    mu = np.tile(2.3, nlayers)
    Mp = 1.0 * pc.mjup
    r0 = 1.0 * pc.rjup
    p0 = 1.0 * pc.bar
    radius = pa.hydro_m(pressure, temperature, mu, Mp, p0, r0) / pc.rjup
    # Radius profile in Jupiter radii:
    np.testing.assert_allclose(radius, radius_m)


def test_hydro_m_ultra_puff():
    nlayers = 15
    pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    tmodel = pa.tmodels.Isothermal(nlayers)
    temperature = tmodel(1500.0)
    mu = np.tile(2.3, nlayers)
    Mp = 0.1 * pc.mjup
    r0 = 2.0 * pc.rjup
    p0 = 1.0 * pc.bar
    radius = pa.hydro_m(pressure, temperature, mu, Mp, p0, r0) / pc.rjup

    puff_radius = np.array([
       23.59979187, 10.78753812,  6.99174271,  5.17191017,  4.10376853,
        3.40130545,  2.90418177,  2.5338435 ,  2.24727349,  2.01893776,
        1.83272274,  1.67795772,  1.54729571])

    assert np.isinf(radius[0])
    assert np.isinf(radius[1])
    np.testing.assert_allclose(radius[2:], puff_radius)


def test_stoich_neutrals():
    species = "H2 He H2O CO CO2 CH4".split()
    elements, stoichs = pa.stoich(species)
    assert elements == 'C H He O'.split()
    np.testing.assert_equal(
        stoichs,
        np.array([[0, 2, 0, 0],
                  [0, 0, 1, 0],
                  [0, 2, 0, 1],
                  [1, 0, 0, 1],
                  [1, 0, 0, 2],
                  [1, 4, 0, 0]]))


def test_stoich_with_electrons():
    species = "H2 He H2O CO CO2 CH4 e-".split()
    elements, stoichs = pa.stoich(species)
    assert elements == 'C H He O e-'.split()
    np.testing.assert_equal(
        stoichs,
        np.array([[0, 2, 0, 0, 0],
                  [0, 0, 1, 0, 0],
                  [0, 2, 0, 1, 0],
                  [1, 0, 0, 1, 0],
                  [1, 0, 0, 2, 0],
                  [1, 4, 0, 0, 0],
                  [0, 0, 0, 0, 1]]))


def test_stoich_with_ions():
    species = "H2 He H2O CO CO2 CH4 H2+ H2- e-".split()
    elements, stoichs = pa.stoich(species)
    assert elements == 'C H He O e-'.split()
    np.testing.assert_equal(
        stoichs,
        np.array([[0, 2, 0, 0, 0],
                  [0, 0, 1, 0, 0],
                  [0, 2, 0, 1, 0],
                  [1, 0, 0, 1, 0],
                  [1, 0, 0, 2, 0],
                  [1, 4, 0, 0, 0],
                  [0, 2, 0, 0, 0],
                  [0, 2, 0, 0, 0],
                  [0, 0, 0, 0, 1]]))


@pytest.mark.parametrize("abundances",
    [ [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4],
     [[0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]]])
def test_mean_weight_molfile(abundances):
    species = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    mu = pa.mean_weight(abundances, species)
    np.testing.assert_allclose(mu, np.array([2.31928918]))


def test_mean_weight_mass():
    abundances = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    species = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    mass = np.array([2.01588, 4.0026020, 18.01528, 28.0101, 44.0095, 16.0425])
    mu = pa.mean_weight(abundances, mass=mass)
    np.testing.assert_allclose(mu, np.array([2.31928918]))


def test_mean_weight_fail():
    abundances = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    with pytest.raises(ValueError,
            match="Either species or mass arguments must be specified"):
        mu = pa.mean_weight(abundances)


def test_ideal_gas_density():
    nlayers = 11
    pressure    = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    temperature = np.tile(1500.0, nlayers)
    species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    qprofiles = pa.uniform(pressure, temperature, species, abundances)
    dens = pa.ideal_gas_density(qprofiles, pressure, temperature)
    for i,density in enumerate(dens):
        np.testing.assert_allclose(density, expected_dens*10**i, rtol=1e-7)


def test_teq():
    tstar  = 6091.0
    rstar  = 1.19 * pc.rsun
    smaxis = 0.04747 * pc.au
    tstar_unc  = 10.0
    rstar_unc  = 0.02 * pc.rsun
    smaxis_unc = 0.00046 * pc.au
    A = 0.3
    f = 1.0
    teq, teq_unc = pa.equilibrium_temp(tstar, rstar, smaxis, A, f,
        tstar_unc, rstar_unc, smaxis_unc)
    np.testing.assert_allclose(teq, 1345.1176526125155)
    np.testing.assert_allclose(teq_unc, 13.2333537885981)


def test_teq_no_uncertainties():
    tstar  = 6091.0
    rstar  = 1.19 * pc.rsun
    smaxis = 0.04747 * pc.au
    A = 0.3
    f = 1.0
    teq, teq_unc = pa.equilibrium_temp(tstar, rstar, smaxis, A, f)
    np.testing.assert_allclose(teq, 1345.1176526125155)
    np.testing.assert_equal(teq_unc, 0.0)


def test_transit_path():
    nlayers = 5
    radius = np.linspace(5.0, 1.0, nlayers)
    path = pa.transit_path(radius)
    assert len(path) == nlayers
    expected_path =[
        np.array([]),
        np.array([3.00000000]),
        np.array([1.35424869, 2.64575131]),
        np.array([1.11847408, 1.22803364, 2.23606798]),
        np.array([1.02599614, 1.04455622, 1.09637632, 1.73205081])
    ]
    for i in range(nlayers):
        np.testing.assert_allclose(path[i], expected_path[i])


def test_transit_path_nskip():
    nlayers = 5
    radius = np.linspace(5.0, 1.0, nlayers)
    path = pa.transit_path(radius, nskip=1)
    assert len(path) == nlayers
    expected_path =[
        np.array([]),
        np.array([]),
        np.array([2.64575131]),
        np.array([1.22803364, 2.23606798]),
        np.array([1.04455622, 1.09637632, 1.73205081])
    ]
    for i in range(nlayers):
        np.testing.assert_allclose(path[i], expected_path[i])


def test_make_atomic_xsolar():
    z, symbol, dex, names, mass = pa.make_atomic(xsolar=0.1)
    np.testing.assert_allclose(dex[symbol=='H'][0],  12.0)
    np.testing.assert_allclose(dex[symbol=='He'][0], 10.93)
    np.testing.assert_allclose(dex[symbol=='C'][0],  7.43)
    np.testing.assert_allclose(dex[symbol=='N'][0],  6.83)
    np.testing.assert_allclose(dex[symbol=='O'][0],  7.69)


def test_make_atomic_escale():
    escale = {'C': 0.1, 'O':10.0}
    z, symbol, dex, names, mass = pa.make_atomic(escale=escale)
    np.testing.assert_allclose(dex[symbol=='H'][0],  12.0)
    np.testing.assert_allclose(dex[symbol=='He'][0], 10.93)
    np.testing.assert_allclose(dex[symbol=='C'][0],  7.43)
    np.testing.assert_allclose(dex[symbol=='N'][0],  7.83)
    np.testing.assert_allclose(dex[symbol=='O'][0],  9.69)


def test_make_atomic_xsolar_escale():
    escale = {'C': 0.1, 'O':10.0}
    z, symbol, dex, names, mass = pa.make_atomic(xsolar=0.1, escale=escale)
    np.testing.assert_allclose(dex[symbol=='H'][0],  12.0)
    np.testing.assert_allclose(dex[symbol=='He'][0], 10.93)
    np.testing.assert_allclose(dex[symbol=='C'][0],  6.43)
    np.testing.assert_allclose(dex[symbol=='N'][0],  6.83)
    np.testing.assert_allclose(dex[symbol=='O'][0],  8.69)


@pytest.mark.skip
def test_make_atomic_file():
    # TBD: generate file in tmp folder, assert it exists
    afile = 'sub_solar_elemental_abundance.txt'
    z, symbol, dex, names, mass = pa.make_atomic(
        xsolar=0.1, atomic_file=afile)


@pytest.mark.parametrize("qcap,qcap_result",
    [(1e-3, False),
     (1e-4, True)])
def test_qcapcheck(qcap, qcap_result):
    nlayers = 11
    pressure    = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    temperature = np.tile(1500.0, nlayers)
    species     = ["H2", "He", "H2O"]
    abundances  = [0.8495, 0.15, 5e-4]
    qprofiles = pa.uniform(pressure, temperature, species, abundances)
    ibulk = [0,1]
    assert pa.qcapcheck(qprofiles, qcap, ibulk) == qcap_result


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


def test_ratio():
    q = np.tile([0.8, 0.2], (5,1))
    q[4] = 0.5, 0.5
    ibulk = [0, 1]
    bratio, invsrat = pa.ratio(q, ibulk)
    np.testing.assert_equal(bratio[:,0], np.tile(1.0, 5))
    np.testing.assert_equal(bratio[:,1], np.array([0.25,0.25,0.25,0.25,1.0]))
    np.testing.assert_equal(invsrat, np.array([0.8, 0.8, 0.8, 0.8, 0.5]))


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


@pytest.mark.parametrize('alkali', ['sodium_vdw', 'potassium_vdw'])
def test_alkali_cutoff_default(alkali):
    na = pa.alkali.get_model('sodium_vdw')
    assert na.cutoff == 4500.0


@pytest.mark.parametrize('cutoff', [4500, 4500.0, 5000.0])
@pytest.mark.parametrize('alkali', ['sodium_vdw', 'potassium_vdw'])
def test_alkali_cutoff(cutoff, alkali):
    na = pa.alkali.get_model('sodium_vdw', cutoff)
    assert na.cutoff == float(cutoff)

