# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os
import pytest

import numpy as np

import pyratbay.atmosphere as pa
import pyratbay.constants as pc

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
    [7.33918917e-01, 9.20325288e-02, 9.04409507e-21, 7.30152372e-05,
     1.73519070e-01, 1.65436843e-04, 2.90924973e-04, 1.07436461e-07,
     3.89720850e-28],
    [4.14510683e-01, 1.15263065e-01, 1.73061626e-20, 5.63391902e-06,
     4.69563070e-01, 2.93053885e-04, 3.64405726e-04, 8.80895544e-08,
     3.93029911e-25],
    [1.72947643e-01, 1.32825259e-01, 2.51143066e-20, 5.27388768e-07,
     6.93462849e-01, 3.43691232e-04, 4.19949733e-04, 8.06172974e-08,
     8.95256339e-23],
    [6.35708181e-02, 1.40776941e-01, 2.87461835e-20, 5.75648620e-08,
     7.94842234e-01, 3.64774280e-04, 4.45096655e-04, 7.91194531e-08,
     9.68868272e-21],
    [2.23443781e-02, 1.43774079e-01, 3.01239077e-20, 6.61311322e-09,
     8.33054288e-01, 3.72594549e-04, 4.54574809e-04, 7.87506649e-08,
     8.02644889e-19],
    [7.73364433e-03, 1.44836270e-01, 3.06131223e-20, 7.72743928e-10,
     8.46596719e-01, 3.75353825e-04, 4.57933858e-04, 7.86416635e-08,
     6.06274321e-17],
    [2.66256822e-03, 1.45204933e-01, 3.07830265e-20, 9.08164972e-11,
     8.51297010e-01, 3.76310163e-04, 4.59099710e-04, 7.86063316e-08,
     4.43633950e-15],
    [9.15014143e-04, 1.45331979e-01, 3.08415900e-20, 1.06941868e-11,
     8.52916787e-01, 3.76639574e-04, 4.59501476e-04, 7.85944482e-08,
     3.21102868e-13],
    [3.14256169e-04, 1.45375654e-01, 3.08617206e-20, 1.26015353e-12,
     8.53473619e-01, 3.76752820e-04, 4.59639568e-04, 7.85903983e-08,
     2.31546006e-11],
    [1.07906288e-04, 1.45390656e-01, 3.08683912e-20, 1.48525643e-13,
     8.53664879e-01, 3.76793354e-04, 4.59685366e-04, 7.85890735e-08,
     1.66751073e-09],
    [3.70491038e-05, 1.45395841e-01, 3.08531298e-20, 1.75125140e-14,
     8.53730403e-01, 3.76925092e-04, 4.59583460e-04, 7.85930900e-08,
     1.19967440e-07],
    [1.27202322e-05, 1.45399989e-01, 2.96606854e-20, 2.10915021e-15,
     8.53742380e-01, 3.85116022e-04, 4.51415838e-04, 7.88727911e-08,
     8.30042617e-06],
    [4.36660565e-06, 1.45460135e-01, 1.05845471e-20, 3.81020567e-16,
     8.53485481e-01, 5.90028815e-04, 2.46876771e-04, 6.61063204e-08,
     2.13045784e-04],
    [1.49892871e-06, 1.45529258e-01, 3.10890062e-22, 6.29938264e-17,
     8.53181783e-01, 8.27255485e-04, 1.01703473e-05, 3.81961453e-09,
     4.50029755e-04],
    [5.14626555e-07, 1.45532249e-01, 4.41627598e-24, 7.51582031e-18,
     8.53169719e-01, 8.37304339e-04, 1.46229494e-07, 5.55864013e-11,
     4.60067093e-04]])

q_H2O = [
    [[3.75844225e-04, 3.76505465e-04, 3.76715031e-04, 3.76781665e-04,
      3.76834538e-04, 3.79982716e-04, 5.29791757e-04, 8.23468187e-04,
      8.37304339e-04],
     [1.47151089e-16, 1.48833605e-14, 1.49374072e-12, 1.49545393e-10,
      1.49595458e-08, 1.49197537e-06, 1.18737317e-04, 7.72554290e-04,
      8.38257597e-04]],
    [[3.73594105e-03, 3.74250184e-03, 3.74458110e-03, 3.74523938e-03,
      3.74547884e-03, 3.74867460e-03, 4.02249162e-03, 7.34889074e-03,
      8.39363322e-03],
     [1.36178455e-16, 1.37736334e-14, 1.38237052e-12, 1.38395812e-10,
      1.38445690e-08, 1.38425073e-06, 1.34907833e-04, 4.65145869e-03,
      8.48049921e-03]]]


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


@pytest.mark.parametrize('metallicity', [0.0, 1.0])
@pytest.mark.parametrize('e_scale', [{}, {'C': 0.7}])
def test_abundances_tea_metallciity_escale(metallicity, e_scale):
    nlayers = 9
    punits  = 'bar'
    pressure = pa.pressure(1e-5, 1e3, nlayers, punits)
    tmodel = pa.tmodels.Isothermal(nlayers)
    temperature = tmodel(1500.0)
    species = 'H He C O H2 H2O CO CO2 CH4'.split()
    elements = 'H He C O'.split()
    q = pa.abundance(
        pressure, temperature, species, elements,
        punits=punits, metallicity=metallicity, e_scale=e_scale)
    np.testing.assert_allclose(q[:,5], q_H2O[metallicity!=0]['C' in e_scale])


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


def test_mean_weight_molfile():
    abundances = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    species = "H2 He H2O CO CO2 CH4".split()
    mu = pa.mean_weight(abundances, species)
    np.testing.assert_allclose(mu, np.array([2.31939114]))


def test_mean_weight_mass():
    # species = "H2 He H2O CO CO2 CH4".split()
    abundances = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
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


def test_chemistry_solar():
    nlayers = 100
    pressure = np.logspace(-8, 3, nlayers) * pc.bar
    temperature = np.tile(900.0, nlayers)
    species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 H2 H He'.split()
    chem_model = 'tea'
    chem_network = pa.chemistry(chem_model, pressure, temperature, species)
    np.testing.assert_equal(
        chem_network.elements,
        ['C', 'H', 'He', 'N', 'O'],
    )
    np.testing.assert_allclose(
        chem_network.element_rel_abundance,
        [2.69153480e-04, 1.0, 8.51138038e-02, 6.76082975e-05, 4.89778819e-04],
    )


def test_chemistry_metallicity():
    nlayers = 100
    pressure = np.logspace(-8, 3, nlayers) * pc.bar
    temperature = np.tile(900.0, nlayers)
    species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 H2 H He'.split()
    chem_model = 'tea'
    chem_network = pa.chemistry(
        chem_model, pressure, temperature, species, metallicity=-1.0)
    np.testing.assert_allclose(
        chem_network.element_rel_abundance,
        [2.69153480e-05, 1.0, 8.51138038e-02, 6.76082975e-06, 4.89778819e-05],
    )


def test_chemistry_escale():
    nlayers = 100
    pressure = np.logspace(-8, 3, nlayers) * pc.bar
    temperature = np.tile(900.0, nlayers)
    species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 H2 H He'.split()
    e_scale = {'C': -1.0, 'O':1.0}
    chem_model = 'tea'
    chem_network = pa.chemistry(
        chem_model, pressure, temperature, species, e_scale=e_scale)
    np.testing.assert_allclose(
        chem_network.element_rel_abundance,
        [2.69153480e-05, 1.0, 8.51138038e-02, 6.76082975e-05, 4.89778819e-03],
    )


def test_chemistry_metallicity_escale():
    nlayers = 100
    pressure = np.logspace(-8, 3, nlayers) * pc.bar
    temperature = np.tile(900.0, nlayers)
    species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 H2 H He'.split()
    e_scale = {'C': -1.0, 'O':1.0}
    chem_model = 'tea'
    chem_network = pa.chemistry(
        chem_model, pressure, temperature, species,
        metallicity=-1.0, e_scale=e_scale)
    np.testing.assert_allclose(
        chem_network.element_rel_abundance,
        [2.69153480e-06, 1.0, 8.51138038e-02, 6.76082975e-06, 4.89778819e-04],
    )


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

