# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import pytest
import re

import numpy as np

import pyratbay.atmosphere as pa
import pyratbay.constants as pc

os.chdir(pc.ROOT+'tests')


expected_pressure = np.logspace(-8, 3, 15)

expected_temp_guillot = np.array(
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

expected_vmr_tea = np.array([
       [7.35957328e-01, 8.90269702e-02, 1.06589477e-20, 6.66282723e-05,
        1.74484284e-01, 1.51805109e-04, 3.12878172e-04, 1.05436507e-07,
        4.64430853e-28],
       [2.19474062e-01, 1.25416170e-01, 2.78559084e-20, 8.51833917e-07,
        6.54361069e-01, 3.06933225e-04, 4.40833548e-04, 8.00915609e-08,
        3.03561249e-23],
       [3.79045613e-02, 1.38205640e-01, 3.49420852e-20, 1.77461950e-08,
        8.23064740e-01, 3.39164796e-04, 4.85798817e-04, 7.75387678e-08,
        1.07129952e-19],
       [5.94139822e-03, 1.40457056e-01, 3.62009876e-20, 4.12812431e-10,
        8.52763045e-01, 3.44709029e-04, 4.93714144e-04, 7.73010740e-08,
        2.11870876e-16],
       [9.17432396e-04, 1.40810933e-01, 3.63990999e-20, 9.76057843e-12,
        8.57431021e-01, 3.45578144e-04, 4.94958268e-04, 7.72682954e-08,
        3.82986201e-13],
       [1.41337243e-04, 1.40865600e-01, 3.64295874e-20, 2.31355566e-13,
        8.58152122e-01, 3.45713033e-04, 4.95149777e-04, 7.72633851e-08,
        6.82773904e-10],
       [2.17662862e-05, 1.40874362e-01, 3.62190415e-20, 5.50509324e-15,
        8.58261672e-01, 3.46941233e-04, 4.93973914e-04, 7.73438678e-08,
        1.20745370e-06],
       [3.35103135e-06, 1.40971226e-01, 5.74741400e-21, 2.58401370e-16,
        8.57843429e-01, 6.86395395e-04, 1.55156394e-04, 4.80863353e-08,
        3.40394582e-04],
       [5.15975068e-07, 1.41015154e-01, 4.70780082e-24, 7.51632062e-18,
        8.57646822e-01, 8.41754217e-04, 1.55892694e-07, 5.92636304e-11,
        4.95597544e-04]])

expected_vmr_tea_H2O = [
      [[3.44846449e-04, 3.45454768e-04, 3.45647563e-04, 3.45708952e-04,
        3.45766288e-04, 3.49506982e-04, 5.15325967e-04, 8.27013608e-04,
        8.41754217e-04],
       [1.34318585e-16, 1.35854196e-14, 1.36347433e-12, 1.36503787e-10,
        1.36549953e-08, 1.36234049e-06, 1.11005412e-04, 7.73010902e-04,
        8.42796258e-04],],
      [[3.42626988e-03, 3.43230384e-03, 3.43421613e-03, 3.43482164e-03,
        3.43505041e-03, 3.43884079e-03, 3.75750334e-03, 7.32956730e-03,
        8.44006218e-03],
       [5.89500056e-03, 5.90536636e-03, 5.90865151e-03, 5.90969120e-03,
        5.91003098e-03, 5.91123064e-03, 6.01443782e-03, 7.84463736e-03,
        8.42759310e-03],]
]


def test_pressure_default_units():
    ptop = 1e-8
    pbottom = 1e3
    nlayers = 15
    pressure = pa.pressure(ptop, pbottom, nlayers)
    np.testing.assert_allclose(pressure, expected_pressure)


def test_pressure_floats():
    ptop = 1e-8
    pbottom = 1e3
    nlayers = 15
    units = "bar"
    pressure = pa.pressure(ptop, pbottom, nlayers, units)
    np.testing.assert_allclose(pressure, expected_pressure)


def test_pressure_with_units():
    ptop = "1e-8 bar"
    pbottom = "1e3 bar"
    nlayers = 15
    pressure = pa.pressure(ptop, pbottom, nlayers)
    np.testing.assert_allclose(pressure, expected_pressure)


@pytest.mark.parametrize(
    "params",
    [
        1500.0,
        [1500.0],
        (1500,),
        np.array([1500.0])
    ],
)
def test_tmodel_isothermal(params):
    nlayers = 100
    pressure = np.logspace(-8, 2, nlayers)
    tmodel = pa.tmodels.Isothermal(pressure)
    np.testing.assert_equal(tmodel(params), np.tile(1500.0, nlayers))


@pytest.mark.parametrize("params",
    [[-4.84, -0.8, -0.8, 0.5, 1200.0, 100.0],
     np.array([-4.84, -0.8, -0.8, 0.5, 1200.0, 100.0])
    ])
def test_tmodel_guillot_floats(params):
    pressure = expected_pressure
    tmodel = pa.tmodels.TCEA(pressure)
    np.testing.assert_allclose(tmodel(params), expected_temp_guillot)


@pytest.mark.parametrize('gravity',
    [None, 2200.0, np.tile(2200.0, len(expected_pressure))])
def test_tmodel_guillot_gravity(gravity):
    params = np.array([-4.84, -0.8, -0.8, 0.5, 1200.0, 100.0])
    pressure = expected_pressure
    if gravity is not None:
        params[0] += np.log10(2200.0)
    tmodel = pa.tmodels.TCEA(pressure, gravity)
    np.testing.assert_allclose(tmodel(params), expected_temp_guillot)


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
    pressure = expected_pressure
    nlayers = len(pressure)
    temp = pa.temperature("isothermal", pressure, params=params)
    np.testing.assert_equal(temp, np.tile(params, nlayers))


def test_temperature_guillot():
    pressure = expected_pressure
    params = [-4.84, -0.8, -0.8, 0.5, 1200.0, 100.0]
    temp = pa.temperature('guillot', pressure, params=params)
    np.testing.assert_allclose(temp, expected_temp_guillot)


def test_temperature_madhu():
    pressure = expected_pressure
    params = -0.77, -3.61, 1.45, 0.85, 0.67, 870.0
    temp = pa.temperature('madhu', pressure, params=params)
    np.testing.assert_allclose(temp, expected_temp_madhu_noinv)


def test_uniform():
    nlayers = 11
    #species = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    vmr = pa.uniform(abundances, nlayers)
    assert np.shape(vmr) == (nlayers, len(abundances))
    for q in vmr:
        np.testing.assert_equal(q, np.array(abundances))


def test_uniform_error():
    nlayers = -11
    abundances = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    match = "The number of layers has to be larger than zero"
    with pytest.raises(ValueError, match=match):
        pa.uniform(abundances, nlayers)


def test_hydro_g():
    nlayers = 11
    pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    tmodel = pa.tmodels.Isothermal(pressure)
    temperature = tmodel(1500.0)
    mu = np.tile(2.3, nlayers)
    g = pc.G * pc.mjup / pc.rjup**2
    r0 = 1.0 * pc.rjup
    p0 = 1.0  # bar
    # Radius profile in Jupiter radii:
    radius = pa.hydro_g(pressure, temperature, mu, g, p0, r0) / pc.rjup
    np.testing.assert_allclose(radius, radius_g)


def test_hydro_m():
    nlayers = 11
    pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    tmodel = pa.tmodels.Isothermal(pressure)
    temperature = tmodel(1500.0)
    mu = np.tile(2.3, nlayers)
    Mp = 1.0 * pc.mjup
    r0 = 1.0 * pc.rjup
    p0 = 1.0  # bar
    radius = pa.hydro_m(pressure, temperature, mu, Mp, p0, r0) / pc.rjup
    # Radius profile in Jupiter radii:
    np.testing.assert_allclose(radius, radius_m)


def test_hydro_m_ultra_puff():
    nlayers = 15
    pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    tmodel = pa.tmodels.Isothermal(pressure)
    temperature = tmodel(1500.0)
    mu = np.tile(2.3, nlayers)
    Mp = 0.1 * pc.mjup
    r0 = 2.0 * pc.rjup
    p0 = 1.0  # bar
    radius = pa.hydro_m(pressure, temperature, mu, Mp, p0, r0) / pc.rjup

    puff_radius = np.array([
       23.59979187, 10.78753812,  6.99174271,  5.17191017,  4.10376853,
        3.40130545,  2.90418177,  2.5338435 ,  2.24727349,  2.01893776,
        1.83272274,  1.67795772,  1.54729571,
    ])

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
    match = "Either species or mass arguments must be specified"
    with pytest.raises(ValueError, match=match):
        pa.mean_weight(abundances)


def test_ideal_gas_density_2D():
    nlayers = 11
    pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    temperature = np.tile(1500.0, nlayers)
    species = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    vmr = pa.uniform(abundances, nlayers)
    dens = pa.ideal_gas_density(vmr, pressure, temperature)
    for i,density in enumerate(dens):
        np.testing.assert_allclose(density, expected_dens*10**i, rtol=1e-7)


def test_ideal_gas_density_1D():
    nlayers = 11
    pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    temperature = np.tile(1500.0, nlayers)
    vmr = np.tile(0.8496, nlayers)
    density = pa.ideal_gas_density(vmr, pressure, temperature)
    expected_density = np.array([
       4.1024185e+10, 4.1024185e+11, 4.1024185e+12, 4.1024185e+13,
       4.1024185e+14, 4.1024185e+15, 4.1024185e+16, 4.1024185e+17,
       4.1024185e+18, 4.1024185e+19, 4.1024185e+20,
    ])
    np.testing.assert_allclose(density, expected_density, rtol=1e-7)


def test_teq():
    tstar = 6091.0
    rstar = 1.19 * pc.rsun
    smaxis = 0.04747 * pc.au
    tstar_unc = 10.0
    rstar_unc = 0.02 * pc.rsun
    smaxis_unc = 0.00046 * pc.au
    A = 0.3
    f = 1.0
    teq, teq_unc = pa.equilibrium_temp(tstar, rstar, smaxis, A, f,
        tstar_unc, rstar_unc, smaxis_unc)
    np.testing.assert_allclose(teq, 1345.1176526125155)
    np.testing.assert_allclose(teq_unc, 13.2333537885981)


def test_teq_no_uncertainties():
    tstar = 6091.0
    rstar = 1.19 * pc.rsun
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


def test_chemistry_uniform():
    nlayers = 11
    pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    tmodel = pa.tmodels.Isothermal(pressure)
    temperature = tmodel(1500.0)
    species = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    chem = pa.chemistry(
        'uniform', pressure, temperature, species, q_uniform=abundances,
    )
    vmr = chem.vmr

    assert np.shape(vmr) == (nlayers, len(species))
    for q in vmr:
        np.testing.assert_equal(q, np.array(abundances))


def test_chemistry_tea_basic():
    nlayers = 9
    pressure = pa.pressure(1e-10, 1e3, nlayers, units='bar')
    tmodel = pa.tmodels.Isothermal(pressure)
    temperature = tmodel(1500.0)
    species = 'H He C O H2 H2O CO CO2 CH4'.split()
    vmr = pa.chemistry('tea', pressure, temperature, species).vmr
    np.testing.assert_allclose(vmr, expected_vmr_tea)


def test_chemistry_tea_solar():
    nlayers = 100
    pressure = pa.pressure(1.0e-08, 1.0e+03, nlayers, units='bar')
    temperature = np.tile(900.0, nlayers)
    species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 H2 H He'.split()
    chem_model = 'tea'
    chem_network = pa.chemistry(chem_model, pressure, temperature, species)

    expected_elements = 'C H He N O'.split()
    np.testing.assert_equal(chem_network.elements, expected_elements)

    expected_rel_abundance = np.array(
        [2.88403150e-04, 1.0, 8.20351544e-02, 6.76082975e-05, 4.89778819e-04])
    np.testing.assert_allclose(
        chem_network.element_rel_abundance, expected_rel_abundance)


def test_chemistry_metallicity():
    nlayers = 100
    pressure = pa.pressure(1.0e-08, 1.0e+03, nlayers, units='bar')
    temperature = np.tile(900.0, nlayers)
    species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 H2 H He'.split()
    chem_model = 'tea'
    chem_network = pa.chemistry(
        chem_model, pressure, temperature, species, metallicity=-1.0,
    )
    expected_rel_abundance = np.array(
        [2.88403150e-05, 1.0, 8.20351544e-02, 6.76082975e-06, 4.89778819e-05],
    )
    e_abundance = chem_network.element_rel_abundance
    np.testing.assert_allclose(e_abundance, expected_rel_abundance)


def test_chemistry_escale():
    nlayers = 100
    pressure = pa.pressure(1.0e-08, 1.0e+03, nlayers, units='bar')
    temperature = np.tile(900.0, nlayers)
    species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 H2 H He'.split()
    e_scale = {'C': -1.0, 'O': 1.0}
    chem_model = 'tea'
    chem_network = pa.chemistry(
        chem_model, pressure, temperature, species, e_scale=e_scale)

    expected_rel_abundance = np.array(
        [2.88403150e-05, 1.0, 8.20351544e-02, 6.76082975e-05, 4.89778819e-03])
    np.testing.assert_allclose(
        chem_network.element_rel_abundance, expected_rel_abundance,
    )


def test_chemistry_eratio():
    nlayers = 100
    pressure = pa.pressure(1.0e-08, 1.0e+03, nlayers, units='bar')
    temperature = pa.tmodels.Isothermal(pressure)(900.0)
    species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 H2 H He'.split()
    e_ratio = {'C_O': 2.0}
    chem_model = 'tea'
    chem_network = pa.chemistry(
        chem_model, pressure, temperature, species, e_ratio=e_ratio)

    expected_rel_abundance = np.array(
        [9.79557639e-04, 1.0, 8.20351544e-02, 6.76082975e-05, 4.89778819e-04])
    np.testing.assert_allclose(
        chem_network.element_rel_abundance, expected_rel_abundance,
    )


def test_chemistry_metallicity_escale():
    nlayers = 100
    pressure = pa.pressure(1.0e-08, 1.0e+03, nlayers, units='bar')
    temperature = np.tile(900.0, nlayers)
    species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 H2 H He'.split()
    metallicity = -1.0
    e_scale = {'C': -1.0, 'O': 1.0}
    chem_model = 'tea'
    chem_network = pa.chemistry(
        chem_model, pressure, temperature, species,
        metallicity=metallicity,
        e_scale=e_scale,
    )

    expected_rel_abundance = np.array(
        [2.88403150e-05, 1.0, 8.20351544e-02, 6.76082975e-06, 4.89778819e-03])
    np.testing.assert_allclose(
        chem_network.element_rel_abundance, expected_rel_abundance,
    )


@pytest.mark.parametrize('e_ratio', [{}, {'C_O': 2.9512092}])
def test_chemistry_tea_metallicity_eratio(e_ratio):
    nlayers = 9
    pressure = pa.pressure(1e-5, 1e3, nlayers, units='bar')
    tmodel = pa.tmodels.Isothermal(pressure)
    temperature = tmodel(1500.0)
    species = 'H He C O H2 H2O CO CO2 CH4'.split()
    i_H2O = species.index('H2O')
    chem_model = pa.chemistry(
        'tea', pressure, temperature, species,
        e_ratio=e_ratio,
    )
    vmr = chem_model.vmr
    # (this C/O ratio leads to same composition as e_scale C=0.7)
    expected_vmr = expected_vmr_tea_H2O[0]['C_O' in e_ratio]
    np.testing.assert_allclose(vmr[:,i_H2O], expected_vmr)


@pytest.mark.parametrize('e_abundances', [{}, {'C': 9.16}])
def test_chemistry_tea_metallicity_eabundances(e_abundances):
    nlayers = 9
    pressure = pa.pressure(1e-5, 1e3, nlayers, units='bar')
    tmodel = pa.tmodels.Isothermal(pressure)
    temperature = tmodel(1500.0)
    species = 'H He C O H2 H2O CO CO2 CH4'.split()
    i_H2O = species.index('H2O')
    chem_model = pa.chemistry(
        'tea', pressure, temperature, species,
        e_abundances=e_abundances,
    )
    vmr = chem_model.vmr
    expected_vmr = expected_vmr_tea_H2O[0]['C' in e_abundances]
    np.testing.assert_allclose(vmr[:,i_H2O], expected_vmr)


def test_chemistry_mismatch_nspecies():
    nlayers = 11
    pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    tmodel = pa.tmodels.Isothermal(pressure)
    temperature = tmodel(1500.0)
    species = ["H2", "He", "H2O", "CO", "CO2"]
    abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    match = "Species (5) and q_uniform array lengths (6) don't match"
    with pytest.raises(ValueError, match=re.escape(match)):
        pa.chemistry(
            'uniform', pressure, temperature, species, q_uniform=abundances,
        )


def test_chemistry_mismatch_nlayers():
    nlayers = 11
    pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    tmodel = pa.tmodels.Isothermal(pressure[:-1])
    temperature = tmodel(1500.0)
    species = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    #abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    match = "pressure (11) and temperature array lengths (10) don't match"
    with pytest.raises(ValueError, match=re.escape(match)):
        pa.chemistry('uniform', pressure, temperature, species)


@pytest.mark.parametrize(
    "qcap,qcap_result",
    [
        (1e-3, False),
        (1e-4, True),
    ],
)
def test_qcapcheck(qcap, qcap_result):
    nlayers = 11
    #species = ["H2", "He", "H2O"]
    abundances = [0.8495, 0.15, 5e-4]
    vmr = pa.uniform(abundances, nlayers)
    ibulk = [0,1]
    assert pa.qcapcheck(vmr, qcap, ibulk) == qcap_result


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


@pytest.mark.skip(reason="TBD: add electrons, should not count as 'metals'")
def test_balance_electrons():
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

