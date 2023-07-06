# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import numpy as np
import pyratbay.opacity as op
import pyratbay.spectrum as ps


# Setup:
wn_min = 1e4/5.0
wn_max = 1e4/0.3
resolution = 5.0
wn = ps.constant_resolution_spectrum(wn_min, wn_max, resolution)


expected_sigma_bf = np.array([
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 7.41444990e-18, 2.72330493e-17,
       3.77481656e-17, 3.98526176e-17, 3.70035432e-17, 3.19222307e-17,
       2.62805159e-17, 2.09559462e-17, 1.63306414e-17,
])

expected_ff_cs_2000 = np.array([
       9.04970217e-38, 6.07480788e-38, 4.08464854e-38, 2.75368934e-38,
       1.86390333e-38, 1.26918335e-38, 8.71515197e-39, 6.05038585e-39,
       4.25417694e-39, 3.02815957e-39, 2.17350879e-39, 1.56255649e-39,
       1.12062985e-39, 8.07423640e-40, 5.95740952e-40,
])

expected_bf_cs_2000 = np.array([
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 6.80626574e-37, 2.51357244e-36,
       3.49134315e-36, 3.68825733e-36, 3.42504489e-36, 2.95477929e-36,
       2.43257678e-36, 1.93972425e-36, 1.51159679e-36,
])


expected_ff_cs_temps = np.array([
       [9.04970217e-38, 6.07480788e-38, 4.08464854e-38, 2.75368934e-38,
        1.86390333e-38, 1.26918335e-38, 8.71515197e-39, 6.05038585e-39,
        4.25417694e-39, 3.02815957e-39, 2.17350879e-39, 1.56255649e-39,
        1.12062985e-39, 8.07423640e-40, 5.95740952e-40],
       [1.09567510e-37, 7.35884215e-38, 4.94828571e-38, 3.33240706e-38,
        2.24840631e-38, 1.52045738e-38, 1.03098517e-38, 7.01417239e-39,
        4.79249602e-39, 3.29302858e-39, 2.27766898e-39, 1.58218090e-39,
        1.09490082e-39, 7.65735318e-40, 5.42652301e-40],
])

expected_bf_cs_temps = np.array([
       [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
        0.00000000e+00, 0.00000000e+00, 6.80626574e-37, 2.51357244e-36,
        3.49134315e-36, 3.68825733e-36, 3.42504489e-36, 2.95477929e-36,
        2.43257678e-36, 1.93972425e-36, 1.51159679e-36],
       [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
        0.00000000e+00, 0.00000000e+00, 1.06889836e-38, 4.16055644e-38,
        6.01534661e-38, 6.53130098e-38, 6.16622007e-38, 5.36525681e-38,
        4.43333985e-38, 3.53960657e-38, 2.75927654e-38],
])

expected_ec_2000 = np.array([
       2.08126499e-12, 1.39709404e-12, 9.39394007e-13, 6.33297879e-13,
       4.28663469e-13, 2.91888817e-13, 1.58535913e-11, 5.79466889e-11,
       8.03923067e-11, 8.48927743e-11, 7.88197170e-11, 6.79904229e-11,
       5.59705624e-11, 4.46286578e-11, 3.47776459e-11,
])

expected_ec_temps = np.array([
       [7.11896048e-10, 4.78128112e-10, 3.21506353e-10, 2.16517417e-10,
        1.46086333e-10, 9.87891033e-11, 1.36436321e-10, 3.15898377e-10,
        4.21975218e-10, 4.45755934e-10, 4.15438306e-10, 3.58878281e-10,
        2.95162561e-10, 2.34955082e-10, 1.82805037e-10],
       [2.08126499e-12, 1.39709404e-12, 9.39394007e-13, 6.33297879e-13,
        4.28663469e-13, 2.91888817e-13, 1.58535913e-11, 5.79466889e-11,
        8.03923067e-11, 8.48927743e-11, 7.88197170e-11, 6.79904229e-11,
        5.59705624e-11, 4.46286578e-11, 3.47776459e-11],
])


def test_h_ion_init():
    h_ion = op.Hydrogen_Ion(wn)

    assert h_ion.name == 'H- bound-free/free-free'
    np.testing.assert_allclose(h_ion.wn, wn)
    assert h_ion.nwave == len(wn)
    np.testing.assert_allclose(h_ion.sigma_bf, expected_sigma_bf)


def test_h_ion_cross_section_single_temp():
    h_ion = op.Hydrogen_Ion(wn)
    temperature = 2000.0
    sigma_ff = h_ion.cross_section_free_free(temperature)
    sigma_bf = h_ion.cross_section_bound_free(temperature)

    np.testing.assert_allclose(sigma_ff, expected_ff_cs_2000)
    np.testing.assert_allclose(sigma_bf, expected_bf_cs_2000)


def test_h_ion_ff_cross_section_multiple_temps():
    h_ion = op.Hydrogen_Ion(wn)
    temperatures = np.array([2000.0, 5000.0])
    sigma_ff = h_ion.cross_section_free_free(temperatures)
    sigma_bf = h_ion.cross_section_bound_free(temperatures)

    np.testing.assert_allclose(sigma_ff, expected_ff_cs_temps)
    np.testing.assert_allclose(sigma_bf, expected_bf_cs_temps)


def test_h_ion_extinction_coefficients():
    h_ion = op.Hydrogen_Ion(wn)
    temperatures = np.array([5000.0, 2000.0])
    number_density = np.array([
        [1.338e+16, 4.856e+11],
        [1.724e+15, 1.334e+10],
    ])
    ec = h_ion.calc_extinction_coefficient(temperatures, number_density)
    np.testing.assert_allclose(ec, expected_ec_temps)


def test_h_ion_extinction_coefficients_single_layer():
    h_ion = op.Hydrogen_Ion(wn)
    temperature = 2000.0
    number_density = np.array([1.724e+15, 1.334e+10])
    ec = h_ion.calc_extinction_coefficient(temperature, number_density)
    np.testing.assert_allclose(ec, expected_ec_2000)


def test_h_ion_wn_order():
    # Order of wn array does not matter:
    wn1 = wn
    wn2 = wn[::-1]
    h_ion1 = op.Hydrogen_Ion(wn1)
    h_ion2 = op.Hydrogen_Ion(wn2)

    temperatures = np.array([5000.0, 2000.0])
    number_density = np.array([
        [1.338e+16, 4.856e+11],
        [1.724e+15, 1.334e+10],
    ])

    ec1 = h_ion1.calc_extinction_coefficient(temperatures, number_density)
    ec2 = h_ion2.calc_extinction_coefficient(temperatures, number_density)

    np.testing.assert_allclose(ec1, np.fliplr(ec2))
