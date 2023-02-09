# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import numpy as np
import pyratbay.opacity as op


expected_bf = np.array([
   5.23741227e-18, 7.39404466e-18, 1.47265849e-17, 2.71485871e-17,
   3.93833584e-17, 2.29099293e-17, 0.00000000e+00, 0.00000000e+00,
   0.00000000e+00, 0.00000000e+00,
])

expected_ff_2000 = np.array([[
    0.00000000e+00, 0.00000000e+00, 5.32586527e-40, 1.17986090e-39,
    2.74411437e-39, 6.63188764e-39, 1.72807983e-38, 4.70045614e-38,
    1.29885863e-37, 3.60666359e-37,
]])

expected_ff_temps = np.array([
    [0.00000000e+00, 0.00000000e+00, 5.47577416e-40, 1.19406230e-39,
     2.72282299e-39, 6.33440250e-39, 1.60636115e-38, 4.33188112e-38,
     1.19526417e-37, 3.31925288e-37],
    [0.00000000e+00, 0.00000000e+00, 5.32586527e-40, 1.17986090e-39,
     2.74411437e-39, 6.63188764e-39, 1.72807983e-38, 4.70045614e-38,
     1.29885863e-37, 3.60666359e-37],
    [0.00000000e+00, 0.00000000e+00, 4.76039649e-40, 1.15975204e-39,
     2.95142126e-39, 7.73420325e-39, 2.08247622e-38, 5.69468078e-38,
     1.57162924e-37, 4.35788070e-37],
])

expected_ec = np.array([
    [4.18992981e-12, 5.91523573e-12, 1.26573918e-11, 2.36293693e-11,
     3.58632035e-11, 2.84629875e-11, 2.57017784e-11, 6.93100979e-11,
     1.91242268e-10, 5.31080461e-10],
    [2.61870613e-15, 3.69702233e-15, 8.42846551e-15, 1.59340153e-14,
     2.51799079e-14, 2.47187399e-14, 3.45615967e-14, 9.40091228e-14,
     2.59771725e-13, 7.21332717e-13],
])

def test_h_ion_init_no_species():
    wl = np.logspace(-1.0, 1.0, 10)
    h_ion = op.Hydrogen_Ion_Opacity(wl)

    assert not h_ion.has_opacity
    assert h_ion.ec is None
    np.testing.assert_allclose(h_ion.wl, wl)
    # BF cross section already exists (depends only on wl)
    np.testing.assert_allclose(h_ion.cross_section_bf, expected_bf)


def test_h_ion_init_with_species():
    wl = np.logspace(-1.0, 1.0, 10)
    species = 'H H2 H- e-'.split()
    h_ion = op.Hydrogen_Ion_Opacity(wl, species)

    # has_opacity is turned on:
    assert h_ion.has_opacity
    assert h_ion.ec is None
    np.testing.assert_allclose(h_ion.wl, wl)
    np.testing.assert_allclose(h_ion.cross_section_bf, expected_bf)


def test_h_ion_init_without_needed_species():
    wl = np.logspace(-1.0, 1.0, 10)
    species = 'H H2 O H2O'.split()
    h_ion = op.Hydrogen_Ion_Opacity(wl, species)
    # has_opacity is turned off (needs all H, H-, and e-):
    assert not h_ion.has_opacity


def test_h_ion_ff_cross_section_single_temp():
    wl = np.logspace(-1.0, 1.0, 10)
    h_ion = op.Hydrogen_Ion_Opacity(wl)
    temperature = 2000.0
    sigma_ff = h_ion.free_free_cross_section(temperature)

    ntemp = 1
    nwave = len(wl)
    # Note, output is 2D even though temperature is a scalar
    assert np.shape(sigma_ff) == (ntemp,nwave)
    np.testing.assert_allclose(sigma_ff, expected_ff_2000)


def test_h_ion_ff_cross_section_multiple_temps():
    wl = np.logspace(-1.0, 1.0, 10)
    h_ion = op.Hydrogen_Ion_Opacity(wl)
    temperatures = np.array([1500.0, 2000.0, 5000.0])
    sigma_ff = h_ion.free_free_cross_section(temperatures)

    ntemp = len(temperatures)
    nwave = len(wl)
    assert np.shape(sigma_ff) == (ntemp,nwave)
    np.testing.assert_allclose(sigma_ff, expected_ff_temps)


def test_h_ion_ff_extinction_coefficients():
    wl = np.logspace(-1.0, 1.0, 10)
    species = 'H H2 H- e-'.split()
    h_ion = op.Hydrogen_Ion_Opacity(wl, species)
    temperatures = np.array([1500.0, 2000.0])
    number_density = np.array([
        [2.0e16, 2.0e17, 8.0e05, 8.0e10],
        [1.0e14, 2.0e16, 5.0e02, 2.0e10],
    ])
    h_ion.absorption(temperatures, number_density)

    assert h_ion.ec is not None
    ntemp = len(temperatures)
    nwave = len(wl)
    assert np.shape(h_ion.ec) == (ntemp,nwave)
    np.testing.assert_allclose(h_ion.ec, expected_ec)


