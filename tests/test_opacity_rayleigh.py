# Copyright (c) 2021-2025 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import pytest
import numpy as np
import pyratbay.spectrum as ps
import pyratbay.opacity as op


# Setup:
wn_min = 1e4/5.0
wn_max = 1e4/0.2
resolution = 5.0
wn = ps.constant_resolution_spectrum(wn_min, wn_max, resolution)


expected_H2_cs = np.array([
       1.30321961e-31, 2.90906397e-31, 6.49464817e-31, 1.45029942e-30,
       3.23972957e-30, 7.24073694e-30, 1.61953527e-29, 3.62658455e-29,
       8.13494717e-29, 1.82952260e-28, 4.13063362e-28, 9.38116916e-28,
       2.14970374e-27, 4.99341283e-27, 1.18406394e-26, 2.89644160e-26,
       7.41827423e-26,
])
expected_H2_ec = np.array([
       [2.34579530e-17, 5.23631514e-17, 1.16903667e-16, 2.61053896e-16,
        5.83151323e-16, 1.30333265e-15, 2.91516349e-15, 6.52785219e-15,
        1.46429049e-14, 3.29314067e-14, 7.43514052e-14, 1.68861045e-13,
        3.86946674e-13, 8.98814310e-13, 2.13131509e-12, 5.21359488e-12,
        1.33528936e-11],
       [1.43354157e-13, 3.19997036e-13, 7.14411298e-13, 1.59532936e-12,
        3.56370253e-12, 7.96481063e-12, 1.78148880e-11, 3.98924301e-11,
        8.94844188e-11, 2.01247486e-10, 4.54369698e-10, 1.03192861e-09,
        2.36467412e-09, 5.49275412e-09, 1.30247033e-08, 3.18608576e-08,
        8.16010165e-08],
])

expected_H_cs = np.array([
       9.28750793e-32, 2.07353044e-31, 4.63047669e-31, 1.03441924e-30,
       2.31206598e-30, 5.17192490e-30, 1.15831334e-29, 2.59885893e-29,
       5.84676534e-29, 1.32074968e-28, 3.00194243e-28, 6.88718557e-28,
       1.60267928e-27, 3.81077623e-27, 9.36004457e-27, 2.41146346e-26,
       6.64300658e-26,
])
expected_H_ec = np.array([
       [1.67175143e-17, 3.73235479e-17, 8.33485804e-17, 1.86195463e-16,
        4.16171876e-16, 9.30946483e-16, 2.08496402e-15, 4.67794607e-15,
        1.05241776e-14, 2.37734942e-14, 5.40349637e-14, 1.23969340e-13,
        2.88482270e-13, 6.85939722e-13, 1.68480802e-12, 4.34063423e-12,
        1.19574118e-11],
       [1.02162587e-13, 2.28088348e-13, 5.09352436e-13, 1.13786116e-12,
        2.54327257e-12, 5.68911739e-12, 1.27414468e-11, 2.85874482e-11,
        6.43144188e-11, 1.45282465e-10, 3.30213667e-10, 7.57590413e-10,
        1.76294720e-09, 4.19185386e-09, 1.02960490e-08, 2.65260981e-08,
        7.30730724e-08],
])

expected_e_ec = np.array([
       [1.19754e-10, 1.19754e-10, 1.19754e-10, 1.19754e-10, 1.19754e-10,
        1.19754e-10, 1.19754e-10, 1.19754e-10, 1.19754e-10, 1.19754e-10,
        1.19754e-10, 1.19754e-10, 1.19754e-10, 1.19754e-10, 1.19754e-10,
        1.19754e-10, 1.19754e-10],
       [7.31830e-07, 7.31830e-07, 7.31830e-07, 7.31830e-07, 7.31830e-07,
        7.31830e-07, 7.31830e-07, 7.31830e-07, 7.31830e-07, 7.31830e-07,
        7.31830e-07, 7.31830e-07, 7.31830e-07, 7.31830e-07, 7.31830e-07,
        7.31830e-07, 7.31830e-07],
])


expected_He_cs = np.array([
       8.77611285e-33, 1.95859554e-32, 4.37127359e-32, 9.75668058e-32,
       2.17792182e-31, 4.86240727e-31, 1.08583333e-30, 2.42565173e-30,
       5.42154014e-30, 1.21271392e-29, 2.71583174e-29, 6.09262745e-29,
       1.37034966e-28, 3.09403355e-28, 7.02551988e-28, 1.60857285e-27,
       3.72779655e-27,
])

expected_He_ec = np.array([
       [1.57970031e-18, 3.52547197e-18, 7.86829247e-18, 1.75620250e-17,
        3.92025928e-17, 8.75233308e-17, 1.95450000e-16, 4.36617311e-16,
        9.75877226e-16, 2.18288506e-15, 4.88849713e-15, 1.09667294e-14,
        2.46662939e-14, 5.56926038e-14, 1.26459358e-13, 2.89543112e-13,
        6.71003379e-13],
       [9.65372413e-15, 2.15445509e-14, 4.80840095e-14, 1.07323486e-13,
        2.39571400e-13, 5.34864799e-13, 1.19441666e-12, 2.66821690e-12,
        5.96369416e-12, 1.33398532e-11, 2.98741491e-11, 6.70189019e-11,
        1.50738463e-10, 3.40343690e-10, 7.72807187e-10, 1.76943013e-09,
        4.10057620e-09],
])


@pytest.mark.parametrize('species', ['H2', 'H', 'He', 'e-'])
def test_rayleigh_init(species):
    rayleigh = op.rayleigh.Kurucz(wn=wn, species=species)

    assert rayleigh.name == 'rayleigh_' + species
    assert rayleigh.species == species
    assert rayleigh.npars == 0
    assert rayleigh.pars == []
    assert rayleigh.pnames == []
    assert rayleigh.texnames == []
    np.testing.assert_allclose(rayleigh.wn, wn)


def test_rayleigh_H2():
    H2_rayleigh = op.rayleigh.Kurucz(wn=wn, species='H2')
    np.testing.assert_allclose(H2_rayleigh.cross_section, expected_H2_cs)

    density = np.array([1.8e+14, 1.1e+18])
    ec = H2_rayleigh.calc_extinction_coefficient(density)
    np.testing.assert_allclose(ec, expected_H2_ec)


def test_rayleigh_H():
    H_rayleigh = op.rayleigh.Kurucz(wn=wn, species='H')
    np.testing.assert_allclose(H_rayleigh.cross_section, expected_H_cs)

    density = np.array([1.8e+14, 1.1e+18])
    ec = H_rayleigh.calc_extinction_coefficient(density)
    np.testing.assert_allclose(ec, expected_H_ec)


def test_rayleigh_He():
    He_rayleigh = op.rayleigh.Kurucz(wn=wn, species='He')
    np.testing.assert_allclose(He_rayleigh.cross_section, expected_He_cs)

    density = np.array([1.8e+14, 1.1e+18])
    ec = He_rayleigh.calc_extinction_coefficient(density)
    np.testing.assert_allclose(ec, expected_He_ec)


def test_rayleigh_e():
    rayleigh = op.rayleigh.Kurucz(wn=wn, species='e-')
    np.testing.assert_allclose(rayleigh.cross_section, 6.653e-25)

    density = np.array([1.8e+14, 1.1e+18])
    ec = rayleigh.calc_extinction_coefficient(density)
    np.testing.assert_allclose(ec, expected_e_ec)
