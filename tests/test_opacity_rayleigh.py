# Copyright (c) 2021-2023 Patricio Cubillos
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

expected_lec_cs = np.array([
       1.27493100e-31, 2.84503350e-31, 6.34874797e-31, 1.41673554e-30,
       3.16147311e-30, 7.05488917e-30, 1.57431234e-29, 3.51310881e-29,
       7.83957112e-29, 1.74941565e-28, 3.90385528e-28, 8.71152952e-28,
       1.94399487e-27, 4.33806264e-27, 9.68047174e-27, 2.16021623e-26,
       4.82056482e-26,
])
expected_lec_ec = np.array([
       [2.29487580e-17, 5.12106029e-17, 1.14277463e-16, 2.55012398e-16,
        5.69065160e-16, 1.26988005e-15, 2.83376220e-15, 6.32359586e-15,
        1.41112280e-14, 3.14894817e-14, 7.02693951e-14, 1.56807531e-13,
        3.49919077e-13, 7.80851275e-13, 1.74248491e-12, 3.88838921e-12,
        8.67701668e-12],
       [1.40242410e-13, 3.12953685e-13, 6.98362277e-13, 1.55840910e-12,
        3.47762042e-12, 7.76037809e-12, 1.73174357e-11, 3.86441969e-11,
        8.62352823e-11, 1.92435721e-10, 4.29424081e-10, 9.58268248e-10,
        2.13839436e-09, 4.77186890e-09, 1.06485189e-08, 2.37623785e-08,
        5.30262130e-08],
])

expected_lec_cs_update_cs = np.array([
       1.27493100e-30, 2.84503350e-30, 6.34874797e-30, 1.41673554e-29,
       3.16147311e-29, 7.05488917e-29, 1.57431234e-28, 3.51310881e-28,
       7.83957112e-28, 1.74941565e-27, 3.90385528e-27, 8.71152952e-27,
       1.94399487e-26, 4.33806264e-26, 9.68047174e-26, 2.16021623e-25,
       4.82056482e-25,
])
expected_lec_cs_update_ec = np.array([
       6.24716190e-34, 2.08249427e-33, 6.94200416e-33, 2.31412026e-32,
       7.71413044e-32, 2.57150891e-31, 8.57213668e-31, 2.85752568e-30,
       9.52557489e-30, 3.17535474e-29, 1.05850595e-28, 3.52853442e-28,
       1.17623856e-27, 3.92099661e-27, 1.30706601e-26, 4.35711052e-26,
       1.45244478e-25,
])
expected_lec_ec_update_ec = np.array([
       [1.12448914e-19, 3.74848969e-19, 1.24956075e-18, 4.16541646e-18,
        1.38854348e-17, 4.62871603e-17, 1.54298460e-16, 5.14354622e-16,
        1.71460348e-15, 5.71563853e-15, 1.90531071e-14, 6.35136196e-14,
        2.11722941e-13, 7.05779390e-13, 2.35271882e-12, 7.84279893e-12,
        2.61440060e-11],
       [6.87187809e-16, 2.29074370e-15, 7.63620458e-15, 2.54553228e-14,
        8.48554349e-14, 2.82865980e-13, 9.42935035e-13, 3.14327824e-12,
        1.04781324e-11, 3.49289021e-11, 1.16435655e-10, 3.88138787e-10,
        1.29386242e-09, 4.31309628e-09, 1.43777261e-08, 4.79282157e-08,
        1.59768926e-07],
])

@pytest.mark.parametrize('species', ['H2', 'H', 'He'])
def test_dalgarno_init(species):
    rayleigh = op.rayleigh.Dalgarno(wn=wn, species=species)

    assert rayleigh.name == 'dalgarno_' + species
    assert rayleigh.species == species
    assert rayleigh.npars == 0
    assert rayleigh.pars == []
    assert rayleigh.pnames == []
    assert rayleigh.texnames == []
    np.testing.assert_allclose(rayleigh.wn, wn)


def test_dalgarno_H2():
    H2_rayleigh = op.rayleigh.Dalgarno(wn=wn, species='H2')
    np.testing.assert_allclose(H2_rayleigh.cross_section, expected_H2_cs)

    density = np.array([1.8e+14, 1.1e+18])
    ec = H2_rayleigh.calc_extinction_coefficient(density)
    np.testing.assert_allclose(ec, expected_H2_ec)


def test_dalgarno_H():
    H_rayleigh = op.rayleigh.Dalgarno(wn=wn, species='H')
    np.testing.assert_allclose(H_rayleigh.cross_section, expected_H_cs)

    density = np.array([1.8e+14, 1.1e+18])
    ec = H_rayleigh.calc_extinction_coefficient(density)
    np.testing.assert_allclose(ec, expected_H_ec)


def test_dalgarno_He():
    He_rayleigh = op.rayleigh.Dalgarno(wn=wn, species='He')
    np.testing.assert_allclose(He_rayleigh.cross_section, expected_He_cs)

    density = np.array([1.8e+14, 1.1e+18])
    ec = He_rayleigh.calc_extinction_coefficient(density)
    np.testing.assert_allclose(ec, expected_He_ec)



def test_lecavelier_default():
    lec_rayleigh = op.rayleigh.Lecavelier(wn=wn)

    assert lec_rayleigh.name == 'lecavelier'
    assert lec_rayleigh.species == 'H2'
    assert lec_rayleigh.npars == 2
    np.testing.assert_allclose(lec_rayleigh.pars, [0.0, -4.0])
    assert lec_rayleigh.pnames == ['log_k_ray', 'alpha_ray']
    assert lec_rayleigh.texnames == [r'$\log\ \kappa_{\rm ray}$', r'$\alpha_{\rm ray}$']
    np.testing.assert_allclose(lec_rayleigh.wn, wn)
    np.testing.assert_allclose(lec_rayleigh.cross_section, expected_lec_cs)

    density = np.array([1.8e+14, 1.1e+18])
    ec = lec_rayleigh.calc_extinction_coefficient(density)
    np.testing.assert_allclose(ec, expected_lec_ec)


def test_lecavelier_update_parameters():
    lec_rayleigh = op.rayleigh.Lecavelier(wn=wn)
    # Update cross sections only:
    pars0 = [1.0, -4.0]
    lec_rayleigh.calc_cross_section(pars0)
    np.testing.assert_allclose(lec_rayleigh.pars, pars0)
    np.testing.assert_allclose(
        lec_rayleigh.cross_section, expected_lec_cs_update_cs,
    )
    # Update cross sections and extinction coefficient:
    density = np.array([1.8e+14, 1.1e+18])
    pars1 = [0.0, -6.0]
    ec = lec_rayleigh.calc_extinction_coefficient(density, pars1)
    np.testing.assert_allclose(
        lec_rayleigh.cross_section, expected_lec_cs_update_ec,
    )
    np.testing.assert_allclose(ec, expected_lec_ec_update_ec)



