# Copyright (c) 2021-2025 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import numpy as np
import pyratbay.spectrum as ps
import pyratbay.opacity as op


# Setup:
wn_min = 1e4/5.0
wn_max = 1e4/0.2
resolution = 5.0
wn = ps.constant_resolution_spectrum(wn_min, wn_max, resolution)


expected_lec_cs = np.array([
       1.27493100e-31, 2.84503350e-31, 6.34874797e-31, 1.41673554e-30,
       3.16147311e-30, 7.05488917e-30, 1.57431234e-29, 3.51310881e-29,
       7.83957112e-29, 1.74941565e-28, 3.90385528e-28, 8.71152952e-28,
       1.94399487e-27, 4.33806264e-27, 9.68047174e-27, 2.16021623e-26,
       4.82056482e-26,
])
expected_lec_ec = np.array([
       [1.62229897e-17, 3.62019193e-17, 8.07852920e-17, 1.80273961e-16,
        4.02284875e-16, 8.97706578e-16, 2.00324981e-15, 4.47029119e-15,
        9.97554234e-15, 2.22606181e-14, 4.96750053e-14, 1.10850747e-13,
        2.47365614e-13, 5.52001211e-13, 1.23180151e-12, 2.74878920e-12,
        6.13397695e-12],
       [1.62229897e-12, 3.62019193e-12, 8.07852920e-12, 1.80273961e-11,
        4.02284875e-11, 8.97706578e-11, 2.00324981e-10, 4.47029119e-10,
        9.97554234e-10, 2.22606181e-09, 4.96750053e-09, 1.10850747e-08,
        2.47365614e-08, 5.52001211e-08, 1.23180151e-07, 2.74878920e-07,
        6.13397695e-07],
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
       [7.94926497e-20, 2.64989111e-19, 8.83342410e-19, 2.94462596e-18,
        9.81592408e-18, 3.27214277e-17, 1.09077029e-16, 3.63608773e-16,
        1.21209151e-15, 4.04051257e-15, 1.34690671e-14, 4.48991967e-14,
        1.49671677e-13, 4.98931219e-13, 1.66318950e-12, 5.54424978e-12,
        1.84817820e-11],
       [7.94926497e-15, 2.64989111e-14, 8.83342410e-14, 2.94462596e-13,
        9.81592408e-13, 3.27214277e-12, 1.09077029e-11, 3.63608773e-11,
        1.21209151e-10, 4.04051257e-10, 1.34690671e-09, 4.48991967e-09,
        1.49671677e-08, 4.98931219e-08, 1.66318950e-07, 5.54424978e-07,
        1.84817820e-06],
])



def test_lecavelier_default():
    pressure = np.logspace(-4.5, 0.5, 2)
    lec_rayleigh = op.clouds.Lecavelier(pressure, wn=wn)

    assert lec_rayleigh.name == 'lecavelier'
    assert lec_rayleigh.npars == 2
    np.testing.assert_allclose(lec_rayleigh.pars, [0.0, -4.0])
    assert lec_rayleigh.pnames == ['log_k_ray', 'alpha_ray']
    texnames = [r'$\log\ \kappa_{\rm ray}$', r'$\alpha_{\rm ray}$']
    assert lec_rayleigh.texnames == texnames

    np.testing.assert_allclose(lec_rayleigh.wn, wn)
    np.testing.assert_allclose(lec_rayleigh.cross_section, expected_lec_cs)

    temperature = np.array(1800.0)
    ec = lec_rayleigh.calc_extinction_coefficient(temperature)
    np.testing.assert_allclose(ec, expected_lec_ec)


def test_lecavelier_update_parameters():
    pressure = np.logspace(-4.5, 0.5, 2)
    lec_rayleigh = op.clouds.Lecavelier(pressure, wn=wn)
    # Update cross sections only:
    pars0 = [1.0, -4.0]
    lec_rayleigh.calc_cross_section(pars0)
    np.testing.assert_allclose(lec_rayleigh.pars, pars0)
    np.testing.assert_allclose(
        lec_rayleigh.cross_section, expected_lec_cs_update_cs,
    )
    # Update cross sections and extinction coefficient:
    temperature = np.array(1800.0)
    ec = lec_rayleigh.calc_extinction_coefficient(temperature)
    pars1 = [0.0, -6.0]
    ec = lec_rayleigh.calc_extinction_coefficient(temperature, pars1)
    np.testing.assert_allclose(
        lec_rayleigh.cross_section, expected_lec_cs_update_ec,
    )
    np.testing.assert_allclose(ec, expected_lec_ec_update_ec)



