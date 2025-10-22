# Copyright (c) 2021-2025 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import pytest

import numpy as np

import pyratbay.opacity.broadening as pb
import pyratbay.constants as pc


os.chdir(pc.ROOT + 'tests')

keys = [
    'lorentz', 'gauss',
    'voigt0.01', 'voigt0.1', 'voigt1.0', 'voigt10.0', 'voigt100.0',
]
expected = {
    key:np.load(f"expected/expected_profile_{key}_test.npz")['arr_0']
    for key in keys
}


def test_Lorentz_hwhm():
    lor = pb.Lorentz(x0=0.0, hwhm=1.0, scale=1.0)
    x = np.linspace(-10.0, 10.0, 100001)
    np.testing.assert_approx_equal(
        0.5*np.ptp(x[lor(x)>0.5*np.amax(lor(x))]), 1.0, 3)


def test_Lorentz_integral():
    lor = pb.Lorentz(x0=0.0, hwhm=1.0, scale=1.0)
    x = np.linspace(-1000.0, 1000.0, 100001)
    np.testing.assert_approx_equal(np.trapezoid(lor(x),x), 1.0, 3)


def test_Lorentz():
    lor = pb.Lorentz(x0=0.0, hwhm=1.0, scale=1.0)
    x = np.linspace(-10.0, 10.0, 1001)
    np.testing.assert_allclose(lor(x), expected['lorentz'], rtol=1e-7)


def test_Gauss_hwhm():
    gauss = pb.Gauss(x0=0.0, hwhm=1.0, scale=1.0)
    x = np.linspace(-10.0, 10.0, 100001)
    np.testing.assert_approx_equal(
        0.5*np.ptp(x[gauss(x)>0.5*np.amax(gauss(x))]), 1.0, 3)


def test_Gauss_integral():
    gauss = pb.Gauss(x0=0.0, hwhm=1.0, scale=1.0)
    x = np.linspace(-100.0, 100.0, 100001)
    np.testing.assert_approx_equal(np.trapezoid(gauss(x),x), 1.0, 7)


def test_Gauss():
    gauss = pb.Gauss(x0=0.0, hwhm=1.0, scale=1.0)
    x = np.linspace(-5.0, 5.0, 1001)
    np.testing.assert_allclose(gauss(x), expected['gauss'], rtol=1e-7)


def test_Voigt_Gauss_hwhm():
    hwhm_L = 1.0e-10
    hwhm_G = 1.0
    voigt = pb.Voigt(x0=0.0, hwhm_L=hwhm_L, hwhm_G=hwhm_G)
    x = np.linspace(-10.0, 10.0, 100001)
    profile = voigt(x)
    voigt_hwhm = 0.5*np.ptp(x[profile>0.5*np.amax(profile)])
    np.testing.assert_approx_equal(voigt_hwhm, 1.0, 3)


def test_Voigt_Lorentz_hwhm():
    hwhm_L = 1.0
    hwhm_G = 1.0e-10
    voigt = pb.Voigt(x0=0.0, hwhm_L=hwhm_L, hwhm_G=hwhm_G)
    x = np.linspace(-10.0, 10.0, 100001)
    profile = voigt(x)
    voigt_hwhm = 0.5*np.ptp(x[profile>0.5*np.amax(profile)])
    np.testing.assert_approx_equal(voigt_hwhm, 1.0, 3)


@pytest.mark.parametrize("hwhm_L, hwhm_G",
    [(1e-10, 1.0),
     (1e-2, 1.0),
     (1.0, 1.0),
     (1.0, 1e-2),
     (1.0, 1e-10)])
def test_Voigt_integral(hwhm_L, hwhm_G):
    voigt = pb.Voigt(x0=0.0, hwhm_L=hwhm_L, hwhm_G=hwhm_G, scale=1.0)
    x = np.linspace(-1000.0, 1000.0, 100001)
    np.testing.assert_approx_equal(np.trapezoid(voigt(x),x), 1.0, 3)


@pytest.mark.parametrize("hwhm_L",
    (0.01, 0.1, 1.0, 10.0, 100.0,))
def test_Voigt(hwhm_L):
    nwidths = 10.0
    hwhm_G = 1.0
    voigt = pb.Voigt(x0=0.0, hwhm_G=hwhm_G, hwhm_L=hwhm_L)
    width = 0.5346*hwhm_L + np.sqrt(0.2166*hwhm_L**2+hwhm_G**2)
    x = np.arange(-nwidths*width, nwidths*width, width/300.0)
    expected_voigt = expected['voigt'+str(hwhm_L)]
    np.testing.assert_allclose(voigt(x), expected_voigt, rtol=1e-7)


def test_doppler_hwhm():
    temperature = 1000.0
    wn = 10000.0
    mass = np.array([18.0, 44.0])
    dop_hw = pb.doppler_hwhm(temperature, mass, wn)
    np.testing.assert_allclose(
        dop_hw, [0.02669241481944169, 0.01707252588229824], rtol=1e-7)


def test_lorentz_hwhm():
    temperature = 1000.0
    pressure = 1.0 # bar
    #                  H2O   CO2   H2    He   
    masses = np.array([18.0, 44.0, 2.0,  4.0])
    radii  = np.array([1.6,  1.9,  1.45, 1.4]) * pc.A
    vmr    = np.array([1e-4, 1e-4, 0.85, 0.15])
    imol = np.array([0, 1])
    lor_hw = pb.lorentz_hwhm(temperature, pressure, masses, radii, vmr, imol)
    np.testing.assert_allclose(
        lor_hw, [0.036911106660883666,0.04308068108378928], rtol=1e-7,
    )


def test_min_widths():
    min_temp = 500.0     # Kelvin
    max_temp = 2500.0  # Kelvin
    min_wn   = 1.0/(0.6*pc.um)
    max_mass = 27.02534  # HCN molecular mass (amu)
    min_rad  = 2.5 * pc.A
    min_press = 1e-5 # bar
    dmin, lmin = pb.min_widths(
        min_temp, max_temp, min_wn, max_mass, min_rad, min_press,
    )
    np.testing.assert_allclose(dmin, 0.025672743788107903)
    np.testing.assert_allclose(lmin, 3.9945391902150206e-07)


def test_max_widths():
    min_temp =  500.0  # Kelvin
    max_temp = 2500.0  # Kelvin
    max_wn   = 1.0/(12.0*pc.um)
    min_mass = 1.007940  # H molecular mass (amu)
    max_rad  = 2.89/2.0 * pc.A
    max_press = 100  # bar
    dmax, lmax = pb.max_widths(
        min_temp, max_temp, max_wn, min_mass, max_rad, max_press,
    )
    np.testing.assert_allclose(dmax, 0.01486262554305687)
    np.testing.assert_allclose(lmax, 8.009255370607491)

