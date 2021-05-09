# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os
import pytest

import numpy as np

import pyratbay.spectrum as ps
import pyratbay.constants as pc
import pyratbay.io as io

os.chdir(pc.ROOT + 'tests')


@pytest.mark.parametrize('wn',
   [[1.0, 10.0, 100.0, 1000.0, 3000.0],
    [1  , 10  , 100  , 1000  , 3000  ],
    (1.0, 10.0, 100.0, 1000.0, 3000.0),
    (1  , 10  , 100  , 1000  , 3000  ),
    np.array([1.0, 10.0, 100.0, 1000.0, 3000.0]),
    np.array([1  , 10  , 100  , 1000  , 3000  ])])
def test_bbflux_type(wn):
    tsun = 5772.0
    flux = ps.bbflux(wn, tsun)
    np.testing.assert_allclose(flux,
        np.array([1.50092461e-01, 1.49924158e+01, 1.48248054e+03,
                  1.32178742e+05, 9.08239148e+05]))

def test_bbflux_sun():
    tsun = 5772.0
    wn = np.logspace(-1, 5, 30000)
    flux = ps.bbflux(wn, tsun)
    # Solar constant:
    s = np.trapz(flux, wn) * (pc.rsun/pc.au)**2
    np.testing.assert_allclose(s, 1361195.40)
    # Wien's displacement law:
    np.testing.assert_allclose(wn[np.argmax(flux)], 5.879e10*tsun/pc.c,
                               rtol=1e-4)


def test_bbflux_error():
    tsun = 5772.0
    wn = 10.0
    with pytest.raises(ValueError, match='Input wn must be an iterable.'):
        dummy = ps.bbflux(wn, tsun)


def test_read_kurucz_sun():
    kfile = 'inputs/fp00k0odfnew.pck'
    tsun = 5772.0
    gsun = 4.44
    flux, wn, ktemp, klogg = ps.read_kurucz(kfile, tsun, gsun)
    s = np.trapz(flux, wn) * (pc.rsun/pc.au)**2
    assert ktemp == 5750.0
    assert klogg == 4.5
    assert len(wn) == len(flux) == 1221
    np.testing.assert_allclose(s, 1339957.11)


def test_read_kurucz_all():
    kfile = 'inputs/fp00k0odfnew.pck'
    fluxes, wn, ktemp, klogg, continua = ps.read_kurucz(kfile)
    assert np.shape(fluxes) == np.shape(continua) == (476, 1221)
    assert len(wn) == 1221
    assert len(ktemp) == len(klogg) == 476
    np.testing.assert_equal(np.unique(ktemp),
        np.array([ 3500.,   3750.,   4000.,   4250.,   4500.,   4750.,   5000.,
                   5250.,   5500.,   5750.,   6000.,   6250.,   6500.,   6750.,
                   7000.,   7250.,   7500.,   7750.,   8000.,   8250.,   8500.,
                   8750.,   9000.,   9250.,   9500.,   9750.,  10000.,  10250.,
                  10500.,  10750.,  11000.,  11250.,  11500.,  11750.,  12000.,
                  12250.,  12500.,  12750.,  13000.,  14000.,  15000.,  16000.,
                  17000.,  18000.,  19000.,  20000.,  21000.,  22000.,  23000.,
                  24000.,  25000.,  26000.,  27000.,  28000.,  29000.,  30000.,
                  31000.,  32000.,  33000.,  34000.,  35000.,  36000.,  37000.,
                  38000.,  39000.,  40000.,  41000.,  42000.,  43000.,  44000.,
                  45000.,  46000.,  47000.,  48000.,  49000.,  50000.]))
    np.testing.assert_equal(np.unique(klogg), np.linspace(0.0, 5.0, 11))
    s = np.trapz(fluxes[108], wn) * (pc.rsun/pc.au)**2
    np.testing.assert_allclose(s, 1339957.11)
    c = np.trapz(continua[108], wn) * (pc.rsun/pc.au)**2
    np.testing.assert_allclose(c, 1618263.50)


def test_tophat_dlambda():
    wl0     = 1.50
    width   = 0.50
    margin  = 0.10
    dlambda = 0.05
    wl, trans = ps.tophat(wl0, width, margin, dlambda)
    np.testing.assert_allclose(wl, np.array(
       [1.15, 1.2 , 1.25, 1.3 , 1.35, 1.4 , 1.45, 1.5 , 1.55, 1.6 , 1.65,
        1.7 , 1.75, 1.8 , 1.85]))
    np.testing.assert_equal(trans, np.array(
       [0., 0., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 0., 0.]))


def test_tophat_resolution():
    wl0     = 1.50
    width   = 0.50
    margin  = 0.10
    resolution = 30.0
    wl, trans = ps.tophat(wl0, width, margin, resolution=resolution)
    np.testing.assert_allclose(wl, np.array(
      [1.14104722, 1.17972679, 1.21971752, 1.26106388, 1.30381181,
       1.34800882, 1.39370403, 1.44094824, 1.48979394, 1.54029543,
       1.59250883, 1.64649218, 1.70230548, 1.76001075, 1.81967213]))
    np.testing.assert_equal(trans, np.array(
      [0., 0., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 0.]))


def test_tophat_savefile(tmpdir):
    ffile = "tophat.dat"
    tmp_file = "{}/{}".format(tmpdir, ffile)
    wl0     = 1.50
    width   = 0.50
    margin  = 0.10
    dlambda = 0.05
    wl, trans = ps.tophat(wl0, width, margin, dlambda, ffile=tmp_file)
    assert ffile in os.listdir(str(tmpdir))
    with open(tmp_file, 'r') as f:
        assert f.readline() == '# Wavelength      transmission\n'
        assert f.readline() == '#         um          unitless\n'
        assert f.readline() == '     1.15000   0.000000000e+00\n'
        assert f.readline() == '     1.20000   0.000000000e+00\n'
        assert f.readline() == '     1.25000   0.000000000e+00\n'
        assert f.readline() == '     1.30000   1.000000000e+00\n'


@pytest.mark.parametrize('wn',
    [np.linspace(1.3, 1.7, 11),
     np.flip(np.linspace(1.3, 1.7, 11), axis=0)])
def test_resample_flip_wn(wn):
    signal = np.array(np.abs(wn-1.5)<0.1, np.double) * wn
    specwn = np.linspace(1, 2, 101)
    resampled, wnidx = ps.resample(signal, wn, specwn)
    np.testing.assert_equal(wnidx, np.arange(31, 70))
    np.testing.assert_allclose(resampled, np.array(
      [0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.355,
       0.71 , 1.065, 1.42 , 1.43 , 1.44 , 1.45 , 1.46 , 1.47 , 1.48 ,
       1.49 , 1.5  , 1.51 , 1.52 , 1.53 , 1.54 , 1.55 , 1.56 , 1.57 ,
       1.58 , 1.185, 0.79 , 0.395, 0.   , 0.   , 0.   , 0.   , 0.   ,
       0.   , 0.   , 0.   ]))


def test_resample_flip_specwn():
    wn = np.linspace(1.3, 1.7, 11)
    signal = np.array(np.abs(wn-1.5)<0.1, np.double) * wn
    specwn = np.flip(np.linspace(1, 2, 101), axis=0)
    resampled, wnidx = ps.resample(signal, wn, specwn)
    np.testing.assert_equal(wnidx, np.arange(31, 70))
    np.testing.assert_allclose(resampled, np.array(
      [0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.395,
       0.79 , 1.185, 1.58 , 1.57 , 1.56 , 1.55 , 1.54 , 1.53 , 1.52 ,
       1.51 , 1.5  , 1.49 , 1.48 , 1.47 , 1.46 , 1.45 , 1.44 , 1.43 ,
       1.42 , 1.065, 0.71 , 0.355, 0.   , 0.   , 0.   , 0.   , 0.   ,
       0.   , 0.   , 0.   ]))


def test_resample_normalize():
    wn = np.linspace(1.3, 1.7, 11)
    signal = np.array(np.abs(wn-1.5)<0.1, np.double)
    specwn = np.linspace(1, 2, 101)
    resampled, wnidx = ps.resample(signal, wn, specwn, normalize=True)
    # For an equi-spaced specwn:
    dx = specwn[1] - specwn[0]
    np.testing.assert_approx_equal(np.sum(resampled)*dx, 1.0)


def test_resample_outbounds():
    wn = np.linspace(1.3, 1.7, 11)
    signal = np.array(np.abs(wn-1.5)<0.1, np.double)
    specwn = np.linspace(1.4, 2, 101)
    with pytest.raises(ValueError,
        match="Resampling signal's wavenumber is not contained in specwn."):
        resampled, wnidx = ps.resample(signal, wn, specwn)


def test_band_integrate_single():
    wn = np.arange(1500, 5000.1, 1.0)
    signal = np.ones_like(wn)
    wn1, irac1 = io.read_spectrum(
        pc.ROOT+"pyratbay/data/filters/spitzer_irac1_sa.dat")
    bandflux = ps.band_integrate(signal, wn, irac1, wn1)
    np.testing.assert_allclose(bandflux, [1.0])


def test_band_integrate_multiple():
    wn = np.arange(1500, 5000.1, 1.0)
    signal = np.ones_like(wn)
    wn1, irac1 = io.read_spectrum(
        pc.ROOT+"pyratbay/data/filters/spitzer_irac1_sa.dat")
    wn2, irac2 = io.read_spectrum(
        pc.ROOT+"pyratbay/data/filters/spitzer_irac2_sa.dat")
    bandflux = ps.band_integrate(signal, wn, [irac1, irac2], [wn1, wn2])
    np.testing.assert_allclose(bandflux, [1.0, 1.0])


def test_band_integrate():
    wn = np.arange(1500, 5000.1, 1.0)
    sflux = ps.bbflux(wn, 1800.0)
    wn1, irac1 = io.read_spectrum(
        pc.ROOT+"pyratbay/data/filters/spitzer_irac1_sa.dat")
    wn2, irac2 = io.read_spectrum(
        pc.ROOT+"pyratbay/data/filters/spitzer_irac2_sa.dat")
    bandfluxes = ps.band_integrate(sflux, wn, [irac1,irac2], [wn1, wn2])
    np.testing.assert_allclose(bandfluxes, [98527.148526, 84171.417692])

