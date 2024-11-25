# Copyright (c) 2021-2022 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import pytest
import re

import numpy as np

import pyratbay.spectrum as ps
import pyratbay.constants as pc
import pyratbay.io as io

os.chdir(pc.ROOT + 'tests')


def test_PassBand_init():
    filter_file = f'{pc.FILTERS}/spitzer_irac2.dat'
    band = ps.PassBand(filter_file)

    # wl0 is passband's wavelength center of mass
    np.testing.assert_allclose(band.wl0, 4.47065351)
    # wn0 is 1/wl0 (which differs from wavenumber center of mass)
    np.testing.assert_allclose(band.wn0, 2236.80944)
    np.testing.assert_equal(band.response, band.input_response)
    np.testing.assert_allclose(band.wl, band.input_wl)
    np.testing.assert_allclose(band.wn, band.input_wn)


@pytest.mark.parametrize('flip', (False, True))
def test_PassBand_wl(flip):
    filter_file = f'{pc.FILTERS}/spitzer_irac2.dat'
    band = ps.PassBand(filter_file)
    wl = np.arange(3.5, 5.5, 0.001)
    if flip:
        wl = np.flip(wl)
    out_wl, out_response = band.set_sampling(wl)

    np.testing.assert_equal(out_response, band.response)
    np.testing.assert_equal(out_wl, band.wl)
    np.testing.assert_allclose(wl[band.idx], band.wl)

    spectrum = np.ones(len(wl))
    wn_integral = band.integrate(spectrum)
    np.testing.assert_allclose(wn_integral, 1.0)


@pytest.mark.parametrize('flip', (False, True))
def test_PassBand_wn(flip):
    filter_file = f'{pc.FILTERS}/spitzer_irac2.dat'
    band = ps.PassBand(filter_file)
    wn = 1e4 / np.flip(np.arange(3.5, 5.5, 0.001))
    if flip:
        wn = np.flip(wn)
    out_wn, out_response = band.set_sampling(wn=wn)

    np.testing.assert_equal(out_response, band.response)
    np.testing.assert_equal(out_wn, band.wn)
    np.testing.assert_allclose(wn[band.idx], band.wn)

    spectrum = np.ones(len(wn))
    wn_integral = band.integrate(spectrum)
    np.testing.assert_allclose(wn_integral, 1.0)


@pytest.mark.parametrize('wave', ('wn', 'wl'))
def test_PassBand_initialize_wave(wave):
    filter_file = f'{pc.FILTERS}/spitzer_irac2.dat'
    wl = np.arange(3.5, 5.5, 0.001)
    if wave == 'wn':
        wn = 1e4 / wl
        band = ps.PassBand(filter_file, wn=wn)
    elif wave == 'wl':
        band = ps.PassBand(filter_file, wl=wl)

    spectrum = np.ones(len(wl))
    wn_integral = band.integrate(spectrum)
    np.testing.assert_allclose(wn_integral, 1.0)


def test_PassBand__call__():
    filter_file = f'{pc.FILTERS}/spitzer_irac2.dat'
    wl = np.arange(3.5, 5.5, 0.001)
    band = ps.PassBand(filter_file, wl)
    # use __call__() method:
    spectrum = np.ones(len(wl))
    wn_integral = band(spectrum)
    np.testing.assert_allclose(wn_integral, 1.0)


def test_PassBand_name():
    filter_file = f'{pc.FILTERS}/spitzer_irac2.dat'
    band = ps.PassBand(filter_file)
    assert band.name == 'spitzer_irac2'
    band_repr = repr(band)
    assert band_repr.startswith("pyratbay.spectrum.PassBand(")
    assert band_repr.endswith("pyratbay/data/filters/spitzer_irac2.dat')")


def test_PassBand_bad_input_both():
    filter_file = f'{pc.FILTERS}/spitzer_irac2.dat'
    band = ps.PassBand(filter_file)
    wl = np.arange(3.5, 5.5, 0.001)
    wn = 1e4 / wl

    error = 'Either provide wavelength or wavenumber array, not both'
    with pytest.raises(ValueError, match=error):
        out_wn, out_response = band.set_sampling(wl, wn=wn)


def test_PassBand_bad_input_none():
    filter_file = f'{pc.FILTERS}/spitzer_irac2.dat'
    band = ps.PassBand(filter_file)
    wl = None
    wn = None

    error = 'Either provide wavelength or wavenumber array, not both'
    with pytest.raises(ValueError, match=re.escape(error)):
        out_wn, out_response = band.set_sampling(wl, wn=wn)


def test_PassBand_bad_spectral_range():
    # Band range not contained in requested wl/wn range
    filter_file = f'{pc.FILTERS}/spitzer_irac2.dat'
    band = ps.PassBand(filter_file)
    wl = np.arange(3.5, 5.5, 0.001)
    wl[10] = wl[11]

    error = (
        'Input wavelength/wavenumber array must be strictly '
        'increasing or decreasing'
    )
    with pytest.raises(ValueError, match=error):
        out_wn, out_response = band.set_sampling(wl)


def test_PassBand_save_filter(tmpdir):
    filter_file = f'{pc.FILTERS}/spitzer_irac2.dat'
    band = ps.PassBand(filter_file)
    wl = np.arange(3.5, 5.5, 0.001)
    out_wl, out_response = band.set_sampling(wl)

    save_file = 'spitzer_irac.dat'
    tmp_file = f'{tmpdir}/{save_file}'
    band.save_filter(tmp_file)

    assert save_file in os.listdir(str(tmpdir))


@pytest.mark.parametrize('flip', (False, True))
def test_Tophat_wl(flip):
    hat = ps.Tophat(4.5, 0.5)
    wl = np.arange(3.5, 5.5, 0.001)
    if flip:
        wl = np.flip(wl)
    out_wl, out_response = hat.set_sampling(wl)

    np.testing.assert_equal(out_response, hat.response)
    np.testing.assert_equal(out_wl, hat.wl)
    np.testing.assert_allclose(wl[hat.idx], hat.wl)

    flux = np.ones(len(wl))
    wn_integral = hat.integrate(flux)
    np.testing.assert_allclose(wn_integral, 1.0)


@pytest.mark.parametrize('flip', (False, True))
def test_Tophat_wn(flip):
    hat = ps.Tophat(4.5, 0.5)
    wn = 1e4 / np.flip(np.arange(3.5, 5.5, 0.001))
    if flip:
        wn = np.flip(wn)
    out_wn, out_response = hat.set_sampling(wn=wn)

    np.testing.assert_equal(out_response, hat.response)
    np.testing.assert_equal(out_wn, hat.wn)
    np.testing.assert_allclose(wn[hat.idx], hat.wn)

    flux = np.ones(len(wn))
    wn_integral = hat.integrate(flux)
    np.testing.assert_allclose(wn_integral, 1.0)


def test_Tophat_gap_error():
    hat = ps.Tophat(4.55, 0.03)
    wl = np.arange(4.0, 5.01, 0.1)
    error = re.escape(
        'Tophat() passband at wl0 = 4.550 um does not cover any spectral point'
    )
    with pytest.raises(ValueError, match=error):
        hat.set_sampling(wl)


def test_Tophat_ignore_gaps():
    hat = ps.Tophat(4.55, 0.03, ignore_gaps=True)
    wl = np.arange(4.0, 5.01, 0.1)
    out_wl, out_response = hat.set_sampling(wl)
    assert out_wl is None
    assert out_response is None
    assert hat.idx is None
    assert hat.wl is None
    assert hat.wn is None
    assert hat.response is None


def test_Tophat_name():
    hat = ps.Tophat(4.5, 0.5)
    assert str(hat) == 'tophat_4.5um'
    assert repr(hat) == 'pyratbay.spectrum.Tophat(4.5, 0.5)'


def test_Tophat_given_name():
    hat = ps.Tophat(1.2, 0.05, name='HST_WFC3')
    assert str(hat) == 'HST_WFC3_1.2um'
    assert repr(hat) == 'pyratbay.spectrum.Tophat(1.2, 0.05)'


@pytest.mark.parametrize('wl_wn', ('both', 'none'))
def test_Tophat_bad_input(wl_wn):
    hat = ps.Tophat(4.5, 0.5)
    wl = None
    wn = None
    if wl_wn == 'both':
        wl = np.arange(3.5, 5.5, 0.001)
        wn = 1e4 / wl

    error = 'Either provide wavelength or wavenumber array, not both'
    with pytest.raises(ValueError, match=error):
        out_wn, out_response = hat.set_sampling(wl, wn=wn)


def test_Tophat_bad_spectral_range():
    # Band range not contained in requested wl/wn range
    hat = ps.Tophat(4.5, 0.5)
    wl = np.arange(3.5, 5.5, 0.001)
    wl[10] = wl[11]

    error = (
        'Input wavelength/wavenumber array must be strictly '
        'increasing or decreasing'
    )
    with pytest.raises(ValueError, match=error):
        out_wn, out_response = hat.set_sampling(wl)


def test_constant_resolution_spectrum():
    wl_min = 0.5
    wl_max = 4.0
    resolution = 5.5
    wl = ps.constant_resolution_spectrum(wl_min, wl_max, resolution)

    expected_wl = np.array([
        0.5,        0.6,        0.72,       0.864,      1.0368,     1.24416,
        1.492992,   1.7915904,  2.14990848, 2.57989018, 3.09586821, 3.71504185,
    ])
    np.testing.assert_allclose(wl, expected_wl)

    # Not exactly the same because wl and res are not centered at same place:
    res = wl[1:] / np.ediff1d(wl)
    expected_res = np.tile(6.0, len(wl)-1)
    np.testing.assert_allclose(res, expected_res)

    res = wl[:-1] / np.ediff1d(wl)
    expected_res = np.tile(5.0, len(wl)-1)
    np.testing.assert_allclose(res, expected_res)


def test_bin_spectrum():
    wl_min = 0.5
    wl_max = 4.0
    resolution = 100
    wl = ps.constant_resolution_spectrum(wl_min, wl_max, resolution)
    spectrum = 3.0 * wl
    spectrum[wl>2.5] = 8.0

    bin_wl = ps.constant_resolution_spectrum(wl_min, wl_max, 10.0)
    bin_spectrum = ps.bin_spectrum(bin_wl, wl, spectrum)

    expected_bin_spectrum = np.array([
        1.54168873, 1.65831041, 1.83271796, 2.02546828, 2.23849051, 2.47391667,
        2.73410304, 3.02165368, 3.33944655, 3.69066229, 4.07881603, 4.50779261,
        4.98188545, 5.50583951, 6.08489877, 6.72485877, 7.60180994, 8.0,
        8.0       , 8.0       , 8.0       ,
    ])
    np.testing.assert_allclose(bin_spectrum, expected_bin_spectrum)


def test_bin_spectrum_gaps_ignore():
    wl_min = 1.0
    wl_max = 3.0
    resolution = 10
    wl = ps.constant_resolution_spectrum(wl_min, wl_max, resolution)
    spectrum = np.ones(len(wl))

    bin_wl = ps.constant_resolution_spectrum(wl_min, wl_max, 12.0)
    bin_spectrum = ps.bin_spectrum(bin_wl, wl, spectrum, gaps='ignore')

    expected_bin_spectrum = np.array([
       1.0, 1.0, 1.0, np.nan, 1.0, 1.0,
       1.0, 1.0, 1.0, np.nan, 1.0, 1.0,
       1.0, np.nan,
    ])
    np.testing.assert_allclose(bin_spectrum, expected_bin_spectrum)


def test_bin_spectrum_gaps_interpolate():
    wl_min = 1.0
    wl_max = 3.0
    resolution = 10
    wl = ps.constant_resolution_spectrum(wl_min, wl_max, resolution)
    spectrum = np.linspace(10.0, 20.0, len(wl))

    bin_wl = ps.constant_resolution_spectrum(wl_min, wl_max, 12.0)
    bin_spectrum = ps.bin_spectrum(bin_wl, wl, spectrum, gaps='interpolate')
    expected_bin_spectrum = np.array([
        10.0, 11.0, 12.0, 12.479167, 13.0, 14.0,
        15.0, 16.0, 17.0, 17.479167, 18.0, 19.0,
        20.0, 20.0,
    ])
    np.testing.assert_allclose(bin_spectrum, expected_bin_spectrum)


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

    expected_flux = np.array([
        1.50092461e-01, 1.49924158e+01, 1.48248054e+03,
        1.32178742e+05, 9.08239148e+05,
    ])
    np.testing.assert_allclose(flux, expected_flux)


def test_bbflux_sun():
    tsun = 5772.0
    wn = np.logspace(-1, 5, 30000)
    flux = ps.bbflux(wn, tsun)
    # Solar constant:
    s = np.trapezoid(flux, wn) * (pc.rsun/pc.au)**2
    np.testing.assert_allclose(s, 1361195.40)
    # Wien's displacement law:
    wn_max = 5.879e10*tsun/pc.c
    np.testing.assert_allclose(wn[np.argmax(flux)], wn_max, rtol=1e-4)


def test_bbflux_error():
    tsun = 5772.0
    wn = 10.0
    with pytest.raises(ValueError, match='Input wn must be an iterable.'):
        dummy = ps.bbflux(wn, tsun)


def test_read_kurucz_sun():
    kfile = 'inputs/mock_fp00k0odfnew.pck'
    tsun = 5772.0
    gsun = 4.44
    flux, wn, ktemp, klogg = ps.read_kurucz(kfile, tsun, gsun)
    # Closest tabulated values:
    assert ktemp == 5750.0
    assert klogg == 4.5
    # Arrays:
    assert len(wn) == len(flux) == 1221
    # Scientific correctness:
    s = np.trapezoid(flux, wn) * (pc.rsun/pc.au)**2
    np.testing.assert_allclose(s, 1339957.11)


def test_read_kurucz_all():
    # This is a mocked Kurucz file with fewer models but same format:
    kfile = 'inputs/mock_fp00k0odfnew.pck'
    fluxes, wn, ktemp, klogg, continua = ps.read_kurucz(kfile)
    assert np.shape(fluxes) == np.shape(continua) == (11, 1221)
    assert len(wn) == 1221
    assert len(ktemp) == 11
    assert len(klogg) == 11

    expected_kurucz_temps = np.array([
       5500.0, 5500.0, 5500.0, 5500.0,
       5750.0, 5750.0, 5750.0, 5750.0,
       6000.0, 6000.0, 6000.0,
    ])
    expected_kurucz_logg = np.array([
        3.5, 4.0, 4.5, 5.0,
        3.5, 4.0, 4.5, 5.0,
        3.5, 4.0, 4.5,
    ])
    np.testing.assert_equal(ktemp, expected_kurucz_temps)
    np.testing.assert_equal(klogg, expected_kurucz_logg)

    s = np.trapezoid(fluxes[6], wn) * (pc.rsun/pc.au)**2
    np.testing.assert_allclose(s, 1339957.11)

    c = np.trapezoid(continua[6], wn) * (pc.rsun/pc.au)**2
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
    error_msg ="Resampling signal's wavenumber is not contained in specwn."
    with pytest.raises(ValueError, match=error_msg):
        resampled, wnidx = ps.resample(signal, wn, specwn)


def test_band_integrate_single():
    wn = np.arange(1500, 5000.1, 1.0)
    signal = np.ones_like(wn)
    filter_file = pc.ROOT+"pyratbay/data/filters/spitzer_irac1.dat"
    wn1, irac1 = io.read_spectrum(filter_file)
    bandflux = ps.band_integrate(signal, wn, irac1, wn1)
    np.testing.assert_allclose(bandflux, [1.0])


def test_band_integrate_multiple_constant():
    wn = np.arange(1500, 5000.1, 1.0)
    signal = np.ones_like(wn)
    filter_file1 = pc.ROOT+"pyratbay/data/filters/spitzer_irac1.dat"
    filter_file2 = pc.ROOT+"pyratbay/data/filters/spitzer_irac2.dat"
    wn1, irac1 = io.read_spectrum(filter_file1)
    wn2, irac2 = io.read_spectrum(filter_file2)

    bandflux = ps.band_integrate(signal, wn, [irac1, irac2], [wn1, wn2])
    np.testing.assert_allclose(bandflux, [1.0, 1.0])


def test_band_integrate_multiple_blackbody():
    wn = np.arange(1500, 5000.1, 1.0)
    sflux = ps.bbflux(wn, 1800.0)
    filter_file1 = pc.ROOT+"pyratbay/data/filters/spitzer_irac1.dat"
    filter_file2 = pc.ROOT+"pyratbay/data/filters/spitzer_irac2.dat"
    wn1, irac1 = io.read_spectrum(filter_file1)
    wn2, irac2 = io.read_spectrum(filter_file2)

    bandfluxes = ps.band_integrate(sflux, wn, [irac1,irac2], [wn1, wn2])
    bandfluxes = ps.band_integrate(sflux, wn, [irac1,irac2], [wn1, wn2])
    np.testing.assert_allclose(bandfluxes, [98527.148526, 84171.417692])

