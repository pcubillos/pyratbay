# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import pytest

import numpy as np

import pyratbay.starspec   as ps
import pyratbay.constants  as pc

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
    kfile = 'fp00k0odfnew.pck'
    tsun = 5772.0
    gsun = 4.44
    flux, wn, ktemp, klogg = ps.read_kurucz(kfile, tsun, gsun)
    s = np.trapz(flux, wn) * (pc.rsun/pc.au)**2
    assert ktemp == 5750.0
    assert klogg == 4.5
    assert len(wn) == len(flux) == 1221
    np.testing.assert_allclose(s, 1339957.11)


def test_read_kurucz_all():
    kfile = 'fp00k0odfnew.pck'
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
