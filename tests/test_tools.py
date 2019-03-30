# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import struct
import pytest

import numpy as np

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)

import pyratbay.tools as pt
import pyratbay.io    as io
import pyratbay.constants as pc
import pyratbay.starspec  as ps

os.chdir(ROOT+'tests')


def mock_pf(epf, temp, pf):
    with open(epf, 'w') as f:
        f.write("\n".join("{:7.1f}  {:.10e}".format(t,z)
                          for t,z in zip(temp,pf)))

def test_binsearch_zero():
    with pytest.raises(ValueError,
    match='Requested binsearch over a zero a zero-sized array.'):
        pt.binsearch('dummy.dat', 1.0, 0, nrec=0, upper=True)


@pytest.mark.parametrize('wn0, upper, result',
    [(0.0, False,  0), (0.0, True, -1),
     (1.0, False,  0), (1.0, True,  0),
     (2.0, False, -1), (2.0, True,  0),])
def test_binsearch_one(wn0, upper, result):
    wn = np.array([1.0])
    with open('tli_test.dat', 'wb') as tli:
        tli.write(struct.pack(str(len(wn))+"d", *list(wn)))
    with open('tli_test.dat', 'rb') as tli:
        assert pt.binsearch(tli, wn0, 0, len(wn), upper) == result


@pytest.mark.parametrize('wn0, upper, result',
    [(0.0, False, 0), (0.0, True,-1),
     (1.0, False, 0), (1.0, True, 0),
     (1.5, False, 1), (1.5, True, 0),
     (2.0, False, 1), (2.0, True, 1),
     (2.5, False,-1), (2.5, True, 1)])
def test_binsearch_two(wn0, upper, result):
    wn = np.array([1., 2.])
    with open('tli_test.dat', 'wb') as tli:
        tli.write(struct.pack(str(len(wn))+"d", *list(wn)))
    with open('tli_test.dat', 'rb') as tli:
        assert pt.binsearch(tli, wn0, 0, len(wn), upper) == result


@pytest.mark.parametrize('wn0, upper, result',
    [(0.0, False, 0), (0.0, True, 0),
     (0.5, False, 1), (0.5, True, 0),
     (1.0, False, 1), (1.0, True, 3),
     (1.5, False, 4), (1.5, True, 3),
     (2.0, False, 4), (2.0, True, 4)])
def test_binsearch_duplicates(wn0, upper, result):
    wn = np.array([0.0, 1.0, 1.0, 1.0, 2.0])
    with open('tli_test.dat', 'wb') as tli:
        tli.write(struct.pack(str(len(wn))+"d", *list(wn)))
    with open('tli_test.dat', 'rb') as tli:
        assert pt.binsearch(tli, wn0, 0, len(wn), upper) == result


@pytest.mark.parametrize('wn0, upper, result',
    [(1.0, False, 0), (1.0, True, 3),
     (1.5, False, 4), (1.5, True, 3),
     (2.0, False, 4), (2.0, True, 5)])
def test_binsearch_duplicates_low_edge(wn0, upper, result):
    wn = np.array([1.0, 1.0, 1.0, 1.0, 2.0, 2.0])
    with open('tli_test.dat', 'wb') as tli:
        tli.write(struct.pack(str(len(wn))+"d", *list(wn)))
    with open('tli_test.dat', 'rb') as tli:
        assert pt.binsearch(tli, wn0, 0, len(wn), upper) == result


@pytest.mark.parametrize('wn0, upper, result',
    [(1.0, False, 0), (1.0, True, 1),
     (1.5, False, 2), (1.5, True, 1),
     (2.0, False, 2), (2.0, True, 5)])
def test_binsearch_duplicates_hi_edge(wn0, upper, result):
    wn = np.array([1.0, 1.0, 2.0, 2.0, 2.0, 2.0])
    with open('tli_test.dat', 'wb') as tli:
        tli.write(struct.pack(str(len(wn))+"d", *list(wn)))
    with open('tli_test.dat', 'rb') as tli:
        assert pt.binsearch(tli, wn0, 0, len(wn), upper) == result


def test_path():
    assert pt.path('file.txt')   == "./file.txt"
    assert pt.path('./file.txt') == "./file.txt"
    assert pt.path('/home/user/file.txt') == "/home/user/file.txt"


@pytest.mark.parametrize('data',
    [[False, True, True, False],
     [0,1,1,0],
     (False, True, True, False),
     np.array([False, True, True, False])])
def test_ifirst_type(data):
    assert pt.ifirst(data) == 1


@pytest.mark.parametrize('data',
    [[False, True, True, False],
     [0,1,1,0],
     (False, True, True, False),
     np.array([False, True, True, False])])
def test_ilast_type(data):
    assert pt.ilast(data) == 2


def test_wrap():
    info = []
    text = "Pyrat atmospheric model\n"
    pt.wrap(info, text)
    assert info == ['Pyrat atmospheric model']
    text = "Pressure = 1.0 bar\nTemperature = 1000.0 K"
    pt.wrap(info, text, indent=2)
    assert info == ['Pyrat atmospheric model',
                    '  Pressure = 1.0 bar',
                    '  Temperature = 1000.0 K']


@pytest.mark.parametrize('db, molecule, isotope',
   [('1H2-16O__POKAZATEL__00400-00500.trans.bz2',   'H2O', '116'),
    ('1H-2H-16O__VTT__00250-00500.trans.bz2',       'H2O', '126'),
    ('12C-16O2__HITEMP.pf',                         'CO2', '266'),
    ('12C-16O-18O__Zak.par',                        'CO2', '268'),
    ('12C-1H4__YT10to10__01100-01200.trans.bz2',    'CH4', '21111'),
    ('12C-1H3-2H__MockName__01100-01200.trans.bz2', 'CH4', '21112')])
def test_get_exomol_mol(db, molecule, isotope):
    mol, iso = pt.get_exomol_mol(db)
    assert mol == molecule
    assert iso == isotope


def test_pf_exomol_single():
    # Mock an Exomol PF:
    epf = '14N-1H4__MockBYTe.pf'
    mock_pf(epf, np.arange(1,6), np.logspace(0,1,5))
    pt.pf_exomol(epf)
    pf, iso, temp = io.read_pf('PF_Exomol_NH4.dat')
    np.testing.assert_allclose(pf[0,:], np.logspace(0,1,5), rtol=1e-5)
    np.testing.assert_allclose(temp, np.arange(1,6), rtol=1e-7)
    assert list(iso) == ['41111']


def test_pf_exomol_listed_single():
    epf = '14N-1H4__MockBYTe.pf'
    mock_pf(epf, np.arange(1,6), np.logspace(0,1,5))
    pt.pf_exomol([epf])
    pf, iso, temp = io.read_pf('PF_Exomol_NH4.dat')
    np.testing.assert_allclose(pf[0,:], np.logspace(0,1,5), rtol=1e-5)
    np.testing.assert_allclose(temp, np.arange(1,6), rtol=1e-7)
    assert list(iso) == ['41111']


def test_pf_exomol_two():
    epf1 = '14N-1H4__MockBYTe.pf'
    epf2 = '15N-1H4__MockBYTe.pf'
    mock_pf(epf1, np.arange(1,6), np.logspace(0,1,5))
    mock_pf(epf2, np.arange(1,9), np.logspace(0,1,8))
    pt.pf_exomol([epf1, epf2])
    pf, iso, temp = io.read_pf('PF_Exomol_NH4.dat')
    np.testing.assert_allclose(pf[0,:5], np.logspace(0,1,5), rtol=1e-5)
    np.testing.assert_allclose(pf[1,:],  np.logspace(0,1,8), rtol=1e-5)
    np.testing.assert_allclose(temp, np.arange(1,9), rtol=1e-7)
    assert list(iso) == ['41111', '51111']


def test_pf_kurucz_H2O():
    with open('Mock_h2opartfn.dat', 'w') as f:
        f.write(
      '''
         s the temperature increases the minor isotopomers become less
         ccurate because of missing levels.
         
           T      1H1H16O      1H1H17O     1H1H18O     1H2H16O
                   170625       30215       30445       42016 levels
           10       1.328       1.330       1.332       1.399
           20       3.349       3.361       3.373       2.945
           30       6.192       6.217       6.240       5.085
           40       9.417       9.457       9.492       7.617
           50      12.962      13.017      13.066      10.476''')
    pt.pf_kurucz('Mock_h2opartfn.dat')
    # Check outputs:
    pf, iso, temp = io.read_pf('PF_kurucz_H2O.dat')
    np.testing.assert_equal(iso,
        np.array(['1H1H16O', '1H1H17O', '1H1H18O', '1H2H16O']))
    np.testing.assert_equal(temp, np.arange(10, 51, 10.0))
    np.testing.assert_allclose(pf,
        np.array([[ 1.328,  3.349,  6.192,  9.417, 12.962],
                  [ 1.33 ,  3.361,  6.217,  9.457, 13.017],
                  [ 1.332,  3.373,  6.24 ,  9.492, 13.066],
                  [ 1.399,  2.945,  5.085,  7.617, 10.476]]))


def test_pf_kurucz_TiO():
    with open('Mock_tiopart.dat', 'w') as f:
        f.write(
      '''  T        46TiO       47TiO       48TiO       49TiO       50TiO
            10      28.829      28.970      29.107      29.240      29.369
            20      54.866      55.151      55.425      55.692      55.950
            30      81.572      82.002      82.417      82.821      83.212
            40     110.039     110.625     111.190     111.741     112.273
            50     141.079     141.834     142.564     143.273     143.960''')
    pt.pf_kurucz('Mock_tiopart.dat')
    # Check outputs:
    pf, iso, temp = io.read_pf('PF_kurucz_TiO.dat')
    np.testing.assert_equal(iso, np.array(['66','76','86','96','06']))
    np.testing.assert_equal(temp, np.arange(10, 51, 10.0))
    np.testing.assert_allclose(pf,
        np.array([[ 28.829,  54.866,  81.572, 110.039, 141.079],
                  [ 28.97 ,  55.151,  82.002, 110.625, 141.834],
                  [ 29.107,  55.425,  82.417, 111.19 , 142.564],
                  [ 29.24 ,  55.692,  82.821, 111.741, 143.273],
                  [ 29.369,  55.95 ,  83.212, 112.273, 143.96 ]]))


@pytest.mark.skip(reason='Do I want to wget this file or mock it?')
def test_cia_hitran():
    ciafile = 'H2-H_2011.cia'
    pt.cia_hitran(ciafile, tstep=1, wstep=1)
    # TBD: implement check


@pytest.mark.skip(reason='Do I want to wget this file or mock it?')
def test_cia_borysow():
    ciafile = 'ciah2he_dh_quantmech'
    pt.cia_borysow(ciafile, 'H2', 'He')
    # TBD: implement check


def test_tophat_dlambda():
    wl0     = 1.50
    width   = 0.50
    margin  = 0.10
    dlambda = 0.05
    wl, trans = pt.tophat(wl0, width, margin, dlambda)
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
    wl, trans = pt.tophat(wl0, width, margin, resolution=resolution)
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
    wl, trans = pt.tophat(wl0, width, margin, dlambda, ffile=tmp_file)
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
    resampled, wnidx = pt.resample(signal, wn, specwn)
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
    resampled, wnidx = pt.resample(signal, wn, specwn)
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
    resampled, wnidx = pt.resample(signal, wn, specwn, normalize=True)
    # For an equi-spaced specwn:
    dx = specwn[1] - specwn[0]
    np.testing.assert_approx_equal(np.sum(resampled)*dx, 1.0)


def test_resample_outbounds():
    wn = np.linspace(1.3, 1.7, 11)
    signal = np.array(np.abs(wn-1.5)<0.1, np.double)
    specwn = np.linspace(1.4, 2, 101)
    with pytest.raises(ValueError,
    match="Resampling signal's wavenumber is not contained in specwn."):
        resampled, wnidx = pt.resample(signal, wn, specwn)


def test_band_integrate_single():
    wn = np.arange(1500, 5000.1, 1.0)
    signal = np.ones_like(wn)
    wn1, irac1 = io.read_spectrum(pc.ROOT+"inputs/filters/spitzer_irac1_sa.dat")
    bandflux = pt.band_integrate(signal, wn, irac1, wn1)
    np.testing.assert_allclose(bandflux, [1.0])


def test_band_integrate_multiple():
    wn = np.arange(1500, 5000.1, 1.0)
    signal = np.ones_like(wn)
    wn1, irac1 = io.read_spectrum(pc.ROOT+"inputs/filters/spitzer_irac1_sa.dat")
    wn2, irac2 = io.read_spectrum(pc.ROOT+"inputs/filters/spitzer_irac2_sa.dat")
    bandflux = pt.band_integrate(signal, wn, [irac1, irac2], [wn1, wn2])
    np.testing.assert_allclose(bandflux, [1.0, 1.0])


def test_band_integrate():
    wn = np.arange(1500, 5000.1, 1.0)
    sflux = ps.bbflux(wn, 1800.0)
    wn1, irac1 = io.read_spectrum(pc.ROOT+"inputs/filters/spitzer_irac1_sa.dat")
    wn2, irac2 = io.read_spectrum(pc.ROOT+"inputs/filters/spitzer_irac2_sa.dat")
    bandfluxes = pt.band_integrate(sflux, wn, [irac1,irac2], [wn1, wn2])
    np.testing.assert_allclose(bandfluxes, [98527.148526, 84171.417692])
