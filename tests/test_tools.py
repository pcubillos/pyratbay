# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import struct
import pytest

import numpy as np

import pyratbay.tools as pt
import pyratbay.io    as io
import pyratbay.constants as pc
import pyratbay.starspec  as ps

os.chdir(pc.ROOT+'tests')


def test_tmp_reset_listed_arguments():
     # All listed arguments are set to None:
     o   = type('obj', (object,), {'x':1.0, 'y':2.0})
     obj = type('obj', (object,), {'z':3.0, 'w':4.0, 'o':o})
     with pt.tmp_reset(obj, 'o.x', 'z'):
         x, y, w, z = obj.o.x, obj.o.y, obj.z, obj.w
     assert (x, y, w, z) == (None, 2.0, None, 4.0)


def test_tmp_reset_keyword_arguments():
     # Keyword arguments can be set to a value, but cannot be recursive:
     o   = type('obj', (object,), {'x':1.0, 'y':2.0})
     obj = type('obj', (object,), {'z':3.0, 'w':4.0, 'o':o})
     with pt.tmp_reset(obj, 'o.x', z=10):
         x, y, w, z = obj.o.x, obj.o.y, obj.z, obj.w
     assert (x,y,w,z) == (None, 2.0, 10, 4.0)


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


@pytest.mark.parametrize('units, value',
    [('cm', 1.0),
     ('m', 100.0),
     ('rearth', 6.3781e8),
     ('kg', 1000.0),
     ('mjup', 1.8982e30),
     ('bar', 1e6),
    ])
def test_u(units, value):
    assert pt.u(units) == value


def test_u_error():
    with pytest.raises(ValueError,
            match="Units 'fake_units' does not exist in pyratbay.constants."):
        pt.u('fake_units')


@pytest.mark.parametrize('value, result',
    [(1.0, 100000.0),
     ('10 cm', 10.0),])
def test_get_param(value, result):
    assert pt.get_param('size', value, 'km') == result


def test_get_param_array():
    value = np.array([10.0, 20.0])
    np.testing.assert_allclose(pt.get_param('size', value, 'km'),
        np.array([1000000.0, 2000000.0]))


def test_get_param_none():
    assert pt.get_param('size', None, 'km') is None


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


def test_isfile_none(tmp_path):
    assert pt.isfile(None) == -1


@pytest.mark.parametrize('inputs',
    ['file1', ['file1'], ['file1', 'file2']])
def test_isfile_exists(tmp_path, inputs):
    path1 = tmp_path / 'file1'
    path1.touch()
    path2 = tmp_path / 'file2'
    path2.touch()
    if isinstance(inputs, str):
        inputs = str(tmp_path / inputs)
    else:
        inputs = [str(tmp_path/path) for path in inputs]
    assert pt.isfile(inputs) == 1


@pytest.mark.parametrize('inputs',
    ['nofile', ['nofile1'], ['file1', 'nofile1'], ['nofile1', 'nofile2']])
def test_isfile_not_exists(tmp_path, inputs):
    path1 = tmp_path / 'file1'
    path1.touch()
    if isinstance(inputs, str):
        inputs = str(tmp_path / inputs)
    else:
        inputs = [str(tmp_path/path) for path in inputs]
    assert pt.isfile(inputs) == 0


def test_file_exists_none():
    pt.file_exists('none', 'None input', None)
    assert True


def test_file_exists_file(tmp_path):
    path = tmp_path / 'new_tmp_file.dat'
    path.touch()
    pt.file_exists('testfile', 'Test', str(path))
    assert True


def test_file_exists_raise(tmp_path):
    with pytest.raises(ValueError,
        match=r"Test file \(testfile\) does not exist: 'no_file.dat'"):
        pt.file_exists('testfile', 'Test', 'no_file.dat')


def test_path():
    assert pt.path('file.txt')   == "./file.txt"
    assert pt.path('./file.txt') == "./file.txt"
    assert pt.path('/home/user/file.txt') == "/home/user/file.txt"


def test_Formatted_Write():
    fmt = pt.Formatted_Write()
    rstar = np.pi/3.14
    default_double_str = str(rstar)
    fmt.write('Stellar radius (rstar, rsun):  {:.2f}', rstar)
    fmt.write('Stellar radius (rstar, rsun):  {:.2f}', None)
    fmt.write('Stellar radius (rstar, rsun):  {}',     rstar)
    fmt.write('Stellar radius (rstar, rsun):  {}',     None)
    assert fmt.text == ('Stellar radius (rstar, rsun):  1.00\n'
                        'Stellar radius (rstar, rsun):  None\n'
                        'Stellar radius (rstar, rsun):  {:s}\n'
                        'Stellar radius (rstar, rsun):  None\n').format(default_double_str)


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


def test_depth_to_radius_scalar():
    depth = 1.44
    depth_err = 0.6
    rprs, rprs_err = pt.depth_to_radius(depth, depth_err)
    assert rprs == 1.2
    assert rprs_err == 0.25


@pytest.mark.parametrize('func',
    [list, tuple, np.array, ])
def test_depth_to_radius_iterable(func):
    depth = 1.44, 2.25
    depth_err = 0.6, 0.9
    depth = func(depth)
    depth_err = func(depth_err)
    rprs, rprs_err = pt.depth_to_radius(depth, depth_err)
    np.testing.assert_allclose(rprs, (1.2,1.5))
    np.testing.assert_allclose(rprs_err, (0.25,0.3))


def test_radius_to_depth_scalar():
    rprs = 1.2
    rprs_err = 0.25
    depth, depth_err = pt.radius_to_depth(rprs, rprs_err)
    assert depth == 1.44
    assert depth_err == 0.6


@pytest.mark.parametrize('func',
    [list, tuple, np.array, ])
def test_radius_to_depth_iterable(func):
    rprs = 1.2, 1.5
    rprs_err = 0.25, 0.3
    rprs = func(rprs)
    rprs_err = func(rprs_err)
    depth, depth_err = pt.radius_to_depth(rprs, rprs_err)
    np.testing.assert_allclose(depth, (1.44, 2.25))
    np.testing.assert_allclose(depth_err, (0.6,0.9))


@pytest.mark.parametrize('flag, output', [(False,1), (True,None)])
def test_ignore_system_exit(flag, output):
    @pt.ignore_system_exit
    def dummy_function(flag):
        if flag:
            sys.exit()
        return 1
    if flag:
        assert dummy_function(flag) is None
    else:
        assert dummy_function(flag) == 1


def test_Namespace_get_path_str():
    ns = pt.Namespace({'path':'file0'})
    assert ns.get_path('path') == os.getcwd() + '/file0'


def test_Namespace_get_path_list():
    ns = pt.Namespace({'path':['file1', 'file2']})
    assert ns.get_path('path')[0] == os.getcwd() + '/file1'
    assert ns.get_path('path')[1] == os.getcwd() + '/file2'


def test_Namespace_get_path_root():
    ns = pt.Namespace({'path':'{ROOT}/rooted_file'})
    assert ns.get_path('path') == pc.ROOT + 'rooted_file'


def test_Namespace_get_path_non_existing():
    ns = pt.Namespace({'path':'file0'})
    with pytest.raises(SystemExit):
        ns.get_path('path', desc='Configuration', exists=True)


@pytest.mark.parametrize('var, val',
    [('10',   10),
     ('-10', -10),
     ('+10',  10),
     ('10.0', 10),
     ('1e5',  100000)])
def test_parse_int(var, val):
    args = {'par':var}
    pt.parse_int(args, 'par')
    assert args['par'] == val


@pytest.mark.parametrize('var',
    ['10.5', 'None', 'True', 'inf', 'nan', '10 20'])
def test_parse_int_fail(var):
    args = {'par':var}
    with pytest.raises(ValueError, match="Invalid data type for par, "
        "could not convert string to integer: '{:s}'".format(var)):
        pt.parse_int(args, 'par')


@pytest.mark.parametrize('var, val',
    [('10',   10.0),
     ('-10', -10.0),
     ('+10',  10.0),
     ('10.0', 10.0),
     ('10.5', 10.5),
     ('1e5',  100000.0),
     ('inf',  np.inf)])
def test_parse_float(var, val):
    args = {'par':var}
    pt.parse_float(args, 'par')
    assert args['par'] == val


def test_parse_float_nan():
    args = {'par':'nan'}
    pt.parse_float(args, 'par')
    assert np.isnan(args['par'])


@pytest.mark.parametrize('var',
    ['None', 'True', '10 20'])
def test_parse_float_fail(var):
    args = {'par':var}
    with pytest.raises(ValueError, match="Invalid data type for par, "
        "could not convert string to float: '{:s}'".format(var)):
        pt.parse_float(args, 'par')


@pytest.mark.parametrize('var, val',
    [('10 20',     np.array([10.0, 20.0])),
     ('10.0 20.0', np.array([10.0, 20.0]))])
def test_parse_array_floats(var, val):
    args = {'par':var}
    pt.parse_array(args, 'par')
    np.testing.assert_equal(args['par'], val)


@pytest.mark.parametrize('var, val',
    [('a b',    ['a', 'b']),
     ('a\nb',   ['a', 'b']),
     ('a b\nc', ['a', 'b', 'c'])])
def test_parse_array_strings(var, val):
    args = {'par':var}
    pt.parse_array(args, 'par')
    assert args['par'] == val


@pytest.mark.parametrize('parser',
   [pt.parse_str, pt.parse_int, pt.parse_float, pt.parse_array])
def test_parse_none(parser):
    args = {}
    pt.parse_array(args, 'par')
    assert args['par'] is None
