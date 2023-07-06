# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import struct
import pytest
import re

import mc3
import numpy as np

import pyratbay.atmosphere as pa
import pyratbay.constants as pc
import pyratbay.io as io
import pyratbay.tools as pt

os.chdir(pc.ROOT+'tests')


def test_tmp_reset_listed_arguments():
     # All listed arguments are set to None:
     o   = type('obj', (object,), {'x':1.0, 'y':2.0})
     obj = type('obj', (object,), {'z':3.0, 'w':4.0, 'o':o})
     with pt.tmp_reset(obj, 'o.x', 'z'):
         x, y, w, z = obj.o.x, obj.o.y, obj.z, obj.w
     assert (x, y, w, z) == (None, 2.0, None, 4.0)


def test_eta_seconds():
    str_time = pt.eta(45.0, 101, 200)
    assert str_time == '44.11 sec'


def test_eta_minutes():
    str_time = pt.eta(145.0, 101, 200)
    assert str_time == '2.37 min'


def test_eta_hours():
    str_time = pt.eta(345.0, 11, 200)
    assert str_time == '1.65 hours'


def test_eta_days():
    str_time = pt.eta(945.0, 11, 2000)
    assert str_time == '1.98 days'


def test_eta_fmt():
    str_time = pt.eta(45.0, 101, 200, fmt='.5e')
    assert str_time == '4.41089e+01 sec'


def test_tmp_reset_keyword_arguments():
     # Keyword arguments can be set to a value, but cannot be recursive:
     o   = type('obj', (object,), {'x':1.0, 'y':2.0})
     obj = type('obj', (object,), {'z':3.0, 'w':4.0, 'o':o})
     with pt.tmp_reset(obj, 'o.x', z=10):
         x, y, w, z = obj.o.x, obj.o.y, obj.z, obj.w
     assert (x,y,w,z) == (None, 2.0, 10, 4.0)


def test_resolve_theme_rgb():
    theme = pt.resolve_theme((1,0,0))
    assert isinstance(theme, mc3.plots.Theme)
    np.testing.assert_allclose(theme.color, np.array([1.0, 0.0, 0.0]))
    np.testing.assert_allclose(theme.dark_color, np.array([0.7, 0.0, 0.0]))
    np.testing.assert_allclose(theme.light_color, np.array([1.00, 0.25, 0.25]))


def test_resolve_theme_str_color():
    theme = pt.resolve_theme('xkcd:green')
    assert isinstance(theme, mc3.plots.Theme)
    assert theme.color == 'xkcd:green'


def test_resolve_theme_str_theme():
    theme = pt.resolve_theme('orange')
    assert isinstance(theme, mc3.plots.Theme)
    assert theme.color == 'darkorange'


def test_resolve_theme_Theme():
    theme = pt.resolve_theme(mc3.plots.THEMES['indigo'])
    assert isinstance(theme, mc3.plots.Theme)
    assert theme.color == 'xkcd:indigo'


def test_resolve_theme_none():
    theme = pt.resolve_theme(None)
    assert theme is None


def test_resolve_theme_not_a_color():
    error = re.escape("Invalid color theme: 'not_a_plt_color'")
    with pytest.raises(ValueError, match=error):
        theme = pt.resolve_theme('not_a_plt_color')


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
    with open('outputs/binsearch.dat', 'wb') as tli:
        tli.write(struct.pack(str(len(wn))+"d", *list(wn)))
    with open('outputs/binsearch.dat', 'rb') as tli:
        assert pt.binsearch(tli, wn0, 0, len(wn), upper) == result


@pytest.mark.parametrize('wn0, upper, result',
    [(0.0, False, 0), (0.0, True,-1),
     (1.0, False, 0), (1.0, True, 0),
     (1.5, False, 1), (1.5, True, 0),
     (2.0, False, 1), (2.0, True, 1),
     (2.5, False,-1), (2.5, True, 1)])
def test_binsearch_two(wn0, upper, result):
    wn = np.array([1., 2.])
    with open('outputs/binsearch.dat', 'wb') as tli:
        tli.write(struct.pack(str(len(wn))+"d", *list(wn)))
    with open('outputs/binsearch.dat', 'rb') as tli:
        assert pt.binsearch(tli, wn0, 0, len(wn), upper) == result


@pytest.mark.parametrize('wn0, upper, result',
    [(0.0, False, 0), (0.0, True, 0),
     (0.5, False, 1), (0.5, True, 0),
     (1.0, False, 1), (1.0, True, 3),
     (1.5, False, 4), (1.5, True, 3),
     (2.0, False, 4), (2.0, True, 4)])
def test_binsearch_duplicates(wn0, upper, result):
    wn = np.array([0.0, 1.0, 1.0, 1.0, 2.0])
    with open('outputs/binsearch.dat', 'wb') as tli:
        tli.write(struct.pack(str(len(wn))+"d", *list(wn)))
    with open('outputs/binsearch.dat', 'rb') as tli:
        assert pt.binsearch(tli, wn0, 0, len(wn), upper) == result


@pytest.mark.parametrize('wn0, upper, result',
    [(1.0, False, 0), (1.0, True, 3),
     (1.5, False, 4), (1.5, True, 3),
     (2.0, False, 4), (2.0, True, 5)])
def test_binsearch_duplicates_low_edge(wn0, upper, result):
    wn = np.array([1.0, 1.0, 1.0, 1.0, 2.0, 2.0])
    with open('outputs/binsearch.dat', 'wb') as tli:
        tli.write(struct.pack(str(len(wn))+"d", *list(wn)))
    with open('outputs/binsearch.dat', 'rb') as tli:
        assert pt.binsearch(tli, wn0, 0, len(wn), upper) == result


@pytest.mark.parametrize('wn0, upper, result',
    [(1.0, False, 0), (1.0, True, 1),
     (1.5, False, 2), (1.5, True, 1),
     (2.0, False, 2), (2.0, True, 5)])
def test_binsearch_duplicates_hi_edge(wn0, upper, result):
    wn = np.array([1.0, 1.0, 2.0, 2.0, 2.0, 2.0])
    with open('outputs/binsearch.dat', 'wb') as tli:
        tli.write(struct.pack(str(len(wn))+"d", *list(wn)))
    with open('outputs/binsearch.dat', 'rb') as tli:
        assert pt.binsearch(tli, wn0, 0, len(wn), upper) == result


def test_unpack_string():
    value = 'H2O'
    with open('outputs/packet.dat', 'wb') as bfile:
        bfile.write(struct.pack('3s', value.encode('utf-8')))
    with open('outputs/packet.dat', 'rb') as bfile:
        output = pt.unpack(bfile, 3, 's')
    assert output == value


def test_unpack_number():
    value = 8
    with open('outputs/packet.dat', 'wb') as bfile:
        bfile.write(struct.pack('h', value))
    with open('outputs/packet.dat', 'rb') as bfile:
        output = pt.unpack(bfile, 1, 'h')
    assert output == value


def test_unpack_tuple():
    value = np.pi, np.e, np.inf
    with open('outputs/packet.dat', 'wb') as bfile:
        bfile.write(struct.pack('3f', *value))
    with open('outputs/packet.dat', 'rb') as bfile:
        output = pt.unpack(bfile, 3, 'f')
    np.testing.assert_allclose(output, value)


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
    assert pt.get_param(value, 'km') == result


def test_get_param_array():
    value = np.array([10.0, 20.0])
    np.testing.assert_allclose(pt.get_param(value, 'm'),
        np.array([1000.0, 2000.0]))


def test_get_param_none():
    assert pt.get_param(None, 'km') is None


def test_get_param_invalid_str_value():
    with pytest.raises(ValueError, match="Invalid value 'val'"):
        pt.get_param('val')


def test_get_param_invalid_str_many_values():
    with pytest.raises(ValueError, match="Invalid value '1 3 km'"):
        pt.get_param('1 3 km')


def test_get_param_invalid_str_units():
    with pytest.raises(ValueError, match="Invalid units for value '1 unit'"):
        pt.get_param('1 unit')


def test_get_param_invalid_units():
    with pytest.raises(ValueError, match="Invalid units 'unit'"):
        pt.get_param('1.0', 'unit')


@pytest.mark.parametrize('value', [-1.0, 0.0])
def test_get_param_invalid_value_gt(value):
    with pytest.raises(ValueError, match=f"Value {value} must be > 0.0"):
        pt.get_param(value, 'kelvin', gt=0.0)


def test_get_param_invalid_value_ge():
    with pytest.raises(ValueError, match="Value -100.0 must be >= 0.0"):
        pt.get_param(-1.0, 'm', ge=0.0)


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


def test_isfile_none():
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


@pytest.mark.parametrize('num_fmt', ['{}', '{:.2f}'])
def test_Formatted_Write_none(num_fmt):
    fmt = pt.Formatted_Write()
    rstar = np.pi/3.14
    default_double_str = str(rstar)
    fmt.write('Stellar radius (rstar, rsun):  ' + num_fmt, None)
    assert fmt.text == 'Stellar radius (rstar, rsun):  None\n'


def test_Formatted_Write_float_fmt():
    fmt = pt.Formatted_Write()
    rstar = np.pi/3.14
    default_double_str = str(rstar)
    fmt.write('Stellar radius (rstar, rsun):  {:.2f}', rstar)
    assert fmt.text == 'Stellar radius (rstar, rsun):  1.00\n'


def test_Formatted_Write_float_default():
    fmt = pt.Formatted_Write()
    rstar = np.pi/3.14
    default_double_str = str(rstar)
    fmt.write('Stellar radius (rstar, rsun):  {}', rstar)
    assert fmt.text == f'Stellar radius (rstar, rsun):  {default_double_str}\n'


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


def test_cia_hitran():
    ciafile = 'inputs/mock_H2-H2_2011.cia'
    pt.cia_hitran(ciafile, tstep=1, wstep=1)
    outfile = 'CIA_HITRAN_H2-H2_333.3-500.0um_0200-0300K.dat'

    assert outfile in os.listdir()
    with open(outfile, 'r') as f:
        text = f.read()
    assert text == (
        '# This file contains the reformated H2-H2 CIA data from\n'
        '# HITRAN file: inputs/mock_H2-H2_2011.cia\n\n'

        '@SPECIES\n'
        'H2  H2\n\n'

        '@TEMPERATURES\n'
        '                 200       225       250       275       300\n\n'

        '# Wavenumber in cm-1, opacity in cm-1 amagat-2:\n'
        '@DATA\n'
        '     20.0  1.926e-08 1.717e-08 1.543e-08 1.396e-08 1.287e-08\n'
        '     21.0  2.116e-08 1.886e-08 1.695e-08 1.534e-08 1.414e-08\n'
        '     22.0  2.312e-08 2.062e-08 1.854e-08 1.678e-08 1.547e-08\n'
        '     23.0  2.516e-08 2.244e-08 2.018e-08 1.827e-08 1.685e-08\n'
        '     24.0  2.725e-08 2.433e-08 2.189e-08 1.983e-08 1.829e-08\n'
        '     25.0  2.941e-08 2.627e-08 2.365e-08 2.143e-08 1.978e-08\n'
        '     26.0  3.163e-08 2.827e-08 2.547e-08 2.309e-08 2.132e-08\n'
        '     27.0  3.390e-08 3.032e-08 2.733e-08 2.479e-08 2.290e-08\n'
        '     28.0  3.622e-08 3.243e-08 2.925e-08 2.655e-08 2.454e-08\n'
        '     29.0  3.859e-08 3.458e-08 3.122e-08 2.835e-08 2.621e-08\n'
        '     30.0  4.100e-08 3.678e-08 3.323e-08 3.019e-08 2.793e-08\n'
        )
    os.remove(outfile)


@pytest.mark.skip(reason='Do I want to wget this file or mock it?')
def test_cia_borysow():
    ciafile = 'ciah2he_dh_quantmech'
    pt.cia_borysow(ciafile, 'H2', 'He')
    # TBD: implement check


def test_interpolate_opacity_no_interp():
    cs_file = 'outputs/exttable_test_300-3000K_1.1-1.7um.npz'
    cs_shape, arrays, cs_data  = io.read_opacity(cs_file, extract='all')
    wn_mask = np.ones(cs_shape[3], bool)
    pressure = arrays[2]
    interp_cs = pt.interpolate_opacity(cs_file, pressure, wn_mask)
    np.testing.assert_allclose(interp_cs, cs_data)


def test_interpolate_opacity_interp():
    nlayers = 30
    pressure = pa.pressure('1e-6 bar', '100 bar', nlayers)

    cs_file = 'outputs/exttable_test_300-3000K_1.1-1.7um.npz'
    cs_shape, arrays, cs_data  = io.read_opacity(cs_file, extract='all')
    wn_mask = np.ones(cs_shape[3], bool)
    interp_cs = pt.interpolate_opacity(cs_file, pressure, wn_mask)

    assert np.shape(interp_cs)[2] == nlayers
    # Test at a couple of temperatures and wavelengths:
    expected_cs = np.array([
        [5.9256602e-26, 5.9256568e-26, 5.9256478e-26, 5.9256355e-26,
         5.9256032e-26, 5.9255582e-26, 5.9254563e-26, 5.9252774e-26,
         5.9249862e-26, 5.9242619e-26, 5.9232402e-26, 5.9205721e-26,
         5.9169112e-26, 5.9073085e-26, 5.8940521e-26, 5.8599102e-26,
         5.8143418e-26, 5.7612622e-26, 5.6636491e-26, 5.5528489e-26,
         5.5917755e-26, 6.0139450e-26, 8.8948915e-26, 1.3198975e-25,
         2.3003751e-25, 3.2905595e-25, 4.9392910e-25, 6.3977799e-25,
         7.7682735e-25, 8.4519621e-25],
        [7.9934423e-29, 1.6021669e-28, 2.8992945e-28, 5.2823678e-28,
         1.0516308e-27, 1.8142800e-27, 3.8143005e-27, 6.5806369e-27,
         1.3834911e-26, 2.3868607e-26, 5.0178817e-26, 8.6575141e-26,
         1.5244410e-25, 3.1401355e-25, 5.4180000e-25, 1.1390651e-24,
         1.9650792e-24, 4.1313340e-24, 7.1248136e-24, 1.4958770e-23,
         2.5692599e-23, 5.2247758e-23, 8.6441724e-23, 1.3506444e-22,
         1.7437439e-22, 1.8541353e-22, 1.9677101e-22, 2.1005127e-22,
         2.1333042e-22, 2.0219771e-22],
        [2.3123446e-22, 2.3123447e-22, 2.3123449e-22, 2.3123454e-22,
         2.3123460e-22, 2.3123477e-22, 2.3123502e-22, 2.3123551e-22,
         2.3123652e-22, 2.3123803e-22, 2.3124199e-22, 2.3124745e-22,
         2.3126165e-22, 2.3128124e-22, 2.3133215e-22, 2.3140087e-22,
         2.3157194e-22, 2.3178886e-22, 2.3220604e-22, 2.3264652e-22,
         2.3264711e-22, 2.2998136e-22, 2.2123921e-22, 1.8936443e-22,
         1.4978443e-22, 9.3358221e-23, 6.4372064e-23, 4.5467870e-23,
         3.5840648e-23, 3.0239903e-23],
        [6.4448852e-23, 6.4448899e-23, 6.4449022e-23, 6.4449191e-23,
         6.4449633e-23, 6.4450249e-23, 6.4451644e-23, 6.4454093e-23,
         6.4458078e-23, 6.4467998e-23, 6.4481996e-23, 6.4518589e-23,
         6.4568896e-23, 6.4701362e-23, 6.4885542e-23, 6.5272418e-23,
         6.6028853e-23, 6.7171766e-23, 7.0189503e-23, 7.4344362e-23,
         8.5163984e-23, 1.0011232e-22, 1.3962457e-22, 1.9414845e-22,
         3.2728892e-22, 4.7195767e-22, 6.4254581e-22, 7.0593098e-22,
         6.5726499e-22, 5.3100702e-22],
    ])
    i_temps = [1,3,7,9]
    i_wave = [2234, 837, 2525, 1091]
    for i in range(4):
        j = i_temps[i]
        k = i_wave[i]
        np.testing.assert_allclose(interp_cs[0,j,:,k], expected_cs[i])


def test_interpolate_opacity_extrapolate():
    nlayers = 30
    pressure = pa.pressure('1e-12 bar', '100 bar', nlayers)

    cs_file = 'outputs/exttable_test_300-3000K_1.1-1.7um.npz'
    cs_shape, arrays, cs_data  = io.read_opacity(cs_file, extract='all')
    press_table = arrays[2]
    wn_mask = np.ones(cs_shape[3], bool)
    interp_cs = pt.interpolate_opacity(cs_file, pressure, wn_mask)

    p_mask = pressure < np.amin(press_table)
    # Everything above min(press_table) is the same:
    assert np.amax(np.abs(1-interp_cs[0,:,p_mask]/cs_data[0,:,0])) < 1e-14


def test_none_div_no_num():
    div = pt.none_div(None, 1.0)
    assert div is None


def test_none_div_no_den():
    div = pt.none_div(1.0, None)
    assert div is None


def test_none_div_all_good():
    div = pt.none_div(1.0, 2.0)
    np.testing.assert_allclose(div, 0.5)


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
    ns = pt.Namespace({'path': 'file0'})
    with pytest.raises(ValueError):
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
