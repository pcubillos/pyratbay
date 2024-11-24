# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import struct
import pytest
import re

import mc3
import numpy as np

import pyratbay as pb
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


@pytest.mark.parametrize(
    'value',
    ['1', '1.0', '-3.14', '+3.14', '1.0e+02', 'inf', 'nan'],
)
def test_is_number_true(value):
    assert pt.is_number(value) is True


@pytest.mark.parametrize(
    'value',
    ['1.0-3.14', '10abcde', '1.0e', 'e1.0', 'true', 'none'],
)
def test_is_number_false(value):
    assert pt.is_number(value) is False


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


@pytest.mark.parametrize(
    'file_path',
    ['file.log', './file.log'],
)
def test_mkdir_in_cwd(tmp_path, file_path):
    with pt.cd(tmp_path):
        pt.mkdir(file_path)
        # No dir made, nothing to test except no error was raised


def test_mkdir_in_existing_absolute_path(tmp_path):
    path_file = tmp_path / 'my_file.log'
    pt.mkdir(path_file)
    # No dir made, nothing to test except no error was raised


@pytest.mark.parametrize(
    'file_path',
    ['NS/file.log', './NS/file.log'],
)
def test_mkdir_make_path_in_cwd(tmp_path, file_path):
    with pt.cd(tmp_path):
        pt.mkdir(file_path)
        assert os.path.exists('NS')


def test_mkdir_make_dir_in_existing_absolute_path(tmp_path):
    file_path = tmp_path / 'NS/my_file.log'
    pt.mkdir(file_path)
    assert os.path.exists(tmp_path / 'NS')


def test_mkdir_nested_new_dirs(tmp_path):
    file_path = tmp_path / 'NS/dir2/my_file.log'
    with pytest.raises(FileNotFoundError):
        pt.mkdir(file_path)


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
    cs_shape, cs_u, arrays, cs_data  = io.read_opacity(cs_file, extract='all')
    temperature = arrays[1]
    pressure = arrays[2]
    interp_cs = pt.interpolate_opacity(cs_file, temperature, pressure)
    np.testing.assert_allclose(interp_cs, cs_data)


def test_interpolate_opacity_wn_mask():
    cs_file = 'outputs/exttable_test_300-3000K_1.1-1.7um.npz'
    cs_shape, cs_u, arrays, cs_data  = io.read_opacity(cs_file, extract='all')
    temperature = arrays[1]
    pressure = arrays[2]
    wn = arrays[3]
    wn_mask = (wn>=6000) & (wn<=8000)
    interp_cs = pt.interpolate_opacity(
        cs_file, temperature, pressure, wn_mask,
    )
    np.testing.assert_allclose(interp_cs, cs_data[:,:,:,wn_mask])


def test_interpolate_opacity_thin_no_interp():
    cs_file = 'outputs/exttable_test_300-3000K_1.1-1.7um.npz'
    cs_shape, cs_u, arrays, cs_data  = io.read_opacity(cs_file, extract='all')
    temperature = arrays[1]
    pressure = arrays[2]
    thin = 10
    interp_cs = pt.interpolate_opacity(
        cs_file, temperature, pressure, wn_thinning=thin,
    )
    np.testing.assert_allclose(interp_cs, cs_data[:,:,:,::thin])


def test_interpolate_opacity_mask_thin():
    cs_file = 'outputs/exttable_test_300-3000K_1.1-1.7um.npz'
    cs_shape, cs_u, arrays, cs_data  = io.read_opacity(cs_file, extract='all')
    temperature = arrays[1]
    pressure = arrays[2]
    wn = arrays[3]
    wn_mask = (wn>=6000) & (wn<=8000)
    thin = 10
    interp_cs = pt.interpolate_opacity(
        cs_file, temperature, pressure, wn_mask, wn_thinning=thin,
    )
    expected_cs = cs_data[:,:,:,wn_mask][:,:,:,::thin]
    np.testing.assert_allclose(interp_cs, expected_cs)


def test_interpolate_opacity_interp_pressure():
    nlayers = 30
    pressure = pa.pressure('1e-6 bar', '100 bar', nlayers)

    cs_file = 'outputs/exttable_test_300-3000K_1.1-1.7um.npz'
    cs_shape, cs_u, arrays, cs_data  = io.read_opacity(cs_file, extract='all')
    interp_cs = pt.interpolate_opacity(cs_file, pressure=pressure)

    assert np.shape(interp_cs)[2] == nlayers
    # Test at a couple of temperatures and wavelengths:
    expected_cs = np.array([
        [5.9256602e-26, 5.9256555e-26, 5.9256496e-26, 5.9256348e-26,
         5.9256057e-26, 5.9255466e-26, 5.9254728e-26, 5.9252899e-26,
         5.9249187e-26, 5.9241789e-26, 5.9232640e-26, 5.9209819e-26,
         5.9162928e-26, 5.9071232e-26, 5.8907561e-26, 5.8675405e-26,
         5.8124616e-26, 5.7650505e-26, 5.6494602e-26, 5.5598046e-26,
         5.5908190e-26, 6.2732416e-26, 9.0794691e-26, 1.2920160e-25,
         2.1522752e-25, 3.3813688e-25, 4.9021336e-25, 6.4348449e-25,
         7.7893727e-25, 8.4519621e-25],
        [7.9934423e-29, 1.5805125e-28, 3.1249628e-28, 5.0706958e-28,
         1.0026088e-27, 1.9824199e-27, 3.9197627e-27, 6.3601197e-27,
         1.2575619e-26, 2.4865286e-26, 4.9165161e-26, 8.0319985e-26,
         1.5773454e-25, 3.1188204e-25, 6.1672524e-25, 1.2193091e-24,
         1.9785237e-24, 3.9120528e-24, 7.7315564e-24, 1.5265081e-23,
         2.4673733e-23, 4.7932903e-23, 8.8928006e-23, 1.4134727e-22,
         1.6978598e-22, 1.8597358e-22, 1.9662256e-22, 2.1173349e-22,
         2.1322711e-22, 2.0219771e-22],
        [2.3123446e-22, 2.3123447e-22, 2.3123449e-22, 2.3123453e-22,
         2.3123460e-22, 2.3123475e-22, 2.3123504e-22, 2.3123561e-22,
         2.3123641e-22, 2.3123824e-22, 2.3124192e-22, 2.3124901e-22,
         2.3125958e-22, 2.3128170e-22, 2.3132630e-22, 2.3141285e-22,
         2.3153971e-22, 2.3176808e-22, 2.3216747e-22, 2.3265329e-22,
         2.3244980e-22, 2.3053653e-22, 2.2034883e-22, 1.9016749e-22,
         1.4033983e-22, 1.0358983e-22, 6.5056047e-23, 4.4927539e-23,
         3.5164807e-23, 3.0239903e-23],
        [6.4448852e-23, 6.4448899e-23, 6.4448998e-23, 6.4449200e-23,
         6.4449598e-23, 6.4450408e-23, 6.4451418e-23, 6.4453922e-23,
         6.4459002e-23, 6.4469135e-23, 6.4481669e-23, 6.4512966e-23,
         6.4577406e-23, 6.4703947e-23, 6.4859073e-23, 6.5258456e-23,
         6.6064676e-23, 6.7637561e-23, 7.0615985e-23, 7.4625355e-23,
         8.4593206e-23, 1.0398544e-22, 1.4225449e-22, 1.9121900e-22,
         3.0693973e-22, 4.8429187e-22, 6.5145866e-22, 7.0003931e-22,
         6.5500481e-22, 5.3100702e-22],
    ])
    i_temps = [1,3,7,9]
    i_wave = [2234, 837, 2525, 1091]
    for i in range(4):
        j = i_temps[i]
        k = i_wave[i]
        np.testing.assert_allclose(interp_cs[0,j,:,k], expected_cs[i])


def test_interpolate_opacity_interp_temperature_pressure():
    cs_file = 'outputs/exttable_test_300-3000K_1.1-1.7um.npz'
    cs_shape, cs_u, arrays, cs_data  = io.read_opacity(cs_file, extract='all')
    nlayers = 30
    pressure = pa.pressure('1e-6 bar', '100 bar', nlayers)
    temperature = arrays[1][1::2]
    ntemps = len(temperature)

    interp_cs = pt.interpolate_opacity(cs_file, temperature, pressure)
    cs_shape = np.shape(interp_cs)
    assert cs_shape[1] == ntemps
    assert cs_shape[2] == nlayers
    # Test at a couple of temperatures and wavelengths:
    expected_cs = np.array([
        [5.9256602e-26, 5.9256555e-26, 5.9256496e-26, 5.9256348e-26,
         5.9256057e-26, 5.9255466e-26, 5.9254728e-26, 5.9252899e-26,
         5.9249187e-26, 5.9241789e-26, 5.9232640e-26, 5.9209819e-26,
         5.9162928e-26, 5.9071232e-26, 5.8907561e-26, 5.8675405e-26,
         5.8124616e-26, 5.7650505e-26, 5.6494602e-26, 5.5598046e-26,
         5.5908190e-26, 6.2732416e-26, 9.0794691e-26, 1.2920160e-25,
         2.1522752e-25, 3.3813688e-25, 4.9021336e-25, 6.4348449e-25,
         7.7893727e-25, 8.4519621e-25],
        [7.9934423e-29, 1.5805125e-28, 3.1249628e-28, 5.0706958e-28,
         1.0026088e-27, 1.9824199e-27, 3.9197627e-27, 6.3601197e-27,
         1.2575619e-26, 2.4865286e-26, 4.9165161e-26, 8.0319985e-26,
         1.5773454e-25, 3.1188204e-25, 6.1672524e-25, 1.2193091e-24,
         1.9785237e-24, 3.9120528e-24, 7.7315564e-24, 1.5265081e-23,
         2.4673733e-23, 4.7932903e-23, 8.8928006e-23, 1.4134727e-22,
         1.6978598e-22, 1.8597358e-22, 1.9662256e-22, 2.1173349e-22,
         2.1322711e-22, 2.0219771e-22],
        [2.3123446e-22, 2.3123447e-22, 2.3123449e-22, 2.3123453e-22,
         2.3123460e-22, 2.3123475e-22, 2.3123504e-22, 2.3123561e-22,
         2.3123641e-22, 2.3123824e-22, 2.3124192e-22, 2.3124901e-22,
         2.3125958e-22, 2.3128170e-22, 2.3132630e-22, 2.3141285e-22,
         2.3153971e-22, 2.3176808e-22, 2.3216747e-22, 2.3265329e-22,
         2.3244980e-22, 2.3053653e-22, 2.2034883e-22, 1.9016749e-22,
         1.4033983e-22, 1.0358983e-22, 6.5056047e-23, 4.4927539e-23,
         3.5164807e-23, 3.0239903e-23],
        [6.4448852e-23, 6.4448899e-23, 6.4448998e-23, 6.4449200e-23,
         6.4449598e-23, 6.4450408e-23, 6.4451418e-23, 6.4453922e-23,
         6.4459002e-23, 6.4469135e-23, 6.4481669e-23, 6.4512966e-23,
         6.4577406e-23, 6.4703947e-23, 6.4859073e-23, 6.5258456e-23,
         6.6064676e-23, 6.7637561e-23, 7.0615985e-23, 7.4625355e-23,
         8.4593206e-23, 1.0398544e-22, 1.4225449e-22, 1.9121900e-22,
         3.0693973e-22, 4.8429187e-22, 6.5145866e-22, 7.0003931e-22,
         6.5500481e-22, 5.3100702e-22],
    ])
    i_temps = [0,1,3,4]
    i_wave = [2234, 837, 2525, 1091]
    for i in range(4):
        j = i_temps[i]
        k = i_wave[i]
        np.testing.assert_allclose(interp_cs[0,j,:,k], expected_cs[i])


def test_interpolate_opacity_interp_and_wn_masking():
    cs_file = 'outputs/exttable_test_300-3000K_1.1-1.7um.npz'
    cs_shape, cs_u, arrays, cs_data  = io.read_opacity(cs_file, extract='all')
    nlayers = 30
    pressure = pa.pressure('1e-6 bar', '100 bar', nlayers)
    temperature = arrays[1][1::2]
    ntemps = len(temperature)
    wn = arrays[3]
    wn_mask = (wn>=6000) & (wn<=8000)
    thin = 10
    expected_wn = wn[wn_mask][::thin]
    nwave = len(expected_wn)

    interp_cs = pt.interpolate_opacity(
        cs_file, temperature, pressure, wn_mask, wn_thinning=thin,
    )
    cs_shape = np.shape(interp_cs)
    assert cs_shape[1] == ntemps
    assert cs_shape[2] == nlayers
    assert cs_shape[3] == nwave
    # Test at a couple of temperatures and wavelengths:
    expected_cs = np.array([
        [7.74448417e-25, 7.74464911e-25, 7.74485639e-25, 7.74537253e-25,
         7.74638920e-25, 7.74845698e-25, 7.75103589e-25, 7.75743079e-25,
         7.77039922e-25, 7.79625009e-25, 7.82820195e-25, 7.90791250e-25,
         8.07133193e-25, 8.39036390e-25, 8.96065534e-25, 9.76402641e-25,
         1.16593626e-24, 1.50914368e-24, 2.08618034e-24, 2.68815554e-24,
         3.68264421e-24, 4.54668842e-24, 4.75500526e-24, 4.26593708e-24,
         3.31417224e-24, 2.77563774e-24, 2.93613872e-24, 3.31473243e-24,
         3.60918512e-24, 3.63204895e-24],
        [8.93169626e-30, 1.76602983e-29, 3.49143743e-29, 5.66588373e-29,
         1.12029303e-28, 2.21511230e-28, 4.37985651e-28, 7.10665765e-28,
         1.40517193e-27, 2.77839227e-27, 5.49360692e-27, 8.97477883e-27,
         1.76254223e-26, 3.48505046e-26, 6.89085295e-26, 1.36217750e-25,
         2.23377669e-25, 4.44885692e-25, 8.85073150e-25, 1.75342161e-24,
         2.87171732e-24, 5.92105094e-24, 1.28200364e-23, 2.58533151e-23,
         4.14267121e-23, 7.62105540e-23, 1.28913002e-22, 1.93630010e-22,
         2.37083156e-22, 2.41601263e-22],
        [4.40350514e-23, 4.40352507e-23, 4.40356427e-23, 4.40363370e-23,
         4.40373532e-23, 4.40398517e-23, 4.40447418e-23, 4.40541591e-23,
         4.40674028e-23, 4.40976390e-23, 4.41584923e-23, 4.42761850e-23,
         4.44518436e-23, 4.48215336e-23, 4.55747585e-23, 4.70685628e-23,
         4.93423961e-23, 5.37457336e-23, 6.28010896e-23, 8.05564734e-23,
         1.13515248e-22, 1.51887332e-22, 2.40781605e-22, 3.97125611e-22,
         6.54022521e-22, 8.85314276e-22, 1.09824343e-21, 1.06357094e-21,
         9.35551345e-22, 8.40289164e-22],
        [5.29568118e-24, 5.29571174e-24, 5.29577759e-24, 5.29591128e-24,
         5.29617462e-24, 5.29671000e-24, 5.29737831e-24, 5.29903506e-24,
         5.30239596e-24, 5.30909906e-24, 5.31739040e-24, 5.33809588e-24,
         5.38071144e-24, 5.46435422e-24, 5.56691213e-24, 5.83095626e-24,
         6.36215042e-24, 7.39628648e-24, 9.36555832e-24, 1.19730996e-23,
         1.84127607e-23, 3.06841533e-23, 5.28102844e-23, 7.53554311e-23,
         1.09576636e-22, 1.30047451e-22, 1.35644410e-22, 1.39477466e-22,
         1.47766760e-22, 1.56840584e-22],
    ])
    i_temps = [0, 1, 3, 4]
    i_wave = [50, 75, 100, 125]
    for i in range(4):
        j = i_temps[i]
        k = i_wave[i]
        np.testing.assert_allclose(interp_cs[0,j,:,k], expected_cs[i])



def test_interpolate_opacity_extrapolate():
    nlayers = 30
    pressure = pa.pressure('1e-12 bar', '100 bar', nlayers)

    cs_file = 'outputs/exttable_test_300-3000K_1.1-1.7um.npz'
    cs_shape, cs_u, arrays, cs_data  = io.read_opacity(cs_file, extract='all')
    press_table = arrays[2]
    interp_cs = pt.interpolate_opacity(cs_file, pressure=pressure)

    p_mask = pressure < np.amin(press_table)
    # Everything above min(press_table) is the same:
    relative_diff = interp_cs[0,:,p_mask]/cs_data[0,:,0]
    expected_diff = np.ones_like(relative_diff)
    np.testing.assert_allclose(relative_diff, expected_diff, rtol=1e-8)
    

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


@pytest.mark.parametrize('var', ['True', 'true', '1', 'yes'])
def test_parse_bool_true(var):
    args = {'par': var}
    pt.parse_bool(args, 'par')
    assert args['par'] is True


@pytest.mark.parametrize('var', ['False', 'false', '0', 'no'])
def test_parse_bool_false(var):
    args = {'par': var}
    pt.parse_bool(args, 'par')
    assert args['par'] is False


def test_parse_bool_default():
    args = {'different_par': 'True'}
    pt.parse_bool(args, 'par')
    assert args['par'] is False


def test_parse_bool_raise():
    args = {'par': 'yes, please!'}
    match = re.escape(
        "Invalid data type for parameter 'par', could not convert string "
        "'yes, please!'to bool"
    )
    with pytest.raises(ValueError, match=match):
        pt.parse_bool(args, 'par')


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



def test_weighted_to_equal_default():
    posterior = pt.weighted_to_equal('inputs/multinest_output.txt')
    assert posterior.shape == (15000,2)
    # Not the highest tolerance because there's a rng in between:
    np.testing.assert_allclose(np.median(posterior), -2.470685, rtol=0.05)
    np.testing.assert_allclose(np.std(posterior), 1.4980580, rtol=0.05)


def test_weighted_to_equal_with_weighted():
    posterior, weighted = pt.weighted_to_equal(
        'inputs/multinest_output.txt',
        get_weighted=True,
    )
    assert posterior.shape == (15000,2)
    assert weighted.shape == (1825,2)


def test_weighted_to_equal_with_size():
    posterior, weighted = pt.weighted_to_equal(
        'inputs/multinest_output.txt',
        get_weighted=True,
        min_size=0,
    )
    assert posterior.shape == (1825,2)
    assert weighted.shape == (1825,2)


def test_get_multinest_map_1mode():
    bestp = pt.get_multinest_map('inputs/multinest_outputstats.txt')
    expected_bestp = np.array([
        -4.45678543, -0.873886357,  1352.33371, -3.35902318, -3.43078850,
    ])
    np.testing.assert_allclose(bestp, expected_bestp)


def test_get_multinest_map_2modes():
    bestp = pt.get_multinest_map('inputs/multinest_outputstats_2modes.txt')
    expected_bestp = np.array([
        472.225158, 2.53499615, -0.102761300, 1.47440028, -1.39857495,
    ])
    np.testing.assert_allclose(bestp, expected_bestp)


def test_loglike():
    pyrat = pb.Pyrat('configs/mcmc_transmission_test.cfg', log=False)
    loglike = pt.Loglike(pyrat)

    # For the log_like, free parameters only
    ifree = pyrat.ret.pstep > 0
    free_pars = pyrat.ret.params[ifree]
    like = loglike(free_pars)

    np.testing.assert_allclose(like, -1627.6484855668375)
    # A non-physical model
    #like = loglike(free_pars)
    #np.testing.assert_allclose(like, -1e+98)

    # Now with better-fitting parameters:
    # map_pars = np.array([
    #     -4.4567854, -0.87388636, 1352.3337, -3.3590232, -3.4307885,
    # ])
    # like = loglike(map_pars)
    # np.testing.assert_allclose(like, -1627.5530504136932)


def test_get_mpi_rank():
    # In reallyty we would like to call this with mpirun, but
    # that needs to be configured in Github Actions, which is work
    os.environ['PBAY_NO_MPI'] = "1"
    rank = pt.get_mpi_rank()
    assert rank == 0


def test_get_mpi_size():
    # In reallyty we would like to call this with mpirun, but
    # that needs to be configured in Github Actions, which is work
    os.environ['PBAY_NO_MPI'] = "1"
    size = pt.get_mpi_size()
    assert size == 1


def test_mpi_barrier():
    # In reallyty we would like to call this with mpirun, but
    # that needs to be configured in Github Actions, which is work
    os.environ['PBAY_NO_MPI'] = "1"
    pt.mpi_barrier()
