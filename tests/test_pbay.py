# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import pytest
import re
import subprocess

import numpy as np

from conftest import make_config

import pyratbay as pb
import pyratbay.constants as pc
import pyratbay.io as io
import pyratbay.atmosphere as pa
from pyratbay.constants import ROOT

os.chdir(ROOT + 'tests')


expected_pressure = pa.pressure('1e-6 bar', '100 bar', nlayers=81)
expected_temperature = np.array([
    1046.89433798, 1046.89534525, 1046.89661946, 1046.89823135,
    1046.90027034, 1046.90284953, 1046.90611195, 1046.91023851,
    1046.91545794, 1046.9220595 , 1046.93040895, 1046.94096875,
    1046.95432364, 1046.97121287, 1046.99257099, 1047.01957934,
    1047.05373107, 1047.09691321, 1047.15151021, 1047.22053447,
    1047.30779086, 1047.41808387, 1047.55747799, 1047.73362488,
    1047.95617334, 1048.23728226, 1048.59226036, 1049.04036124,
    1049.60576707, 1050.31879834, 1051.21739055, 1052.34887907,
    1053.77212894, 1055.56003335, 1057.80237786, 1060.60902021,
    1064.11326045, 1068.47516244, 1073.88443242, 1080.56226117,
    1088.76131307, 1098.76283701, 1110.86975475, 1125.39465039,
    1142.64193973, 1162.88419029, 1186.3335257 , 1213.11007182,
    1243.21017731, 1276.47740788, 1312.57904127, 1350.99026501,
    1390.98801603, 1431.65689115, 1471.91099843, 1510.53770339,
    1546.27096545, 1577.9014716 , 1604.42506087, 1625.21631993,
    1640.18896334, 1649.87550785, 1655.34982468, 1657.96419233,
    1658.9854222 , 1659.31327079, 1659.42226561, 1659.49173767,
    1659.56937861, 1659.66627037, 1659.78818839, 1659.94163515,
    1660.13475269, 1660.37777747, 1660.68357589, 1661.06831324,
    1661.55228907, 1662.16097781, 1662.92632193, 1663.88833303,
    1665.09706555,
])
expected_temperature_madhu = np.array([
    1300.39928128, 1300.63210647, 1300.95565962, 1301.38579756,
    1301.93577227, 1302.61570817, 1303.43267611, 1304.39119013,
    1305.49388442, 1306.74216399, 1308.13671376, 1309.6778389 ,
    1311.36566367, 1313.20023427, 1315.18156644, 1317.30966395,
    1319.58452681, 1322.00615501, 1324.57454856, 1327.28970745,
    1330.15163169, 1333.16032128, 1336.31577621, 1339.61799648,
    1343.0669821 , 1346.66273307, 1350.40524939, 1354.29453104,
    1358.33057805, 1362.5133904 , 1366.8429681 , 1371.31931114,
    1375.94241952, 1380.71229326, 1385.62893233, 1390.69233676,
    1395.90250653, 1401.25944164, 1406.7631421 , 1412.41360791,
    1418.21080526, 1424.15465089, 1430.24492698, 1436.4810888 ,
    1442.86190456, 1449.38486209, 1456.04535703, 1462.83584552,
    1469.74537081, 1476.76002857, 1483.86482936, 1491.04693861,
    1498.29954996, 1505.62455979, 1513.03363807, 1520.54541108,
    1528.17920864, 1535.94559525, 1543.83431159, 1551.80067553,
    1559.75278637, 1567.54405188, 1574.97722796, 1581.8250059 ,
    1587.8668003 , 1592.93330178, 1596.94409851, 1599.92401619,
    1601.99180122, 1603.32648976, 1604.1252013 , 1604.56713048,
    1604.79269689, 1604.89868562, 1604.94442906, 1604.96249712,
    1604.96897056, 1604.97101225, 1604.971507  , 1604.971507  ,
    1604.971507,
])

expected_vmr = np.load(
    f'{ROOT}/tests/expected/expected_tea_profile.npz')['arr_0']

expected_vmr_sub_solar = np.load(
    f'{ROOT}/tests/expected/expected_tea_sub_solar_profile.npz')['arr_0']

expected_radius = np.array([
    7.33194831e+09, 7.32830466e+09, 7.32466463e+09, 7.32102821e+09,
    7.31739539e+09, 7.31376616e+09, 7.31014053e+09, 7.30651847e+09,
    7.30289998e+09, 7.29928506e+09, 7.29567369e+09, 7.29206586e+09,
    7.28846155e+09, 7.28486075e+09, 7.28126344e+09, 7.27766960e+09,
    7.27407920e+09, 7.27049221e+09, 7.26690859e+09, 7.26332829e+09,
    7.25975125e+09, 7.25617739e+09, 7.25260662e+09, 7.24903883e+09,
    7.24547387e+09, 7.24191156e+09, 7.23835167e+09, 7.23479392e+09,
    7.23123794e+09, 7.22768330e+09, 7.22412943e+09, 7.22057562e+09,
    7.21702100e+09, 7.21346446e+09, 7.20990464e+09, 7.20633984e+09,
    7.20276796e+09, 7.19918642e+09, 7.19559205e+09, 7.19198102e+09,
    7.18834871e+09, 7.18468964e+09, 7.18099737e+09, 7.17726447e+09,
    7.17348248e+09, 7.16964207e+09, 7.16573310e+09, 7.16174496e+09,
    7.15766688e+09, 7.15348842e+09, 7.14920000e+09, 7.14479352e+09,
    7.14026300e+09, 7.13560519e+09, 7.13082018e+09, 7.12591185e+09,
    7.12088816e+09, 7.11576107e+09, 7.11054613e+09, 7.10526146e+09,
    7.09992625e+09, 7.09455885e+09, 7.08917488e+09, 7.08378593e+09,
    7.07839926e+09, 7.07301860e+09, 7.06764540e+09, 7.06228008e+09,
    7.05692265e+09, 7.05157307e+09, 7.04623124e+09, 7.04089706e+09,
    7.03557039e+09, 7.03025107e+09, 7.02493891e+09, 7.01963368e+09,
    7.01433507e+09, 7.00904271e+09, 7.00375615e+09, 6.99847482e+09,
    6.99319800e+09])


@pytest.mark.parametrize('call', ['-v', '--version'])
def test_command_line_version(capfd, call):
    subprocess.call(['pbay', call])
    captured = capfd.readouterr()
    assert f'Pyrat Bay version {pb.__version__}' in captured.out


def test_command_line_root(capfd):
    subprocess.call('pbay --root'.split())
    captured = capfd.readouterr()
    assert pb.constants.ROOT in captured.out


# Warm up, check when units are well or wrongly set:
def test_units_variable_not_needed(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg')
    atmosphere = pb.run(cfg)


def test_units_separate(tmp_path):
    reset = {
        'mplanet':'1.0',
        'mpunits':'mjup',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset=reset,
    )
    atmosphere = pb.run(cfg)
    # TBD: assert atmosphere.munits == 'mjup'


def test_units_in_value(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0 rjup'})
    atmosphere = pb.run(cfg)
    # TBD: assert atmosphere.munits == 'mjup'


def test_units_missing(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0'},
    )
    error = re.escape("Invalid units 'None' for parameter mplanet")
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_units_invalid(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0', 'mpunits':'nope'},
    )
    error = re.escape("Invalid planet mass units (mpunits): nope")
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_units_in_value_invalid(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0 nope'},
    )
    error = re.escape(
        "Invalid units for value '1.0 nope' for parameter mplanet"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.sort(order=1)
def test_tli_hitran_wfc3():
    pb.run(ROOT+'tests/configs/tli_hitran_1.1-1.7um_test.cfg')
    assert 'HITRAN_H2O_1.1-1.7um_test.tli' in os.listdir('outputs/')


@pytest.mark.skip(reason="TBD")
def test_tli_repack():
    pb.run('configs/tli_repack_test.cfg')
    # TBD: asserts on output file


@pytest.mark.skip(reason="Outdated dataset")
def test_tli_tio_schwenke():
    pb.run('configs/tli_tio_schwenke_test.cfg')
    # TBD: asserts on output file


def test_pt_isothermal(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg')
    atmosphere = pb.run(cfg)
    expected_temp = np.tile(1500.0, 81)
    np.testing.assert_allclose(atmosphere.press, expected_pressure)
    np.testing.assert_allclose(atmosphere.temp, expected_temp)


def test_pt_guillot(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_guillot.cfg')
    atmosphere = pb.run(cfg)
    np.testing.assert_allclose(atmosphere.press, expected_pressure)
    np.testing.assert_allclose(atmosphere.temp, expected_temperature)


def test_pt_tcea(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_guillot.cfg',
        reset={'tmodel': 'tcea'},
    )
    # Allowed, but discouraged:
    with pytest.warns(DeprecationWarning):
        atmosphere = pb.run(cfg)
    np.testing.assert_allclose(atmosphere.press, expected_pressure)
    np.testing.assert_allclose(atmosphere.temp, expected_temperature)


def test_pt_madhu(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_madhu.cfg')
    atmosphere = pb.run(cfg)
    np.testing.assert_allclose(atmosphere.press, expected_pressure)
    np.testing.assert_allclose(atmosphere.temp, expected_temperature_madhu)


def test_atmosphere_uniform(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/atmosphere_uniform_test.cfg',
        reset={'output_atmfile':atmfile},
    )

    atmosphere = pb.run(cfg)
    expected_vmr = np.tile(
        [0.85, 0.149, 3.0e-6, 4.0e-4, 1.0e-4, 5.0e-4, 1.0e-7],
        (81,1),
    )
    np.testing.assert_allclose(atmosphere.press, expected_pressure)
    np.testing.assert_allclose(atmosphere.temp, expected_temperature)
    np.testing.assert_equal(atmosphere.vmr, expected_vmr)
    # Compare against the atmospheric file now:
    atm = io.read_atm(atmfile)
    assert atm[0] == ('bar', 'kelvin', 'volume', None)
    np.testing.assert_equal(atm[1], np.array('H2 He Na H2O CH4 CO CO2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atm[2], expected_pressure, rtol=3e-5)
    np.testing.assert_allclose(atm[3], expected_temperature, rtol=1e-6)
    np.testing.assert_equal(atm[4], expected_vmr)


def test_atmosphere_tea_no_vmr_models(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/atmosphere_tea_test.cfg',
        reset={'output_atmfile':atmfile},
    )
    expected_species = np.array('H2 He Na K H2O CH4 CO CO2 NH3 HCN N2'.split())

    atmosphere = pb.run(cfg)
    np.testing.assert_allclose(atmosphere.press, expected_pressure)
    np.testing.assert_allclose(atmosphere.temp, expected_temperature)
    np.testing.assert_allclose(atmosphere.vmr, expected_vmr)
    # Compare against the atmospheric file now:
    atmf = io.read_atm(atmfile)
    assert atmf[0] == ('bar', 'kelvin', 'volume', None)
    np.testing.assert_equal(atmf[1], expected_species)
    # File read-write loses precision:
    np.testing.assert_allclose(atmf[2], expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atmf[3], expected_temperature, rtol=1e-6)
    np.testing.assert_allclose(atmf[4], expected_vmr, rtol=1e-6)


def test_atmosphere_tea_no_molpars_but_with_vmr_models(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    reset = {
        'output_atmfile': atmfile,
        'mol_vars': '[M/H]',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/atmosphere_tea_test.cfg',
        reset=reset,
    )
    expected_species = np.array('H2 He Na K H2O CH4 CO CO2 NH3 HCN N2'.split())

    atmosphere = pb.run(cfg)
    np.testing.assert_allclose(atmosphere.press, expected_pressure)
    np.testing.assert_allclose(atmosphere.temp, expected_temperature)
    np.testing.assert_allclose(atmosphere.vmr, expected_vmr)
    # Compare against the atmospheric file now:
    atmf = io.read_atm(atmfile)
    assert atmf[0] == ('bar', 'kelvin', 'volume', None)
    np.testing.assert_equal(atmf[1], expected_species)
    # File read-write loses precision:
    np.testing.assert_allclose(atmf[2], expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atmf[3], expected_temperature, rtol=1e-6)
    np.testing.assert_allclose(atmf[4], expected_vmr, rtol=1e-6)


def test_atmosphere_tea_with_vmr_models(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    reset = {
        'output_atmfile': atmfile,
        'vmr_vars': '[M/H] -1.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/atmosphere_tea_test.cfg',
        reset=reset,
    )
    expected_species = np.array('H2 He Na K H2O CH4 CO CO2 NH3 HCN N2'.split())

    atmosphere = pb.run(cfg)
    np.testing.assert_allclose(atmosphere.press, expected_pressure)
    np.testing.assert_allclose(atmosphere.temp, expected_temperature)
    np.testing.assert_allclose(atmosphere.vmr, expected_vmr_sub_solar)
    # Compare against the atmospheric file now:
    atmf = io.read_atm(atmfile)
    assert atmf[0] == ('bar', 'kelvin', 'volume', None)
    np.testing.assert_equal(atmf[1], expected_species)
    # File read-write loses precision:
    np.testing.assert_allclose(atmf[2], expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atmf[3], expected_temperature, rtol=1e-6)
    np.testing.assert_allclose(atmf[4], expected_vmr_sub_solar, rtol=1e-6)


@pytest.mark.parametrize(
    'arg',
    ['molmodel', 'molfree', 'molvars', 'molpars'],
)
def test_deprecated_vmr_arguments(tmp_path, arg):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/atmosphere_tea_test.cfg',
        reset={arg:'tcea'},
    )

    error = re.escape(
        f"The '{arg}' argument is deprecated, use 'vmr_vars' instead"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_atmosphere_hydro(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/atmosphere_hydro_test.cfg',
        reset={'output_atmfile':atmfile},
    )
    atmosphere = pb.run(cfg)
    expected_vmr = np.tile(
        [0.85, 0.149, 3.0e-6, 4.0e-4, 1.0e-4, 5.0e-4, 1.0e-7],
        (81,1),
    )
    #print(atmosphere.runits)
    np.testing.assert_allclose(atmosphere.press, expected_pressure)
    np.testing.assert_allclose(atmosphere.temp, expected_temperature)
    np.testing.assert_allclose(atmosphere.vmr, expected_vmr)
    np.testing.assert_allclose(atmosphere.radius, expected_radius)
    # Compare against the atmospheric file now:
    atmf = io.read_atm(atmfile)
    assert atmf[0] == ('bar', 'kelvin', 'volume', 'rjup')
    np.testing.assert_equal(atmf[1],
        np.array('H2 He Na H2O CH4 CO CO2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atmf[2], expected_pressure, rtol=1e-6)
    np.testing.assert_allclose(atmf[3], expected_temperature, rtol=1e-6)
    np.testing.assert_allclose(atmf[4], expected_vmr)
    np.testing.assert_allclose(atmf[5]*pc.rjup, expected_radius, rtol=1e-6)


# See test_emission.py and test_transmission.py for spectrum tests


@pytest.mark.sort(order=10)
def test_opacity_pbay(capfd):
    pyrat = pb.run(ROOT+'tests/configs/opacity_test.cfg')
    captured = capfd.readouterr()
    assert "Cross-section table written to file:" in captured.out
    assert "exttable_test_300-3000K_1.1-1.7um.npz'." in captured.out
    assert "exttable_test_300-3000K_1.1-1.7um.npz" in os.listdir('outputs/')


@pytest.mark.skip
def test_mcmc():
    pyrat = pb.run(ROOT+'tests/configs/mcmc_transmission_test.cfg')
