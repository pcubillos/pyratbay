# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os
import pytest
import subprocess

import numpy as np

from conftest import make_config

import pyratbay as pb
import pyratbay.constants as pc
import pyratbay.io as io
from pyratbay.constants import ROOT

os.chdir(ROOT + 'tests')


expected_pressure = np.logspace(0, 8, 81)
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
    1665.09706555])

expected_abundance = np.load(
    f'{ROOT}/tests/expected/expected_tea_profile.npz')['arr_0']

expected_radius = np.array([
    7.33195650e+09, 7.32831268e+09, 7.32467248e+09, 7.32103589e+09,
    7.31740291e+09, 7.31377352e+09, 7.31014771e+09, 7.30652549e+09,
    7.30290684e+09, 7.29929175e+09, 7.29568022e+09, 7.29207222e+09,
    7.28846775e+09, 7.28486679e+09, 7.28126932e+09, 7.27767531e+09,
    7.27408475e+09, 7.27049760e+09, 7.26691382e+09, 7.26333335e+09,
    7.25975615e+09, 7.25618213e+09, 7.25261120e+09, 7.24904325e+09,
    7.24547813e+09, 7.24191566e+09, 7.23835561e+09, 7.23479770e+09,
    7.23124157e+09, 7.22768677e+09, 7.22413274e+09, 7.22057877e+09,
    7.21702399e+09, 7.21346729e+09, 7.20990732e+09, 7.20634236e+09,
    7.20277032e+09, 7.19918862e+09, 7.19559409e+09, 7.19198290e+09,
    7.18835043e+09, 7.18469120e+09, 7.18099877e+09, 7.17726570e+09,
    7.17348355e+09, 7.16964296e+09, 7.16573383e+09, 7.16174551e+09,
    7.15766725e+09, 7.15348861e+09, 7.14920000e+09, 7.14479333e+09,
    7.14026261e+09, 7.13560460e+09, 7.13081938e+09, 7.12591084e+09,
    7.12088693e+09, 7.11575961e+09, 7.11054445e+09, 7.10525955e+09,
    7.09992411e+09, 7.09455648e+09, 7.08917228e+09, 7.08378310e+09,
    7.07839620e+09, 7.07301531e+09, 7.06764188e+09, 7.06227633e+09,
    7.05691868e+09, 7.05156887e+09, 7.04622681e+09, 7.04089240e+09,
    7.03556550e+09, 7.03024596e+09, 7.02493358e+09, 7.01962812e+09,
    7.01432929e+09, 7.00903671e+09, 7.00374992e+09, 6.99846837e+09,
    6.99319134e+09])


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
    pressure, temperature, abundances, species, radius = pb.run(cfg)


def test_units_separate(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0',
               'mpunits':'mjup'})
    pressure, temperature, abundances, species, radius = pb.run(cfg)


def test_units_in_value(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0 rjup'})
    pressure, temperature, abundances, species, radius = pb.run(cfg)


def test_units_missing(tmp_path, capfd):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Invalid units 'None' for parameter mplanet." in captured.out


def test_units_invalid(tmp_path, capfd):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0', 'mpunits':'nope'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Invalid planet mass units (mpunits): nope" in captured.out


def test_units_in_value_invalid(tmp_path, capfd):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'mplanet':'1.0 nope'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Invalid units for value '1.0 nope' for parameter mplanet." \
        in captured.out


@pytest.mark.sort(order=1)
def test_tli_hitran_wfc3():
    pb.run(ROOT+'tests/configs/tli_hitran_1.1-1.7um_test.cfg')
    assert 'HITRAN_H2O_1.1-1.7um_test.tli' in os.listdir('outputs/')


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_repack():
    pb.run('configs/tli_repack_test.cfg')
    # TBD: asserts on output file


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_tio_schwenke():
    pb.run('configs/tli_tio_schwenke_test.cfg')
    # TBD: asserts on output file


def test_pt_isothermal(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg')

    pressure, temperature, abundances, species, radius = pb.run(cfg)
    np.testing.assert_allclose(pressure, expected_pressure, rtol=1e-7)
    np.testing.assert_equal(temperature, np.tile(1500.0, 81))


def test_pt_TCEA(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_tcea.cfg')
    pressure, temperature, abundances, species, radius = pb.run(cfg)
    np.testing.assert_allclose(pressure, expected_pressure, rtol=1e-7)
    np.testing.assert_allclose(temperature, expected_temperature, rtol=1e-7)


def test_atmosphere_uniform(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_uniform_test.cfg',
        reset={'atmfile':atmfile})

    press, temp, abund, species, radius = pb.run(cfg)
    np.testing.assert_allclose(press, expected_pressure, rtol=1e-7)
    np.testing.assert_allclose(temp, expected_temperature, rtol=1e-7)
    q = np.tile([0.85, 0.149, 3.0e-6, 4.0e-4, 1.0e-4, 5.0e-4, 1.0e-7], (81,1))
    np.testing.assert_equal(abund, q)
    # Compare against the atmospheric file now:
    atm = io.read_atm(atmfile)
    assert atm[0] == ('bar', 'kelvin', 'volume', None)
    np.testing.assert_equal(atm[1], np.array('H2 He Na H2O CH4 CO CO2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atm[2]*pc.bar, expected_pressure, rtol=3e-5)
    np.testing.assert_allclose(atm[3], expected_temperature, rtol=1e-6)
    np.testing.assert_equal(atm[4], q)


def test_atmosphere_tea(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/atmosphere_tea_test.cfg',
        reset={'atmfile':atmfile},
    )

    press, temp, abundance, species, radius = pb.run(cfg)
    np.testing.assert_allclose(press, expected_pressure, rtol=1e-7)
    np.testing.assert_allclose(temp,  expected_temperature, atol=1e-7)
    np.testing.assert_allclose(abundance, expected_abundance, rtol=1e-7)
    # Compare against the atmospheric file now:
    atmf = io.read_atm(atmfile)
    assert atmf[0] == ('bar', 'kelvin', 'volume', None)
    np.testing.assert_equal(atmf[1],
        np.array('H2 He Na K H2O CH4 CO CO2 NH3 HCN N2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atmf[2]*pc.bar, expected_pressure, rtol=3e-5)
    np.testing.assert_allclose(atmf[3], expected_temperature, rtol=5e-6)
    np.testing.assert_allclose(atmf[4], expected_abundance, rtol=1e-6)


def test_atmosphere_hydro(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_hydro_test.cfg',
        reset={'atmfile':atmfile})

    press, temp, abund, species, radius = pb.run(cfg)
    np.testing.assert_allclose(press, expected_pressure, rtol=1e-7)
    np.testing.assert_allclose(temp, expected_temperature, rtol=1e-7)
    q = np.tile([0.85, 0.149, 3.0e-6, 4.0e-4, 1.0e-4, 5.0e-4, 1.0e-7], (81,1))
    np.testing.assert_equal(abund, q)
    np.testing.assert_allclose(radius, expected_radius, rtol=1e-7)
    # Compare against the atmospheric file now:
    atmf = io.read_atm(atmfile)
    assert atmf[0] == ('bar', 'kelvin', 'volume', 'rjup')
    np.testing.assert_equal(atmf[1],
        np.array('H2 He Na H2O CH4 CO CO2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atmf[2]*pc.bar, expected_pressure, rtol=3e-5)
    np.testing.assert_allclose(atmf[3], expected_temperature, rtol=5e-6)
    np.testing.assert_equal(atmf[4], q)
    np.testing.assert_allclose(atmf[5]*pc.rjup, expected_radius, rtol=5e-5)


def test_atmosphere_hydro_default_runits(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_hydro_test.cfg',
        reset={'atmfile':atmfile, 'gplanet':'2478.7504116251885'},
        remove=['rplanet'])

    press, temp, abund, species, radius = pb.run(cfg)
    np.testing.assert_allclose(radius, expected_radius, rtol=1e-7)
    atmf = io.read_atm(atmfile)
    assert atmf[0] == ('bar', 'kelvin', 'volume', 'rjup')
    np.testing.assert_allclose(atmf[5]*pc.rjup, expected_radius, rtol=5e-5)



# See tests/test_spectrum.py for spectrum tests


def test_spectrum_emission(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'rt_path':'emission', 'cpars':'-0.5'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    # TBD: implement asserts


@pytest.mark.sort(order=10)
def test_opacity_pbay(capfd):
    pyrat = pb.run(ROOT+'tests/configs/opacity_test.cfg')
    captured = capfd.readouterr()
    assert "Extinction-coefficient table written to file:" in captured.out
    assert "exttable_test_300-3000K_1.1-1.7um.npz'." in captured.out
    assert "exttable_test_300-3000K_1.1-1.7um.npz" in os.listdir('outputs/')


@pytest.mark.skip
def test_mcmc():
    pyrat = pb.run(ROOT+'tests/configs/mcmc_transmission_test.cfg')
