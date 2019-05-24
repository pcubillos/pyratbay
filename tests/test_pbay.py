import os
import sys
import subprocess
import pytest

import numpy as np

from conftest import make_config

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay as pb
import pyratbay.constants  as pc
import pyratbay.atmosphere as pa

os.chdir(ROOT+'tests')

expected_pressure    = np.logspace(0, 8, 81)
expected_temperature = np.array(
      [1047.04535531, 1047.04636729, 1047.04764749, 1047.04926694,
       1047.05131548, 1047.05390677, 1047.05718449, 1047.06133039,
       1047.06657429, 1047.07320679, 1047.08159536, 1047.09220465,
       1047.10562211, 1047.12259045, 1047.1440486 , 1047.17118342,
       1047.20549505, 1047.24887932, 1047.30373179, 1047.37307893,
       1047.46074334, 1047.57155184, 1047.7115971 , 1047.88856622,
       1048.11215262, 1048.39457122, 1048.75120098, 1049.2013834 ,
       1049.76941034, 1050.48573869, 1051.38847286, 1052.52515618,
       1053.95490802, 1055.75092978, 1058.00337628, 1060.82254126,
       1064.34222961, 1068.72307524, 1074.15540631, 1080.8610604 ,
       1089.09332826, 1099.13399769, 1111.28635212, 1125.86305112,
       1143.16818192, 1163.4734689 , 1186.98959518, 1213.83461222,
       1244.00218101, 1277.3326441 , 1313.48964651, 1351.94449886,
       1391.97022283, 1432.64772742, 1472.88802251, 1511.47646601,
       1547.14676037, 1578.69185734, 1605.11307825, 1625.79393924,
       1640.65971886, 1650.25479478, 1655.66159741, 1658.23445306,
       1659.23535547, 1659.55574759, 1659.66292568, 1659.73229419,
       1659.81020803, 1659.90748667, 1660.02989399, 1660.1839565 ,
       1660.37784873, 1660.62184804, 1660.92887213, 1661.31515062,
       1661.80106362, 1662.4121864 , 1663.18058734, 1664.14643497,
       1665.35997887])


@pytest.mark.sort(order=1)
def test_tli_hitran_wfc3():
    pb.run(ROOT+'tests/tli_hitran_1.1-1.7um_test.cfg')
    # TBD: asserts on output file


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_pands():
    pb.run('tli_pands_test.cfg')
    # TBD: asserts on output file


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_repack():
    pb.run('tli_repack_test.cfg')
    # TBD: asserts on output file


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_tio_schwenke():
    pb.run('tli_tio_schwenke_test.cfg')
    # TBD: asserts on output file


def test_pt_isothermal(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg')

    pressure, temperature = pb.run(cfg)
    np.testing.assert_allclose(pressure, expected_pressure, rtol=1e-7)
    np.testing.assert_equal(temperature, np.tile(1500.0, 81))


def test_pt_TCEA(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/pt_tcea.cfg')

    pressure, temperature = pb.run(cfg)
    np.testing.assert_allclose(pressure, expected_pressure, rtol=1e-7)
    np.testing.assert_allclose(temperature, expected_temperature, atol=1e-10)


def test_atmosphere_uniform(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(tmp_path, ROOT+'tests/atmosphere_uniform_test.cfg',
        reset={'atmfile':atmfile})

    atm = pb.run(cfg)
    np.testing.assert_allclose(atm[0], expected_pressure,    rtol=1e-7)
    np.testing.assert_allclose(atm[1], expected_temperature, atol=1e-10)
    q = np.tile([0.85, 0.149, 3.0e-6, 4.0e-4, 1.0e-4, 5.0e-4, 1.0e-7], (81,1))
    np.testing.assert_equal(atm[2], q)
    # Compare against the atmospheric file now:
    atm = pa.readatm(atmfile)
    assert atm[0] == ('bar', 'kelvin', 'number', None)
    np.testing.assert_equal(atm[1], np.array('H2 He Na H2O CH4 CO CO2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atm[2]*pc.bar, expected_pressure, rtol=3e-5)
    np.testing.assert_allclose(atm[3], expected_temperature,     atol=5e-4)
    np.testing.assert_equal(atm[4], q)


@pytest.mark.skip(reason="Skip until implementing the fast TEA")
def test_atmosphere_tea(tmp_path):
    atmfile = str(tmp_path / 'test.atm')
    cfg = make_config(tmp_path, ROOT+'tests/atmosphere_tea_test.cfg',
        reset={'atmfile':atmfile})

    atm = pb.run(cfg)
    np.testing.assert_allclose(atm[0], expected_pressure,    rtol=1e-7)
    np.testing.assert_allclose(atm[1], expected_temperature, atol=1e-10)
    np.testing.assert_allclose(atm[2], expected_abundances,  rtol=1e-7)
    # Compare against the atmospheric file now:
    atm = pa.readatm(atmfile)
    assert atm[0] == ('bar', 'kelvin', 'number', None)
    np.testing.assert_equal(atm[1],
        np.array('H2 He Na K H2O CH4 CO CO2 NH3 HCN N2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atm[2]*pc.bar, expected_pressure, rtol=3e-5)
    np.testing.assert_allclose(atm[3], expected_temperature,     atol=5e-4)
    np.testing.assert_allclose(atm[4], expected_abundances,      rtol=1e-4)


# See tests/test_spectrum.py for spectrum tests


def test_spectrum_emission(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'path':'eclipse', 'cpars':'-0.5'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    # TBD: implement asserts


@pytest.mark.sort(order=10)
def test_opacity(capfd):
    pyrat = pb.run(ROOT+'tests/opacity_test.cfg')
    captured = capfd.readouterr()
    assert "Extinction-coefficient table written to file:" in captured.out
    assert "exttable_test_300-3000K_1.1-1.7um.dat'." in captured.out
    assert "exttable_test_300-3000K_1.1-1.7um.dat" in os.listdir('.')


@pytest.mark.skip
def test_mcmc():
    pyrat = pb.run(ROOT+'tests/mcmc_transmission_test.cfg')
