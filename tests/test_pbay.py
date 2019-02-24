import os
import sys
import subprocess
import pytest

import numpy as np

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay as pb
import pyratbay.constants  as pc
import pyratbay.atmosphere as pa

os.chdir(ROOT+'tests')

expected_pressure = np.array(
      [1.00000000e+00, 1.25892541e+00, 1.58489319e+00, 1.99526231e+00,
       2.51188643e+00, 3.16227766e+00, 3.98107171e+00, 5.01187234e+00,
       6.30957344e+00, 7.94328235e+00, 1.00000000e+01, 1.25892541e+01,
       1.58489319e+01, 1.99526231e+01, 2.51188643e+01, 3.16227766e+01,
       3.98107171e+01, 5.01187234e+01, 6.30957344e+01, 7.94328235e+01,
       1.00000000e+02, 1.25892541e+02, 1.58489319e+02, 1.99526231e+02,
       2.51188643e+02, 3.16227766e+02, 3.98107171e+02, 5.01187234e+02,
       6.30957344e+02, 7.94328235e+02, 1.00000000e+03, 1.25892541e+03,
       1.58489319e+03, 1.99526231e+03, 2.51188643e+03, 3.16227766e+03,
       3.98107171e+03, 5.01187234e+03, 6.30957344e+03, 7.94328235e+03,
       1.00000000e+04, 1.25892541e+04, 1.58489319e+04, 1.99526231e+04,
       2.51188643e+04, 3.16227766e+04, 3.98107171e+04, 5.01187234e+04,
       6.30957344e+04, 7.94328235e+04, 1.00000000e+05, 1.25892541e+05,
       1.58489319e+05, 1.99526231e+05, 2.51188643e+05, 3.16227766e+05,
       3.98107171e+05, 5.01187234e+05, 6.30957344e+05, 7.94328235e+05,
       1.00000000e+06, 1.25892541e+06, 1.58489319e+06, 1.99526231e+06,
       2.51188643e+06, 3.16227766e+06, 3.98107171e+06, 5.01187234e+06,
       6.30957344e+06, 7.94328235e+06, 1.00000000e+07, 1.25892541e+07,
       1.58489319e+07, 1.99526231e+07, 2.51188643e+07, 3.16227766e+07,
       3.98107171e+07, 5.01187234e+07, 6.30957344e+07, 7.94328235e+07,
       1.00000000e+08])

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


def test_tli_hitran():
    #subprocess.call(['wget', '--user=HITRAN', '--password=getdata', '-N',
    #                 'https://www.cfa.harvard.edu/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Compressed-files/01_hit12.zip'])
    #subprocess.call(['unzip', '01_hit12.zip'])
    pb.pbay.run(ROOT+'tests/tli_hitran_test.cfg')
    # asserts on output file


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_pands():
    pb.pbay.run('tli_pands_test.cfg')
    # asserts on output file


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_exomol():
    pb.pbay.run('tli_exomol_test.cfg')
    # asserts on output file


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_repack():
    pb.pbay.run('tli_repack_test.cfg')
    # asserts on output file


@pytest.mark.skip(reason="Skip until implementing in Python3")
def test_tli_tio_schwenke():
    pb.pbay.run('tli_tio_schwenke_test.cfg')
    # asserts on output file



def test_pt_isothermal():
    pressure, temperature = pb.pbay.run(ROOT+'tests/pt_isothermal_test.cfg')
    np.testing.assert_allclose(pressure, expected_pressure, rtol=1e-7)
    np.testing.assert_equal(temperature, np.tile(1500.0, 81))


def test_pt_TCEA():
    pressure, temperature = pb.pbay.run(ROOT+'tests/pt_tcea_test.cfg')
    np.testing.assert_allclose(pressure, expected_pressure, rtol=1e-7)
    np.testing.assert_allclose(temperature, expected_temperature, atol=1e-10)


def test_atmosphere_uniform():
    atm = pb.pbay.run(ROOT+"tests/atmosphere_uniform_test.cfg")
    np.testing.assert_allclose(atm[0], expected_pressure,    rtol=1e-7)
    np.testing.assert_allclose(atm[1], expected_temperature, atol=1e-10)
    q = np.tile([0.85, 0.149, 3.0e-6, 4.0e-4, 1.0e-4, 5.0e-4, 1.0e-7], (81,1))
    np.testing.assert_equal(atm[2], q)
    # Compare against the atmospheric file now:
    atm = pa.readatm(ROOT+'tests/atmosphere_uniform_test.atm')
    assert atm[0] == ('bar', 'kelvin', 'number', None)
    np.testing.assert_equal(atm[1], np.array('H2 He Na H2O CH4 CO CO2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atm[2]*pc.bar, expected_pressure, rtol=3e-5)
    np.testing.assert_allclose(atm[3], expected_temperature,     atol=5e-4)
    np.testing.assert_equal(atm[4], q)


@pytest.mark.skip(reason="Skip until implementing the fast TEA")
def test_atmosphere_tea():
    atm = pb.pbay.run(ROOT+"tests/atmosphere_tea_test.cfg")
    np.testing.assert_allclose(atm[0], expected_pressure,    rtol=1e-7)
    np.testing.assert_allclose(atm[1], expected_temperature, atol=1e-10)
    np.testing.assert_allclose(atm[2], expected_abundances,  rtol=1e-7)
    # Compare against the atmospheric file now:
    atm = pa.readatm(ROOT+'tests/atmosphere_tea_test.atm')
    assert atm[0] == ('bar', 'kelvin', 'number', None)
    np.testing.assert_equal(atm[1],
        np.array('H2 He Na K H2O CH4 CO CO2 NH3 HCN N2'.split()))
    # File read-write loses precision:
    np.testing.assert_allclose(atm[2]*pc.bar, expected_pressure, rtol=3e-5)
    np.testing.assert_allclose(atm[3], expected_temperature,     atol=5e-4)
    np.testing.assert_allclose(atm[4], expected_abundances,      rtol=1e-4)


@pytest.mark.skip(reason="See tests/test_spectrum.py")
def test_spectrum_transmission():
    pyrat = pb.pbay.run(ROOT+'tests/spectrum_transmission_test.cfg')
    # implement asserts


def test_spectrum_emission():
    pyrat = pb.pbay.run(ROOT+'tests/spectrum_emission_test.cfg')
    # implement asserts


def test_opacity():
    pyrat = pb.pbay.run(ROOT+'tests/opacity_test.cfg')
    # implement asserts


@pytest.mark.skip
def test_mcmc():
    pass
