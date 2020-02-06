# Copyright (c) 2016-2020 Patricio Cubillos.
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE).

import os
import pytest
import pickle

import numpy as np

import pyratbay.io as io
import pyratbay.constants as pc
import pyratbay.part_func as pf

os.chdir(pc.ROOT+'tests')


def mock_pf(epf, temp, pf):
    with open(epf, 'w') as f:
        f.write("\n".join("{:7.1f}  {:.10e}".format(t,z)
                          for t,z in zip(temp,pf)))


def test_get_tips_molname():
    assert pf.get_tips_molname(1) == 'H2O'
    assert pf.get_tips_molname(6) == 'CH4'
    assert pf.get_tips_molname(45) == 'H2'


def test_get_tips_molname_error():
    with pytest.raises(ValueError,
            match='TIPS 2017 database does not contain molecule ID: 0'):
        dummy = pf.get_tips_molname(0)


def test_pf_tips():
    expected_temp = np.array([1] + [i for i in range(5,19,5)]
                                 + [i for i in range(20,5001,10)], float)
    with open(pc.ROOT+'inputs/tips_2017.pkl', 'rb') as p:
        expected_pf = pickle.load(p)['H2O']
    pf_data, isotopes, temp = pf.tips('H2O', outfile='default')
    np.testing.assert_equal(temp, expected_temp)
    np.testing.assert_equal(pf_data[0], expected_pf['161'])
    assert isotopes == list(expected_pf.keys())

    pf_read, iso, temp_read = io.read_pf('PF_tips_H2O.dat')
    np.testing.assert_allclose(pf_read[0,:], expected_pf['161'], rtol=1e-7)
    np.testing.assert_allclose(temp_read, expected_temp, rtol=1e-7)
    assert list(iso) == list(expected_pf.keys())


def test_pf_exomol_single():
    epf = '14N-1H4__MockBYTe.pf'
    expected_temp = np.arange(1,6)
    expected_pf   = np.logspace(0,1,5)
    mock_pf(epf, expected_temp, expected_pf)

    pf_data, isotopes, temp = pf.exomol(epf, outfile='default')
    np.testing.assert_allclose(pf_data[0,:], expected_pf)
    np.testing.assert_allclose(temp, expected_temp)
    assert list(isotopes) == ['41111']

    pf_read, iso, temp_read = io.read_pf('PF_exomol_NH4.dat')
    np.testing.assert_allclose(pf_read[0,:], expected_pf, rtol=1e-5)
    np.testing.assert_allclose(temp_read, expected_temp, rtol=1e-7)
    assert list(iso) == ['41111']


def test_pf_exomol_listed_single():
    epf = '14N-1H4__MockBYTe.pf'
    expected_temp = np.arange(1,6)
    expected_pf   = np.logspace(0,1,5)
    mock_pf(epf, expected_temp, expected_pf)

    pf_data, isotopes, temp = pf.exomol([epf], outfile='default')
    np.testing.assert_allclose(pf_data[0,:], expected_pf)
    np.testing.assert_allclose(temp, expected_temp)
    assert list(isotopes) == ['41111']

    pf_read, iso, temp_read = io.read_pf('PF_exomol_NH4.dat')
    np.testing.assert_allclose(pf_read[0,:], expected_pf, rtol=1e-5)
    np.testing.assert_allclose(temp_read, expected_temp, rtol=1e-7)
    assert list(iso) == ['41111']


def test_pf_exomol_two():
    epf1 = '14N-1H4__MockBYTe.pf'
    epf2 = '15N-1H4__MockBYTe.pf'
    expected_temp1 = np.arange(1,6)
    expected_temp2 = np.arange(1,9)
    expected_pf1   = np.logspace(0,1,5)
    expected_pf2   = np.logspace(0,1,8)
    mock_pf(epf1, expected_temp1, expected_pf1)
    mock_pf(epf2, expected_temp2, expected_pf2)

    pf_data, isotopes, temp = pf.exomol([epf1, epf2], outfile='default')
    np.testing.assert_allclose(pf_data[0,:], expected_pf1,     rtol=1e-5)
    np.testing.assert_allclose(pf_data[1,:], expected_pf2[:5], rtol=1e-5)
    np.testing.assert_allclose(temp, expected_temp1)
    assert list(isotopes) == ['41111', '51111']

    pf_read, iso, temp_read = io.read_pf('PF_exomol_NH4.dat')
    np.testing.assert_allclose(pf_read[0,:], expected_pf1,     rtol=1e-5)
    np.testing.assert_allclose(pf_read[1,:], expected_pf2[:5], rtol=1e-5)
    np.testing.assert_allclose(temp_read, expected_temp1)
    assert list(iso) == ['41111', '51111']


@pytest.mark.sort(order=1)
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

    pf_data, isotopes, temp = pf.kurucz('Mock_h2opartfn.dat', outfile='default')
    np.testing.assert_equal(temp, np.arange(10, 51, 10.0))
    np.testing.assert_allclose(pf_data,
        np.array([[ 1.328,  3.349,  6.192,  9.417, 12.962],
                  [ 1.33 ,  3.361,  6.217,  9.457, 13.017],
                  [ 1.332,  3.373,  6.24 ,  9.492, 13.066],
                  [ 1.399,  2.945,  5.085,  7.617, 10.476]]))
    np.testing.assert_equal(isotopes,
        np.array(['1H1H16O', '1H1H17O', '1H1H18O', '1H2H16O']))

    pf_read, iso, temp_read = io.read_pf('PF_kurucz_H2O.dat')
    np.testing.assert_allclose(pf_data,
        np.array([[ 1.328,  3.349,  6.192,  9.417, 12.962],
                  [ 1.33 ,  3.361,  6.217,  9.457, 13.017],
                  [ 1.332,  3.373,  6.24 ,  9.492, 13.066],
                  [ 1.399,  2.945,  5.085,  7.617, 10.476]]))
    np.testing.assert_equal(temp_read, np.arange(10, 51, 10.0))
    np.testing.assert_equal(iso,
        np.array(['1H1H16O', '1H1H17O', '1H1H18O', '1H2H16O']))


def test_pf_kurucz_TiO():
    with open('Mock_tiopart.dat', 'w') as f:
        f.write(
      '''  T        46TiO       47TiO       48TiO       49TiO       50TiO
            10      28.829      28.970      29.107      29.240      29.369
            20      54.866      55.151      55.425      55.692      55.950
            30      81.572      82.002      82.417      82.821      83.212
            40     110.039     110.625     111.190     111.741     112.273
            50     141.079     141.834     142.564     143.273     143.960''')

    pf_data, isotopes, temp = pf.kurucz('Mock_tiopart.dat', outfile='default')
    np.testing.assert_allclose(pf_data,
        np.array([[ 28.829,  54.866,  81.572, 110.039, 141.079],
                  [ 28.97 ,  55.151,  82.002, 110.625, 141.834],
                  [ 29.107,  55.425,  82.417, 111.19 , 142.564],
                  [ 29.24 ,  55.692,  82.821, 111.741, 143.273],
                  [ 29.369,  55.95 ,  83.212, 112.273, 143.96 ]]))
    np.testing.assert_equal(temp, np.arange(10, 51, 10.0))
    np.testing.assert_equal(isotopes, np.array(['66','76','86','96','06']))

    pf_read, iso, temp_read = io.read_pf('PF_kurucz_TiO.dat')
    np.testing.assert_allclose(pf_read,
        np.array([[ 28.829,  54.866,  81.572, 110.039, 141.079],
                  [ 28.97 ,  55.151,  82.002, 110.625, 141.834],
                  [ 29.107,  55.425,  82.417, 111.19 , 142.564],
                  [ 29.24 ,  55.692,  82.821, 111.741, 143.273],
                  [ 29.369,  55.95 ,  83.212, 112.273, 143.96 ]]))
    np.testing.assert_equal(temp_read, np.arange(10, 51, 10.0))
    np.testing.assert_equal(iso, np.array(['66','76','86','96','06']))
