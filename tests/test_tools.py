import os
import sys
import pytest

import numpy as np

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)

import pyratbay.tools as pt
import pyratbay.io    as io

os.chdir(ROOT+'tests')


def mock_pf(epf, temp, pf):
    with open(epf, 'w') as f:
        f.write("\n".join("{:7.1f}  {:.10e}".format(t,z)
                          for t,z in zip(temp,pf)))


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
