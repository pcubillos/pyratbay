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
    pf, temp, iso = io.read_pf('PF_Exomol_NH4.dat')
    np.testing.assert_allclose(pf[0,:], np.logspace(0,1,5), rtol=1e-5)
    np.testing.assert_allclose(temp, np.arange(1,6), rtol=1e-7)
    assert list(iso) == ['41111']


def test_pf_exomol_listed_single():
    epf = '14N-1H4__MockBYTe.pf'
    mock_pf(epf, np.arange(1,6), np.logspace(0,1,5))
    pt.pf_exomol([epf])
    pf, temp, iso = io.read_pf('PF_Exomol_NH4.dat')
    np.testing.assert_allclose(pf[0,:], np.logspace(0,1,5), rtol=1e-5)
    np.testing.assert_allclose(temp, np.arange(1,6), rtol=1e-7)
    assert list(iso) == ['41111']


def test_pf_exomol_two():
    epf1 = '14N-1H4__MockBYTe.pf'
    epf2 = '15N-1H4__MockBYTe.pf'
    mock_pf(epf1, np.arange(1,6), np.logspace(0,1,5))
    mock_pf(epf2, np.arange(1,9), np.logspace(0,1,8))
    pt.pf_exomol([epf1, epf2])
    pf, temp, iso = io.read_pf('PF_Exomol_NH4.dat')
    np.testing.assert_allclose(pf[0,:5], np.logspace(0,1,5), rtol=1e-5)
    np.testing.assert_allclose(pf[1,:], np.logspace(0,1,8), rtol=1e-5)
    np.testing.assert_allclose(temp, np.arange(1,9), rtol=1e-7)
    assert list(iso) == ['41111', '51111']

