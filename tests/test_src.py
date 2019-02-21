import os
import sys
import pytest

import numpy as np
import scipy.integrate as si

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'

sys.path.append(ROOT + 'pyratbay/lib')
import simpson as s
import trapz   as t
import cutils  as cu
import indices

os.chdir(ROOT+'tests')


@pytest.mark.skip
def test_simps():
    nx = 999
    x = np.linspace(-3, 3, nx)
    sigma = 0.5
    y = np.exp(-0.5*(x/sigma)**2) / np.sqrt(2*np.pi*sigma**2)

    I = si.simps(y, x, even='first')

    h = np.ediff1d(x)
    hsum, hrat, hfac = s.geth(h)
    integ = s.simps(y, h, hsum, hrat, hfac)


def test_ifirst_default():
    assert indices.ifirst(np.array([0])) == -1


@pytest.mark.parametrize("def_ret", [-1, 0, 1])
def test_ifirst_with_default(def_ret):
    assert indices.ifirst(np.array([0]), def_ret) == def_ret


def test_ifirst():
    assert indices.ifirst(np.array([1,0,0])) == 0
    assert indices.ifirst(np.array([1,1,0])) == 0
    assert indices.ifirst(np.array([0,1,0])) == 1
    assert indices.ifirst(np.array([0,1,1])) == 1
    assert indices.ifirst(np.array([0,0,1])) == 2
    assert indices.ifirst(np.array([0,2,1])) == 2


def test_ilast_default():
    assert indices.ilast(np.array([0])) == -1


@pytest.mark.parametrize("def_ret", [-1, 0, 1])
def test_ilast_with_default(def_ret):
    assert indices.ilast(np.array([0]), def_ret) == def_ret


def test_ilast():
    assert indices.ilast(np.array([1,0,0])) == 0
    assert indices.ilast(np.array([0,1,0])) == 1
    assert indices.ilast(np.array([1,1,0])) == 1
    assert indices.ilast(np.array([0,1,1])) == 2
    assert indices.ilast(np.array([0,0,1])) == 2
    assert indices.ilast(np.array([0,2,1])) == 2
