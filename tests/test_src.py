import sys
import pytest

import numpy as np
import scipy.integrate as si


from pyratbay.constants import ROOT
sys.path.append(ROOT + 'pyratbay/lib')
import simpson as s
import trapz   as t
import cutils  as cu
import _indices


@pytest.mark.parametrize('odd', [0,1])
def test_simps(odd):
    nx = 98 + odd
    x = np.linspace(-3, 3, nx)
    sigma = 0.5
    y = np.exp(-0.5*(x/sigma)**2) / np.sqrt(2*np.pi*sigma**2)

    h = np.ediff1d(x)
    hsum, hrat, hfac = s.geth(h)
    integ = s.simps(y, h, hsum, hrat, hfac)

    np_integ = si.simps(y, x, even='last')
    np.testing.assert_approx_equal(integ, np_integ)


def test_ifirst_default():
    assert _indices.ifirst(np.array([0])) == -1


@pytest.mark.parametrize("def_ret", [-1, 0, 1])
def test_ifirst_with_default(def_ret):
    assert _indices.ifirst(np.array([0]), def_ret) == def_ret


def test_ifirst():
    assert _indices.ifirst(np.array([1,0,0])) == 0
    assert _indices.ifirst(np.array([1,1,0])) == 0
    assert _indices.ifirst(np.array([0,1,0])) == 1
    assert _indices.ifirst(np.array([0,1,1])) == 1
    assert _indices.ifirst(np.array([0,0,1])) == 2
    assert _indices.ifirst(np.array([0,2,1])) == 2


def test_ilast_default():
    assert _indices.ilast(np.array([0])) == -1


@pytest.mark.parametrize("def_ret", [-1, 0, 1])
def test_ilast_with_default(def_ret):
    assert _indices.ilast(np.array([0]), def_ret) == def_ret


def test_ilast():
    assert _indices.ilast(np.array([1,0,0])) == 0
    assert _indices.ilast(np.array([0,1,0])) == 1
    assert _indices.ilast(np.array([1,1,0])) == 1
    assert _indices.ilast(np.array([0,1,1])) == 2
    assert _indices.ilast(np.array([0,0,1])) == 2
    assert _indices.ilast(np.array([0,2,1])) == 2
