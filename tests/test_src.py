# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import sys
import pytest

import numpy as np
import scipy.integrate as si


from pyratbay.constants import ROOT
sys.path.append(ROOT + 'pyratbay/lib')
import _pt as pt
import _simpson as s
import _trapz as t
import cutils as cu
import _indices


tcea_temp = np.array([
    1247.34007597, 1247.25420737, 1247.05317072, 1246.58981369,
    1245.54307503, 1243.2402402 , 1238.35469967, 1228.5250307 ,
    1210.37410313, 1181.94657437, 1152.76433729, 1158.74965957,
    1231.41559729, 1343.09186968, 1431.70753501, 1456.94915967,
    1458.09017515, 1458.86603858, 1460.90529873, 1466.24155921])

@pytest.mark.parametrize('g', [2200.0, None])
def test_src_tcea(g):
    nlayers = 20
    press = np.logspace(-6, 2, nlayers) * 1e6
    grav = np.tile(2200.0, nlayers)
    params = np.array([-1.5, -0.8, 0.4, 0.5, 1200.0, 100.0])
    if g is not None:
        temp = pt.tcea(params, press, grav)
    else:
        params[0] -= np.log10(2200.0)
        temp = pt.tcea(params, press)
    np.testing.assert_allclose(temp, tcea_temp)


@pytest.mark.parametrize('odd', [0,1])
def test_simps(odd):
    nx = 98 + odd
    x = np.linspace(-3, 3, nx)
    sigma = 0.5
    y = np.exp(-0.5*(x/sigma)**2) / np.sqrt(2*np.pi*sigma**2)

    h = np.ediff1d(x)
    hsum, hrat, hfac = s.geth(h)
    integ = s.simps(y, h, hsum, hrat, hfac)

    np_integ = si.simps(y, x, even='first')
    np.testing.assert_approx_equal(integ, np_integ)


def test_ediff():
    array = np.array([2.0, 1.0, 0.5, 0.25, 0.0])
    np.testing.assert_allclose(
        cu.ediff(array),
        np.array([-1.0, -0.5, -0.25, -0.25])
    )


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
