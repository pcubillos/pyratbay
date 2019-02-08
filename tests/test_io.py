import os
import sys

import numpy as np

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay.io as io


def test_read_write_opacity():
    ofile = 'opacity_test.dat'
    molID = np.array([101, 105])
    temp  = np.linspace(300, 3000, 28)
    press = np.logspace(-6, 2, 21)
    wn    = np.linspace(1000, 2000, 1001)

    nmol    = len(molID)
    ntemp   = len(temp)
    nlayers = len(press)
    nwave   = len(wn)

    etable = np.linspace(0.0, 1.0, nmol*ntemp*nlayers*nwave).reshape(
                         (nmol, ntemp, nlayers, nwave))

    io.write_opacity(ofile, molID, temp, press, wn, etable)
    edata = io.read_opacity(ofile)

    assert nmol    ==    2
    assert ntemp   ==   28
    assert nlayers ==   21
    assert nwave   == 1001
    assert edata[0] == (nmol, ntemp, nlayers, nwave)

    np.testing.assert_almost_equal(molID, edata[1][0], decimal=7)
    np.testing.assert_almost_equal(temp,  edata[1][1], decimal=7)
    np.testing.assert_almost_equal(press, edata[1][2], decimal=7)
    np.testing.assert_almost_equal(wn,    edata[1][3], decimal=7)

    np.testing.assert_almost_equal(etable, edata[2],   decimal=7)
