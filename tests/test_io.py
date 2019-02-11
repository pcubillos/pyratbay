import os
import sys

import numpy as np

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay.io as io
import pyratbay.atmosphere as pa
import pyratbay.constants  as pc


def test_read_write_opacity(tmpdir):
    ofile = "{}/opacity_test.dat".format(tmpdir)
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


def test_read_write_atm(tmpdir):
    atmfile = "WASP-99b.atm"
    atm = "{}/{}".format(tmpdir,atmfile)
    nlayers = 11
    pressure    = np.logspace(-8, 2, nlayers)
    temperature = np.tile(1500.0, nlayers)
    species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    qprofiles = pa.uniform(pressure, temperature, species, abundances)
    pa.writeatm(atm, pressure, temperature, species, qprofiles,
                punits='bar', header='# Test write atm\n')
    assert atmfile in os.listdir(tmpdir)

    atm_input = pa.readatm(atm)
    assert atm_input[0] == ('bar', 'kelvin', 'number', None)
    np.testing.assert_equal(atm_input[1], np.array(species))
    np.testing.assert_almost_equal(atm_input[2], pressure/pc.bar)
    np.testing.assert_almost_equal(atm_input[3], temperature)
    np.testing.assert_almost_equal(atm_input[4], qprofiles)
    assert atm_input[5] is None
