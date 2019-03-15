import os
import sys
import pytest

import numpy as np

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay.io as io
import pyratbay.atmosphere as pa
import pyratbay.constants  as pc

os.chdir(ROOT+'tests')


def test_read_write_spectrum(tmpdir):
    sfile = "{}/spectrum_test.dat".format(tmpdir)
    wl = np.linspace(1.1, 1.7, 7) * 1e-4
    spectrum = np.ones(7)
    io.write_spectrum(wl, spectrum, filename=sfile, path='transit',
                      wlunits='um')
    # Take a look at the output file:
    with open(sfile, 'r') as f:
        content = "".join(f.readlines())
    assert content == ('# Wavelength        (Rp/Rs)**2\n'
                       '#         um          unitless\n'
                       '     1.10000   1.000000000e+00\n'
                       '     1.20000   1.000000000e+00\n'
                       '     1.30000   1.000000000e+00\n'
                       '     1.40000   1.000000000e+00\n'
                       '     1.50000   1.000000000e+00\n'
                       '     1.60000   1.000000000e+00\n'
                       '     1.70000   1.000000000e+00\n')

    wn, flux = io.read_spectrum(sfile)
    np.testing.assert_allclose(wn, np.array(
      [9090.90909091, 8333.33333333, 7692.30769231, 7142.85714286,
       6666.66666667, 6250.        , 5882.35294118]), rtol=1e-7)
    np.testing.assert_equal(flux, np.ones(7))
    wl, flux = io.read_spectrum(sfile, wn=False)
    np.testing.assert_almost_equal(wl, np.linspace(1.1, 1.7, 7))
    np.testing.assert_equal(flux, np.ones(7))


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
    atm = "{}/{}".format(tmpdir, atmfile)
    nlayers = 11
    pressure    = np.logspace(-8, 2, nlayers)
    temperature = np.tile(1500.0, nlayers)
    species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    qprofiles = pa.uniform(pressure, temperature, species, abundances)
    pa.writeatm(atm, pressure, temperature, species, qprofiles,
                punits='bar', header='# Test write atm\n')
    assert atmfile in os.listdir(str(tmpdir))

    atm_input = pa.readatm(atm)
    assert atm_input[0] == ('bar', 'kelvin', 'number', None)
    np.testing.assert_equal(atm_input[1], np.array(species))
    np.testing.assert_almost_equal(atm_input[2], pressure/pc.bar)
    np.testing.assert_almost_equal(atm_input[3], temperature)
    np.testing.assert_almost_equal(atm_input[4], qprofiles)
    assert atm_input[5] is None


def test_read_atm_no_temp():
    with pytest.raises(ValueError, match="Atmospheric file does not have "
                       "'@TEMPERATURE' header"):
        atm_input = pa.readatm('uniform_notemp_test.atm')


def test_read_write_pf(tmpdir):
    pffile = 'PF_Exomol_NH3_test.dat'
    pff = "{}/{}".format(tmpdir, pffile)
    isotopes = ['4111', '5111']
    temp   = np.linspace(10,100,4)
    pf     = np.array([np.logspace(0,3,4), np.logspace(1,4,4)])
    header = '# Mock partition function for NH3.\n'
    io.write_pf(pff, pf, isotopes, temp, header)

    assert pffile in os.listdir(str(tmpdir))
    PF, iso, T = io.read_pf(pff)
    np.testing.assert_allclose(PF, pf, rtol=1e-5)
    np.testing.assert_allclose(temp, T, rtol=1e-7)
    assert list(iso) == isotopes


def test_write_pf_mismatch_iso():
    pffile = 'PF_Exomol_NH3_test.dat'
    isotopes = ['4111']
    temp   = np.linspace(10,100,4)
    pf     = np.array([np.logspace(0,3,4), np.logspace(1,4,4)])
    with pytest.raises(ValueError, match='Shape of the partition-function '
                       'array does not match with the number of isotopes.'):
        io.write_pf(pffile, pf, isotopes, temp)


def test_write_pf_mismatch_temp():
    pffile = 'PF_Exomol_NH3_test.dat'
    isotopes = ['4111', '5111']
    temp   = np.linspace(10,100,5)
    pf     = np.array([np.logspace(0,3,4), np.logspace(1,4,4)])
    with pytest.raises(ValueError, match='Shape of the partition-function '
            'array does not match with the number of temperature samples.'):
        io.write_pf(pffile, pf, isotopes, temp)

@pytest.mark.parametrize('species',
    [['CH4'], ['H2', 'H2']])
def test_read_write_cs(species, tmpdir):
    csfile = 'CS_Mock.dat'
    csf = "{}/{}".format(tmpdir, csfile)
    temp = np.linspace(100, 1000, 3)
    wn   = np.arange(10, 15, 1.0)
    cs   = np.array([np.logspace( 0,-4,5),
                     np.logspace(-1,-5,5),
                     np.logspace(-2,-6,5)])
    header = '# Mock cross-section.\n'
    io.write_cs(csf, cs, species, temp, wn, header)
    assert csfile in os.listdir(str(tmpdir))
    # TBD: assert file mentions the right opacity units

    cross_sec, specs, t, wave = io.read_cs(csf)
    np.testing.assert_allclose(cs, cross_sec, rtol=1e-5)
    np.testing.assert_allclose(temp, t,  rtol=1e-7)
    np.testing.assert_allclose(wn, wave, rtol=1e-7)
    assert specs == species


def test_write_cs_mismatch_temp():
    csfile = 'CS_Mock.dat'
    species = ['H2', 'H2']
    temp = np.linspace(100, 1000, 10)
    wn   = np.arange(10, 15, 1.0)
    cs   = np.array([np.logspace( 0,-4,5),
                     np.logspace(-1,-5,5),
                     np.logspace(-2,-6,5)])
    with pytest.raises(ValueError, match='Shape of the cross-section array '
                       'does not match the number of temperature samples.'):
        io.write_cs(csfile, cs, species, temp, wn)


def test_write_cs_mismatch_wn():
    csfile = 'CS_Mock.dat'
    species = ['H2', 'H2']
    temp = np.linspace(100, 1000, 3)
    wn   = np.arange(10, 15, 0.5)
    cs   = np.array([np.logspace( 0,-4,5),
                     np.logspace(-1,-5,5),
                     np.logspace(-2,-6,5)])
    with pytest.raises(ValueError, match='Shape of the cross-section array '
                       'does not match the number of wavenumber samples.'):
        io.write_cs(csfile, cs, species, temp, wn)


def test_read_write_filter(tmpdir):
    ffile = 'filter.dat'
    ff = "{}/{}".format(tmpdir, ffile)
    wl = np.linspace(1.4, 1.5, 20)
    transmission = np.array(np.abs(wl-1.45) < 0.035, np.double)
    io.write_filter(ff, wl, transmission)
    assert ffile in os.listdir(str(tmpdir))

    wn, trans = io.read_filter(ff)
    np.testing.assert_equal(trans, transmission)
    np.testing.assert_allclose(1e4/wn, np.array(
      [1.4    , 1.40526, 1.41053, 1.41579, 1.42105, 1.42632, 1.43158,
       1.43684, 1.44211, 1.44737, 1.45263, 1.45789, 1.46316, 1.46842,
       1.47368, 1.47895, 1.48421, 1.48947, 1.49474, 1.5    ]), atol=1e-7)
