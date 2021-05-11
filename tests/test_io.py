# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os
import pytest

import numpy as np

from conftest import make_config

import pyratbay as pb
import pyratbay.io as io
import pyratbay.atmosphere as pa
import pyratbay.constants as pc
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


def test_load_save_pyrat(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'cpars': '1.0'})
    pyrat = pb.run(cfg)
    spectrum = np.copy(pyrat.spec.spectrum)
    io.save_pyrat(pyrat)

    pfile = pyrat.log.logname.replace('.log', '.pickle')
    new_pyrat = io.load_pyrat(pfile)
    # Check previous spectrum value stand:
    np.testing.assert_equal(new_pyrat.spec.spectrum, spectrum)
    # Check re-run reproduces previous spectrum:
    new_pyrat.run()
    np.testing.assert_equal(new_pyrat.spec.spectrum, spectrum)



def test_read_write_spectrum(tmpdir):
    sfile = "spectrum_test.dat"
    tmp_file = "{}/{}".format(tmpdir, sfile)

    wl = np.linspace(1.1, 1.7, 7) * 1e-4
    spectrum = np.ones(7)
    io.write_spectrum(
        wl, spectrum, filename=tmp_file, type='transit', wlunits='um')
    # Take a look at the output file:
    assert sfile in os.listdir(str(tmpdir))
    with open(tmp_file, 'r') as f:
        assert f.readline() == '# Wavelength        (Rp/Rs)**2\n'
        assert f.readline() == '#         um          unitless\n'

    wn, flux = io.read_spectrum(tmp_file)
    np.testing.assert_allclose(wn, np.array(
      [9090.90909091, 8333.33333333, 7692.30769231, 7142.85714286,
       6666.66666667, 6250.        , 5882.35294118]), rtol=1e-7)
    np.testing.assert_equal(flux, np.ones(7))
    wl, flux = io.read_spectrum(tmp_file, wn=False)
    np.testing.assert_allclose(wl, np.linspace(1.1, 1.7, 7))
    np.testing.assert_equal(flux, np.ones(7))


def test_read_write_spectrum_filter(tmpdir):
    ffile = 'filter_test.dat'
    tmp_file = "{}/{}".format(tmpdir, ffile)

    # Assert write:
    wl = np.linspace(1.4, 1.5, 20)
    transmission = np.array(np.abs(wl-1.45) < 0.035, np.double)
    io.write_spectrum(wl*pc.um, transmission, tmp_file, 'filter')
    assert ffile in os.listdir(str(tmpdir))
    with open(tmp_file, 'r') as f:
        assert f.readline() == '# Wavelength      transmission\n'
        assert f.readline() == '#         um          unitless\n'

    # Assert read:
    wn, trans = io.read_spectrum(tmp_file)
    np.testing.assert_equal(trans, transmission)
    np.testing.assert_allclose(1e4/wn, np.array(
      [1.4    , 1.40526, 1.41053, 1.41579, 1.42105, 1.42632, 1.43158,
       1.43684, 1.44211, 1.44737, 1.45263, 1.45789, 1.46316, 1.46842,
       1.47368, 1.47895, 1.48421, 1.48947, 1.49474, 1.5    ]), atol=1e-7)


def test_write_spectrum_bad_type():
    wl = np.linspace(1.4, 1.5, 20)
    trans = np.ones(20)
    match = "Input 'type' argument must be 'transit', 'emission', or 'filter'."
    with pytest.raises(ValueError, match=match):
        io.write_spectrum(wl, trans, "tophat_filter.dat", 'bad_type')


@pytest.mark.parametrize('header',
    ['# Wavelength flux\n',
     '# Something something\n# Wavelength flux\n# um funits\n',
     '',
     '\n',
     '# Wavelength flux\n\n# um funits\n',
     '# um funits\n',
     '# um\n'
    ])
def test_read_spectrum_custom_header(tmpdir, header):
    ffile = 'filter_test.dat'
    tmp_file = "{}/{}".format(tmpdir, ffile)
    data = '1.0 1.0\n2.0 1.5\n4.0 1.0\n'
    with open(tmp_file, 'w') as f:
        f.write(header + data)
    wn, spec = io.read_spectrum(tmp_file)
    np.testing.assert_allclose(wn, np.array([10000.0, 5000.0, 2500.0]))


def test_read_write_opacity(tmpdir):
    ofile = "{}/opacity_test.npz".format(tmpdir)
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

    np.testing.assert_allclose(molID, edata[1][0])
    np.testing.assert_allclose(temp,  edata[1][1])
    np.testing.assert_allclose(press, edata[1][2])
    np.testing.assert_allclose(wn,    edata[1][3])
    np.testing.assert_allclose(etable, edata[2])


def test_read_write_atm_pt(tmpdir):
    atmfile = "WASP-00b.atm"
    atm = "{}/{}".format(tmpdir, atmfile)
    nlayers = 11
    pressure    = pa.pressure('1e-8 bar', '1e2 bar', nlayers)
    temperature = pa.tmodels.Isothermal(nlayers)(1500.0)
    io.write_atm(atm, pressure, temperature, punits='bar')
    assert atmfile in os.listdir(str(tmpdir))

    atm_input = io.read_atm(atm)
    assert atm_input[0] == ('bar', 'kelvin', None, None)
    assert atm_input[1] is None
    np.testing.assert_allclose(atm_input[2], pressure/pc.bar)
    np.testing.assert_allclose(atm_input[3], temperature)
    assert atm_input[4] is None
    assert atm_input[5] is None


def test_read_write_atm_ptq(tmpdir):
    atmfile = "WASP-00b.atm"
    atm = "{}/{}".format(tmpdir, atmfile)
    nlayers = 11
    pressure    = pa.pressure('1e-8 bar', '1e2 bar', nlayers)
    temperature = pa.tmodels.Isothermal(nlayers)(1500.0)
    species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    qprofiles = pa.uniform(pressure, temperature, species, abundances)
    io.write_atm(atm, pressure, temperature, species, qprofiles,
                 punits='bar', header='# Test write atm\n')
    assert atmfile in os.listdir(str(tmpdir))

    atm_input = io.read_atm(atm)
    assert atm_input[0] == ('bar', 'kelvin', 'volume', None)
    np.testing.assert_equal(atm_input[1], np.array(species))
    np.testing.assert_allclose(atm_input[2], pressure/pc.bar)
    np.testing.assert_allclose(atm_input[3], temperature)
    np.testing.assert_allclose(atm_input[4], qprofiles)
    assert atm_input[5] is None


def test_read_write_atm_ptqr(tmpdir):
    atmfile = "WASP-00b.atm"
    atm = "{}/{}".format(tmpdir, atmfile)
    nlayers = 11
    pressure    = pa.pressure('1e-8 bar', '1e2 bar', nlayers)
    temperature = pa.tmodels.Isothermal(nlayers)(1500.0)
    species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    qprofiles = pa.uniform(pressure, temperature, species, abundances)
    radius = pa.hydro_g(pressure, temperature, 2.3, 2479.0, pc.bar, pc.rjup)
    io.write_atm(atm, pressure, temperature, species, qprofiles,
        radius=radius, punits='bar', runits='km',
        header='# Test write atm\n')
    assert atmfile in os.listdir(str(tmpdir))

    atm_input = io.read_atm(atm)
    assert atm_input[0] == ('bar', 'kelvin', 'volume', 'km')
    np.testing.assert_equal(atm_input[1], np.array(species))
    np.testing.assert_allclose(atm_input[2], pressure/pc.bar)
    np.testing.assert_allclose(atm_input[3], temperature)
    np.testing.assert_allclose(atm_input[4], qprofiles)
    np.testing.assert_allclose(atm_input[5], radius/pc.km, rtol=1e-5)


def test_read_atm_no_temp():
    match = "Atmospheric file does not have '@TEMPERATURE' header"
    with pytest.raises(ValueError, match=match):
        dummy = io.read_atm('inputs/uniform_notemp_test.atm')


def test_read_write_pf(tmpdir):
    pffile = 'PF_Exomol_NH3_test.dat'
    pff = "{}/{}".format(tmpdir, pffile)
    isotopes = ['4111', '5111']
    temp   = np.linspace(10, 100, 4)
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


def test_read_pt(tmpdir):
    ptfile = 'mock_pt.dat'
    ptf = "{}/{}".format(tmpdir, ptfile)
    temp  = np.array([100.0, 150.0, 200.0, 175.0, 150.0])
    press = np.array([1e-6,  1e-4,  1e-2,  1e0,   1e2])
    with open(ptf, 'w') as f:
        for p,t in zip(press, temp):
            f.write('{:.3e}  {:5.1f}\n'.format(p, t))
    pressure, temperature = io.read_pt(ptf)
    np.testing.assert_allclose(pressure, press*pc.bar,  rtol=1e-7)
    np.testing.assert_allclose(temperature, temp,  rtol=1e-7)


def test_read_molecs():
    names, mass, diam = io.read_molecs(ROOT+'pyratbay/data/molecules.dat')
    assert 'H2O' in names
    assert 'CH4' in names
    assert 'CO' in names
    assert 'CO2' in names
    assert 'H2' in names
    np.testing.assert_allclose(mass[names == 'H2O'], 18.01528)
    np.testing.assert_allclose(diam[names == 'H2'], 2.89)


def test_read_isotopes():
    ID, mol, hit_iso, exo_iso, ratio, mass = io.read_isotopes(
        pc.ROOT+'pyratbay/data/isotopes.dat')
    assert 'H2O' in mol
    assert np.all(ID[mol=='H2O'] == 1)
    assert list(hit_iso[mol=='H2O']) == \
        ['161', '181', '171', '162', '182', '172', '262', '282', '272']
    assert list(exo_iso[mol=='H2O']) == \
        ['116', '118', '117', '126', '128', '127', '226', '228', '227']
    np.testing.assert_allclose(ratio[mol=='H2O'],
        np.array([9.973e-01, 1.999e-03, 3.719e-04, 3.107e-04, 6.230e-07,
                  1.158e-07, 2.420e-08, 0.000e+00, 0.000e+00]))
    np.testing.assert_allclose(mass[mol=='H2O'],
        np.array([18.010565, 20.014811, 19.014781, 19.016841, 21.021088,
                  20.021058, 20.021,    22.0000, 21.0000]))

@pytest.mark.skip(
    reason='This requires either to download a huge file or mock it up.')
def test_import_xs():
    # wget this file first: http://www.exomol.com/db/H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL__R15000_0.3-50mu.xsec.TauREx.h5
    filename = '1H2-16O__POKAZATEL__R15000_0.3-50mu.xsec.TauREx.h5'
    xs_H2O, press, temp, wn, species = io.import_xs(filename, 'exomol')

