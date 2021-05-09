# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os
import pytest

import numpy as np

from conftest import make_config

import pyratbay as pb
import pyratbay.constants as pc
import pyratbay.io as io
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


# Expected spectra:
keys = [
    'lec', 'cia', 'alkali', 'deck', 'tli', 'patchy', 'patchy_clear',
    'patchy_cloudy', 'all', 'etable',
    'tmodel', 'vert', 'scale', 'fit1', 'fit2', 'fit3', 'fit4',
    'bandflux4', 'resolution', 'odd_even']
expected = {
    key:np.load(f"{ROOT}tests/expected/"
                f"expected_spectrum_transmission_{key}_test.npz")['arr_0']
    for key in keys}
#np.savez('expected/expected_spectrum_transmission_fit_test.npz', model1[0])

INPUTS = f'{ROOT}tests/inputs/'
OUTPUTS = f'{ROOT}tests/outputs/'

# Relative tolerance of less than 0.01% difference:
rtol = 0.01 / 100.0

def test_transmission_clear(tmp_path):
    # No opacity whatsoever:
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'alkali', 'clouds'])
    pyrat = pb.run(cfg)
    depth_bottom = (pyrat.atm.radius[-1] / pyrat.phy.rstar)**2
    np.testing.assert_allclose(pyrat.spec.spectrum, depth_bottom, rtol=rtol)


def test_transmission_lecavelier(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'alkali', 'clouds'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['lec'], rtol=rtol)


def test_transmission_CIA(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'rayleigh', 'alkali', 'clouds'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['cia'], rtol=rtol)


def test_transmission_alkali(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'clouds'],
        reset={'wllow':'0.45 um', 'wlhigh':'1.0 um'})
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['alkali'], rtol=rtol)


def test_transmission_deck(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'alkali'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['deck'], rtol=rtol)


def test_transmission_tli(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['csfile', 'rayleigh', 'clouds', 'alkali'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['tli'], rtol=rtol)


def test_transmission_all(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=rtol)


def test_transmission_patchy(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'fpatchy': '0.5',
               'rpars': '10.0 -15.0'})
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['patchy'], rtol=rtol)
    np.testing.assert_allclose(
        pyrat.spec.clear, expected['patchy_clear'], rtol=rtol)
    np.testing.assert_allclose(
        pyrat.spec.cloudy, expected['patchy_cloudy'], rtol=rtol)

    pyrat.cloud.fpatchy = 0.0
    pyrat.run()
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['patchy_clear'], rtol=rtol)

    pyrat.cloud.fpatchy = 1.0
    pyrat.run()
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['patchy_cloudy'], rtol=rtol)


def test_transmission_resolution(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'resolution':'5000.0'},
        remove=['clouds'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['resolution'], rtol=rtol)


# Optical-depth integration is home made, which depends on whether there is
# an odd or even number of layers. Thus, the need for this test.
def test_transmission_odd_even(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'rpars':'1.0 -4.0'},
        remove=['tlifile', 'csfile', 'alkali', 'clouds'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['odd_even'], rtol=rtol)

    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'atmfile':f'{INPUTS}atmosphere_uniform_even_layers.atm',
               'rpars':'1.0 -4.0'},
        remove=['tlifile', 'csfile', 'alkali', 'clouds'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['odd_even'], rtol=rtol)


def test_transmission_etable(tmp_path):
    # LBL from extinction table:
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'clouds'],
        reset={'extfile':f'{OUTPUTS}exttable_test_300-3000K_1.1-1.7um.npz'})
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['etable'], rtol=rtol)


def test_transmission_input_radius(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'atmfile':f'{INPUTS}atmosphere_uniform_radius.atm'},
        remove=['radmodel'])
    pyrat = pb.run(cfg)
    atm = io.read_atm('inputs/atmosphere_uniform_radius.atm')
    np.testing.assert_allclose(pyrat.atm.radius, atm[5]*pc.km, rtol=rtol)


def test_transmission_input_radius_overwrite(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'atmfile':f'{INPUTS}/atmosphere_uniform_radius.atm'})
    pyrat = pb.run(cfg)
    atm = io.read_atm('inputs/atmosphere_uniform_radius.atm')
    np.testing.assert_raises(AssertionError, np.testing.assert_array_equal,
        pyrat.atm.radius, atm[5]*pc.km)


@pytest.mark.skip(reason="TBI")
def test_transmission_qmass_input():
    # This is the gist of it, prepare a qmass atmospheric file:
    atm = io.read_atm('inputs/atmosphere_uniform_test.atm')
    units, species, press, temp, q = atm
    symbol, mass, diam = io.read_molecs(pc.ROOT+"pyratbay/data/molecules.dat")
    mm = pa.mean_weight(q, species)
    qmass = qprofiles * molmass / mm
    io.write_atm(atmfile)
    # Then run spectrum, results must be the same as qnumber run:
    pyrat = pb.run(ROOT+'tests/configs/spectrum_transmission_qmass_test.cfg')
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=rtol)


# Now try some forward models that modify the atmospheric profile:
def test_transmission_tmodel_none(tmp_path):
    # include tmodel, but tpars is None
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'tmodel':'tcea'})
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=1e-7)
    # Now, re-run with user-input tpars:
    pyrat.atm.tpars = np.array([-4.67, -0.8, -0.8, 0.5, 1486.0, 100.0])
    pyrat.run()
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['tmodel'], rtol=rtol)


def test_transmission_tmodel(tmp_path):
    # Include tmodel and tpars in input config file:
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'tmodel':'tcea', 'tpars':'-4.67 -0.8 -0.8 0.5 1486.0 100.0'})
    pyrat = pb.run(cfg)
    tmodel2 = pyrat.spec.spectrum
    np.testing.assert_allclose(tmodel2, expected['tmodel'], rtol=rtol)


def test_transmission_vert_none_model(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'molmodel':'vert', 'molfree':'H2O', 'bulk':'H2 He'})
    pyrat = pb.run(cfg)
    vmodel0 = pyrat.spec.spectrum
    np.testing.assert_allclose(vmodel0, expected['all'], rtol=rtol)
    pyrat.atm.molpars = [-5]
    pyrat.run()
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['vert'], rtol=rtol)
    np.testing.assert_equal(pyrat.atm.q[:,3], 1e-5)


def test_transmission_vert_model(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'molmodel':'vert', 'molfree':'H2O', 'molpars':'-5',
               'bulk':'H2 He'})
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['vert'], rtol=rtol)


def test_transmission_scale_model(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'molmodel':'scale', 'molfree':'H2O', 'molpars':'-1',
               'bulk':'H2 He'})
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['scale'], rtol=rtol)
    np.testing.assert_equal(pyrat.atm.q[:,3], 0.1*pyrat.atm.qbase[:,3])


def test_transmission_fit(tmp_path):
    # Without evaulating params:
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'tmodel':'tcea', 'cpars':'2.0',
               'molmodel':'vert', 'molfree':'H2O', 'bulk':'H2 He',
               'retflag':'temp mol ray cloud',
               'params':'-4.67 -0.8 -0.8 0.5 1486.0 100.0 -4.0 0.0 -4.0 2.0'})
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=rtol)
    # Eval default params:
    model1 = pyrat.eval(pyrat.ret.params, retmodel=True)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['fit1'], rtol=rtol)
    np.testing.assert_equal(model1[0], pyrat.spec.spectrum)
    assert model1[1] is None
    # Cloud deck:
    params = [-4.67, -0.8, -0.8, 0.5, 1486.0, 100.0, -4.0, 0.0, -4.0, -3.0]
    model2 = pyrat.eval(params, retmodel=True)
    rmin = np.amin(np.sqrt(pyrat.spec.spectrum)) * pyrat.phy.rstar
    rexpected = pyrat.cloud.models[0].rsurf
    np.testing.assert_allclose(rmin, rexpected, rtol=rtol)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['fit2'], rtol=rtol)
    # Check pyrat.ret.params has been updated:
    np.testing.assert_equal(pyrat.ret.params, params)
    # Depleted H2O:
    params = [-4.67, -0.8, -0.8, 0.5, 1486.0, 100.0, -8.0, 0.0, -4.0, 2.0]
    model3 = pyrat.eval(params, retmodel=True)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['fit3'], rtol=rtol)


def test_transmission_fit_filters():
    pyrat = pb.run(ROOT+'tests/configs/spectrum_transmission_filters_test.cfg')
    model4 = pyrat.eval(pyrat.ret.params, retmodel=True)
    np.testing.assert_allclose(model4[0], expected['fit4'], rtol=rtol)
    np.testing.assert_allclose(model4[1], expected['bandflux4'], rtol=rtol)


def test_multiple_opacities(tmp_path):
    # Generate TLI files:
    cfg = make_config(tmp_path, f'{ROOT}tests/configs/tli_multiple_opacity.cfg')
    pyrat = pb.run(cfg)
    cfg = make_config(tmp_path,
        f'{ROOT}tests/configs/tli_multiple_opacity.cfg',
        reset={'dblist':f'{INPUTS}02_hit12.par',
               'tlifile':f'{OUTPUTS}HITRAN_CO2_1.5-1.6um_test.tli'})
    pyrat = pb.run(cfg)
    cfg = make_config(tmp_path,
        f'{ROOT}tests/configs/tli_multiple_opacity.cfg',
        reset={'dblist':f'{INPUTS}06_hit12.par',
               'tlifile':f'{OUTPUTS}HITRAN_CH4_1.5-1.6um_test.tli'})
    pyrat = pb.run(cfg)

    # Generate opacity files:
    cfg = make_config(tmp_path, ROOT+'tests/configs/opacity_multiple.cfg')
    pyrat = pb.run(cfg)
    assert pyrat is not None
    cfg = make_config(tmp_path,
        f'{ROOT}tests/configs/opacity_multiple.cfg',
        reset={'tlifile':f'{OUTPUTS}HITRAN_CO2_1.5-1.6um_test.tli',
               'extfile':f'{OUTPUTS}exttable_CO2_300-3000K_1.5-1.6um.npz'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    cfg = make_config(tmp_path,
        f'{ROOT}tests/configs/opacity_multiple.cfg',
        reset={'tlifile':f'{OUTPUTS}HITRAN_CH4_1.5-1.6um_test.tli',
               'extfile':f'{OUTPUTS}exttable_CH4_300-3000K_1.5-1.6um.npz'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    cfg = make_config(tmp_path,
        f'{ROOT}tests/configs/opacity_multiple.cfg',
        reset={'tlifile':f'{OUTPUTS}HITRAN_CO2_1.5-1.6um_test.tli'
                   f'\n    {OUTPUTS}HITRAN_CH4_1.5-1.6um_test.tli',
               'extfile':f'{OUTPUTS}exttable_CO2-CH4_300-3000K_1.5-1.6um.npz'})
    pyrat = pb.run(cfg)
    assert pyrat is not None

    # Compute spectra from opacities:
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'clouds'],
        reset={'extfile':f'{OUTPUTS}exttable_H2O_300-3000K_1.5-1.6um.npz'
                     f'\n  {OUTPUTS}exttable_CO2_300-3000K_1.5-1.6um.npz'
                     f'\n  {OUTPUTS}exttable_CH4_300-3000K_1.5-1.6um.npz',
               'wllow':'1.5 um', 'wlhigh':'1.6 um'})
    pyrat1 = pb.run(cfg)
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'clouds'],
        reset={'extfile':f'{OUTPUTS}exttable_H2O_300-3000K_1.5-1.6um.npz'
                     f'\n  {OUTPUTS}exttable_CO2-CH4_300-3000K_1.5-1.6um.npz',
               'wllow':'1.5 um', 'wlhigh':'1.6 um'})
    pyrat2 = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat1.spec.spectrum, pyrat2.spec.spectrum, rtol=rtol)


@pytest.mark.parametrize('wllow,wlhigh',
    [('1.1 um', '1.6 um'),
     ('1.2 um', '1.7 um'),
     ('1.2 um', '1.6 um')])
def test_opacity_reset_wn(tmp_path, wllow, wlhigh):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'clouds'],
        reset={'extfile':f'{OUTPUTS}exttable_test_300-3000K_1.1-1.7um.npz',
               'wllow':wllow, 'wlhigh':wlhigh})
    pyrat = pb.run(cfg)
    wn = np.arange(1/1.7e-4, 1/1.1e-4, 1.0)
    etab = expected['etable'][(wn>= pyrat.spec.wnlow)
                            & (wn <=pyrat.spec.wnhigh)]
    np.testing.assert_allclose(pyrat.spec.spectrum, etab, rtol=rtol)


# These are extra bits for testing the tests before testing:
def spectrum_fm():
    from scipy.ndimage.filters import gaussian_filter1d as gaussf
    pyrat = pb.run(ROOT+'tests/configs/spectrum_transmission_filters_test.cfg')
    params = [-1.5, -0.8, 0.0,  1.0,  1.0,  71500.0, -3.4, 2.0]
    model = pyrat.eval(params, retmodel=True)
    bandflux = pyrat.obs.bandflux

    plt.figure(1)
    plt.clf()
    plt.plot(1e4/pyrat.spec.wn, pyrat.spec.spectrum, 'blue')
    plt.plot(1e4/pyrat.spec.wn, gaussf(pyrat.spec.spectrum, 5), 'navy')
    plt.plot(1e4/pyrat.obs.bandwn, bandflux, "o", c='orange')

    # Now, add noise:
    SNR = 350.0
    std = 35.0
    np.random.seed(209458)
    snr = np.random.normal(SNR, std, len(bandflux))
    uncert = bandflux/snr

    data = bandflux + np.random.normal(0.0, uncert)
    plt.errorbar(1e4/pyrat.obs.bandwn, data, uncert, fmt='or', zorder=10)


def spectrum_exomol():
    pyrat = pb.run(ROOT+'tests/configs/spectrum_transmission_exomol.cfg')

    plt.figure(1)
    plt.clf()
    plt.plot(1e4/hcn[0], hcn[1], 'b')
    plt.plot(1e4/pyrat.spec.wn, pyrat.spec.spectrum, 'orange', alpha=0.5)

