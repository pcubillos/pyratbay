# Copyright (c) 2021-2022 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

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
    'bandflux4', 'resolution', 'odd_even',
]

root = f"{ROOT}tests/expected/expected_spectrum_transmission"
expected = {
    key: np.load(f"{root}_{key}_test.npz")['arr_0']
    for key in keys
}
#np.savez('expected/expected_spectrum_transmission_fit_test.npz', model1[0])

INPUTS = f'{ROOT}tests/inputs/'
OUTPUTS = f'{ROOT}tests/outputs/'

# Relative tolerance of less than 0.01% difference:
rtol = 0.01 / 100.0

def test_transmission_clear(tmp_path):
    # No opacity whatsoever:
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'alkali', 'clouds'])
    pyrat = pb.run(cfg)
    depth_bottom = (pyrat.atm.radius[-1] / pyrat.phy.rstar)**2
    np.testing.assert_allclose(pyrat.spec.spectrum, depth_bottom, rtol=rtol)


def test_transmission_lecavelier(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'alkali', 'clouds'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['lec'], rtol=rtol)


def test_transmission_CIA(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'rayleigh', 'alkali', 'clouds'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['cia'], rtol=rtol)


def test_transmission_alkali(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'clouds'],
        reset={'wllow':'0.45 um', 'wlhigh':'1.0 um'})
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['alkali'], rtol=rtol)


def test_transmission_deck(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'alkali'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['deck'], rtol=rtol)


def test_transmission_tli(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['csfile', 'rayleigh', 'clouds', 'alkali'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['tli'], rtol=rtol)


def test_transmission_all(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=rtol)


def test_transmission_patchy(tmp_path):
    reset = {
        'fpatchy': '0.5',
        'rpars': '10.0 -15.0',
    }
    cfg_file = ROOT+'tests/configs/spectrum_transmission_test.cfg'
    cfg = make_config(tmp_path, cfg_file, reset=reset)
    pyrat = pb.run(cfg)

    spectrum = pyrat.spec.spectrum
    clear = pyrat.spec.clear
    cloudy = pyrat.spec.cloudy
    np.testing.assert_allclose(spectrum, expected['patchy'], rtol=rtol)
    np.testing.assert_allclose(clear, expected['patchy_clear'], rtol=rtol)
    np.testing.assert_allclose(cloudy, expected['patchy_cloudy'], rtol=rtol)

    pyrat.cloud.fpatchy = 0.0
    pyrat.run()
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['patchy_clear'], rtol=rtol)

    pyrat.cloud.fpatchy = 1.0
    pyrat.run()
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['patchy_cloudy'], rtol=rtol)


def test_transmission_resolution(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'resolution':'5000.0'},
        remove=['clouds'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['resolution'], rtol=rtol)


@pytest.mark.skip(reason="TBD")
def test_transmission_wl_step(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'wl_step':'1e-5 um'},
        remove=['clouds'],
    )
    pyrat = pb.run(cfg)
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['wl_step'], rtol=rtol)


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
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['odd_even'], rtol=rtol)


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
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'atmfile':f'{INPUTS}/atmosphere_uniform_radius.atm'},
    )
    pyrat = pb.run(cfg)
    atm = io.read_atm('inputs/atmosphere_uniform_radius.atm')
    np.testing.assert_raises(
        AssertionError,
        np.testing.assert_array_equal,
        pyrat.atm.radius,
        atm[5]*pc.km,
    )


# Now try some forward models that modify the atmospheric profile:
def test_transmission_tmodel_none(tmp_path):
    # include tmodel, but tpars is None
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'tmodel': 'guillot'},
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=rtol)
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
        reset={'tmodel':'guillot', 'tpars':'-4.67 -0.8 -0.8 0.5 1486.0 100.0'})
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
    reset = {
        'molmodel': 'vert',
        'molfree': 'H2O',
        'molpars': '-5',
        'bulk': 'H2 He',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset=reset,
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['vert'], rtol=rtol)


def test_transmission_scale_model(tmp_path):
    reset = {
        'molmodel': 'scale',
        'molfree': 'H2O',
        'molpars': '-1',
        'bulk': 'H2 He',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset=reset,
    )
    pyrat = pb.run(cfg)
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['scale'], rtol=rtol)
    np.testing.assert_equal(pyrat.atm.q[:,3], 0.1*pyrat.atm.qbase[:,3])


def test_transmission_fit(tmp_path):
    # Without evaulating params:
    reset = {
        'tmodel': 'guillot',
        'cpars': '2.0',
        'molmodel': 'vert',
        'molfree': 'H2O',
        'bulk': 'H2 He',
        'retflag': 'temp mol ray cloud',
        'params': '-4.67 -0.8 -0.8 0.5 1486.0 100.0 -4.0 0.0 -4.0 2.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset=reset,
    )
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
    tli_cfg = f'{ROOT}tests/configs/tli_multiple_opacity.cfg'
    cfg = make_config(tmp_path, tli_cfg)
    pyrat = pb.run(cfg)

    reset = {
        'dblist': f'{INPUTS}02_hit12.par',
        'tlifile': f'{OUTPUTS}HITRAN_CO2_1.5-1.6um_test.tli',
    }
    cfg = make_config(tmp_path, tli_cfg, reset=reset)
    pyrat = pb.run(cfg)

    reset = {
        'dblist': f'{INPUTS}06_hit12.par',
        'tlifile': f'{OUTPUTS}HITRAN_CH4_1.5-1.6um_test.tli',
    }
    cfg = make_config(tmp_path, tli_cfg, reset=reset)
    pyrat = pb.run(cfg)

    # Generate opacity files:
    opac_cfg = f'{ROOT}tests/configs/opacity_multiple.cfg'
    cfg = make_config(tmp_path, opac_cfg)
    pyrat = pb.run(cfg)
    assert pyrat is not None

    reset = {
        'tlifile': f'{OUTPUTS}HITRAN_CO2_1.5-1.6um_test.tli',
        'extfile': f'{OUTPUTS}exttable_CO2_300-3000K_1.5-1.6um.npz',
    }
    cfg = make_config(tmp_path, opac_cfg, reset=reset)
    pyrat = pb.run(cfg)
    assert pyrat is not None

    reset = {
        'tlifile':f'{OUTPUTS}HITRAN_CH4_1.5-1.6um_test.tli',
        'extfile': f'{OUTPUTS}exttable_CH4_300-3000K_1.5-1.6um.npz',
    }
    cfg = make_config(tmp_path, opac_cfg, reset=reset)
    pyrat = pb.run(cfg)
    assert pyrat is not None

    reset = {
        'tlifile':
            f'{OUTPUTS}HITRAN_CO2_1.5-1.6um_test.tli'
            f'\n    {OUTPUTS}HITRAN_CH4_1.5-1.6um_test.tli',
        'extfile': f'{OUTPUTS}exttable_CO2-CH4_300-3000K_1.5-1.6um.npz',
    }
    cfg = make_config(tmp_path, opac_cfg, reset=reset)
    pyrat = pb.run(cfg)
    assert pyrat is not None

    # Compute spectra from opacities:
    spec_cfg = f'{ROOT}tests/configs/spectrum_transmission_test.cfg'
    reset = {
        'extfile':
            f'{OUTPUTS}exttable_H2O_300-3000K_1.5-1.6um.npz'
            f'\n  {OUTPUTS}exttable_CO2_300-3000K_1.5-1.6um.npz'
            f'\n  {OUTPUTS}exttable_CH4_300-3000K_1.5-1.6um.npz',
        'wllow':'1.5 um',
        'wlhigh':'1.6 um',
    }
    cfg = make_config(
        tmp_path, spec_cfg, reset=reset, remove=['tlifile','clouds'],
    )
    pyrat = pb.run(cfg)
    spectrum1 = pyrat.spec.spectrum

    reset = {
        'extfile':
            f'{OUTPUTS}exttable_H2O_300-3000K_1.5-1.6um.npz'
            f'\n  {OUTPUTS}exttable_CO2-CH4_300-3000K_1.5-1.6um.npz',
        'wllow':'1.5 um',
        'wlhigh':'1.6 um',
    }
    cfg = make_config(
        tmp_path, spec_cfg, reset=reset, remove=['tlifile','clouds'],
    )
    pyrat = pb.run(cfg)
    spectrum2 = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum1, spectrum2, rtol=rtol)


@pytest.mark.parametrize('wllow,wlhigh',
    [('1.1 um', '1.6 um'),
     ('1.2 um', '1.7 um'),
     ('1.2 um', '1.6 um')])
def test_opacity_reset_wn(tmp_path, wllow, wlhigh):
    reset = {
        'extfile': f'{OUTPUTS}exttable_test_300-3000K_1.1-1.7um.npz',
        'wllow': wllow,
        'wlhigh': wlhigh,
    }
    cfg = make_config(
        tmp_path,
        f'{ROOT}tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'clouds'],
        reset=reset,
    )
    pyrat = pb.run(cfg)
    wn = np.arange(1/1.7e-4, 1/1.1e-4, 1.0)
    wn_range = (wn>= pyrat.spec.wnlow) & (wn <=pyrat.spec.wnhigh)
    etab = expected['etable'][wn_range]
    np.testing.assert_allclose(pyrat.spec.spectrum, etab, rtol=rtol)

