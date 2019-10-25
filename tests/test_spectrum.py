# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

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
keys = ['lec', 'cia', 'alkali', 'deck', 'tli', 'all', 'etable',
        'tmodel', 'vert', 'scale', 'fit1', 'fit2', 'fit3', 'fit4',
        'bandflux4', 'resolution', 'odd_even']
expected = {key:np.load("expected_spectrum_transmission_{:s}_test.npz".
                        format(key))['arr_0']
            for key in keys}
#np.savez('expected_spectrum_transmission_fit_test.npz', model1[0])


def test_transmission_clear(tmp_path):
    # No opacity whatsoever:
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'alkali', 'clouds'])
    clear = pb.run(cfg)
    depth_bottom = (clear.atm.radius[-1] / clear.phy.rstar)**2
    np.testing.assert_allclose(clear.spec.spectrum, depth_bottom, rtol=1e-7)


def test_transmission_lecavelier(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'alkali', 'clouds'])
    ray = pb.run(cfg)
    np.testing.assert_allclose(ray.spec.spectrum, expected['lec'], rtol=1e-7)


def test_transmission_CIA(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['tlifile', 'rayleigh', 'alkali', 'clouds'])
    cia = pb.run(cfg)
    np.testing.assert_allclose(cia.spec.spectrum, expected['cia'], rtol=1e-7)


def test_transmission_alkali(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'clouds'],
        reset={'wllow':'0.45 um', 'wlhigh':'1.0 um'})
    alkali = pb.run(cfg)
    np.testing.assert_allclose(alkali.spec.spectrum, expected['alkali'],
        rtol=1e-7)


def test_transmission_deck(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'alkali'])
    deck = pb.run(cfg)
    np.testing.assert_allclose(deck.spec.spectrum, expected['deck'], rtol=1e-7)


def test_transmission_tli(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['csfile', 'rayleigh', 'clouds', 'alkali'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['tli'], rtol=1e-7)


def test_transmission_all(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['clouds'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=1e-7)


def test_transmission_resolution(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'resolution':'5000.0'},
        remove=['clouds'])
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['resolution'],
        rtol=1e-7)


# Optical-depth integration is home made, which depends on whether there is
# an odd or even number of layers. Thus, the need for this test.
def test_transmission_odd_even(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'rpars':'1.0 -4.0'},
        remove=['tlifile', 'csfile', 'alkali', 'clouds'])
    odd = pb.run(cfg)

    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'atmfile':'atmosphere_uniform_even_layers.atm',
               'rpars':'1.0 -4.0'},
        remove=['tlifile', 'csfile', 'alkali', 'clouds'])
    even = pb.run(cfg)
    np.testing.assert_allclose(odd.spec.spectrum, expected['odd_even'],
        rtol=1e-7)
    np.testing.assert_allclose(even.spec.spectrum, expected['odd_even'],
        rtol=1e-7)


def test_transmission_etable(tmp_path):
    # LBL from extinction table:
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['tlifile', 'clouds'],
        reset={'extfile':'exttable_test_300-3000K_1.1-1.7um.dat'})
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['etable'],
                               rtol=1e-7)


@pytest.mark.parametrize('param', ['refpressure', 'rplanet'])
def test_transmission_input_radius(tmp_path, param):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'atmfile':'atmosphere_uniform_radius.atm'},
        remove=[param, 'mplanet'])
    pyrat = pb.run(cfg)
    atm = io.read_atm('atmosphere_uniform_radius.atm')
    np.testing.assert_allclose(pyrat.atm.radius, atm[5]*pc.km, rtol=1e-7)


@pytest.mark.skip(reason="TBI")
def test_transmission_qmass_input():
    # This is the gist of it, prepare a qmass atmospheric file:
    units, species, press, temp, q = io.read_atm('atmosphere_uniform_test.atm')
    molID, symbol, mass, diam = pa.readmol(pc.ROOT+"inputs/molecules.dat")
    mm = pa.meanweight(q, species)
    qmass = qprofiles * molmass / mm
    io.write_atm(atmfile)
    # Then run spectrum, results must be the same as qnumber run:
    pyrat = pb.run(ROOT+'tests/spectrum_transmission_qmass_test.cfg')
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=1e-7)


def plot_transmission():
    import matplotlib.pyplot as plt
    plt.figure(10)
    plt.clf()
    plt.plot(1e4/clear.spec.wn,  clear.spec.spectrum,  'k')
    plt.plot(1e4/alkali.spec.wn, alkali.spec.spectrum, 'orange')
    plt.plot(1e4/ray.spec.wn,    ray.spec.spectrum,    'b')
    plt.plot(1e4/cia.spec.wn,    cia.spec.spectrum,    'r')
    plt.plot(1e4/deck.spec.wn,   deck.spec.spectrum,   'limegreen')
    plt.plot(1e4/tli.spec.wn,    tli.spec.spectrum,    '0.5')
    plt.plot(1e4/pyrat.spec.wn,  pyrat.spec.spectrum,  'navy', alpha=0.6)
    plt.plot(1e4/epyrat.spec.wn, epyrat.spec.spectrum, 'deepskyblue', alpha=0.6)


# Now try some forward models that modify the atmospheric profile:
def test_transmission_tmodel_none(tmp_path):
    # include tmodel, but tpars is None
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'tmodel':'tcea'})
    pyrat = pb.run(cfg)
    tmodel0 = pyrat.spec.spectrum
    np.testing.assert_allclose(tmodel0, expected['all'], rtol=1e-7)
    # Now, re-run with user-input tpars:
    pyrat.atm.tpars = np.array([-1.5, -0.8, -0.8,  0.5,  1.0])
    pyrat.run()
    tmodel1 = pyrat.spec.spectrum
    np.testing.assert_allclose(tmodel1, expected['tmodel'], rtol=1e-7)
    # np.testing.assert_allclose(tpyrat.atm.temp, tmodel, rtol=1e-7)


def test_transmission_tmodel(tmp_path):
    # Include tmodel and tpars in input config file:
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'tmodel':'tcea', 'tpars':'-1.5 -0.8 -0.8 0.5 1.0'})
    pyrat = pb.run(cfg)
    tmodel2 = pyrat.spec.spectrum
    np.testing.assert_allclose(tmodel2, expected['tmodel'], rtol=1e-7)


def test_transmission_vert_none_model(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'molmodel':'vert', 'molfree':'H2O', 'bulk':'H2 He'})
    pyrat = pb.run(cfg)
    vmodel0 = pyrat.spec.spectrum
    np.testing.assert_allclose(vmodel0, expected['all'], rtol=1e-7)
    pyrat.atm.molpars = [-5]
    pyrat.run()
    vmodel1 = pyrat.spec.spectrum
    np.testing.assert_allclose(vmodel1, expected['vert'], rtol=1e-7)
    np.testing.assert_equal(pyrat.atm.q[:,3], 1e-5)


def test_transmission_vert_model(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'molmodel':'vert', 'molfree':'H2O', 'molpars':'-5',
               'bulk':'H2 He'})
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['vert'], rtol=1e-7)


def test_transmission_scale_model(tmp_path):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'molmodel':'scale', 'molfree':'H2O', 'molpars':'-1',
               'bulk':'H2 He'})
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['scale'],
        rtol=1e-7)
    np.testing.assert_equal(pyrat.atm.q[:,3], 0.1*pyrat.atm.qbase[:,3])


def test_fit(tmp_path):
    # Without evaulating params:
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'tmodel':'tcea', 'cpars':'2.0',
               'molmodel':'vert', 'molfree':'H2O', 'bulk':'H2 He',
               'retflag':'temp mol ray cloud',
               'params':'-1.5 -0.8 -0.8 0.5 1.0 -4.0 0.0 -4.0 2.0'})
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=1e-7)
    # Eval default params:
    model1 = pyrat.eval(pyrat.ret.params, retmodel=True)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['fit1'], rtol=1e-7)
    np.testing.assert_equal(model1[0], pyrat.spec.spectrum)
    assert model1[1] is None
    # Cloud deck:
    params = [-1.5, -0.8, -0.8,  0.5,  1.0, -4.0,  0.0, -4.0,  -3.0]
    model2 = pyrat.eval(params, retmodel=True)
    rmin = np.amin(np.sqrt(pyrat.spec.spectrum)) * pyrat.phy.rstar
    assert np.amax(pyrat.atm.press[pyrat.atm.radius > rmin]) == 1e-3*pc.bar
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['fit2'], rtol=1e-7)
    # Depleted H2O:
    params = [-1.5, -0.8, -0.8,  0.5,  1.0, -8.0,  0.0, -4.0,  2.0]
    model3 = pyrat.eval(params, retmodel=True)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['fit3'], rtol=1e-7)


def test_fit_filters():
    pyrat = pb.run(ROOT+'tests/spectrum_transmission_filters_test.cfg')
    model4 = pyrat.eval(pyrat.ret.params, retmodel=True)
    np.testing.assert_allclose(model4[0], expected['fit4'],      rtol=1e-7)
    np.testing.assert_allclose(model4[1], expected['bandflux4'], rtol=1e-7)


def test_multiple_opacities(tmp_path):
    # Generate TLI files:
    cfg = make_config(tmp_path, ROOT+'tests/tli_multiple_opacity.cfg')
    pb.run(cfg)
    cfg = make_config(tmp_path, ROOT+'tests/tli_multiple_opacity.cfg',
        reset={'dblist':'02_hit12.par',
               'tlifile':'HITRAN_CO2_1.5-1.6um_test.tli'})
    pb.run(cfg)
    cfg = make_config(tmp_path, ROOT+'tests/tli_multiple_opacity.cfg',
        reset={'dblist':'06_hit12.par',
               'tlifile':'HITRAN_CH4_1.5-1.6um_test.tli'})
    pb.run(cfg)

    # Generate opacity files:
    cfg = make_config(tmp_path, ROOT+'tests/opacity_multiple.cfg')
    pb.run(cfg)
    cfg = make_config(tmp_path, ROOT+'tests/opacity_multiple.cfg',
        reset={'tlifile':'HITRAN_CO2_1.5-1.6um_test.tli',
               'extfile':'exttable_CO2_300-3000K_1.5-1.6um.dat'})
    pb.run(cfg)
    cfg = make_config(tmp_path, ROOT+'tests/opacity_multiple.cfg',
        reset={'tlifile':'HITRAN_CH4_1.5-1.6um_test.tli',
               'extfile':'exttable_CH4_300-3000K_1.5-1.6um.dat'})
    pb.run(cfg)
    cfg = make_config(tmp_path, ROOT+'tests/opacity_multiple.cfg',
        reset={'tlifile':'HITRAN_CO2_1.5-1.6um_test.tli'
                   '\n    HITRAN_CH4_1.5-1.6um_test.tli',
               'extfile':'exttable_CO2-CH4_300-3000K_1.5-1.6um.dat'})
    pb.run(cfg)

    # Compute spectra from opacities:
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['tlifile', 'clouds'],
        reset={'extfile':'exttable_H2O_300-3000K_1.5-1.6um.dat'
                     '\n  exttable_CO2_300-3000K_1.5-1.6um.dat'
                     '\n  exttable_CH4_300-3000K_1.5-1.6um.dat',
               'wllow':'1.5', 'wlhigh':'1.6'})
    pyrat1 = pb.run(cfg)
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['tlifile', 'clouds'],
        reset={'extfile':'exttable_H2O_300-3000K_1.5-1.6um.dat'
                     '\n  exttable_CO2-CH4_300-3000K_1.5-1.6um.dat',
               'wllow':'1.5', 'wlhigh':'1.6'})
    pyrat2 = pb.run(cfg)
    np.testing.assert_allclose(pyrat1.spec.spectrum,
                               pyrat2.spec.spectrum, rtol=1e-5)


# These are extra bits for testing the tests before testing:
def plot_fit():
    import matplotlib.pyplot as plt
    plt.figure(0)
    plt.clf()
    plt.plot(1e4/pyrat.spec.wn, model0, 'blue')
    plt.plot(1e4/pyrat.spec.wn, model1[0], 'orange',    alpha=0.6)
    plt.plot(1e4/pyrat.spec.wn, model2[0], 'limegreen', alpha=0.6)
    plt.plot(1e4/pyrat.spec.wn, model3[0], 'red',       alpha=0.6)
    plt.plot(1e4/pyrat.spec.wn, model4[0], 'black',     alpha=0.6)


def spectrum_fm():
    from scipy.ndimage.filters import gaussian_filter1d as gaussf
    pyrat = pb.run(ROOT+'tests/spectrum_transmission_filters_test.cfg')
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
    pyrat = pb.run(ROOT+'tests/spectrum_transmission_exomol.cfg')

    plt.figure(1)
    plt.clf()
    plt.plot(1e4/hcn[0], hcn[1], 'b')
    plt.plot(1e4/pyrat.spec.wn, pyrat.spec.spectrum, 'orange', alpha=0.5)

