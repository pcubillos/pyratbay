import os
import sys
import pytest

import numpy as np

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay as pb
import pyratbay.io as io
import pyratbay.constants  as pc
import pyratbay.atmosphere as pa

os.chdir(ROOT+'tests')


# Expected spectra:
keys = ['lec', 'cia', 'alkali', 'deck', 'tli', 'all', 'etable',
        'tmodel', 'vert', 'scale', 'fit1', 'fit2', 'fit3', 'fit4',
        'bandflux4']
expected = {key:np.load("expected_spectrum_transmission_{:s}_test.npz".
                        format(key))['arr_0']
            for key in keys}
#np.savez('expected_spectrum_transmission_fit_test.npz', model1[0])


# TBD: Check output files
def test_transmission_clear():
    # No opacity whatsoever:
    clear = pb.pbay.run(ROOT+'tests/spectrum_transmission_clear_test.cfg')
    depth_bottom = (clear.atm.radius[-1] / clear.phy.rstar)**2
    np.testing.assert_allclose(clear.spec.spectrum, depth_bottom, rtol=1e-7)


def test_transmission_lecavelier():
    ray = pb.pbay.run(ROOT+'tests/spectrum_transmission_lecavelier_test.cfg')
    np.testing.assert_allclose(ray.spec.spectrum, expected['lec'], rtol=1e-7)


def test_transmission_CIA():
    cia = pb.pbay.run(ROOT+'tests/spectrum_transmission_CIA_test.cfg')
    np.testing.assert_allclose(cia.spec.spectrum, expected['cia'], rtol=1e-7)


def test_transmission_alkali():
    alkali = pb.pbay.run(ROOT+'tests/spectrum_transmission_alkali_test.cfg')
    np.testing.assert_allclose(alkali.spec.spectrum, expected['alkali'],
        rtol=1e-7)


def test_transmission_deck():
    deck = pb.pbay.run(ROOT+'tests/spectrum_transmission_deck_test.cfg')
    np.testing.assert_allclose(deck.spec.spectrum, expected['deck'], rtol=1e-7)


def test_transmission_tli():
    tli = pb.pbay.run(ROOT+'tests/spectrum_transmission_tli_test.cfg')
    np.testing.assert_allclose(tli.spec.spectrum, expected['tli'], rtol=1e-7)


def test_transmission():
    # Transmission including all types of opacity:
    pyrat = pb.pbay.run(ROOT+'tests/spectrum_transmission_test.cfg')
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=1e-7)


def test_transmission_etable():
    # LBL from extinction table:
    epyrat = pb.pbay.run(ROOT+'tests/spectrum_transmission_etable_test.cfg')
    np.testing.assert_allclose(epyrat.spec.spectrum, expected['etable'],
    rtol=1e-7)


def plot_transmission():
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
def test_transmission_tmodel_none():
    # include tmodel, but tparams is None
    tpyrat = pb.pbay.run(ROOT+'tests/spectrum_transmission_tmodel_none_test.cfg')
    tmodel0 = tpyrat.spec.spectrum
    np.testing.assert_allclose(tmodel0, expected['all'], rtol=1e-7)
    # Now, re-run with user-input tparams:
    tpyrat.atm.tpars = np.array([-1.5, -0.8, -0.8,  0.5,  1.0])
    tpyrat.run()
    tmodel1 = tpyrat.spec.spectrum
    np.testing.assert_allclose(tmodel1, expected['tmodel'], rtol=1e-7)
    # np.testing.assert_allclose(tpyrat.atm.temp, tmodel, rtol=1e-7)


def test_transmission_tmodel():
    # Include tmodel and tpars  in input config file:
    tpyrat = pb.pbay.run(ROOT+'tests/spectrum_transmission_tmodel_test.cfg')
    tmodel2 = tpyrat.spec.spectrum
    np.testing.assert_allclose(tmodel2, expected['tmodel'], rtol=1e-7)


def test_transmission_vert_none_model():
    vpyrat = pb.pbay.run(ROOT+'tests/spectrum_transmission_vert_none_test.cfg')
    vmodel0 = vpyrat.spec.spectrum
    np.testing.assert_allclose(vmodel0, expected['all'], rtol=1e-7)
    vpyrat.atm.molpars = [-5]
    vpyrat.run()
    vmodel1 = vpyrat.spec.spectrum
    np.testing.assert_allclose(vmodel1, expected['vert'], rtol=1e-7)
    np.testing.assert_equal(vpyrat.atm.q[:,3], 1e-5)


def test_transmission_vert_model():
    vpyrat = pb.pbay.run(ROOT+'tests/spectrum_transmission_vert_test.cfg')
    vmodel2 = vpyrat.spec.spectrum
    np.testing.assert_allclose(vmodel2, expected['vert'], rtol=1e-7)


def test_transmission_scale_model():
    spyrat = pb.pbay.run(ROOT+'tests/spectrum_transmission_scale_test.cfg')
    smodel1 = spyrat.spec.spectrum
    np.testing.assert_allclose(smodel1, expected['scale'], rtol=1e-7)
    np.testing.assert_equal(spyrat.atm.q[:,3], 0.1*spyrat.atm.qbase[:,3])


def test_fit():
    # Without evaulating params:
    pyrat = pb.pbay.run(ROOT+'tests/spectrum_transmission_fit_test.cfg')
    model0 = np.copy(pyrat.spec.spectrum)
    np.testing.assert_allclose(model0, expected['all'])
    # Eval default params:
    model1 = pyrat.eval(pyrat.ret.params, retmodel=True)
    np.testing.assert_allclose(model1[0], expected['fit1'], rtol=1e-7)
    assert model1[1] is None
    # Cloud deck:
    params = [-1.5, -0.8, -0.8,  0.5,  1.0, -4.0,  0.0, -4.0,  -3.0]
    model2 = pyrat.eval(params, retmodel=True)
    rmin = np.amin(np.sqrt(pyrat.spec.spectrum)) * pyrat.phy.rstar
    assert np.amax(pyrat.atm.press[pyrat.atm.radius > rmin]) == 1e-3*pc.bar
    np.testing.assert_allclose(model2[0], expected['fit2'], rtol=1e-7)
    assert model2[1] is None
    # Depleted H2O:
    params = [-1.5, -0.8, -0.8,  0.5,  1.0, -8.0,  0.0, -4.0,  2.0]
    model3 = pyrat.eval(params, retmodel=True)
    np.testing.assert_allclose(model3[0], expected['fit3'], rtol=1e-7)
    assert model3[1] is None


def test_fit_filters():
    pyrat = pb.pbay.run(ROOT+'tests/spectrum_transmission_filters_test.cfg')
    model4 = pyrat.eval(pyrat.ret.params, retmodel=True)
    np.testing.assert_allclose(model4[0], expected['fit4'],      rtol=1e-7)
    np.testing.assert_allclose(model4[1], expected['bandflux4'], rtol=1e-7)


def plot_fit():
    plt.figure(0)
    plt.clf()
    plt.plot(1e4/pyrat.spec.wn, model0, 'blue')
    plt.plot(1e4/pyrat.spec.wn, model1[0], 'orange',    alpha=0.6)
    plt.plot(1e4/pyrat.spec.wn, model2[0], 'limegreen', alpha=0.6)
    plt.plot(1e4/pyrat.spec.wn, model3[0], 'red',       alpha=0.6)
    plt.plot(1e4/pyrat.spec.wn, model4[0], 'black',     alpha=0.6)


def spectrum_fm():
    from scipy.ndimage.filters import gaussian_filter1d as gaussf
    import pyratbay.wine as pw
    pyrat = pb.pbay.run(ROOT+'tests/spectrum_transmission_filters_test.cfg')
    params = [-1.5, -0.8, 0.0,  1.0,  1.0,  71500.0, -3.4, 2.0]
    model = pyrat.eval(params, retmodel=True)
    #bandflux = pw.bandintegrate(pyrat=pyrat)
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

