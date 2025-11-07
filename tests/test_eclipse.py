# Copyright (c) 2021-2025 Cubillos & Blecic
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os
import pytest
import re

import numpy as np

from conftest import make_config

import pyratbay as pb
import pyratbay.constants as pc
import pyratbay.spectrum as ps
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


# Expected spectra:
keys = [
    'sampled_cs', 'tli', 'lec', 'cia', 'alkali', 'deck', 'all',
    'patchy', 'patchy_clear', 'patchy_cloudy',
    'quadrature', 'resolution', 'two_stream',
    'tmodel', 'vert', 'scale',
]
expected = {
    key:np.load(
        f"{ROOT}tests/expected/expected_spectrum_eclipse_{key}_test.npz"
    )['arr_0']
    for key in keys
}

INPUTS = f'{ROOT}tests/inputs/'
OUTPUTS = f'{ROOT}tests/outputs/'

# Relative tolerance of less than 0.01% difference:
rtol = 0.01 / 100.0


def test_eclipse_clear(tmp_path):
    # No opacity whatsoever:
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        remove=['sampled_cross_sec', 'continuum_cross_sec', 'alkali', 'clouds'],
    )
    pyrat = pb.run(cfg)
    spectrum = (
        ps.bbflux(pyrat.spec.wn, pyrat.atm.temp[-1]) / pyrat.spec.starflux
        * (pyrat.atm.rplanet/pyrat.atm.rstar)**2
    )
    # TBD: Should be last layer, check ideep calculation
    np.testing.assert_allclose(pyrat.spec.spectrum, spectrum, rtol=rtol)


def test_eclipse_sampled_cs(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        remove=['continuum_cross_sec', 'clouds', 'alkali'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['sampled_cs'], rtol=rtol)


def test_eclipse_tli(tmp_path):
    reset = {
        'tlifile': f'{OUTPUTS}HITRAN_H2O_1.1-1.7um_test.tli',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        remove=['sampled_cross_sec', 'continuum_cross_sec', 'clouds', 'alkali'],
        reset=reset,
    )
    pyrat = pb.run(cfg)
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['tli'], rtol=rtol)


def test_eclipse_lecavelier(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        remove=['sampled_cross_sec', 'continuum_cross_sec', 'alkali'],
        reset={'clouds': 'lecavelier 2.0 -4.0'}
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['lec'], rtol=rtol)


def test_eclipse_CIA(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        remove=['sampled_cross_sec', 'alkali', 'clouds'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['cia'], rtol=rtol)


def test_eclipse_alkali(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        remove=['sampled_cross_sec', 'continuum_cross_sec', 'clouds'],
        reset={'wl_low':'0.45 um', 'wl_high':'1.0 um'},
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['alkali'], rtol=rtol,
    )


def test_eclipse_deck(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        remove=['sampled_cross_sec', 'continuum_cross_sec', 'alkali'],
        reset={'clouds': 'deck -1.0'},
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['deck'], rtol=rtol)

    tsurf_cloud = pyrat.opacity.models[0].tsurf
    spectrum = ps.bbflux(pyrat.spec.wn, tsurf_cloud)
    np.testing.assert_allclose(pyrat.spec.fplanet, spectrum, rtol=rtol)


def test_eclipse_all(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=rtol)


def test_eclipse_patchy(tmp_path):
    reset = {
        'fpatchy': '0.5',
        'clouds': 'deck -3.0\nlecavelier 10.0 -15.0',
    }
    cfg_file = ROOT+'tests/configs/spectrum_eclipse_test.cfg'
    cfg = make_config(tmp_path, cfg_file, reset=reset)
    pyrat = pb.run(cfg)

    spectrum = pyrat.spec.spectrum
    clear = pyrat.spec.clear
    cloudy = pyrat.spec.cloudy
    np.testing.assert_allclose(spectrum, expected['patchy'], rtol=rtol)
    np.testing.assert_allclose(clear, expected['patchy_clear'], rtol=rtol)
    np.testing.assert_allclose(cloudy, expected['patchy_cloudy'], rtol=rtol)

    pyrat.opacity.fpatchy = 0.0
    pyrat.run()
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['patchy_clear'], rtol=rtol)

    pyrat.opacity.fpatchy = 1.0
    pyrat.run()
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['patchy_cloudy'], rtol=rtol)


def test_eclipse_quadrature(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        reset={'quadrature': '5'},
    )
    pyrat = pb.run(cfg)
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['quadrature'], rtol=rtol)


def test_eclipse_two_stream(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        reset={'rt_path': 'eclipse_two_stream'},
    )
    pyrat = pb.run(cfg)
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['two_stream'], rtol=rtol)


def test_eclipse_resolution(tmp_path):
    reset = {
        'tlifile': f'{OUTPUTS}HITRAN_H2O_1.1-1.7um_test.tli',
        'resolution': '5000.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        remove=['sampled_cross_sec'],
        reset=reset,
    )
    pyrat = pb.run(cfg)
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['resolution'], rtol=rtol)


def test_eclipse_dilution(tmp_path):
    reset = {
        'f_dilution': '0.75',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum,
        0.75*expected['all'],
        rtol=rtol,
    )


# Optical-depth integration is home made, which depends on whether there is
# an odd or even number of layers. Thus, the need for this test.
def test_eclipse_odd_even(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        remove=['sampled_cross_sec', 'alkali'],
        reset={'clouds':'lecavelier  1.0 -4.0'},
    )
    pyrat = pb.run(cfg)
    odd_spectrum = pyrat.spec.spectrum

    reset = {
        'atmfile': f'{INPUTS}atmosphere_uniform_even_layers.atm',
        'clouds': 'lecavelier 1.0 -4.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        remove=['sampled_cross_sec', 'alkali'],
        reset=reset,
    )
    pyrat = pb.run(cfg)
    even_spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(odd_spectrum, even_spectrum, rtol=rtol)


def test_eclipse_tmodel(tmp_path):
    # Include tmodel and tpars in input config file:
    reset = {
        'tmodel':'guillot',
        'tpars':'-4.67 -0.8 -0.8 0.5 1486.0 100.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['tmodel'], rtol=rtol,
    )


def test_eclipse_tmodel_no_tpars(tmp_path):
    # include tmodel, but tpars is None
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        reset={'tmodel': 'guillot'},
    )
    error = re.escape('Not all temperature parameters were defined (tpars)')
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_eclipse_vert_model(tmp_path):
    reset = {
        'vmr_vars': 'log_H2O -5',
        'bulk': 'H2 He',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['vert'], rtol=rtol)


def test_eclipse_vert_model_no_molpars(tmp_path):
    reset = {
        'vmr_vars': 'log_H2O',
        'bulk': 'H2 He',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        reset=reset,
    )
    error = re.escape('Not all vmr parameter values were defined (vmr_vars)')
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_eclipse_scale_model(tmp_path):
    reset={
        'vmr_vars': 'scale_H2O -1',
        'bulk': 'H2 He',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_eclipse_test.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['scale'], rtol=rtol)

    i_H2O = list(pyrat.atm.species).index('H2O')
    np.testing.assert_allclose(
        pyrat.atm.vmr[:,i_H2O],
        0.1*pyrat.atm.base_vmr[:,i_H2O],
        rtol=rtol,
    )


def test_eclipse_band_integrate_no_data():
    pyrat = pb.run(ROOT+'tests/configs/spectrum_eclipse_filters_test.cfg')
    ev_spectrum, ev_bandflux = pyrat.eval(pyrat.ret.params, retmodel=True)
    spectrum = np.copy(pyrat.spec.spectrum)
    bandflux = pyrat.band_integrate()

    expected_bandflux = [
        2.0600788319e-04, 2.2919658941e-04, 2.4065562405e-04, 2.8029874812e-04,
        3.0155157632e-04, 3.1963498262e-04, 3.3379724117e-04, 3.3539207179e-04,
        3.2066247922e-04, 2.4138123734e-04, 2.6113192476e-04, 2.6163135874e-04,
        2.6030040874e-04, 3.0940015836e-04, 3.3163344643e-04, 3.9102537977e-04,
        4.5804719858e-04, 4.9061725664e-04, 5.3090560619e-04, 5.5291950762e-04,
        5.6982335860e-04,
    ]
    #print(' '.join([f'{flux:.10e},' for flux in bandflux]))
    np.testing.assert_allclose(pyrat.spec.spectrum, spectrum)
    np.testing.assert_allclose(bandflux, expected_bandflux)
    np.testing.assert_allclose(ev_bandflux, expected_bandflux)


def test_eclipse_plot_spectrum_band_integrate():
    pyrat = pb.run(ROOT+'tests/configs/spectrum_eclipse_filters_test.cfg')
    ev_spectrum, ev_bandflux = pyrat.eval(pyrat.ret.params, retmodel=True)
    filename = f'{pc.ROOT}tests/outputs/eclipse_spectrum_filters.png'
    ax = pyrat.plot_spectrum(filename=filename)

