# Copyright (c) 2021-2022 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os
import pytest
import re

import numpy as np

from conftest import make_config

import pyratbay as pb
import pyratbay.constants as pc
import pyratbay.io as io
import pyratbay.spectrum as ps
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


# Expected spectra:
keys = [
    'lec', 'cia', 'alkali', 'deck', 'tli', 'all', 'quadrature', 'etable',
    'resolution', 'two_stream',
    'tmodel', 'vert', 'scale',
    #'fit1', 'fit2', 'fit3', 'fit4',
    #'bandflux4',
    ]
expected = {
    key:np.load(f"{ROOT}tests/expected/"
                f"expected_spectrum_emission_{key}_test.npz")['arr_0']
    for key in keys}
#np.savez('expected/expected_spectrum_emission_fit_test.npz', model1[0])

INPUTS = f'{ROOT}tests/inputs/'
OUTPUTS = f'{ROOT}tests/outputs/'

# Relative tolerance of less than 0.01% difference:
rtol = 0.01 / 100.0


def test_emission_clear(tmp_path):
    # No opacity whatsoever:
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'alkali', 'clouds'],
    )
    pyrat = pb.run(cfg)
    spectrum = ps.bbflux(pyrat.spec.wn, pyrat.atm.temp[-2])
    # TBD: Should be last layer, check ideep calculation
    np.testing.assert_allclose(pyrat.spec.spectrum, spectrum, rtol=rtol)


def test_emission_lecavelier(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['tlifile', 'csfile', 'alkali', 'clouds'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['lec'], rtol=rtol)


def test_emission_CIA(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['tlifile', 'rayleigh', 'alkali', 'clouds'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['cia'], rtol=rtol)


def test_emission_alkali(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'clouds'],
        reset={'wllow':'0.45 um', 'wlhigh':'1.0 um'},
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['alkali'], rtol=rtol,
    )


def test_emission_deck(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'alkali'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['deck'], rtol=rtol)

    tsurf_cloud = pyrat.opacity.models[0].tsurf
    spectrum = ps.bbflux(pyrat.spec.wn, tsurf_cloud)
    np.testing.assert_allclose(pyrat.spec.spectrum, spectrum, rtol=rtol)


def test_emission_tli(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['csfile', 'rayleigh', 'clouds', 'alkali'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['tli'], rtol=rtol)


def test_emission_all(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['clouds'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=rtol)


def test_emission_quadrature(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        reset={'quadrature': '5'},
        remove=['clouds'],
    )
    pyrat = pb.run(cfg)
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['quadrature'], rtol=rtol)


def test_emission_two_stream(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        reset={'rt_path': 'emission_two_stream'},
        remove=['clouds'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['two_stream'], rtol=rtol)


def test_emission_resolution(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        reset={'resolution':'5000.0'},
        remove=['clouds'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['resolution'], rtol=rtol)


def test_emission_etable(tmp_path):
    # LBL from extinction table:
    reset = {
        'extfile': f'{OUTPUTS}exttable_test_300-3000K_1.1-1.7um.npz',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['tlifile', 'clouds'],
        reset=reset,
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['etable'], rtol=rtol,
    )


def test_emission_dilution(tmp_path):
    reset = {
        'f_dilution': '0.75',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        reset=reset,
        remove=['clouds'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum,
        0.75*expected['all'],
        rtol=rtol,
    )


# Optical-depth integration is home made, which depends on whether there is
# an odd or even number of layers. Thus, the need for this test.
def test_emission_odd_even(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        reset={'rpars':'1.0 -4.0'},
        remove=['tlifile', 'alkali', 'clouds'],
    )
    pyrat = pb.run(cfg)
    odd_spectrum = pyrat.spec.spectrum

    reset = {
        'atmfile': f'{INPUTS}atmosphere_uniform_even_layers.atm',
        'rpars': '1.0 -4.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['tlifile', 'alkali', 'clouds'],
        reset=reset,
    )
    pyrat = pb.run(cfg)
    even_spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(odd_spectrum, even_spectrum, rtol=rtol)


@pytest.mark.skip(reason="")
def test_emission_input_radius(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        reset={'atmfile':f'{INPUTS}atmosphere_uniform_radius.atm'},
        remove=['radmodel'],
    )
    pyrat = pb.run(cfg)
    atm = io.read_atm('inputs/atmosphere_uniform_radius.atm')
    np.testing.assert_allclose(pyrat.atm.radius, atm[5]*pc.km, rtol=rtol)


@pytest.mark.skip(reason="")
def test_emission_input_radius_overwrite(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        reset={'atmfile': f'{INPUTS}/atmosphere_uniform_radius.atm'},
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
def test_emission_tmodel(tmp_path):
    # Include tmodel and tpars in input config file:
    reset = {
        'tmodel':'guillot',
        'tpars':'-4.67 -0.8 -0.8 0.5 1486.0 100.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['clouds', 'cpars'],
        reset=reset,
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['tmodel'], rtol=rtol,
    )


def test_emission_tmodel_no_tpars(tmp_path):
    # include tmodel, but tpars is None
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'tmodel': 'guillot'},
    )
    error = re.escape('Not all temperature parameters were defined (tpars)')
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_emission_vert_model(tmp_path):
    mol_vars = "log_H2O -5"
    reset = {
        'vmr_vars': 'log_H2O -5',
        'bulk': 'H2 He',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['clouds', 'cpars'],
        reset=reset,
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['vert'], rtol=rtol)


def test_emission_vert_model_no_molpars(tmp_path):
    reset = {
        'vmr_vars': 'log_H2O',
        'bulk': 'H2 He',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['clouds', 'cpars'],
        reset=reset,
    )
    error = re.escape('Not all vmr parameter values were defined (vmr_vars)')
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_emission_scale_model(tmp_path):
    reset={
        'vmr_vars': 'scale_H2O -1',
        'bulk': 'H2 He',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['clouds', 'cpars'],
        reset=reset,
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['scale'], rtol=rtol,
    )
    i_H2O = list(pyrat.atm.species).index('H2O')
    np.testing.assert_allclose(
        pyrat.atm.vmr[:,i_H2O],
        0.1*pyrat.atm.base_vmr[:,i_H2O],
        rtol=rtol,
    )


@pytest.mark.skip(reason="")
def test_emission_fit(tmp_path):
    # Without evaulating params:
    reset = {
        'tmodel': 'guillot',
        'cpars': '2.0',
        'vmr_vars': 'log_H2O',
        'bulk': 'H2 He',
        'retflag': 'temp mol ray cloud',
        'params': '-4.67 -0.8 -0.8 0.5 1486.0 100.0 -4.0 0.0 -4.0 2.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=rtol)
    # Eval default params:
    model1 = pyrat.eval(pyrat.ret.params, retmodel=True)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['fit1'], rtol=rtol)
    np.testing.assert_allclose(model1[0], pyrat.spec.spectrum, rtol=rtol)
    assert model1[1] is None
    # Cloud deck:
    params = [-4.67, -0.8, -0.8, 0.5, 1486.0, 100.0, -4.0, 0.0, -4.0, -3.0]
    model2 = pyrat.eval(params, retmodel=True)
    rmin = np.amin(np.sqrt(pyrat.spec.spectrum)) * pyrat.phy.rstar
    rexpected = pyrat.cloud.models[0].rsurf
    np.testing.assert_allclose(rmin, rexpected, rtol=rtol)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['fit2'], rtol=rtol)
    # Check pyrat.ret.params has been updated:
    np.testing.assert_allclose(pyrat.ret.params, params, rtol=rtol)
    # Depleted H2O:
    params = [-4.67, -0.8, -0.8, 0.5, 1486.0, 100.0, -8.0, 0.0, -4.0, 2.0]
    model3 = pyrat.eval(params, retmodel=True)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['fit3'], rtol=rtol)


def test_emission_band_integrate_no_data():
    pyrat = pb.run(ROOT+'tests/configs/spectrum_emission_filters_test.cfg')
    ev_spectrum, ev_bandflux = pyrat.eval(pyrat.ret.params, retmodel=True)
    spectrum = np.copy(pyrat.spec.spectrum)
    bandflux = pyrat.band_integrate()

    expected_bandflux = [
        2.0600785289e-04, 2.2919655505e-04, 2.4065559114e-04, 2.8029860209e-04,
        3.0155143043e-04, 3.1963485203e-04, 3.3379717207e-04, 3.3539204122e-04,
        3.2066224764e-04, 2.4138103230e-04, 2.6113180837e-04, 2.6163118428e-04,
        2.6030023197e-04, 3.0939984882e-04, 3.3163315401e-04, 3.9102530827e-04,
        4.5804716173e-04, 4.9061722081e-04, 5.3090552207e-04, 5.5291946981e-04,
        5.6982334683e-04
    ]
    #print(' '.join([f'{flux:.10e},' for flux in bandflux]))
    np.testing.assert_allclose(pyrat.spec.spectrum, spectrum)
    np.testing.assert_allclose(bandflux, expected_bandflux)
    np.testing.assert_allclose(ev_bandflux, expected_bandflux)


def test_emission_plot_spectrum_band_integrate():
    pyrat = pb.run(ROOT+'tests/configs/spectrum_emission_filters_test.cfg')
    ev_spectrum, ev_bandflux = pyrat.eval(pyrat.ret.params, retmodel=True)
    filename = f'{pc.ROOT}tests/outputs/emission_spectrum_filters.png'
    ax = pyrat.plot_spectrum(filename=filename)


@pytest.mark.skip(reason="")
def test_multiple_opacities(tmp_path):
    # Generate TLI files:
    cfg = make_config(tmp_path, f'{ROOT}tests/configs/tli_multiple_opacity.cfg')
    pyrat = pb.run(cfg)
    cfg = make_config(
        tmp_path,
        f'{ROOT}tests/configs/tli_multiple_opacity.cfg',
        reset={'dblist':f'{INPUTS}02_hit12.par',
               'tlifile':f'{OUTPUTS}HITRAN_CO2_1.5-1.6um_test.tli'})
    pyrat = pb.run(cfg)
    cfg = make_config(
        tmp_path,
        f'{ROOT}tests/configs/tli_multiple_opacity.cfg',
        reset={'dblist':f'{INPUTS}06_hit12.par',
               'tlifile':f'{OUTPUTS}HITRAN_CH4_1.5-1.6um_test.tli'})
    pyrat = pb.run(cfg)

    # Generate opacity files:
    cfg = make_config(tmp_path, ROOT+'tests/configs/opacity_multiple.cfg')
    pyrat = pb.run(cfg)
    assert pyrat is not None
    cfg = make_config(
        tmp_path,
        f'{ROOT}tests/configs/opacity_multiple.cfg',
        reset={'tlifile':f'{OUTPUTS}HITRAN_CO2_1.5-1.6um_test.tli',
               'extfile':f'{OUTPUTS}exttable_CO2_300-3000K_1.5-1.6um.npz'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    cfg = make_config(
        tmp_path,
        f'{ROOT}tests/configs/opacity_multiple.cfg',
        reset={'tlifile':f'{OUTPUTS}HITRAN_CH4_1.5-1.6um_test.tli',
               'extfile':f'{OUTPUTS}exttable_CH4_300-3000K_1.5-1.6um.npz'})
    pyrat = pb.run(cfg)
    assert pyrat is not None
    cfg = make_config(
        tmp_path,
        f'{ROOT}tests/configs/opacity_multiple.cfg',
        reset={'tlifile':f'{OUTPUTS}HITRAN_CO2_1.5-1.6um_test.tli'
                   f'\n    {OUTPUTS}HITRAN_CH4_1.5-1.6um_test.tli',
               'extfile':f'{OUTPUTS}exttable_CO2-CH4_300-3000K_1.5-1.6um.npz'})
    pyrat = pb.run(cfg)
    assert pyrat is not None

    # Compute spectra from opacities:
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['tlifile', 'clouds'],
        reset={'extfile':f'{OUTPUTS}exttable_H2O_300-3000K_1.5-1.6um.npz'
                     f'\n  {OUTPUTS}exttable_CO2_300-3000K_1.5-1.6um.npz'
                     f'\n  {OUTPUTS}exttable_CH4_300-3000K_1.5-1.6um.npz',
               'wllow':'1.5 um', 'wlhigh':'1.6 um'})
    pyrat1 = pb.run(cfg)
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['tlifile', 'clouds'],
        reset={'extfile':f'{OUTPUTS}exttable_H2O_300-3000K_1.5-1.6um.npz'
                     f'\n  {OUTPUTS}exttable_CO2-CH4_300-3000K_1.5-1.6um.npz',
               'wllow':'1.5 um', 'wlhigh':'1.6 um'})
    pyrat2 = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat1.spec.spectrum, pyrat2.spec.spectrum, rtol=rtol)


@pytest.mark.skip(reason="")
@pytest.mark.parametrize('wllow,wlhigh',
    [('1.1 um', '1.6 um'),
     ('1.2 um', '1.7 um'),
     ('1.2 um', '1.6 um')])
def test_opacity_reset_wn(tmp_path, wllow, wlhigh):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_emission_test.cfg',
        remove=['tlifile', 'clouds'],
        reset={'extfile':f'{OUTPUTS}exttable_test_300-3000K_1.1-1.7um.npz',
               'wllow':wllow, 'wlhigh':wlhigh},
    )
    pyrat = pb.run(cfg)
    wn = np.arange(1/1.7e-4, 1/1.1e-4, 1.0)
    wn_range = (wn>= pyrat.spec.wnlow) & (wn <=pyrat.spec.wnhigh)
    etab = expected['etable'][wn_range]
    np.testing.assert_allclose(pyrat.spec.spectrum, etab, rtol=rtol)

