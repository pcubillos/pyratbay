# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import pytest
import re

import numpy as np
import pyratbay as pb
import pyratbay.constants as pc
import pyratbay.io as io
from pyratbay.constants import ROOT

from conftest import make_config
os.chdir(ROOT+'tests')


# Expected spectra:
keys = [
    'lec', 'cia', 'alkali', 'deck', 'tli',
    'patchy', 'patchy_clear', 'patchy_cloudy', 'h_ion', 'all', 'etable',
    'tmodel', 'vert', 'scale', 'fit1', 'fit2', 'fit3', 'fit4',
    'bandflux4', 'resolution', 'wl_step',
    'skip_ls', 'skip_lbl', 'skip_cia', 'skip_H2_H2_cia',
    'skip_alkali', 'skip_sodium', 'skip_rayleigh', 'skip_dalgarno',
    'skip_cloud', 'skip_deck', 'all_ls'
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
        remove=['tlifile', 'csfile', 'rayleigh', 'alkali', 'clouds'],
    )
    pyrat = pb.run(cfg)
    depth_bottom = (pyrat.atm.radius[-1] / pyrat.phy.rstar)**2
    np.testing.assert_allclose(pyrat.spec.spectrum, depth_bottom, rtol=rtol)


def test_transmission_lecavelier(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'alkali', 'clouds'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['lec'], rtol=rtol)


def test_transmission_CIA(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'rayleigh', 'alkali', 'clouds'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['cia'], rtol=rtol)


def test_transmission_alkali(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'clouds'],
        reset={'wllow':'0.45 um', 'wlhigh':'1.0 um'},
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['alkali'], rtol=rtol)


def test_transmission_deck(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'csfile', 'rayleigh', 'alkali'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['deck'], rtol=rtol)


def test_transmission_tli(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['csfile', 'rayleigh', 'clouds', 'alkali'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['tli'], rtol=rtol)


def test_transmission_h_ion(tmp_path):
    cfg = f'{ROOT}tests/configs/spectrum_transmission_h_ion.cfg'
    pyrat = pb.run(cfg)
    desired = expected['h_ion']
    np.testing.assert_allclose(pyrat.spec.spectrum, desired, rtol=rtol)


def test_transmission_all(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds'],
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['all'], rtol=rtol)


def test_transmission_patchy_from_fpatchy_param(tmp_path):
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

    pyrat.opacity.fpatchy = 0.0
    pyrat.run()
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['patchy_clear'], rtol=rtol)

    pyrat.opacity.fpatchy = 1.0
    pyrat.run()
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['patchy_cloudy'], rtol=rtol)


def test_transmission_f_patchy_from_retrieval_params(tmp_path):
    reset = {
        'rpars': '10.0 -15.0',
        'retrieval_params': 'f_patchy  0.5  0.0  1.0 0.1',
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

    pyrat.eval([0.0])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['patchy_clear'], rtol=rtol)

    pyrat.eval([1.0])
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
        pyrat.spec.spectrum, expected['resolution'], rtol=rtol,
    )


def test_transmission_wl_step(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'wlstep': '1e-4 um'},
        remove=['clouds'],
    )
    pyrat = pb.run(cfg)
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['wl_step'], rtol=rtol)


# Optical-depth integration is home made, which depends on whether there is
# an odd or even number of layers. Thus, the need for this test.
def test_transmission_odd_even(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'rpars':'1.0 -4.0'},
        remove=['tlifile', 'csfile', 'alkali', 'clouds'],
    )
    pyrat = pb.run(cfg)
    odd_spectrum = pyrat.spec.spectrum

    reset = {
        'atmfile': f'{INPUTS}atmosphere_uniform_even_layers.atm',
        'rpars': '1.0 -4.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset=reset,
        remove=['tlifile', 'csfile', 'alkali', 'clouds'],
    )
    pyrat = pb.run(cfg)
    even_spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(odd_spectrum, even_spectrum, rtol=rtol)


def test_transmission_etable(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['tlifile', 'clouds'],
        reset={'extfile':f'{OUTPUTS}exttable_test_300-3000K_1.1-1.7um.npz'},
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(
        pyrat.spec.spectrum, expected['etable'], rtol=rtol)


def test_transmission_interpolate_temp_etable(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_interp_temp.cfg',
    )
    pyrat = pb.run(cfg)

    expected_temp = np.array([ 500., 1000., 1500., 2000.])
    ntemp = len(expected_temp)
    assert pyrat.opacity.models[0].ntemp == ntemp
    np.testing.assert_allclose(
        pyrat.opacity.models[0].temp,
        expected_temp,
        rtol=rtol,
    )


def test_transmission_interpolate_press_etable(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_interp_press.cfg',
    )
    pyrat = pb.run(cfg)

    nlayers = 7
    expected_press = np.logspace(-5, 1, nlayers)
    assert pyrat.opacity.models[0].nlayers == nlayers
    np.testing.assert_allclose(
        pyrat.opacity.models[0].press,
        expected_press,
        rtol=rtol,
    )


def test_transmission_input_radius(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'atmfile': f'{INPUTS}atmosphere_uniform_radius.atm'},
        remove=['radmodel'],
    )
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
def test_transmission_tmodel(tmp_path):
    # Include tmodel and tpars in input config file:
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'tmodel':'guillot', 'tpars':'-4.67 -0.8 -0.8 0.5 1486.0 100.0'},
    )
    pyrat = pb.run(cfg)
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['tmodel'], rtol=rtol)


def test_transmission_tmodel_no_tpars(tmp_path):
    # include tmodel, but tpars is None
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'tmodel': 'guillot'},
    )
    error = re.escape('Not all temperature parameters were defined (tpars)')
    with pytest.raises(ValueError, match=error):
        pb.run(cfg)


def test_transmission_vert_model(tmp_path):
    reset = {
        'vmr_vars': 'log_H2O -5',
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


def test_transmission_vert_model_no_molpars(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds', 'cpars'],
        reset={'vmr_vars':'log_H2O', 'bulk':'H2 He'},
    )
    error = re.escape('Not all vmr parameter values were defined (vmr_vars)')
    with pytest.raises(ValueError, match=error):
        pb.run(cfg)


def test_transmission_scale_model(tmp_path):
    reset = {
        'vmr_vars': 'scale_H2O -1.0',
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
    i_H2O = list(pyrat.atm.species).index('H2O')
    np.testing.assert_allclose(
        pyrat.atm.vmr[:,i_H2O],
        0.1*pyrat.atm.base_vmr[:,i_H2O],
    )
    np.testing.assert_allclose(spectrum, expected['scale'], rtol=rtol)


def test_transmission_fit(tmp_path):
    # Without evaulating params:
    retrieval_params = """
        log_kappa'  -4.67
        log_gamma1  -0.8
        log_gamma2  -0.8
        alpha        0.5
        T_irr        1486.0
        T_int        100.0
        log_H2O      -4.0
        log_k_ray     0.0
        alpha_ray    -4.0
        log_p_cl      2.0
    """
    reset = {
        'tmodel': 'guillot',
        'cpars': '2.0',
        'vmr_vars': 'log_H2O',
        'bulk': 'H2 He',
        'retrieval_params': retrieval_params,
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['fit1'], rtol=rtol)
    # Cloud deck:
    params = [-4.67, -0.8, -0.8, 0.5, 1486.0, 100.0, -4.0, 0.0, -4.0, -3.0]
    pyrat.eval(params, retmodel=True)
    rmin = np.amin(np.sqrt(pyrat.spec.spectrum)) * pyrat.phy.rstar
    cloud_deck = pyrat.opacity.models[5]
    rexpected = cloud_deck.rsurf
    np.testing.assert_allclose(rmin, rexpected, rtol=rtol)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['fit2'], rtol=rtol)
    # Check pyrat.ret.params has been updated:
    np.testing.assert_allclose(pyrat.ret.params, params, rtol=rtol)
    # Depleted H2O:
    params = [-4.67, -0.8, -0.8, 0.5, 1486.0, 100.0, -8.0, 0.0, -4.0, 2.0]
    pyrat.eval(params, retmodel=True)
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['fit3'], rtol=rtol)


def test_transmission_fit_filters():
    pyrat = pb.run(ROOT+'tests/configs/spectrum_transmission_filters_test.cfg')
    model4 = pyrat.eval(pyrat.ret.params, retmodel=True)
    np.testing.assert_allclose(model4[0], expected['fit4'], rtol=rtol)
    np.testing.assert_allclose(model4[1], expected['bandflux4'], rtol=rtol)


def test_multiple_opacities(tmp_path):
    # Generate TLI files:
    cfg = f'{ROOT}tests/configs/tli_multiple_opacity_H2O.cfg'
    pyrat = pb.run(cfg)

    cfg = f'{ROOT}tests/configs/tli_multiple_opacity_CH4.cfg'
    pyrat = pb.run(cfg)

    cfg = f'{ROOT}tests/configs/tli_multiple_opacity_CO2.cfg'
    pyrat = pb.run(cfg)

    # Generate opacity files:
    opac_cfg = f'{ROOT}tests/configs/opacity_multiple.cfg'
    cfg = make_config(tmp_path, opac_cfg)
    pyrat = pb.run(cfg)
    assert pyrat is not None

    reset = {
        'tlifile': f'{OUTPUTS}HITRAN_CO2_1.5-1.51um_test.tli',
        'extfile': f'{OUTPUTS}exttable_CO2_300-3000K_1.5-1.51um.npz',
    }
    cfg = make_config(tmp_path, opac_cfg, reset=reset)
    pyrat = pb.run(cfg)
    assert pyrat is not None

    reset = {
        'tlifile': f'{OUTPUTS}HITRAN_CH4_1.5-1.51um_test.tli',
        'extfile': f'{OUTPUTS}exttable_CH4_300-3000K_1.5-1.51um.npz',
    }
    cfg = make_config(tmp_path, opac_cfg, reset=reset)
    pyrat = pb.run(cfg)
    assert pyrat is not None

    # Two species at a time:
    reset = {
        'tlifile':
            f'{OUTPUTS}HITRAN_CO2_1.5-1.51um_test.tli'
            f'\n    {OUTPUTS}HITRAN_CH4_1.5-1.51um_test.tli',
        'extfile': f'{OUTPUTS}exttable_CO2-CH4_300-3000K_1.5-1.51um.npz',
    }
    cfg = make_config(tmp_path, opac_cfg, reset=reset)
    pyrat = pb.run(cfg)
    assert pyrat is not None

    # Compute spectra from opacities:
    cfg = f'{ROOT}tests/configs/spectrum_transmission_multiple_opacities.cfg'
    pyrat = pb.run(cfg)
    spectrum1 = pyrat.spec.spectrum

    reset = {
        'extfile':
            f'{OUTPUTS}exttable_H2O_300-3000K_1.5-1.51um.npz'
            f'\n  {OUTPUTS}exttable_CO2-CH4_300-3000K_1.5-1.51um.npz',
    }
    cfg = make_config(tmp_path, cfg, reset=reset)
    pyrat = pb.run(cfg)
    spectrum2 = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum1, spectrum2, rtol=rtol)


@pytest.mark.parametrize(
    'wllow,wlhigh',
    [('1.1 um', '1.6 um'),
     ('1.2 um', '1.7 um'),
     ('1.2 um', '1.6 um'),]
)
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
    wn_range = (wn >= pyrat.spec.wnlow) & (wn <= pyrat.spec.wnhigh)
    etab = expected['etable'][wn_range]
    np.testing.assert_allclose(pyrat.spec.spectrum, etab, rtol=rtol)


# Now, the skips
@pytest.mark.parametrize(
    'model',
    ['line_sample', 'line sampling'],
)
def test_transmission_skip_sampling(tmp_path, model):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        remove=['clouds'],
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=[model])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['skip_ls'], rtol=rtol)


@pytest.mark.parametrize(
    'model',
    ['lbl', 'line by line'],
)
def test_transmission_skip_lbl(tmp_path, model):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds'],
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=[model])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['skip_lbl'], rtol=rtol)



def test_transmission_skip_alkali(tmp_path):
    reset = {
        'alkali': 'sodium_vdw potassium_vdw',
        'wllow': '0.4 um',
        'wlhigh': '1.2 um',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset=reset,
        remove=['clouds'],
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=['alkali'])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['skip_alkali'], rtol=rtol)


def test_transmission_skip_sodium_vdw(tmp_path):
    reset = {
        'alkali': 'sodium_vdw potassium_vdw',
        'wllow': '0.4 um',
        'wlhigh': '1.2 um',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset=reset,
        remove=['clouds'],
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=['sodium_vdw'])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['skip_sodium'], rtol=rtol)



def test_transmission_skip_cia(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        remove=['clouds'],
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=['cia'])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['skip_cia'], rtol=rtol)


def test_transmission_skip_H2_H2_cia(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        remove=['clouds'],
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=['CIA H2-H2'])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['skip_H2_H2_cia'], rtol=rtol)


def test_transmission_skip_rayleigh(tmp_path):
    reset = {
        'rayleigh': 'dalgarno_H2 dalgarno_He',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
        remove=['clouds', 'rpars'],
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=['rayleigh'])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['skip_rayleigh'], rtol=rtol)



def test_transmission_skip_dalgarno(tmp_path):
    reset = {
        'rayleigh': 'dalgarno_H2 dalgarno_He',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
        remove=['clouds', 'rpars'],
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=['dalgarno_H2'])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['skip_dalgarno'], rtol=rtol)


def test_transmission_skip_lecavelier(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        remove=['clouds'],
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=['lecavelier'])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['skip_rayleigh'], rtol=rtol)


def test_transmission_skip_cloud(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=['cloud'])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['skip_cloud'], rtol=rtol)


def test_transmission_skip_deck(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=['deck'])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['skip_deck'], rtol=rtol)


def test_transmission_skip_H2O_line_sample(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        remove=['clouds'],
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=['H2O'])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['skip_ls'], rtol=rtol)


def test_transmission_skip_H2O_lbl(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['clouds'],
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=['H2O'])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['skip_lbl'], rtol=rtol)


def test_transmission_skip_nothing(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        remove=['clouds'],
    )
    pyrat = pb.Pyrat(cfg)
    pyrat.run(skip=['this_wont_match'])
    spectrum = pyrat.spec.spectrum
    np.testing.assert_allclose(spectrum, expected['all_ls'], rtol=rtol)


