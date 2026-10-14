# Copyright (c) 2021-2026 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import pytest
import re

from conftest import make_config

import numpy as np
import pyratbay as pb
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# evals
keys = [
    'sampled_res', 'conv_800', 'conv_400',
    'rv_0000', 'rv_2000', 'rv_1000',
]

expected = {
    key:np.load(f"{ROOT}tests/expected/"
                f"expected_spectrum_eval_{key}_test.npz")['arr_0']
    for key in keys
}



# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
def test_retrieval_texnames(tmp_path):
    pyrat = pb.run(ROOT+'tests/configs/spectrum_emission_texnames.cfg')

    expected_texnames = [
        '$T_{\\rm irr} (K)$',
        '$M_{\\rm p}$ ($M_{\\rm Jup}$)',
        '$\\alpha$',
        '$\\log\\ \\kappa_{\\rm ray}$',
        '$\\alpha_{\\rm ray}$',
        '[O/H]',
        '[M/H]',
        'C/O',
    ]
    assert pyrat.ret.texnames == expected_texnames


def test_eval_offset_data_defaults(tmp_path):
    offset_inst = """
        offset_STIS
        offset_WFC3
    """
    reset = {
        'offset_inst': offset_inst,
        'obsfile': '{ROOT}tests/inputs/obs_file_err_scaling.dat',
        'dunits': 'ppm',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    expected_data = np.array([
        0.006647, 0.006626, 0.006578, 0.006549, 0.00655,
        0.006591, 0.00671,  0.006721, 0.006721, 0.006683,
        0.006646, 0.006605, 0.00655,  0.006537, 0.006541,
    ])
    np.testing.assert_allclose(pyrat.obs.data, expected_data)


def test_eval_offset_data_fixed(tmp_path):
    offset_inst = """
        offset_STIS 5000
        offset_WFC3 2000
    """
    reset = {
        'offset_inst': offset_inst,
        'obsfile': '{ROOT}/tests/inputs/obs_file_err_scaling.dat',
        'dunits': 'ppm',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    expected_data = np.array([
        0.011647, 0.011626, 0.011578, 0.011549, 0.01155,
        0.008591, 0.00871,  0.008721, 0.008721, 0.008683,
        0.008646, 0.008605, 0.00855,  0.008537, 0.008541,
    ])
    np.testing.assert_allclose(pyrat.obs.data, expected_data)


def test_eval_offset_data_fit_and_fix(tmp_path):
    offset_inst = """
        offset_STIS 3000
        offset_WFC3
    """
    ret_params = """
        offset_WFC3  2000 -3.0  3.0 1.0
    """
    reset = {
        'offset_inst': offset_inst,
        'obsfile': '{ROOT}/tests/inputs/obs_file_err_scaling.dat',
        'dunits': 'ppm',
        'retrieval_params': ret_params,
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    # Only STIS uncertainties have changed
    expected_data = np.array([
        0.009647, 0.009626, 0.009578, 0.009549, 0.00955,
        0.006591, 0.00671,  0.006721, 0.006721, 0.006683,
        0.006646, 0.006605, 0.00655,  0.006537, 0.006541,
    ])
    np.testing.assert_allclose(pyrat.obs.data, expected_data)

    # Now modify WFC3 uncertainties
    pyrat.eval(pyrat.ret.params)
    expected_data = np.array([
        0.009647, 0.009626, 0.009578, 0.009549, 0.00955,
        0.008591, 0.00871,  0.008721, 0.008721, 0.008683,
        0.008646, 0.008605, 0.00855,  0.008537, 0.008541,
    ])
    np.testing.assert_allclose(pyrat.obs.data, expected_data)


def test_eval_offset_data_no_data(tmp_path):
    offset_inst = "offset_WFC3 2.0"
    reset = {
        'offset_inst': offset_inst,
        'dunits': 'ppm',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        "Invalid instrumental offset parameter 'offset_WFC3'. There is no "
        "instrument matching the name 'WFC3'"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)



def test_eval_uncert_scaling_defaults(tmp_path):
    uncert_scaling = """
        err_scale_STIS
        err_quad_WFC3
    """
    reset = {
        'uncert_scaling': uncert_scaling,
        'obsfile': '{ROOT}/tests/inputs/obs_file_err_scaling.dat',
        'dunits': 'ppm',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    expected_uncert = np.array([
        2.e-05, 2.e-05, 2.e-05, 2.e-05, 2.e-05, 2.e-05, 2.e-05, 2.e-05,
        2.e-05, 2.e-05, 2.e-05, 2.e-05, 2.e-05, 2.e-05, 2.e-05,
    ])
    np.testing.assert_allclose(pyrat.obs.uncert, expected_uncert)


def test_eval_uncert_scaling_fixed(tmp_path):
    uncert_scaling = """
        err_scale_STIS 1.5
        err_quad_WFC3 2.0
    """
    reset = {
        'uncert_scaling': uncert_scaling,
        'obsfile': '{ROOT}/tests/inputs/obs_file_err_scaling.dat',
        'dunits': 'ppm',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    expected_uncert = np.array([
       0.00063245553, 0.00063245553, 0.00063245553, 0.00063245553,
       0.00063245553, 0.00010198039, 0.00010198039, 0.00010198039,
       0.00010198039, 0.00010198039, 0.00010198039, 0.00010198039,
       0.00010198039, 0.00010198039, 0.00010198039,
    ])
    np.testing.assert_allclose(pyrat.obs.uncert, expected_uncert)


def test_eval_uncert_scaling_fit_and_fix(tmp_path):
    uncert_scaling = """
        err_scale_STIS 1.5
        err_quad_WFC3
    """
    ret_params = """
        err_quad_WFC3  2.5 -3.0  3.0 1.0
    """
    reset = {
        'uncert_scaling': uncert_scaling,
        'obsfile': '{ROOT}/tests/inputs/obs_file_err_scaling.dat',
        'dunits': 'ppm',
        'retrieval_params': ret_params,
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    # Only STIS uncertainties have changed
    expected_uncert = np.array([
       0.00063245553, 0.00063245553, 0.00063245553, 0.00063245553,
       0.00063245553, 0.00002, 0.00002, 0.00002, 0.00002, 0.00002,
       0.00002, 0.00002, 0.00002, 0.00002, 0.00002,
    ])
    np.testing.assert_allclose(pyrat.obs.uncert, expected_uncert)

    # Now modify WFC3 uncertainties
    pyrat.eval(pyrat.ret.params)
    expected_uncert = np.array([
       0.00063245553, 0.00063245553, 0.00063245553, 0.00063245553,
       0.00063245553, 0.00031685959, 0.00031685959, 0.00031685959,
       0.00031685959, 0.00031685959, 0.00031685959, 0.00031685959,
       0.00031685959, 0.00031685959, 0.00031685959,
    ])
    np.testing.assert_allclose(pyrat.obs.uncert, expected_uncert)


def test_eval_uncert_scaling_no_data(tmp_path):
    uncert_scaling = "err_quad_WFC3 2.0"
    reset = {
        'uncert_scaling': uncert_scaling,
        'dunits': 'ppm',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        "Invalid retrieval parameter 'err_quad_WFC3'. There is no "
        "instrument matching the name 'WFC3'"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_eval_inst_res_no_obs(tmp_path):
    cfg = ROOT+'tests/configs/spectrum_inst_res_no_obs.cfg'
    pyrat = pb.run(cfg)
    pyrat.eval(pyrat.ret.params)
    # full-resolution
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['sampled_res'])
    # convolved
    assert len(pyrat.obs.convolved_spectrum) == 1
    np.testing.assert_allclose(pyrat.obs.convolved_spectrum[0], expected['conv_800'])
    assert len(pyrat.obs.bandflux) == 0


@pytest.mark.parametrize(
    'inst_res',
    ['all', 'named'],
)
def test_eval_inst_res_single_obs(inst_res):
    cfg = f'{ROOT}tests/configs/spectrum_inst_res_sinfoni_{inst_res}.cfg'
    pyrat = pb.run(cfg)
    pyrat.eval(pyrat.ret.params)

    expected_band_flux = [
        6.58482163e-03, 6.60168415e-03, 6.54978523e-03, 6.47581211e-03,
        6.48394974e-03,
    ]

    # full-resolution
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['sampled_res'])
    # convolved
    assert len(pyrat.obs.convolved_spectrum) == 1
    np.testing.assert_allclose(pyrat.obs.convolved_spectrum[0], expected['conv_800'])
    np.testing.assert_allclose(pyrat.obs.bandflux, expected_band_flux)


def test_eval_inst_res_single_obs_missing():
    cfg = ROOT+'tests/configs/spectrum_inst_res_sinfoni_missing.cfg'
    error = re.escape(
        "Missing instrumental resolution for ['sinfoni'] data, all "
        "high-resolution data must have an instrumental resolution"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_eval_inst_res_multi_obs_all():
    cfg = f'{ROOT}tests/configs/spectrum_inst_res_gravity_sinfoni_all.cfg'
    pyrat = pb.run(cfg)
    pyrat.eval(pyrat.ret.params)

    expected_band_flux = [
        6.58482163e-03, 6.60168415e-03, 6.54978523e-03, 6.47581211e-03,
        6.48394974e-03, 6.53280975e-03, 6.73603680e-03, 6.69654605e-03,
        6.71894321e-03, 6.71866380e-03, 6.61848780e-03, 6.57120856e-03,
        6.50441699e-03, 6.53682492e-03, 6.51843522e-03,
    ]

    # full-resolution
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['sampled_res'])
    # convolved
    assert len(pyrat.obs.convolved_spectrum) == 1
    np.testing.assert_allclose(pyrat.obs.convolved_spectrum[0], expected['conv_800'])
    np.testing.assert_allclose(pyrat.obs.bandflux, expected_band_flux)


def test_eval_inst_res_multi_obs_named():
    cfg = f'{ROOT}tests/configs/spectrum_inst_res_gravity_sinfoni_named.cfg'
    pyrat = pb.run(cfg)
    pyrat.eval(pyrat.ret.params)

    expected_band_flux = [
        6.58482163e-03, 6.60168415e-03, 6.54978523e-03, 6.47581211e-03,
        6.48394974e-03, 6.53731146e-03, 6.71874367e-03, 6.71524846e-03,
        6.72491824e-03, 6.70641183e-03, 6.62658782e-03, 6.58294680e-03,
        6.51023138e-03, 6.52724592e-03, 6.52589355e-03,
    ]

    # full-resolution
    np.testing.assert_allclose(pyrat.spec.spectrum, expected['sampled_res'])
    # convolved
    assert len(pyrat.obs.convolved_spectrum) == 2
    np.testing.assert_allclose(pyrat.obs.convolved_spectrum[0], expected['conv_800'])
    np.testing.assert_allclose(pyrat.obs.convolved_spectrum[1], expected['conv_400'])
    np.testing.assert_allclose(pyrat.obs.bandflux, expected_band_flux)


def test_eval_inst_res_multi_obs_missing():
    cfg = ROOT+'tests/configs/spectrum_inst_res_gravity_sinfoni_missing.cfg'
    error = re.escape(
        "Missing instrumental resolution for ['gravity' 'sinfoni'] data, all "
        "high-resolution data must have an instrumental resolution"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_eval_inst_res_multi_obs_misnamed():
    cfg = ROOT+'tests/configs/spectrum_inst_res_gravity_sinfoni_misnamed.cfg'
    error = re.escape(
        "Invalid inst_resolution input, there is no instrument matching the "
        "name 'soss'"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_eval_high_resolution():
    cfg = f'{ROOT}tests/configs/spectrum_inst_res_gravity_sinfoni_named.cfg'
    pyrat = pb.run(cfg)
    pyrat.eval(pyrat.ret.params)

    # half-width set by inst_resolution
    expected_half_widths = [
       7.05625e-04, 7.28750e-04, 7.51875e-04, 7.75000e-04, 7.98125e-04,
       1.64125e-03, 1.68750e-03, 1.73375e-03, 1.78000e-03, 1.82625e-03,
       1.87125e-03, 1.91750e-03, 1.96375e-03, 2.01000e-03, 2.05625e-03,
    ]
    np.testing.assert_allclose(pyrat.obs.half_widths, expected_half_widths)


def test_eval_rv_shift_zero():
    cfg = f'{ROOT}tests/configs/spectrum_inst_res_gravity_sinfoni_named.cfg'
    pyrat = pb.run(cfg)
    pyrat.eval([0.0])

    expected_band_flux = [
        6.58482163e-03, 6.60168415e-03, 6.54978523e-03, 6.47581211e-03,
        6.48394974e-03, 6.53731146e-03, 6.71874367e-03, 6.71524846e-03,
        6.72491824e-03, 6.70641183e-03, 6.62658782e-03, 6.58294680e-03,
        6.51023138e-03, 6.52724592e-03, 6.52589355e-03,
    ]
    assert len(pyrat.obs.rv_wn) == 2
    np.testing.assert_allclose(pyrat.obs.rv_wn[0], expected['rv_0000'])
    np.testing.assert_allclose(pyrat.obs.rv_wn[0], pyrat.spec.wn)
    np.testing.assert_allclose(pyrat.obs.bandflux, expected_band_flux)


def test_eval_rv_shift_redshift():
    cfg = f'{ROOT}tests/configs/spectrum_inst_res_gravity_sinfoni_named.cfg'
    pyrat = pb.run(cfg)
    pyrat.eval([2000.0])

    expected_band_flux = [
        6.59412279e-03, 6.57537475e-03, 6.54950339e-03, 6.47174060e-03,
        6.45805756e-03, 6.53410677e-03, 6.69524995e-03, 6.72523065e-03,
        6.69117916e-03, 6.66708765e-03, 6.63316833e-03, 6.58372580e-03,
        6.52332997e-03, 6.50911580e-03, 6.53868303e-03,
    ]
    assert len(pyrat.obs.rv_wn) == 2
    np.testing.assert_allclose(pyrat.obs.rv_wn[0], expected['rv_2000'])
    np.testing.assert_allclose(pyrat.obs.bandflux, expected_band_flux)


def test_eval_rv_shift_named():
    cfg = f'{ROOT}tests/configs/spectrum_inst_res_gravity_sinfoni_rv_named.cfg'
    pyrat = pb.run(cfg)
    pyrat.eval([2000.0, 1000.0])

    expected_band_flux = [
        6.59412279e-03, 6.57537475e-03, 6.54950339e-03, 6.47174060e-03,
        6.45805756e-03, 6.55074306e-03, 6.71978349e-03, 6.73919255e-03,
        6.72296712e-03, 6.66128073e-03, 6.67522315e-03, 6.60662684e-03,
        6.52594817e-03, 6.51660535e-03, 6.52494796e-03,
    ]
    assert len(pyrat.obs.rv_wn) == 3
    np.testing.assert_allclose(pyrat.obs.rv_wn[0], expected['rv_2000'])
    np.testing.assert_allclose(pyrat.obs.rv_wn[1], expected['rv_1000'])
    np.testing.assert_allclose(pyrat.obs.bandflux, expected_band_flux)


def test_eval_rv_shift_misnamed():
    cfg = f'{ROOT}tests/configs/spectrum_inst_res_gravity_sinfoni_rv_misnamed.cfg'
    error = re.escape(
        "Invalid RV parameter 'rv_soss'. There is no instrument matching "
        "the name 'soss'"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_hybrid_eval_reject_VMR_OOBs():
    cfg = ROOT+'tests/configs/retrieval_transmission_hybrid.cfg'
    pyrat = pb.Pyrat(cfg)
    pyrat.eval(pyrat.ret.params)

    # Spectrum exists
    expected_spec = [0.00679258, 0.00695102, 0.0067942 , 0.00682667, 0.00679244]
    np.testing.assert_allclose(pyrat.spec.spectrum[0:5], expected_spec, rtol=1e-6)
    # Evaluation at data points is np.inf (rejected)
    np.testing.assert_array_equal(pyrat.obs.bandflux, np.inf)


@pytest.mark.skip(reason='Too long to execute')
@pytest.mark.filterwarnings("ignore: divide by zero encountered")
@pytest.mark.filterwarnings("ignore: Attempting to set identical low and high")
def test_mcmc_transmission(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/retrieval_transmission_test.cfg',
    )
    pyrat = pb.run(cfg)


@pytest.mark.skip(reason='TBD')
def test_mcmc_emission(tmp_path):
    reset = {
        'rt_path': 'emission',
        'kurucz': f'{ROOT}tests/inputs/mock_fp00k0odfnew.pck',
        'log_gstar': '4.5',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/mcmc_transmission_test.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    #error = re.escape(undefined_mcmc[param])
    #with pytest.raises(ValueError, match=error):
    #    pyrat = pb.run(cfg)

