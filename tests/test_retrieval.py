# Copyright (c) 2021-2024 Patricio Cubillos
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
def test_retrieval_parameters(tmp_path):
    reset = {
        'rt_path': 'emission',
        'kurucz': f'{ROOT}tests/inputs/mock_fp00k0odfnew.pck',
        'log_gstar': '4.5',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/retrieval_transmission_tea.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)


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


# Short run will trigger warning when evaluating GR tests, ignore it.
@pytest.mark.filterwarnings("ignore: divide by zero encountered")
def test_mcmc_transmission(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/mcmc_transmission_test.cfg',
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

