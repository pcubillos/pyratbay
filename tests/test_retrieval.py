# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import pytest
import re

from conftest import make_config

import numpy as np
import pyratbay as pb
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


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
    print(repr(pyrat.obs.uncert*1000))
    expected_uncert = np.array([
       0.00063245553, 0.00063245553, 0.00063245553, 0.00063245553,
       0.00063245553, 0.00031685959, 0.00031685959, 0.00031685959,
       0.00031685959, 0.00031685959, 0.00031685959, 0.00031685959,
       0.00031685959, 0.00031685959, 0.00031685959,
    ])
    np.testing.assert_allclose(pyrat.obs.uncert, expected_uncert)


def test_eval_uncert_scaling_no_data(tmp_path):
    uncert_scaling = """
        err_quad_WFC3 2.0
    """
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

