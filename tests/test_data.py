# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import struct
import pytest
import re

import mc3
import numpy as np

import pyratbay.atmosphere as pa
import pyratbay.constants as pc
import pyratbay.io as io
import pyratbay.tools as pt

os.chdir(pc.ROOT+'tests')


# Sample data
base_data = np.array([
    0.0221654, 0.0219554, 0.0218184, 0.0214919, 0.0213856,
    0.0209675, 0.0210958, 0.0213826, 0.0216047, 0.0226758,
    0.0220769, 0.0213577, 0.0214181, 0.0215911, 0.0213866,
    0.0210104, 0.021041 , 0.0199637, 0.0189196,
])

base_uncert = np.array([
    3.746e-04, 8.370e-05, 6.950e-05, 7.710e-05, 6.450e-05,
    7.390e-05, 6.770e-05, 7.290e-05, 7.740e-05, 9.270e-05,
    1.241e-04, 1.230e-04, 1.747e-04, 1.936e-04, 2.233e-04,
    2.885e-04, 2.954e-04, 4.220e-04, 8.458e-04,
])

base_names = np.array([
    'nirspec_g395h_nrs1', 'nirspec_g395h_nrs1', 'nirspec_g395h_nrs1',
    'nirspec_g395h_nrs1', 'nirspec_g395h_nrs1',
    'nirspec_g395h_nrs2', 'nirspec_g395h_nrs2', 'nirspec_g395h_nrs2',
    'nirspec_g395h_nrs2', 'nirspec_g395h_nrs2',
    'miri_lrs', 'miri_lrs', 'miri_lrs', 'miri_lrs', 'miri_lrs',
    'miri_lrs', 'miri_lrs', 'miri_lrs', 'miri_lrs',
])



def test_parse_error_param_scale():
    inst, texname, scaling = pt.parse_error_param('err_scale_WFC3')
    assert inst == 'WFC3'
    assert texname == '$\\log\\ S^\\sigma_{\\rm WFC3}$'
    assert scaling == 'scale'


def test_parse_error_param_quad():
    inst, texname, scaling = pt.parse_error_param('err_quad_WFC3')
    assert inst == 'WFC3'
    assert texname == '$\\log\\ \\sigma_{\\rm WFC3}$'
    assert scaling == 'quadrature'


def test_parse_error_param_underscores():
    inst, texname, scaling = pt.parse_error_param('err_quad_HST_WFC3')
    assert inst == 'HST WFC3'
    assert texname == '$\\log\\ \\sigma_{\\rm HST WFC3}$'
    assert scaling == 'quadrature'


def test_parse_error_param_fail():
    match = re.escape(
        "Invalid error scaling parameter 'err_fudging_IRAC1'. Valid "
        "options begin with: ['err_scale_', 'err_quad_']"
    )
    with pytest.raises(ValueError, match=match):
        pt.parse_error_param('err_fudging_IRAC1')


def test_Data_no_models():
    data = pt.Data(base_data, base_uncert, base_names)

    assert data.n_epars == 0
    np.testing.assert_allclose(data.data, base_data)
    np.testing.assert_allclose(data.uncert, base_uncert)


def test_Data_offset_str_input():
    data = pt.Data(base_data, base_uncert, base_names, offset_models='offset_nrs1')
    assert data.n_offsets == 1
    assert data.offset_models == ['offset_nrs1']

    offset_data = data.offset_data([4000.0], 'ppm')
    expected_offset = np.array([
        0.0261654, 0.0259554, 0.0258184, 0.0254919, 0.0253856,
        0.0209675, 0.0210958, 0.0213826, 0.0216047, 0.0226758,
        0.0220769, 0.0213577, 0.0214181, 0.0215911, 0.0213866,
        0.0210104, 0.021041 , 0.0199637, 0.0189196,
    ])
    np.testing.assert_allclose(offset_data, expected_offset)


def test_Data_offset_iterable_input():
    offset_models = ['nirspec', 'miri']
    data = pt.Data(base_data, base_uncert, base_names, offset_models)
    assert data.n_offsets == 2

    offset_data = data.offset_data([400.0, 800.0], 'ppm')
    expected_offset = np.array([
        0.0225654, 0.0223554, 0.0222184, 0.0218919, 0.0217856,
        0.0213675, 0.0214958, 0.0217826, 0.0220047, 0.0230758,
        0.0228769, 0.0221577, 0.0222181, 0.0223911, 0.0221866,
        0.0218104, 0.021841 , 0.0207637, 0.0197196,
    ])
    np.testing.assert_allclose(offset_data, expected_offset)


def test_Data_offset_not_found_band_name():
    match = re.escape(
        "Invalid instrumental offset parameter 'offset_prism'. There is no "
        "instrument matching the name 'prism'"
    )
    with pytest.raises(ValueError, match=match):
        data = pt.Data(
            base_data, base_uncert, base_names, offset_models='offset_prism',
        )


def test_Data_offset_overlapping_bands():
    err_models = ['nirspec', 'nrs1']
    match = "Multiple instrumental offsets apply to a same data point"
    with pytest.raises(ValueError, match=re.escape(match)):
        data = pt.Data(base_data, base_uncert, base_names, err_models)


def test_Data_error_multiplicative_str():
    err_model = 'err_scale_nrs1'
    data = pt.Data(base_data, base_uncert, base_names, err_models=err_model)
    assert data.n_epars == 1
    assert data.scaling_modes == ['scale']

    inflated_err = data.scale_errors([1.0])
    expected_err = np.array([
        3.746e-03, 8.370e-04, 6.950e-04, 7.710e-04, 6.450e-04,
        7.390e-05, 6.770e-05, 7.290e-05, 7.740e-05, 9.270e-05,
        1.241e-04, 1.230e-04, 1.747e-04, 1.936e-04, 2.233e-04,
        2.885e-04, 2.954e-04, 4.220e-04, 8.458e-04,
    ])
    np.testing.assert_allclose(inflated_err, expected_err)


def test_Data_error_quadrature_list():
    err_models = ['err_quad_nirspec', 'err_quad_miri']
    data = pt.Data(base_data, base_uncert, base_names, err_models=err_models)
    assert data.n_epars == 2
    assert data.scaling_modes == ['quadrature', 'quadrature']

    inflated_err = data.scale_errors([2.3, 2.5], 'ppm')
    expected_err = np.array([
       0.00042442417, 0.00021637099, 0.00021128409, 0.00021390448,
       0.00020969255, 0.00021277201, 0.00021069885, 0.00021242676,
       0.0002140128 , 0.00022000911, 0.00033970695, 0.00033930665,
       0.00036127564, 0.00037078425, 0.00038712129, 0.00042805636,
       0.00043273683, 0.0005273367 , 0.00090298264,
    ])
    np.testing.assert_allclose(inflated_err, expected_err)


def test_Data_error_not_found_band_name():
    match = re.escape(
        "Invalid retrieval parameter 'err_scale_prism'. There is no "
        "instrument matching the name 'prism'"
    )
    with pytest.raises(ValueError, match=match):
        e_models = 'err_scale_prism'
        data = pt.Data(base_data, base_uncert, base_names, err_models=e_models)


def test_Data_error_overlapping_bands():
    e_models = ['err_quad_nirspec', 'err_quad_nrs1']
    match = re.escape("Multiple uncertainty scaling apply to a same data point")
    with pytest.raises(ValueError, match=match):
        data = pt.Data(base_data, base_uncert, base_names, err_models=e_models)

