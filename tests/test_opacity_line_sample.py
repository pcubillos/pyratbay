# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import pytest
import re

import numpy as np
import pyratbay.atmosphere as pa
import pyratbay.constants as pc
import pyratbay.opacity as op



@pytest.mark.parametrize('cs_input', ['str', 'list'])
def test_line_sample_init(cs_input):
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    if cs_input == 'list':
        ls = op.Line_Sample([cs_files])
    else:
        ls = op.Line_Sample(cs_files)

    assert ls.cs_files[0] == cs_files
    assert ls.ntemp == 10
    assert ls.nwave == 3209
    assert ls.nspec == 1
    assert ls.nlayers == 51

    expected_temp = np.array([
         300.,  600.,  900., 1200., 1500., 1800., 2100., 2400., 2700., 3000.
    ])
    expected_wn = np.arange(5882.35294118, 9091.0, 1.0)
    expected_pressure = pa.pressure('1e-6 bar', '1e2 bar', nlayers=51)

    assert ls.species == ['H2O']
    np.testing.assert_allclose(ls.wn, expected_wn)
    np.testing.assert_allclose(ls.temp, expected_temp)
    np.testing.assert_allclose(ls.press, expected_pressure, rtol=3e-5)

    np.testing.assert_allclose(ls.tmin, 300.0)
    np.testing.assert_allclose(ls.tmax, 3000.0)
    assert hasattr(ls, 'cs_table')


def test_line_sample_trim_wl():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files, min_wl=1.4, max_wl=1.6)

    assert ls.ntemp == 10
    assert ls.nwave == 893
    assert ls.nspec == 1
    assert ls.nlayers == 51

    expected_wn = np.arange(6250.35294118, 7143.0, 1.0)
    np.testing.assert_allclose(ls.wn, expected_wn)


def test_line_sample_trim_wn():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files, min_wn=6000, max_wn=7000)

    assert ls.ntemp == 10
    assert ls.nwave == 1000
    assert ls.nspec == 1
    assert ls.nlayers == 51

    expected_wn = np.arange(6000.35294118, 7000.0, 1.0)
    np.testing.assert_allclose(ls.wn, expected_wn)


def test_line_sample_thin_wn():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    thinning = 10
    ls = op.Line_Sample(cs_files, wn_thinning=thinning)

    assert ls.ntemp == 10
    assert ls.nwave == 321
    assert ls.nspec == 1
    assert ls.nlayers == 51

    expected_wn = np.arange(5882.35294118, 9091.0, 1.0)[::thinning]
    np.testing.assert_allclose(ls.wn, expected_wn)


def test_line_sample_duplicate_wn_min():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"

    error = re.escape("Either define min_wn or max_wl, not both")
    with pytest.raises(ValueError, match=error):
        ls = op.Line_Sample(cs_files, min_wn=6000, max_wl=1.6)


def test_line_sample_duplicate_wn_max():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"

    error = re.escape("Either define min_wl or max_wn, not both")
    with pytest.raises(ValueError, match=error):
        ls = op.Line_Sample(cs_files, max_wn=7000, min_wl=1.4)


def test_line_sample_trim_get_wl():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files, min_wl=1.4, max_wl=1.6)

    expected_wl = 1e4/np.arange(6250.35294118, 7143.0, 1.0)
    wl = ls.get_wl()
    np.testing.assert_allclose(wl, expected_wl)


def test_line_sample_trim_get_units_wl():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files, min_wl=1.4, max_wl=1.6)

    expected_wl = 1e8/np.arange(6250.35294118, 7143.0, 1.0)
    wl = ls.get_wl('A')
    np.testing.assert_allclose(wl, expected_wl)


def test_line_sample_trim_bad_get_wl():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files, min_wl=1.4, max_wl=1.6)

    error = re.escape(
        "Invalid wavelength units 'rjup', select one from "
        "['A', 'nm', 'um', 'mm', 'cm', 'm', 'km']"
    )
    with pytest.raises(ValueError, match=error):
        ls.get_wl('rjup')


def test_line_sample_resample_pressure():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    pressure = pa.pressure('1e-8 bar', '1e-3 bar', nlayers=90)
    ls = op.Line_Sample(cs_files, pressure=pressure)

    assert ls.ntemp == 10
    assert ls.nwave == 3209
    assert ls.nspec == 1
    assert ls.nlayers == 90

    np.testing.assert_allclose(ls.press, pressure, rtol=3e-5)


def test_line_sample_extrapolate_low_pressure():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    pressure = pa.pressure('1e-8 bar', '10 bar', nlayers=90)
    ls = op.Line_Sample(cs_files, pressure=pressure)

    assert ls.ntemp == 10
    assert ls.nwave == 3209
    assert ls.nspec == 1
    assert ls.nlayers == 90

    np.testing.assert_allclose(ls.press, pressure, rtol=3e-5)


# Cross sections
def test_line_sample_cross_section():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files, min_wn=9000.0, max_wn=9040.0)

    with np.load(f"{pc.ROOT}tests/expected/expected_ls_H2O_opacity.npz") as d:
        expected_cs = d['expected_cs']

    temp = np.tile(1200.0, ls.nlayers)
    cross_section = ls.calc_cross_section(temp)
    np.testing.assert_allclose(cross_section, expected_cs)


def test_line_sample_single_cross_section():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files, min_wn=9000.0, max_wn=9040.0)

    with np.load(f"{pc.ROOT}tests/expected/expected_ls_H2O_opacity.npz") as d:
        expected_cs = d['expected_cs']

    temp = np.tile(1200.0, ls.nlayers)
    layer = 40
    cross_section = ls.calc_cross_section(temp, layer=layer)
    np.testing.assert_allclose(cross_section, expected_cs[layer])


def test_line_sample_per_mol_cross_section():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files, min_wn=9000.0, max_wn=9040.0)

    with np.load(f"{pc.ROOT}tests/expected/expected_ls_H2O_opacity.npz") as d:
        expected_cs = d['expected_cs']

    temp = np.tile(1200.0, ls.nlayers)
    cross_section = ls.calc_cross_section(temp, per_mol=True)
    assert np.shape(cross_section) == (1, 51, 40)
    np.testing.assert_allclose(cross_section[0], expected_cs)


def test_line_sample_per_mol_single_cross_section():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files, min_wn=9000.0, max_wn=9040.0)

    with np.load(f"{pc.ROOT}tests/expected/expected_ls_H2O_opacity.npz") as d:
        expected_cs = d['expected_cs']

    temp = np.tile(1200.0, ls.nlayers)
    cross_section = ls.calc_cross_section(temp, per_mol=True)
    assert np.shape(cross_section) == (1, 51, 40)
    np.testing.assert_allclose(cross_section[0], expected_cs)


# Extinctions
def test_line_sample_extinction_coefficient():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files, min_wn=9000.0, max_wn=9040.0)

    with np.load(f"{pc.ROOT}tests/expected/expected_ls_H2O_opacity.npz") as d:
        expected_ec = d['expected_ec']

    temp = np.tile(1200.0, ls.nlayers)
    vmr = np.tile(3e-4, (ls.nlayers,ls.nspec))
    densities = pa.ideal_gas_density(vmr, ls.press, temp)
    extinction = ls.calc_extinction_coefficient(temp, densities)
    np.testing.assert_allclose(extinction, expected_ec)


def test_line_sample_single_extinction_coefficient():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files, min_wn=9000.0, max_wn=9040.0)

    with np.load(f"{pc.ROOT}tests/expected/expected_ls_H2O_opacity.npz") as d:
        expected_ec = d['expected_ec']

    temp = np.tile(1200.0, ls.nlayers)
    vmr = np.tile(3e-4, (ls.nlayers,ls.nspec))
    densities = pa.ideal_gas_density(vmr, ls.press, temp)

    layer = 40
    extinction = ls.calc_extinction_coefficient(temp, densities, layer=layer)
    np.testing.assert_allclose(extinction, expected_ec[layer])


def test_line_sample_per_mol_extinction_coefficient():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files, min_wn=9000.0, max_wn=9040.0)

    with np.load(f"{pc.ROOT}tests/expected/expected_ls_H2O_opacity.npz") as d:
        expected_ec = d['expected_ec']

    temp = np.tile(1200.0, ls.nlayers)
    vmr = np.tile(3e-4, (ls.nlayers,ls.nspec))
    densities = pa.ideal_gas_density(vmr, ls.press, temp)
    extinction = ls.calc_extinction_coefficient(temp, densities, per_mol=True)

    assert np.shape(extinction) == (1, 51, 40)
    np.testing.assert_allclose(extinction[0], expected_ec)


def test_line_sample_per_mol_single_extinction_coefficient():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files, min_wn=9000.0, max_wn=9040.0)

    with np.load(f"{pc.ROOT}tests/expected/expected_ls_H2O_opacity.npz") as d:
        expected_ec = d['expected_ec']

    temp = np.tile(1200.0, ls.nlayers)
    vmr = np.tile(3e-4, (ls.nlayers,ls.nspec))
    densities = pa.ideal_gas_density(vmr, ls.press, temp)
    layer = 40
    extinction = ls.calc_extinction_coefficient(
        temp, densities, layer=layer, per_mol=True,
    )

    assert np.shape(extinction) == (1, 40)
    np.testing.assert_allclose(extinction[0], expected_ec[layer])


# Errors
def test_line_sample_missing_cs_file():
    cs_files = "non_exttable_test_300-3000K_1.1-1.7um.npz"
    error = re.escape(
        "Missing opacity files: ['non_exttable_test_300-3000K_1.1-1.7um.npz']"
    )
    with pytest.raises(ValueError, match=error):
        ls = op.Line_Sample(cs_files)


@pytest.mark.skip(reason='TBD')
@pytest.mark.parametrize(
    'mismatch',
    ['wavenumber', 'temperature'],
)
def test_line_sample_mismatch_sizes(mismatch):
    cs_files = [
        f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz",
        # TBD: a second CS with different shape
    ]
    error = re.escape(
        f"Tabulated {mismatch} values in file '{cs_files}' "
        "do not match with previous arrays"
    )
    with pytest.raises(ValueError, match=error):
        ls = op.Line_Sample(cs_files)


@pytest.mark.parametrize('call', ['cross_section', 'extinction'])
def test_line_sample_temp_outbounds(call):
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    ls = op.Line_Sample(cs_files)

    temp = np.tile(9000.0, ls.nlayers)
    vmr = np.tile(3e-4, (ls.nlayers,ls.nspec))
    densities = pa.ideal_gas_density(vmr, ls.press, temp)

    error = re.escape('Temperatures are out of line-sample bounds')
    if call == 'cross_section':
        with pytest.raises(ValueError, match=error):
            cs = ls.calc_cross_section(temp)
    if call == 'extinction':
        with pytest.raises(ValueError, match=error):
            extinction = ls.calc_extinction_coefficient(temp, densities)


def test_line_sample_extrapolate_high_pressure():
    cs_files = f"{pc.ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz"
    pressure = pa.pressure('1e-6 bar', '1e3 bar', nlayers=90)
    error = re.escape(
        "Pressure profile extends beyond the maximum tabulated pressure"
    )
    with pytest.raises(ValueError, match=error):
        op.Line_Sample(cs_files, pressure=pressure)

