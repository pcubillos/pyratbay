# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import pytest
import re
import numpy as np
import pyratbay.atmosphere as pa
import pyratbay.constants as pc
import pyratbay.opacity as op
import pyratbay.spectrum as ps


# Setup:
wn_min = 1e4/1.0
wn_max = 1e4/0.5
resolution = 15000.0
wn = ps.constant_resolution_spectrum(wn_min, wn_max, resolution)

# For the atmoshere, consider a simple solar-abundance isothermal atmosphere
nlayers = 6
pressure = pa.pressure('1e-8 bar', '1e2 bar', nlayers)
temperature = np.tile(1000.0, nlayers)
species = ['Na', 'K', 'H2', 'H', 'He']

vmr = pa.chemistry('tea', pressure, temperature, species).vmr
number_densities = pa.ideal_gas_density(vmr, pressure, temperature)

Na_density = number_densities[:,0]
K_density = number_densities[:,1]

root = f"{pc.ROOT}tests/expected/"


@pytest.mark.parametrize(
    'model',
    ['sodium', 'potassium'],
)
def test_alkali_init(model):
    if model == 'sodium':
        alkali = op.alkali.SodiumVdW(pressure, wn=wn)
        mol = 'Na'
    elif model == 'potassium':
        alkali = op.alkali.PotassiumVdW(pressure, wn=wn)
        mol = 'K'

    assert alkali.name == f'{model}_vdw'
    assert alkali.species == mol
    assert alkali.nlayers == len(pressure)
    assert alkali.nwave == len(wn)
    assert alkali.cross_section is None
    np.testing.assert_allclose(alkali.wn, wn)
    np.testing.assert_allclose(alkali.pressure, pressure)
    np.testing.assert_allclose(alkali.cutoff, 4500.0)


@pytest.mark.parametrize(
    'model',
    ['sodium', 'potassium'],
)
def test_alkali_wl_init(model):
    wl = 1/(wn*pc.um)
    if model == 'sodium':
        alkali = op.alkali.SodiumVdW(pressure, wl=wl)
        mol = 'Na'
    elif model == 'potassium':
        alkali = op.alkali.PotassiumVdW(pressure, wl=wl)
        mol = 'K'

    assert alkali.name == f'{model}_vdw'
    assert alkali.species == mol
    assert alkali.nlayers == len(pressure)
    assert alkali.nwave == len(wn)
    assert alkali.cross_section is None
    np.testing.assert_allclose(alkali.wn, wn)
    np.testing.assert_allclose(alkali.pressure, pressure)
    np.testing.assert_allclose(alkali.cutoff, 4500.0)


@pytest.mark.parametrize(
    'model',
    ['sodium', 'potassium'],
)
def test_alkali_duplicated_wave(model):
    wl = 1/(wn*pc.um)
    if model == 'sodium':
        alkali_model = op.alkali.SodiumVdW
    elif model == 'potassium':
        alkali_model = op.alkali.PotassiumVdW

    error = 'Either provide wavelength or wavenumber array, not both'
    with pytest.raises(ValueError, match=re.escape(error)):
        alkali_model(pressure, wl=wl, wn=wn)


@pytest.mark.parametrize(
    'model',
    ['sodium', 'potassium'],
)
def test_alkali_missing_wave(model):
    if model == 'sodium':
        alkali_model = op.alkali.SodiumVdW
    elif model == 'potassium':
        alkali_model = op.alkali.PotassiumVdW

    error = 'Neither of wavelength (wl) nor wavenumber (wn) were provided'
    with pytest.raises(ValueError, match=re.escape(error)):
        alkali_model(pressure)


def test_alkali_wn_order():
    # Order of wn array does not matter:
    wn1 = wn
    wn2 = wn[::-1]
    sodium1 = op.alkali.SodiumVdW(pressure, wn=wn1)
    sodium2 = op.alkali.SodiumVdW(pressure, wn=wn2)

    Na_ext1 = sodium1.calc_extinction_coefficient(temperature, Na_density)
    Na_ext2 = sodium2.calc_extinction_coefficient(temperature, Na_density)

    np.testing.assert_allclose(Na_ext1, np.fliplr(Na_ext2))


def test_sodium_cross_section():
    wn_min = 1e4/0.65
    wn_max = 1e4/0.55
    resolution = 15000.0
    wn = ps.constant_resolution_spectrum(wn_min, wn_max, resolution)
    alkali = op.alkali.SodiumVdW(pressure, wn=wn, cutoff=1000.0)

    with np.load(f'{root}expected_alkali_Na_opacity.npz') as d:
        expected_cs1 = d['expected_cs1']
        expected_cs2 = d['expected_cs2']

    temperature1 = np.tile(1000.0, nlayers)
    alkali.calc_cross_section(temperature1)
    np.testing.assert_allclose(alkali.cross_section, expected_cs1)

    temperature2 = np.tile(2500.0, nlayers)
    alkali.calc_cross_section(temperature2)
    np.testing.assert_allclose(alkali.cross_section, expected_cs2)


def test_sodium_extinction_coefficient():
    wn_min = 1e4/0.65
    wn_max = 1e4/0.55
    resolution = 15000.0
    wn = ps.constant_resolution_spectrum(wn_min, wn_max, resolution)
    alkali = op.alkali.SodiumVdW(pressure, wn=wn, cutoff=1000.0)

    with np.load(f'{root}expected_alkali_Na_opacity.npz') as d:
        expected_cs1 = d['expected_cs1']
        expected_cs2 = d['expected_cs2']
        expected_ec1 = d['expected_ec1']
        expected_ec2 = d['expected_ec2']

    temperature1 = np.tile(1000.0, nlayers)
    ec1 = alkali.calc_extinction_coefficient(temperature1, Na_density)
    np.testing.assert_allclose(alkali.cross_section, expected_cs1)
    np.testing.assert_allclose(ec1, expected_ec1)

    temperature2 = np.tile(2500.0, nlayers)
    ec2 = alkali.calc_extinction_coefficient(temperature2, Na_density)
    np.testing.assert_allclose(alkali.cross_section, expected_cs2)
    np.testing.assert_allclose(ec2, expected_ec2)


def test_potassium_cross_section():
    wn_min = 1e4/0.84
    wn_max = 1e4/0.70
    resolution = 15000.0
    wn = ps.constant_resolution_spectrum(wn_min, wn_max, resolution)
    alkali = op.alkali.PotassiumVdW(pressure, wn=wn, cutoff=1000.0)

    with np.load(f'{root}expected_alkali_K_opacity.npz') as d:
        expected_cs1 = d['expected_cs1']
        expected_cs2 = d['expected_cs2']

    temperature1 = np.tile(1000.0, nlayers)
    alkali.calc_cross_section(temperature1)
    np.testing.assert_allclose(alkali.cross_section, expected_cs1)

    temperature2 = np.tile(2500.0, nlayers)
    alkali.calc_cross_section(temperature2)
    np.testing.assert_allclose(alkali.cross_section, expected_cs2)


def test_potassium_extinction_coefficient():
    wn_min = 1e4/0.84
    wn_max = 1e4/0.70
    resolution = 15000.0
    wn = ps.constant_resolution_spectrum(wn_min, wn_max, resolution)
    alkali = op.alkali.PotassiumVdW(pressure, wn=wn, cutoff=1000.0)

    with np.load(f'{root}expected_alkali_K_opacity.npz') as d:
        expected_cs1 = d['expected_cs1']
        expected_cs2 = d['expected_cs2']
        expected_ec1 = d['expected_ec1']
        expected_ec2 = d['expected_ec2']

    temperature1 = np.tile(1000.0, nlayers)
    ec1 = alkali.calc_extinction_coefficient(temperature1, K_density)
    np.testing.assert_allclose(alkali.cross_section, expected_cs1)
    np.testing.assert_allclose(ec1, expected_ec1)

    temperature2 = np.tile(2500.0, nlayers)
    ec2 = alkali.calc_extinction_coefficient(temperature2, K_density)
    np.testing.assert_allclose(alkali.cross_section, expected_cs2)
    np.testing.assert_allclose(ec2, expected_ec2)



def test_sodium_layer_extinction_coefficient():
    wn_min = 1e4/0.65
    wn_max = 1e4/0.55
    resolution = 15000.0
    wn = ps.constant_resolution_spectrum(wn_min, wn_max, resolution)
    alkali = op.alkali.SodiumVdW(pressure, wn=wn, cutoff=1000.0)

    with np.load(f'{root}expected_alkali_Na_opacity.npz') as d:
        expected_cs1 = d['expected_cs1']
        expected_cs2 = d['expected_cs2']
        expected_ec1 = d['expected_ec1']
        expected_ec2 = d['expected_ec2']

    layer = 1
    temperature1 = np.tile(1000.0, nlayers)
    ec1 = alkali.calc_extinction_coefficient(temperature1, Na_density, layer)
    np.testing.assert_allclose(ec1, expected_ec1[layer])

