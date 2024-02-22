# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import pytest
import re

import numpy as np
import pyratbay.atmosphere as pa
import pyratbay.constants as pc
import pyratbay.opacity as op
import pyratbay.spectrum as ps


# Setup:
wn_min = 1e4/10.0
wn_max = 1e4/0.5
resolution = 15.0
wn = ps.constant_resolution_spectrum(wn_min, wn_max, resolution)

# For the atmoshere, consider a simple solar-abundance isothermal atmosphere
nlayers = 6
pressure = pa.pressure('1e-8 bar', '1e2 bar', nlayers)
temperature = np.tile(1200.0, nlayers)
species = ['H2', 'H', 'He']

vmr = pa.chemistry('tea', pressure, temperature, species).vmr
number_densities = pa.ideal_gas_density(vmr, pressure, temperature)

H2_density = number_densities[:,0]

c_root = f"{pc.ROOT}pyratbay/data/CIA/"
e_root = f"{pc.ROOT}tests/expected/"


def test_cia_init():
    cs_file = 'CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
    cia = op.Collision_Induced(f'{c_root}{cs_file}')

    default_temps = np.array([
         60.,  100.,  150.,  200.,  250.,  300.,  350.,  400.,  500.,
        600.,  700.,  800.,  900., 1000., 2000., 3000., 4000., 5000.,
       6000., 7000.,
    ])

    assert os.path.split(cia.cia_file)[1] == cs_file
    assert cia.species == ['H2', 'H2']
    assert cia.ntemp == 20
    assert cia.nwave == 824
    assert cia.nspec == 2
    np.testing.assert_allclose(cia.tmin, 60.0)
    np.testing.assert_allclose(cia.tmax, 7000.0)
    assert hasattr(cia, 'tab_cross_section')
    assert hasattr(cia, '_dcs_dt')
    np.testing.assert_allclose(cia.wn, np.arange(20, 16500, 20))
    np.testing.assert_allclose(cia.temps, default_temps)


def test_cia_wn_init():
    # Ensure that order of wn array does not matter:
    cs_file = 'CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
    wn1 = wn
    wn2 = wn[::-1]
    cia1 = op.Collision_Induced(f'{c_root}{cs_file}', wn=wn1)
    cia2 = op.Collision_Induced(f'{c_root}{cs_file}', wn=wn2)

    temperature = np.tile(1200.0, nlayers)
    cia_cs1 = cia1.calc_cross_section(temperature)
    cia_cs2 = cia2.calc_cross_section(temperature)

    np.testing.assert_allclose(cia_cs1, np.fliplr(cia_cs2))


def test_cia_wl_init():
    # Ensure that order of wl array does not matter:
    cs_file = 'CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
    wl1 = 1/wn/pc.um
    wl2 = 1/wn[::-1]/pc.um
    cia1 = op.Collision_Induced(f'{c_root}{cs_file}', wl=wl1)
    cia2 = op.Collision_Induced(f'{c_root}{cs_file}', wl=wl2)

    temperature = np.tile(1200.0, nlayers)
    cia_cs1 = cia1.calc_cross_section(temperature)
    cia_cs2 = cia2.calc_cross_section(temperature)

    np.testing.assert_allclose(cia_cs1, np.fliplr(cia_cs2))


def test_cia_both_wl_wn():
    # Ensure that order of wl array does not matter:
    cs_file = 'CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
    wl = 1/wn/pc.um
    error = 'Either provide wl or wn array for CIA, not both'
    with pytest.raises(ValueError, match=re.escape(error)):
        op.Collision_Induced(f'{c_root}{cs_file}', wl=wl, wn=wn)


def test_cia_cross_section():
    cs_file = 'CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
    cia = op.Collision_Induced(f'{c_root}{cs_file}', wn=wn)

    with np.load(f'{e_root}expected_cia_H2H2_opacity.npz') as d:
        expected_cs1 = d['expected_cs1']
        expected_cs2 = d['expected_cs2']

    temp1 = np.tile(1200.0, nlayers)
    cross_section1 = cia.calc_cross_section(temp1)
    np.testing.assert_allclose(cia.cross_section, expected_cs1)
    np.testing.assert_allclose(cia.cross_section, cross_section1)

    temp2 = np.tile(3050.0, nlayers)
    cross_section2 = cia.calc_cross_section(temp2)
    np.testing.assert_allclose(cia.cross_section, expected_cs2)
    np.testing.assert_allclose(cia.cross_section, cross_section2)


def test_cia_single_cross_section():
    cs_file = 'CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
    cia = op.Collision_Induced(f'{c_root}{cs_file}', wn=wn)

    with np.load(f'{e_root}expected_cia_H2H2_opacity.npz') as d:
        expected_cs3 = d['expected_cs3']

    temp3 = 1200.0
    cross_section = cia.calc_cross_section(temp3)
    np.testing.assert_allclose(cross_section, expected_cs3)
    # Single-layer runs don't update
    assert not hasattr(cia, 'cross_section')



def test_cia_extinction_coefficient():
    cs_file = 'CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
    cia = op.Collision_Induced(f'{c_root}{cs_file}', wn=wn)

    cia_indices = [species.index(mol) for mol in cia.species]
    densities = number_densities[:,cia_indices]

    with np.load(f'{e_root}expected_cia_H2H2_opacity.npz') as d:
        expected_cs1 = d['expected_cs1']
        expected_ec1 = d['expected_ec1']

    temp1 = np.tile(1200.0, nlayers)
    extinction = cia.calc_extinction_coefficient(temp1, densities)
    np.testing.assert_allclose(extinction, expected_ec1)
    # Extinction coeff call on atmosphere also updates cross section
    np.testing.assert_allclose(cia.cross_section, expected_cs1)


def test_cia_single_extinction_coefficient():
    cs_file = 'CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
    cia = op.Collision_Induced(f'{c_root}{cs_file}', wn=wn)

    cia_indices = [species.index(mol) for mol in cia.species]
    densities = number_densities[:,cia_indices]

    with np.load(f'{e_root}expected_cia_H2H2_opacity.npz') as d:
        expected_ec3 = d['expected_ec3']

    temp = 1200.0
    extinction = cia.calc_extinction_coefficient(temp, densities[3])
    np.testing.assert_allclose(extinction, expected_ec3)
    # Single-layer runs don't update cross_section
    assert not hasattr(cia, 'cross_section')


def test_cia_single_mismatch():
    cs_file = 'CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
    cia = op.Collision_Induced(f'{c_root}{cs_file}', wn=wn)

    densities = [5.18507926e+16]
    temp = 1200.0
    error = re.escape(
        "Incompatible dimensions, if temperature is scalar density "
        "must have self.nspec elements"
    )
    with pytest.raises(ValueError, match=error):
        cia.calc_extinction_coefficient(temp, densities)


def test_cia_array_mismatch():
    cs_file = 'CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
    cia = op.Collision_Induced(f'{c_root}{cs_file}', wn=wn)

    densities = [5.18507926e+16]
    temp = np.tile(1200.0, nlayers)
    error = re.escape(
        "Incompatible dimensions, density must be a 2D array of "
        "shape [6, 2], i.e., [ntemp, nspec]"
    )
    with pytest.raises(ValueError, match=error):
        cia.calc_extinction_coefficient(temp, densities)


@pytest.mark.parametrize('call', ['cross_section', 'extinction'])
def test_cia_temp_outbounds(call):
    cs_file = 'CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
    cia = op.Collision_Induced(f'{c_root}{cs_file}', wn=wn)

    cia_indices = [species.index(mol) for mol in cia.species]
    densities = number_densities[:,cia_indices]

    temp = np.tile(9000.0, nlayers)
    error = re.escape(
        "Invalid temperature, values must be in the 60.0-7000.0 K range"
    )
    with pytest.raises(ValueError, match=error):
        cia.calc_extinction_coefficient(temp, densities)

