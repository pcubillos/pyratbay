# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os

import pytest
import numpy as np
from cycler import cycler
import matplotlib.pyplot as plt

import pyratbay as pb
import pyratbay.plots as pp
import pyratbay.atmosphere as pa


os.chdir(pb.constants.ROOT+'tests')

nlayers = 51
pressure = pa.pressure('1e-6 bar', '1e2 bar', nlayers)
temperature = pa.temperature('isothermal', pressure,  params=1000.0)
species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 TiO VO H2 H He Na K'.split()
vmr = pa.chemistry('tea', pressure, temperature, species).vmr

# Templates, I don't know how to truly automatize this since I'd
# need to see the output plots to be OK, but a not-breaking-code
# is the minimum testing I can guarantee at the moment.

@pytest.mark.skip
def test_spectrum():
    import mc3
    d = np.load('multinest_MCMC_WASP18b_equil_group02_dilut_yes_posteriors_info.npz')
    spectra = d['depth_posterior']
    wl = d['wl']
    ibounds = [3,1,2,4]
    bounds = spectra[ibounds] - spectra[0]

    bands_wl0 = np.linspace(1.0, 2.75, 12)
    data = np.random.normal(0.0, 3e-5, 12)
    uncert = np.tile(3e-5, 12)
    bands_flux = np.random.normal(0.0, 1e-5, 12)

    theme = mc3.plots.Theme('darkorange')
    theme.dark_color = 'maroon'
    theme.light_color = 'gold'
    pp.spectrum(
        spectra[0]*0, wl, 'emission', bounds=bounds, theme=theme, gaussbin=50,
        data=data, uncert=uncert, bands_wl0=bands_wl0, bands_flux=bands_flux,
    )


@pytest.mark.parametrize("ndim", ['1d', '2d'])
def test_temperature_minimal(ndim):
    tmodel = pa.tmodels.Guillot(pressure)
    profile = tmodel([-4.0, -1.0, 0.0, 0.0, 1000.0, 0.0])
    if ndim == '2d':
        profile = [profile]
    ax = pp.temperature(pressure, profiles=profile)


def test_temperature_labels():
    tmodel = pa.tmodels.Guillot(pressure)
    profile = tmodel([-4.0, -1.0, 0.0, 0.0, 1000.0, 0.0])
    ax = pp.temperature(pressure, profiles=profile, labels=['model'])


def test_temperature_custom_colors():
    tmodel = pa.tmodels.Guillot(pressure)
    profile = tmodel([-4.0, -1.0, 0.0, 0.0, 1000.0, 0.0])
    ax = pp.temperature(pressure, profiles=profile, colors=['r'])


def test_temperature_bounds_1sigma():
    tmodel = pa.tmodels.Guillot(pressure)
    profile = tmodel([-4.0, -1.0, 0.0, 0.0, 1000.0, 0.0])
    bounds = [
        tmodel([-4.0, -0.9, 0.0, 0.0, 1000.0, 0.0]),
        tmodel([-4.0, -1.1, 0.0, 0.0, 1000.0, 0.0]),
    ]
    ax = pp.temperature(pressure, profiles=profile, bounds=bounds)


def test_temperature_bounds_2sigma():
    tmodel = pa.tmodels.Guillot(pressure)
    profile = tmodel([-4.0, -1.0, 0.0, 0.0, 1000.0, 0.0])
    bounds = [
        tmodel([-4.0, -0.9, 0.0, 0.0, 1000.0, 0.0]),
        tmodel([-4.0, -1.1, 0.0, 0.0, 1000.0, 0.0]),
        tmodel([-4.1, -0.88, 0.0, 0.0, 990.0, 0.0]),
        tmodel([-3.9, -1.11, 0.0, 0.0, 1010.0, 0.0]),
    ]
    ax = pp.temperature(pressure, profiles=profile, bounds=bounds)


def test_temperature_ax():
    tmodel = pa.tmodels.Guillot(pressure)
    profile = tmodel([-4.0, -1.0, 0.0, 0.0, 1000.0, 0.0])
    plt.figure(0)
    plt.clf()
    ax = plt.subplot(221)
    ax = pp.temperature(pressure, profiles=profile, ax=ax)


def test_abundance_default_molecs():
    ax = pp.abundance(vmr, pressure, species, colors='default',
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())


def test_abundance_defaults_extra_molecs():
    spec2 = np.copy(species)
    spec2[1] = 'CH3'
    spec2[2] = 'O'
    ax = pp.abundance(vmr, pressure, spec2, colors='default',
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())


def test_abundance_custom_colors():
    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(vmr, pressure, species, colors=colors,
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())


def test_abundance_custom_dashes_and_colors():
    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(
        vmr, pressure, species, colors=colors,
        dashes=cycler(dashes=[(8,1), (3,1)]),
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split(),
    )


def test_abundance_custom_colors_a():
    # Different highlight changes color assignments:
    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(vmr, pressure, species, colors=colors,
        highlight='CH4 CO CO2 NH3 HCN H2 H He'.split())


def test_abundance_matplotlib_default_colors():
    # Matplotlib colors:
    ax = pp.abundance(vmr, pressure, species, colors=None,
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())


def test_abundance_one_to_one_colors():
    # Custom colors, ncolors == nspecies:
    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(vmr, pressure, species[:7], colors=colors,
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())


def test_abundance_one_to_one_dashes_and_colors():
    # Custom colors/dashes, ncolors == nspecies:
    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(vmr, pressure, species[:7], colors=colors,
        dashes=[(),(3,1),(),(),(),(),()],
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())


def test_abundance_one_to_one_colors_alternative_highlight():
    # ncolors >= nspecies, different highlight preserves coloring:
    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(vmr, pressure, species[:7], colors=colors,
        highlight='CH4 CO CO2 NH3 HCN H2 H He'.split())


def test_abundance_one_to_one_dashes_and_colors_alternative_highlight():
    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(vmr, pressure, species[:7], colors=colors,
        dashes=[(),(3,1),(),(),(),(),()],
        highlight='CH4 CO CO2 NH3 HCN H2 H He'.split())


