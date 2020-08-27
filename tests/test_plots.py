# Copyright (c) 2016-2020 Patricio Cubillos.
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE).

import os

import numpy as np

import pyratbay as pb
import pyratbay.plots as pp

os.chdir(pb.constants.ROOT+'tests')


def posterior_pt():
    pyrat = pb.run('configs/mcmc_transmission_test.cfg', init=True)
    with np.load('MCMC_transmission_test.npz') as d:
        Z  = d['Z']
    itemp = np.arange(np.sum(pyrat.ret.pstep[pyrat.ret.itemp] > 0))
    tpars = pyrat.ret.params[pyrat.ret.itemp]

    posterior = Z[:,itemp]
    tmodel = pyrat.atm.tmodel
    pressure = pyrat.atm.press
    pp.pt(posterior, pyrat.atm.tmodel, tpars, itemp, pressure)


def abundance():
    # Template, I don't know how to truly automatize this since I'd
    # need to see the output plots to be OK, but a not-breaking-code
    # is the minimum testing I can guarantee at the moment.
    import pyratbay.atmosphere as pa
    import pyratbay.plots as pp

    nlayers = 51
    pressure = pa.pressure('1e-6 bar', '1e2 bar', nlayers)
    temperature = pa.temperature('isothermal', pressure,  params=1000.0)
    species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 TiO VO H2 H He Na K'.split()
    Q = pa.abundance(pressure, temperature, species, ncpu=3)
    ax = pp.abundance(Q, pressure, species, colors='default',
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())

    spec2 = np.copy(species)
    spec2[1] = 'CH3'
    spec2[2] = 'O'
    ax = pp.abundance(Q, pressure, spec2, colors='default',
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())

    # Custom colors:
    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(Q, pressure, species, colors=colors,
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())

    # Cycler dashes:
    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(Q, pressure, species, colors=colors,
        dashes=cycler(dashes=[(8,1), (3,1)]),
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())

    # Custom colors (different highlighting):
    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(Q, pressure, species, colors=colors,
        highlight='CH4 CO CO2 NH3 HCN H2 H He'.split())

    # Matplotlib colors:
    ax = pp.abundance(Q, pressure, species, colors=None,
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())

    # Custom colors, ncolors == nspecies:
    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(Q, pressure, species[:7], colors=colors,
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())

    # Custom colors/dashes, ncolors == nspecies:
    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(Q, pressure, species[:7], colors=colors,
        dashes=[(),(3,1),(),(),(),(),()],
        highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())

    # Custom colors, ncolors == nspecies, preserve coloring:
    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(Q, pressure, species[:7], colors=colors,
        highlight='CH4 CO CO2 NH3 HCN H2 H He'.split())

    colors = 'b r g y m k 0.5'.split()
    ax = pp.abundance(Q, pressure, species[:7], colors=colors,
        dashes=[(),(3,1),(),(),(),(),()],
        highlight='CH4 CO CO2 NH3 HCN H2 H He'.split())


