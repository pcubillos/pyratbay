# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'CCSgray',
    'Deck',
]

import numpy as np
import scipy.interpolate as si

from ... import constants as pc
from ... import tools     as pt


class CCSgray():
    """
    Constant cross-section gray cloud model.
    """
    def __init__(self, pressure, wn):
        self.name = 'ccsgray' # Model name is lowercased class name
        self.pressure = pressure / pc.bar
        self.wn = wn
        self.nwave = len(wn)
        self.nlayers = len(pressure)
        # log10 of cross-section scale factor, top, and bottom pressure (bar)
        self.pars = [0.0, -4.0, 2.0]
        self.npars = len(self.pars)  # Number of model fitting parameters
        self.ec = np.zeros((self.nlayers, self.nwave))
        self.mol = 'H2'            # Species causing the extinction
        # Fitting-parameter names (plain text and figure labels):
        self.pnames = ['log_k_gray', 'log_p_top', 'log_p_bot']
        self.texnames = [
            r'$\log_{10}(f_{\rm gray})$',
            r'$\log_{10}(p_{\rm top})\ ({\rm bar})$',
            r'$\log_{10}(p_{\rm bot})\ ({\rm bar})$',
        ]
        self.s0 = 5.31e-27  # Default coss-section (cm-2 molec-1)

    def calc_cross_section(self):
        """
        Calculate a uniform gray-cloud cross section in cm2 molec-1:
           cross section = s0 * 10**pars[0],
        between layers with pressure 10**pars[1] -- 10**pars[2] bar
        (top and bottom layers, respectively).
        s0 is the H2 Rayleigh cross section at 0.35 um.

        Parameters
        ----------
        wn:  1D float ndarray
           Wavenumber array in cm-1.
        """
        # Get indices for cloud layer boundaries:
        p_top = 10**self.pars[1]*pc.bar
        p_bottom = 10**self.pars[1]*pc.bar
        p_mask = (self.pressure >= p_bottom) & (self.pressure <= p_top)

        # Gray opacity cross section in cm2 molec-1
        self.ec[p_mask,:] = 10**self.pars[0] * self.s0


    def calc_extinction_coefficient(self, density, pars=None, layer=None):
        if pars is not None:
            self.pars[:] = pars
        # Densities in molecules cm-3:
        # Cross section (in cm2 molecule-1):
        cs = self.calc_cross_section()
        # Cloud absorption (cm-1):
        if layer is not None:
            return cs[layer] * density[layer]
        return cs * np.expand_dims(density, axis=1)


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write("Model name (name): '{}'", self.name)
        fw.write('Model species (mol): {}', self.mol)
        fw.write('Number of model parameters (npars): {}', self.npars)
        fw.write('Parameter name     Value\n'
                 '  (pnames)         (pars)\n')
        for pname, param in zip(self.pnames, self.pars):
            fw.write('  {:15s}  {: .3e}', pname, param)
        fw.write('Extinction-coefficient (ec, cm2 molec-1):\n{}', self.ec,
            fmt={'float':'{: .3e}'.format}, edge=3)
        return fw.text


class Deck():
    """
    Instantly opaque gray cloud deck at given pressure.
    """
    def __init__(self, pressure, wn):
        self.name = 'deck'
        self.pressure = pressure / pc.bar
        self.wn = wn
        self.nwave = len(wn)
        self.nlayers = len(pressure)
        self.pars = [-1.0]          # log10(Pressure[bar]) of cloud top
        self.npars = len(self.pars)  # Number of model fitting parameters
        self.ec = np.zeros((self.nlayers, self.nwave))
        # Fitting-parameter names (plain text and figure labels):
        self.pnames = ['log_p_cl']
        self.texnames = [r'$\log\ p_{\rm cl}$']
        self.itop = None
        self.rsurf = 0.0
        self.tsurf = 0.0

    def calc_extinction_coefficient(
        self, radius, temperature, pars=None, layer=None,
    ):
        """
        Calculate gray-cloud deck that's transparent above ptop,
        and becomes instantly opaque at ptop, with
        ptop (bar) = 10**pars[0].

        Parameters
        ----------
        radius: 1D float ndarray
            Atmospheric radius profile (in cm).
        temp: 1D float ndarray
            Atmospheric temperature profile (in Kelvin degree).
        """
        if pars is not None:
            self.pars[:] = pars

        ptop = 10**self.pars[0]
        # Index of layer directly below cloud top:
        if ptop >= self.pressure[-1]:  # Atmosphere boundary cases
            self.itop = self.nlayers-1
        elif ptop < self.pressure[0]:
            self.itop = 1
        else:
            self.itop = np.where(self.pressure>=ptop)[0][0]

        if layer is not None:
            # Just a boolean indicating whether the cloud top is above layer:
            return np.zeros(self.nwave) + int(layer > self.itop)

        # Radius and temperature at the cloud top:
        self.tsurf = float(si.interp1d(self.pressure, temperature)(ptop))
        self.rsurf = float(si.interp1d(self.pressure, radius)(ptop))
        return self.ec


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write("Model name (name): '{}'", self.name)
        fw.write('Number of model parameters (npars): {}', self.npars)
        fw.write('Parameter name     Value\n'
                 '  (pnames)         (pars)\n')
        for pname, param in zip(self.pnames, self.pars):
            fw.write('  {:15s}  {: .3e}', pname, param)
        fw.write('Index of atmospheric layer at or directly below cloud '
                 f'top: {self.itop}')
        fw.write('Cloud-top pressure: {:.4e} bar', 10**self.pars[0])
        fw.write('Cloud-top altitude: {:.2f} km', self.rsurf/pc.km)
        fw.write('Cloud-top temperature: {:.2f} K', self.tsurf)
        return fw.text

