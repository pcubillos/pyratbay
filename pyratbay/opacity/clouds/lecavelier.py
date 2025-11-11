# Copyright (c) 2021-2025 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Lecavelier',
]

import numpy as np

from ... import constants as pc
from ... import tools as pt


class Lecavelier():
    """
    Rayleigh-like scattering model from Lecavelier des Etangs et al. (2008).
    AA, 485, 865.

    Examples
    --------
    >>> import pyratbay.opacity as op
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.spectrum as ps
    >>> import numpy as np

    >>> wl = ps.constant_resolution_spectrum(0.2, 6.0, resolution=1000)
    >>> nlayers = 81
    >>> pressure = pa.pressure('1e-8 bar', '1e2 bar', nlayers)
    >>> lecavelier = op.clouds.Lecavelier(pressure, wl=wl)

    >>> temperature = np.tile(1800.0, nlayers)
    >>> extinction = lecavelier.calc_extinction_coefficient(temperature)
    """
    def __init__(self, pressure, wl=None, wn=None):
        self.name = 'lecavelier'
        self.pressure = pressure

        if wl is not None and wn is not None:
            raise ValueError('Either provide wl or wn array, not both')
        if wl is None and wn is None:
            msg = 'Missing spectral array input, either wl (micron) or wn (cm-1)'
            raise ValueError(msg)
        if wl is not None:
            wn = 1.0 / (wl*pc.um)
        self.wn = wn

        # Model parameters: cross-section scale factor and power-law exponent
        self.pars = [0.0, -4.0]
        self.npars = len(self.pars)
        self.pnames = ['log_k_ray', 'alpha_ray']
        self.texnames = [r'$\log\ \kappa_{\rm ray}$', r'$\alpha_{\rm ray}$']
        self.cross_section = None
        self.s0 = 5.31e-27  # Cross section (cm2 molec-1) at l0
        self.l0 = 3.5e-5    # Nominal wavelength (cm)
        self.calc_cross_section()
        self.ec = np.zeros((len(self.pressure), len(self.wn)))


    def calc_cross_section(self, pars=None):
        """
        Calculate the Rayleigh cross section in cm2 molec-1:
            cross section = k_ray * s0 * (lambda/l0)**alpha_ray,
        parameterized as params = [log10(k_ray), alpha_ray].

        Parameters
        ----------
        pars: 1D iterable
            If not None, update the model parameters with given values.

        Returns
        -------
        cross_section: 1D float array
            Cross section (cm2 molecule-1) as function of wavelength
        """
        if pars is not None:
            self.pars[:] = pars

        # Rayleigh-like cross section in cm2 molec-1:
        self.cross_section = \
            10.0**self.pars[0] * self.s0 * (self.wn*self.l0)**(-self.pars[1])
        return self.cross_section


    def calc_extinction_coefficient(self, temperature, pars=None, layer=None):
        """
        Calculate extinction-coefficient (cm-1 units) over wavelength
        and pressure arrays.
        The nominal density profile of the absorber is assumed as
        density = pressure / (k*temperature)

        temperature: 1D float array
            Temperature profile (K)
        pars: 1D iterable
            If not None, update the model parameters with given values
            and re-calculate the cross sections.
        layer: Integer
            If not None, compute the extinction coefficient only
            at the given index in density array.

        Returns
        -------
        extinction_coefficient: 2D float array
            The Rayleigh extinction coefficient (cm-1 units).
        """
        if pars is not None:
            self.calc_cross_section(pars)

        # Densities in molecules cm-3:
        density = self.pressure*pc.bar / temperature / pc.k

        # Extinction coefficient (cm-1):
        if layer is not None:
            return self.cross_section * density[layer]

        self.ec = self.cross_section * np.expand_dims(density, axis=1)
        return self.ec


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write("Model name (name): '{}'", self.name)
        p_top, p_bot = self.pressure[0], self.pressure[-1]
        fw.write('Pressure array (pressure, bar): {} ... {}', p_top, p_bot)
        fw.write('Number of model parameters (npars): {}', self.npars)
        fw.write('Parameter name     Value\n'
                 '  (pnames)         (pars)\n')
        for pname, param in zip(self.pnames, self.pars):
            fw.write('  {:15s}  {: .3e}', pname, param)
        fw.write(
            'Wavenumber (wn, cm-1):\n   {}',
            self.wn,
            fmt={'float':'{:.2f}'.format},
        )
        fw.write(
            'Cross section (cross_section, cm2 molec-1):\n   {}',
            self.cross_section,
            fmt={'float':'{: .3e}'.format}, edge=3)
        return fw.text
