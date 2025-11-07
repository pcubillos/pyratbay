# Copyright (c) 2021-2025 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Lecavelier',
]

import numpy as np

from ... import tools as pt


class Lecavelier():
    """
    Rayleigh-scattering model from Lecavelier des Etangs et al. (2008).
    AA, 485, 865.
    """
    def __init__(self, wn, species='H2'):
        self.name = 'lecavelier'
        self.species = species  # Species causing the extinction
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

    def calc_cross_section(self, pars=None):
        """
        Calculate the Rayleigh cross section in cm2 molec-1:
            cross section = k_ray * s0 * (lambda/l0)**alpha_ray,
        parameterized as params = [log10(k_ray), alpha_ray].

        Parameters
        ----------
        pars: 1D iterable
            If not None, update the model parameters with given values.
        """
        if pars is not None:
            self.pars[:] = pars
        # Rayleigh cross section opacity in cm2 molec-1:
        self.cross_section = \
            10.0**self.pars[0] * self.s0 * (self.wn*self.l0)**(-self.pars[1])


    def calc_extinction_coefficient(self, density, pars=None, layer=None):
        """
        Calculate extinction-coefficient (cm-1 units) over wavelength
        and layers grid.

        density: 1D float array
            Number density of this species over an atmospheric profile
            (molecules cm-3 units)
        pars: 1D iterable
            If not None, update the model parameters with given values
            and re-calculate the cross sections.
        layer: Integer
            If not None, compute the extinction coefficient only
            at the given index in density array.

        Returns
        -------
        extinction_coefficient: 2D float array
            The Rayleigh extinction for this species for the given
            atmosphere (cm-1 units).
        """
        if pars is not None:
            self.calc_cross_section(pars)

        if layer is not None:
            return self.cross_section * density[layer]

        # Extinction coefficient profile (cm-1):
        return self.cross_section * np.expand_dims(density, axis=1)


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write("Model name (name): '{}'", self.name)
        fw.write('Model species (species): {}', self.species)
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
