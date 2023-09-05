# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Dalgarno',
    'Lecavelier',
]

import numpy as np

from ... import tools as pt


class Dalgarno():
    """
    Rayleigh-scattering models from Dalgarno (1962), Kurucz (1970), and
    Dalgarno & Williams (1962).
    """
    def __init__(self, wn, species):
        """
        Parameters
        ----------
        wn: 1D float ndarray
           Wavenumber in cm-1.
        species: String
           The species, which can be H, He, or H2.
        """
        self.name = f'dalgarno_{species}'
        self.species = species
        self.set_wn(wn)
        self.npars = 0
        self.pars = []
        self.pnames = []
        self.texnames = []

    def set_wn(self, wn):
        """
        When wn is updated the cross-sections must be re-calculated.
        """
        self.wn = wn
        if self.species == 'H':
            self.coef = np.array([5.799e-45, 1.422e-54, 2.784e-64])
            self._calc_H_cross_section()
        elif self.species == 'He':
            self.coef = np.array([5.484e-46, 2.440e-11, 5.940e-42, 2.900e-11])
            self._calc_He_cross_section()
        elif self.species == 'H2':
            self.coef = np.array([8.140e-45, 1.280e-54, 1.610e-64])
            self._calc_H_cross_section()

    def _calc_H_cross_section(self):
        """
        Rayleigh H/H2 cross section in cm2 molec-1 units.
        Sections 5.4 and 5.13, Kurucz (1970).
        """
        self.cross_section = (
            self.coef[0]*self.wn**4.0 +
            self.coef[1]*self.wn**6.0 +
            self.coef[2]*self.wn**8.0
        )

    def _calc_He_cross_section(self):
        """
        Rayleigh He cross section in cm2 molec-1 units.
        Section 5.8, Kurucz (1970)
        """
        self.cross_section = self.coef[0]*self.wn**4 * (
            1.0 +
            self.coef[1]*self.wn**2 +
            self.coef[2]*self.wn**4/(1 - self.coef[3]*self.wn**2)
        )**2.0

    def calc_extinction_coefficient(self, density, layer=None):
        """
        Calculate extinction-coefficient (cm-1 units) over wavelength
        and layers grid.

        Parameters
        ----------
        density: 1D float array
            Number density of this species over an atmospheric profile
            (molecules cm-3 units)
        layer: Integer
            If not None, compute the extinction coefficient only
            at the given index in density array.

        Returns
        -------
        extinction_coefficient: 2D float array
            The Rayleigh extinction for this species for the given
            atmosphere (cm-1 units).
        """
        # Extinction coefficient at single layer (cm-1):
        if layer is not None:
            return self.cross_section * density[layer]

        # Extinction coefficient profile (cm-1):
        return self.cross_section * np.expand_dims(density, axis=1)


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write(f"Model name (name): '{self.name}'")
        fw.write(f'Model species (species): {self.species}')
        fw.write(f'Number of model parameters (npars): {self.npars}')
        fw.write(
            'Wavenumber (wn, cm-1):\n   {}',
            self.wn,
            fmt={'float':'{:.2f}'.format},
        )
        fw.write(
            'Cross section (cross_section, cm2 molec-1):\n   {}',
            self.cross_section,
            fmt={'float':'{:.3e}'.format},
        )
        return fw.text


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
