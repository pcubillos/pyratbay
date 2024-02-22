# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Collision_Induced',
]

import numpy as np
import mc3.utils as mu
from .. import io as io
from .. import constants as pc
from .. import tools as pt
from ..lib import _spline as sp


class Collision_Induced():
    """Collision Induced Absorption opacities"""
    def __init__(self, cia_file, *, wn=None, wl=None, log=None):
        """
        Parameters
        ----------
        cia_file: String
            A CIA cross section file.
        wn: 1D float array
            Wavenumber array (cm-1 units) where to sample the CIA
            (only one of wl or wn should be provided).
        wl: 1D float array
            Wavelength array (micron units) where to sample the CIA
            (only one of wl or wn should be provided).

        Examples
        --------
        >>> import pyratbay.spectrum as ps
        >>> import pyratbay.constants as pc
        >>> import pyratbay.opacity as op

        >>> wl = ps.constant_resolution_spectrum(0.61, 10.01, 15000.0)
        >>> cs_file = f'{pc.ROOT}/pyratbay/data/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
        >>> cia = op.Collision_Induced(cs_file, wl=wl)
        """
        if log is None:
            log = mu.Log()

        self.cia_file = cia_file
        absorption, species, temps, tab_wn = io.read_cs(cia_file)

        self.species = species
        self.nspec = len(self.species)
        self.name = 'CIA ' + '-'.join(self.species)

        # Temperature must be sorted in increasing order
        t_sort = np.argsort(temps)
        absorption = absorption[t_sort]

        self.temps = temps[t_sort]
        self.ntemp = len(self.temps)
        self.tmin = np.amin(self.temps)
        self.tmax = np.amax(self.temps)

        if wl is not None and wn is not None:
            raise ValueError('Either provide wl or wn array for CIA, not both')
        if wl is not None:
            wn = 1.0 / (wl*pc.um)

        if wn is None:
            self.wn = tab_wn
            self.tab_cross_section = absorption
            self.nwave = len(self.wn)
        else:
            self.wn = wn
            self.nwave = len(self.wn)
            # Wavenumber range check:
            if np.amin(tab_wn) > np.amin(wn) or np.amax(tab_wn) < np.amax(wn):
                log.warning(
                    f"The tabulated wavenumber range [{np.amin(tab_wn):.2f}, "
                    f"{np.amax(tab_wn):.2f}] cm-1 "
                    "does not cover the whole requested wavenumber range: "
                    f"[{np.amin(wn):.2f}, {np.amax(wn):.2f}] cm-1 "
                    f"for cross-section file: '{cia_file}'"
                )
            # Make sure wavenumbers are in increasing order:
            if self.wn[1] < self.wn[0]:
                sorted_self_wn = self.wn[::-1]
            else:
                sorted_self_wn = self.wn
            if tab_wn[1] < tab_wn[0]:
                sorted_tab_wn = tab_wn[::-1]
            else:
                sorted_tab_wn = tab_wn

            # Interpolate along wavenumber axis:
            cross_section = np.zeros((self.ntemp, self.nwave), np.double)
            extrapolate_value = 0.0
            for j in range(self.ntemp):
                ddev = sp.second_deriv(absorption[j], sorted_tab_wn)
                cross_section[j] = sp.splinterp_1D(
                    absorption[j], sorted_tab_wn, ddev,
                    sorted_self_wn,
                    extrapolate_value,
                )
            if self.wn[1] < self.wn[0]:
                cross_section = np.fliplr(cross_section)
            self.tab_cross_section = cross_section

        # per amagat**N to per (molec cm-3)**N
        self.tab_cross_section /= pc.amagat**self.nspec

        # Wavenumber indices where data exists:
        wn_good = (self.wn >= np.amin(tab_wn)) & (self.wn <= np.amax(tab_wn))
        self._wn_lo_idx = np.where(wn_good)[0][0]
        self._wn_hi_idx = np.where(wn_good)[0][-1] + 1
        # Delta-cross section over delta-temperature
        self._dcs_dt = (
            np.diff(self.tab_cross_section, axis=0) /
            np.expand_dims(np.ediff1d(self.temps), axis=1)
        )


    def calc_cross_section(self, temperature):
        """
        Calculate cross-section spectra (cm-1) over temperature
        profiles by interpolating from tabulated values.

        Parameters
        ----------
        temperature: float or 1D float array
            Temperature array (Kelvin) at which to interpolate the
            cross section.

        Returns
        -------
        cross_section: 1D or 2D float array
            Cross section spectra (cm-1 /(molec cm-3)**N)
            where N is the number of species.
            If temperature is scalar, cross_section is 1D of length self.nwave.
            Otherwise, cross_section is a 2D array of shape [nlayers,nwave]
        """
        # Enforce np.array()
        scalar_temp = np.isscalar(temperature)
        if scalar_temp:
            temperature = np.array([temperature])

        if np.any(temperature < self.tmin) or np.any(temperature > self.tmax):
            raise ValueError(
                'Invalid temperature, values must be in the '
                f'{self.tmin:.1f}-{self.tmax:.1f} K range'
            )

        nlayers = len(temperature)
        cross_section = np.zeros((nlayers, self.nwave))
        sp.lin_interp_2D(
            self.tab_cross_section, self.temps, self._dcs_dt,
            np.array(temperature), cross_section,
            self._wn_lo_idx, self._wn_hi_idx,
        )

        if scalar_temp:
            return cross_section[0]

        self.cross_section = cross_section
        return cross_section


    def calc_extinction_coefficient(self, temperature, density):
        """
        Calculate extinction-coefficient spectra (cm-1) over temperature
        and density profiles by interpolating in temperature from
        tabulated values.

        Parameters
        ----------
        temperature: 1D float array
            Temperature array (Kelvin) at which to interpolate the
            cross section.
        density: 2D float array
            Number density array (molecules cm-3) of self.species
            over the temperature array.
            If temperature is 1D array, must be of shape [self.nspec, ntemp]
            If temperature is scalar, must be of length self.nspec.

        Returns
        -------
        extinction_coefficient: 1D or 2D float array
            Extinction coefficient spectra (cm-1)
            If temperature is scalar, output is a 1D array of length nwave.
            Otherwise, outout is a 2D array of shape [nlayers,nwave]
        """
        # Dimentional checks:
        is_scalar = np.isscalar(temperature)
        dens_shape = np.shape(density)
        if is_scalar:
            if np.size(density) != self.nspec:
                raise ValueError(
                    'Incompatible dimensions, if temperature is scalar '
                    'density must have self.nspec elements'
                )
        elif dens_shape != (len(temperature), self.nspec):
            ntemp = len(temperature)
            nspec = self.nspec
            raise ValueError(
                'Incompatible dimensions, density must be a 2D array '
                f'of shape [{ntemp}, {nspec}], i.e., [ntemp, nspec]'
            )

        cross_section = self.calc_cross_section(temperature)

        if is_scalar:
            return cross_section * np.prod(density)

        extinction_coefficient = (
            cross_section *
            np.prod(density, axis=1, keepdims=True)
        )
        return extinction_coefficient


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write(f"CIA file name (cia_file): '{self.cia_file}'")
        fw.write(f'CIA species (species): {self.species}')
        fw.write(f'Number of temperature samples (ntemp): {self.ntemp}')
        fw.write(f'Number of wavenumber samples (nwave): {self.nwave}')
        fw.write(f'Temperature array (temps, K):\n{self.temps}')
        fw.write(
            'Wavenumber array (wn, cm-1):\n{}',
            self.wn,
            fmt={'float': '{:.3f}'.format},
        )
        fw.write(
            'Wavelength array (um):\n{}',
            1/self.wn/pc.um,
            fmt={'float': '{:.5f}'.format},
        )
        cm_pow = 3*self.nspec - 1
        mol_pow = -self.nspec
        fw.write(
            'Tabulated cross section (tab_cross_section, cm{} '
            'molec{}):\n{}', cm_pow, mol_pow, self.tab_cross_section,
            fmt={'float': '{:.2e}'.format})
        fw.write(
            'Minimum and maximum temperatures (tmin, tmax) in K: '
            f'[{self.tmin:.1f}, {self.tmax:.1f}]'
        )
        return fw.text

