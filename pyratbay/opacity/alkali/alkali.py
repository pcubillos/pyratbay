# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'SodiumVdW',
    'PotassiumVdW',
]

import numpy as np

from ... import constants as pc
from ... import tools as pt
from ... import io as io
from .. import broadening

from ...lib import _alkali


class VanderWaals():
    """
    Base class for Van der Waals plus statistical theory model
    Burrows et al. (2000), ApJ, 531, 438.
    """
    def __init__(self, pressure, wn, cutoff):
        self.pressure = pressure
        self.wn = wn
        self.nwave = len(wn)
        self.nlayers = len(pressure)
        # Hard cutoff from line center (cm-1)
        self.cutoff = cutoff
        self.cross_section = None
        self.voigt = broadening.Voigt()

        names, masses, _ = io.read_molecs(pc.ROOT+'pyratbay/data/molecules.dat')
        self.mass = masses[list(names).index(self.species)]

        # Spectral sampling rate at alkali wn0:
        self._dwave = np.zeros(self.nlines)
        for i,wn0 in enumerate(self.wn0):
            i_wn = np.argmin(np.abs(wn-wn0))
            is_last = i_wn == self.nwave-1
            is_over = i_wn > 0 and wn[i_wn] > wn0
            if is_over or is_last:
                i_wn -= 1
            self._dwave[i] = np.abs(wn[i_wn+1]-wn[i_wn])


    def voigt_det(self, temperature):
        """
        Calculate Voigt profile value at the detuning wavelength
        from the center of the lines (at each layer).

        Parameters
        ----------
        temperature: 1D float array
            A temperature profile (kelvin).

        Returns
        -------
        voigt_det: 2D float array
            Voigt-profile values at the detuning wavelength
            of shape [nlayers,nlines].
        """
        # Detuning frequency (cm-1):
        dsigma = self.detuning * (temperature/500.0)**0.6
        # Lorentz half width (cm-1):
        lor = (
            self.lpar * (temperature/2000.0)**(-0.7) *
            self.pressure * pc.bar/pc.atm
        )
        voigt_det = np.zeros((self.nlayers,self.nlines))
        for j in range(self.nlines):
            wn0 = self.wn0[j]
            # Doppler half width (cm-1):
            dop = np.sqrt(2*pc.k*temperature / (self.mass*pc.amu)) * wn0 / pc.c
            self.voigt.x0 = wn0
            for i in range(self.nlayers):
                self.voigt.hwhm_L = lor[i]
                self.voigt.hwhm_G = dop[i]
                # EC at the detuning boundary:
                voigt_det[i,j] = self.voigt(wn0+dsigma[i])
        return voigt_det


    def calc_cross_section(self, temperature, layer=None):
        """
        Calculate cross section (cm2 molecule-1) at given temperatures
        Saves value into self.cross_section.

        Parameters
        ----------
        temperature: 1D float array
            Temperature profile in Kelvin (must match self.pressure array)
        layer: Interger
            If not None, index at which to calculate the cross section.

        Returns
        -------
        cross_section: 2D/1D float array
            If layer is None, cross-section spectra (cm2 molec-1) at each
            pressure-temperature point (also sets self.cross_section).
            If layer is not None, cross-section spectrum at a single layer.
        """
        voigt_det = self.voigt_det(temperature)
        if layer is None:
            nlayers = self.nlayers
            pressure = self.pressure*pc.bar
        else:
            nlayers = 1
            pressure = self.pressure[layer:layer+1]*pc.bar
            temperature = temperature[layer:layer+1]
            voigt_det = voigt_det[layer:layer+1]

        i_nearest_wn0 = np.argmin(
            np.abs(np.expand_dims(self.wn0,1)-self.wn),
            axis=1,
        )
        cross_section = np.zeros((nlayers, self.nwave), np.double)
        _alkali.alkali_cross_section(
            pressure, self.wn, temperature, voigt_det,
            cross_section,
            self.detuning, self.mass, self.lpar, self.Z, self.cutoff,
            np.array(self.wn0), np.array(self.gf), np.array(self._dwave),
            i_nearest_wn0,
        )
        if layer is None:
            self.cross_section = cross_section
        else:
            cross_section = cross_section[0]
        return cross_section


    def _calc_cross_section(self, temperature):
        """
        Calculate cross section (cm2 molecule-1)

        You are not supposed to be here. Use calc_cross_section() instead.
        """
        self.cross_section = np.zeros((self.nlayers, self.nwave), np.double)

        # Detuning frequency (cm-1):
        dsigma = self.detuning * (temperature/500.0)**0.6
        # Doppler half width (cm-1):
        doppler = (
            np.sqrt(2*pc.k*temperature / (self.mass*pc.amu))
            * np.expand_dims(self.wn0, axis=1) / pc.c
        )
        # Lorentz half width (cm-1):
        lor = (
            self.lpar * (temperature/2000.0)**(-0.7) *
            self.pressure * pc.bar/pc.atm
        )
        # Calculate cross section:
        for wn0, gf, dwave, dop in zip(self.wn0, self.gf, self._dwave, doppler):
            iwave = np.abs(self.wn-wn0) < self.cutoff
            wn = self.wn[iwave]
            cs = np.zeros((self.nlayers, len(wn)), np.double)
            # Update Voigt model:
            self.voigt.x0 = wn0
            fwidth = 2*(0.5346*lor + np.sqrt(0.2166*lor**2 + dop**2))
            for j in range(self.nlayers):
                # EC at the detuning boundary:
                self.voigt.hwhm_L = lor[j]
                self.voigt.hwhm_G = dop[j]
                edet = self.voigt(wn0+dsigma[j])
                # Extinction outside the detuning region (power law):
                cs[j] += edet * (np.abs(wn-wn0)/dsigma[j])**-1.5
                # Profile ranges:
                det = np.abs(wn - wn0) < dsigma[j]
                wlo = pt.ifirst(det, default_ret=-1)
                whi = pt.ilast( det, default_ret=-2) + 1
                wndet = wn[wlo:whi]
                # Extinction in the detuning region (Voigt profile):
                profile = self.voigt(wndet)
                # Correction for undersampled line:
                if whi > wlo and fwidth[j] < 2.0*dwave:
                    i0 = np.argmin(np.abs(wn0-wndet))
                    profile[i0] = 0.0
                    profile[i0] = 1.0 - np.trapezoid(profile, wndet)
                cs[j, wlo:whi] = profile
            # Add up contribution (include exponential cutoff):
            exp = np.exp(np.outer(1/temperature, -pc.C2*np.abs(wn-wn0)))
            self.cross_section[:,iwave] += pc.C3 * cs * gf/self.Z * exp
            # Note this equation neglects the exp(-Elow/T)*(1-exp(-wn0/T))
            # terms because they are approximately 1.0 at T <= 4000K


    def calc_extinction_coefficient(self, temperature, density, layer=None):
        """
        Calculate extinction coefficient (cm-1) for given temperature
        and number-density profiles.

        Parameters
        ----------
        temperature: 1D float array
            Temperature profile in Kelvin (must match self.pressure array)
        density: 1D float array
            Number density profile (molecules per cm3) of self.species.
            (must match self.pressure array)
        layer: Interger
            If not None, calculate the extinction coefficient at a
            single layer pointed by this index. In this case the output
            is a 1D array.

        Returns
        -------
        extinction_coefficient: 2D float array
            Alkali extinction coefficient spectrum (units of cm-1)
            of shape [nlayers, nwave].
        """
        cross_section = self.calc_cross_section(temperature, layer)
        if layer is not None:
            return cross_section * density[layer]
        return cross_section * np.expand_dims(density, axis=1)


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write("Model name (name): '{}'", self.name)
        fw.write('Model species (species): {}', self.species)
        fw.write('Species mass (mass, amu): {}', self.mass)
        fw.write(
            'Profile hard cutoff from line center (cutoff, cm-1): {}',
            self.cutoff,
        )
        fw.write('Detuning parameter (detuning): {}', self.detuning)
        fw.write('Lorentz-width parameter (lpar): {}', self.lpar)
        fw.write('Partition function (Z): {}', self.Z)
        fw.write(
            'Wavenumber  Wavelength          gf   Lower-state energy\n'
            '      cm-1          um               cm-1\n'
            '      (wn)                    (gf)   (elow)',
        )
        for wn, gf, elow in zip(self.wn0, self.gf, self.elow):
            fw.write(
                '  {:8.2f}  {:10.6f}   {:.3e}   {:.3e}',
                wn, 1.0/(wn*pc.um), gf, elow,
            )
        fw.write(
            'Wavenumber (wn, cm-1):\n   {}',
            self.wn,
            fmt={'float':'{:.2f}'.format},
        )
        fw.write(
            'Pressure (pressure, bar):\n{}',
            self.pressure,
            fmt={'float':'{:.3e}'.format},
        )
        fw.write(
            'Cross section (cross_section, cm2 molecule-1):\n{}',
            self.cross_section,
            fmt={'float': '{:.3e}'.format}, edge=2,
        )
        return fw.text


class SodiumVdW(VanderWaals):
    """
    Sodium Van der Waals model (Burrows et al. 2000, ApJ, 531).

    Examples
    --------
    >>> import pyratbay.opacity as op
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.spectrum as ps

    >>> pressure = pa.pressure('1e-8 bar', '1e2 bar', nlayers=81)
    >>> wl = ps.constant_resolution_spectrum(0.5, 1.0, resolution=15000.0)
    >>> sodium = op.alkali.SodiumVdW(pressure, 1e4/wl)

    >>> temperature = np.tile(1500.0, 81)
    >>> vmr = np.tile(1.6e-6, 81)
    >>> density = pa.ideal_gas_density(vmr, pressure, temperature)
    >>> ec = sodium.calc_extinction_coefficient(temperature, density)
    """
    def __init__(self, pressure, *, wn=None, wl=None, cutoff=4500.0):
        """
        Parameters
        ----------
        pressure: 1D float array
            Pressure profile (bars) over which the opacities will be
            evalulated.
        wn: 1D float array
            Wavenumber array (cm-1 units) over which the opacities
            will be sampled (only one of wl or wn should be provided).
        wl: 1D float array
            Wavelength array (micron units) over which the opacities
            will be sampled (only one of wl or wn should be provided).
        cutoff: Float
            Maximum wavenumber extent (cm-1) of the line-profiles
            from the center of each line.
        """
        self.name = 'sodium_vdw'
        self.species = 'Na'

        if wl is None and wn is None:
            raise ValueError(
                'Neither of wavelength (wl) nor wavenumber (wn) were provided'
            )
        if wl is not None and wn is not None:
            raise ValueError(
                'Either provide wavelength or wavenumber array, not both'
            )
        if wn is None:
            wn = 1.0 / (wl*pc.um)

        # Line-transition properties (from VALD, Piskunov 1995):
        # Wavenumber (cm-1), lower-state energy (cm-1), gf (unitless)
        self.wn0 = [16960.87, 16978.07]
        self.elow = [0.0, 0.0]
        self.gf = [0.65464, 1.30918]
        self.nlines = len(self.wn0)

        # Lorentz width parameter (Iro et al. 2005)
        self.lpar = 0.071
        # Partition function (temp < 4000 K, Barklem 2016)
        self.Z = 2.0
        self.detuning = 30.0
        # Default cutoff of 4500 cm-1 following Baudino et al. (2017).
        super().__init__(pressure, wn, cutoff)


class PotassiumVdW(VanderWaals):
    """
    Potassium Van der Waals model (Burrows et al. 2000, ApJ, 531)

    Examples
    --------
    >>> import pyratbay.opacity as op
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.spectrum as ps

    >>> pressure = pa.pressure('1e-8 bar', '1e2 bar', nlayers=81)
    >>> wl = ps.constant_resolution_spectrum(0.5, 1.0, resolution=15000.0)
    >>> potassium = op.alkali.PotassiumVdW(pressure, 1e4/wl)

    >>> temperature = np.tile(1500.0, 81)
    >>> vmr = np.tile(1.1e-7, 81)
    >>> density = pa.ideal_gas_density(vmr, pressure, temperature)
    >>> ec = potassium.calc_extinction_coefficient(temperature, density)
    """
    def __init__(self, pressure, *, wn=None, wl=None, cutoff=4500.0):
        """
        Parameters
        ----------
        pressure: 1D float array
            Pressure profile (bar) over which the opacities will be
            evalulated.
        wn: 1D float array
            Wavenumber array (cm-1 units) over which the opacities
            will be sampled (only one of wl or wn should be provided).
        wl: 1D float array
            Wavelength array (micron units) over which the opacities
            will be sampled (only one of wl or wn should be provided).
        cutoff: Float
            Maximum wavenumber extent (cm-1) of the line-profiles
            from the center of each line.
        """
        self.name = 'potassium_vdw'
        self.species = 'K'

        if wl is None and wn is None:
            raise ValueError(
                'Neither of wavelength (wl) nor wavenumber (wn) were provided'
            )
        if wl is not None and wn is not None:
            raise ValueError(
                'Either provide wavelength or wavenumber array, not both'
            )
        if wn is None:
            wn = 1.0 / (wl*pc.um)

        # Line-transition properties (from VALD, Piskunov 1995):
        # Wavenumber (cm-1), lower-state energy (cm-1), gf (unitless)
        self.wn0 = [12988.76, 13046.486]
        self.elow = [0.0, 0.0]
        self.gf = [0.701455, 1.40929]
        self.nlines = len(self.wn0)

        self.lpar = 0.14
        self.Z = 2.0
        self.detuning = 20.0
        # Default cutoff of 4500 cm-1 following Baudino et al. (2017).
        super().__init__(pressure, wn, cutoff)

