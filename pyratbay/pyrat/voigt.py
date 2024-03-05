# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import numpy as np

from .. import tools as pt
from .. import opacity as op
from ..lib import vprofile as vp


class Voigt():
    """Grid of Voigt profiles"""
    def __init__(self, inputs, lbl, ex, spec, atm, log):
        self.profile = None  # Voigt profile [sum(2*size+1)]
        self.dlratio = None  # Doppler-Lorentz ratio threshold
        self.size = None  # Profile wavenumber half-size [nlor, ndop]
        self.index = None  # Index where each profile starts [nlor, ndop]

        # Profiles extent (in HWHMs) and cutoff (in cm-1)
        self.extent = inputs.voigt_extent
        self.cutoff = inputs.voigt_cutoff
        self.dlratio = inputs.voigt_dlratio

        if inputs.tlifile is None:
            log.head('\nSkip LBL Voigt-profile calculation.')
            return

        log.head('\nCalculate LBL Voigt profiles:')
        # Estimate the boundaries for the Doppler and Lorentz widths.
        min_wn = np.amin(spec.wn)
        max_wn = np.amax(spec.wn)
        min_pressure = np.amin(atm.press)
        max_pressure = np.amax(atm.press)

        min_temp = 100.0 if ex.tmin is None else ex.tmin
        max_temp = 3000.0 if ex.tmax is None else ex.tmax

        # Get mass of line-transition molecules:
        mol_indices = np.unique(lbl.iso_atm_index)
        min_mass = np.amin(atm.mol_mass[mol_indices])
        max_mass = np.amax(atm.mol_mass[mol_indices])
        min_rad = np.amin(atm.mol_radius[mol_indices])
        max_rad = np.amax(atm.mol_radius[mol_indices])

        # Estimate min/max Doppler/Lorentz HWHMs from atmospheric properties:
        dmin, lmin = op.broadening.min_widths(
            min_temp, max_temp, min_wn,
            max_mass, min_rad, min_pressure,
        )
        dmax, lmax = op.broadening.max_widths(
            min_temp, max_temp, max_wn,
            min_mass, max_rad, max_pressure,
        )

        # Doppler-width boundaries:
        self.dmin = inputs.voigt_dmin
        if self.dmin is None:
            self.dmin = dmin
        self.dmax = inputs.voigt_dmax
        if self.dmax is None:
            self.dmax = dmax
        if self.dmax <= self.dmin:
            log.error(
                f'Voigt dmax ({self.dmax:.4e} cm-1) must be > dmin '
                f'({self.dmin:.4e} cm-1)'
            )

        self.ndop = inputs.voigt_ndop
        log.msg(
            f'Doppler width limits: {self.dmin:.1e} -- {self.dmax:.1e} '
            f'cm-1 ({self.ndop:d} samples).',
            indent=2,
        )
        # Grid of Doppler HWHMs
        self.doppler = np.logspace(
            np.log10(self.dmin),
            np.log10(self.dmax),
            self.ndop,
        )

        # Lorentz-width boundaries:
        self.lmin = inputs.voigt_lmin
        if self.lmin is None:
            self.lmin = lmin
        self.lmax = inputs.voigt_lmax
        if self.lmax is None:
            self.lmax = lmax
        if self.lmax <= self.lmin:
            log.error(
                f'Voigt lmax ({self.lmax:.4e} cm-1) must be > lmin '
                f'({self.lmin:.4e} cm-1)'
            )

        self.nlor = inputs.voigt_nlor
        log.msg(
            f'Lorentz width limits: {self.lmin:.1e} -- {self.lmax:.1e} '
            f'cm-1 ({self.nlor:d} samples).',
            indent=2,
        )
        # Grid of Lorentz HWHMs
        self.lorentz = np.logspace(
            np.log10(self.lmin),
            np.log10(self.lmax),
            self.nlor,
        )


        # Calculate grid of Voigt profiles:
        self.size = np.zeros((self.nlor, self.ndop), int)
        self.index = np.zeros((self.nlor, self.ndop), int)
        # Calculate the half-size of the profiles:
        for i in range(self.nlor):
            # Profile half-width in cm-1:
            pwidth = self.extent * (
                0.5346*self.lorentz[i]
                + np.sqrt(0.2166*self.lorentz[i]**2 + self.doppler**2)
            )
            # Apply fixed cutoff:
            if self.cutoff > 0:
                pwidth = np.minimum(pwidth, self.cutoff)
            # Width in number of spectral samples:
            psize = 1 + 2*np.asarray(pwidth/spec.ownstep + 0.5, int)
            # Clip to max and min values:
            psize = np.clip(psize, 3, 1+2*spec.onwave)
            # Temporarily set the size to 0 for not calculated profiles:
            # (sizes will be reset in vp.grid())
            skip = self.doppler/self.lorentz[i] < self.dlratio
            skip[0] = False
            psize[skip] = 0
            self.size[i] = psize//2
        log.debug(f'Voigt half-sizes:\n{self.size}', indent=2)

        cutoff_text = ''
        if self.cutoff > 0:
            cutoff_text = f'\nand fixed cutoff: {self.cutoff} cm-1'
        log.msg(
            'Calculating Voigt profiles with max extent: '
            f'{self.extent:.1f} HWHM{cutoff_text}.',
            indent=2,
        )
        # Allocate profile arrays (concatenated in a 1D array):
        self.profile = np.zeros(np.sum(2*self.size+1), np.double)

        # Calculate the Voigt profiles in C:
        vp.grid(
            self.profile, self.size, self.index,
            self.lorentz, self.doppler,
            spec.ownstep, log.verb,
        )
        log.debug(f'Voigt indices:\n{self.index}', indent=2)
        log.head('Voigt grid pre-calculation done.')


    def __str__(self):
        fw = pt.Formatted_Write(fmt={'float':'{:.3e}'.format}, edge=3)
        fw.write('Voigt-profile information:')
        fw.write('\nNumber of Doppler-width samples (ndop): {:d}', self.ndop)
        fw.write('Number of Lorentz-width samples (nlor): {:d}', self.nlor)
        fw.write('Doppler HWHM (doppler, cm-1):\n    {}', self.doppler)
        fw.write('Lorentz HWMH (lorentz, cm-1):\n    {}', self.lorentz)
        fw.write(
            f'Doppler--Lorentz ratio threshold (dlratio): {self.dlratio:.3e}'
        )
        fw.write(
            f"\nVoigt-profiles extent (extent, in HWHMs): {self.extent:.1f}",
        )
        fw.write(
            f"Voigt-profiles cutoff extent (cutoff in cm-1): {self.cutoff:.1f}",
        )
        fw.write(
            'Voigt-profile half-sizes (size) of shape [nlor, ndop]:\n{}',
            self.size,
            edge=2,
        )
        fw.write(
            'Voigt-profile indices (index) of shape [nlor, ndop]:\n{}',
            self.index,
            edge=2,
        )

        index, size = self.index[0,0], 2*self.size[0,0]+1
        fw.write(
            '\nVoigt profiles:\n  profile[ 0, 0]: {}',
            self.profile[index:index+size],
            fmt={'float':'{:.5e}'.format}, edge=2)
        index = self.index[self.nlor-1,self.ndop-1]
        size = 2*self.size[self.nlor-1,self.ndop-1] + 1
        fw.write(
            '  ...\n  profile[{:2d},{:2d}]: {}',
            self.nlor-1, self.ndop-1,
            self.profile[index:index+size],
            fmt={'float':'{:.5e}'.format}, edge=2)
        return fw.text

