# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import numpy as np

from .. import tools as pt
from ..opacity import broadening
from ..lib import vprofile as vp


def voigt(pyrat):
    """
    Driver to calculate a grid of Voigt profiles.
    """
    # Check if reading extinction-coefficient table or no TLI files:
    if pt.isfile(pyrat.ex.extfile) == 1 or pyrat.lt.tlifile is None:
        pyrat.log.head('\nSkip LBL Voigt-profile calculation.')
        return

    pyrat.log.head('\nCalculate LBL Voigt profiles:')
    # Calculate Doppler and Lorentz-width boundaries:
    width_limits(pyrat)

    # Make Voigt-width arrays:
    voigt = pyrat.voigt
    voigt.doppler = np.logspace(
        np.log10(voigt.dmin), np.log10(voigt.dmax), voigt.ndop)
    voigt.lorentz = np.logspace(
        np.log10(voigt.lmin), np.log10(voigt.lmax), voigt.nlor)

    # Calculate profiles:
    calc_voigt(pyrat)
    pyrat.log.head('Voigt grid pre-calculation done.')


def width_limits(pyrat):
    """
    Calculate the boundaries for the Doppler and Lorentz widths.
    """
    voigt = pyrat.voigt
    # Get minimum and maximum temperatures:
    tmin =  100.0 if pyrat.ex.tmin is None else pyrat.ex.tmin
    tmax = 3000.0 if pyrat.ex.tmax is None else pyrat.ex.tmax

    # Get mass of line-transition molecules:
    mols = np.unique(pyrat.iso.imol) # Molecules with transitions
    mols = mols[np.where(mols>=0)]   # Remove -1's

    # Estimate min/max Doppler/Lorentz HWHMs from atmospheric properties:
    dmin, lmin = broadening.min_widths(
        tmin, tmax, np.amin(pyrat.spec.wn), np.amax(pyrat.mol.mass[mols]),
        np.amin(pyrat.mol.radius[mols]), np.amin(pyrat.atm.press))

    dmax, lmax = broadening.max_widths(
        tmin, tmax, np.amax(pyrat.spec.wn), np.amin(pyrat.mol.mass[mols]),
        np.amax(pyrat.mol.radius[mols]), np.amax(pyrat.atm.press))

    # Doppler-width boundaries:
    if voigt.dmin is None:
        voigt.dmin = dmin
    if voigt.dmax is None:
        voigt.dmax = dmax
    pyrat.log.msg(
        f'Doppler width limits: {voigt.dmin:.1e} -- {voigt.dmax:.1e} '
        f'cm-1 ({voigt.ndop:d} samples).',
        indent=2)

    # Lorentz-width boundaries:
    if voigt.lmin is None:
        voigt.lmin = lmin
    if voigt.lmax is None:
        voigt.lmax = lmax
    pyrat.log.msg(
        f'Lorentz width limits: {voigt.lmin:.1e} -- {voigt.lmax:.1e} '
        f'cm-1 ({voigt.nlor:d} samples).',
        indent=2)


def calc_voigt(pyrat):
    """
    Wrapper to the Voigt-profile calculator.

    Determine the size of each voigt profile, find the ones that don't need
    to be recalculated (small Doppler/Lorentz width ratio) and get the profiles.
    """
    # Voigt object from pyrat:
    voigt = pyrat.voigt
    voigt.size  = np.zeros((voigt.nlor, voigt.ndop), np.int)
    voigt.index = np.zeros((voigt.nlor, voigt.ndop), np.int)
    # Calculate the half-size of the profiles:
    for i in range(voigt.nlor):
        # Profile half-width in cm-1:
        pwidth = voigt.extent * (
            0.5346*voigt.lorentz[i]
            + np.sqrt(0.2166*voigt.lorentz[i]**2 + voigt.doppler**2))
        # Apply fixed cutoff:
        if voigt.cutoff > 0:
            pwidth = np.minimum(pwidth, voigt.cutoff)
        # Width in number of spectral samples:
        psize = 1 + 2*np.asarray(pwidth/pyrat.spec.ownstep + 0.5, np.int)
        # Clip to max and min values:
        psize = np.clip(psize, 3, 1+2*pyrat.spec.onwave)
        # Temporarily set the size to 0 for not calculated profiles:
        # (sizes will be set in vp.grid())
        skip_voigt = np.where(voigt.doppler/voigt.lorentz[i] < voigt.dlratio)
        psize[skip_voigt[0][1:]] = 0
        voigt.size[i] = psize//2
    pyrat.log.debug(f'Voigt half-sizes:\n{voigt.size}', indent=2)

    cutoff_text = ''
    if voigt.cutoff > 0:
        cutoff_text = f'\nand fixed cutoff: {voigt.cutoff} cm-1'
    pyrat.log.msg(
        f'Calculating Voigt profiles with max extent: {voigt.extent:.1f} HWHM'
        f'{cutoff_text}.', indent=2)
    # Allocate profile arrays (concatenated in a 1D array):
    voigt.profile = np.zeros(np.sum(2*voigt.size+1), np.double)

    # Calculate the Voigt profiles in C:
    vp.grid(
        voigt.profile, voigt.size, voigt.index,
        voigt.lorentz, voigt.doppler,
        pyrat.spec.ownstep, pyrat.verb)
    pyrat.log.debug(f'Voigt indices:\n{voigt.index}', indent=2)
