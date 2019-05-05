# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import numpy as np

from .. import constants  as pc
from .. import broadening as broad

sys.path.append(pc.ROOT+'/lib')
import vprofile as vp


def voigt(pyrat):
  """
  Driver to calculate a grid of Voigt profiles.
  """
  # Check if reading extinction-coefficient table or no TLI files:
  if ((pyrat.ex.extfile is not None and os.path.isfile(pyrat.ex.extfile))
      or pyrat.lt.nTLI == 0):
      pyrat.log.msg('\nSkip LBL Voigt-profile calculation.')
      return

  pyrat.log.msg('\nCalculate LBL Voigt profiles:')
  # Calculate Doppler and Lorentz-width boundaries:
  width_limits(pyrat)

  # Make Voigt-width arrays:
  voigt = pyrat.voigt
  voigt.doppler = np.logspace(np.log10(voigt.Dmin), np.log10(voigt.Dmax),
                              voigt.nDop)
  voigt.lorentz = np.logspace(np.log10(voigt.Lmin), np.log10(voigt.Lmax),
                              voigt.nLor)

  # Calculate profiles:
  calc_voigt(pyrat)
  pyrat.log.msg('Voigt grid pre-calculation done.')


def width_limits(pyrat):
  """
  Calculate the boundaries for the Doppler and Lorentz widths.
  """
  voigt = pyrat.voigt
  # Get minimum temperature:
  tmin = pyrat.ex.tmin
  if tmin is None:
      tmin = np.amin(pyrat.atm.temp)
  # Get maximum temperature:
  tmax = pyrat.ex.tmax
  if tmax is None:
      tmax = np.amax(pyrat.atm.temp)

  # Get mass of line-transition molecules:
  mols = np.unique(pyrat.iso.imol) # Molecules with transitions
  mols = mols[np.where(mols>=0)]   # Remove -1's

  # Estimate min/max Doppler/Lorentz HWHMs from atmospheric properties:
  Dmin, Lmin = broad.min_widths(tmin, np.amin(pyrat.spec.wn),
      np.amax(pyrat.mol.mass[mols]), voigt.DLratio)

  Dmax, Lmax = broad.max_widths(tmin, tmax, np.amax(pyrat.spec.wn),
      np.amin(pyrat.mol.mass[mols]), np.amax(pyrat.mol.radius[mols]),
      np.amax(pyrat.atm.press))

  # Doppler-width boundaries:
  if voigt.Dmin is None:
      voigt.Dmin = Dmin
  if voigt.Dmax is None:
      voigt.Dmax = Dmax
  pyrat.log.msg('Doppler width limits: {:.1e} -- {:.1e} cm-1 ({:d} samples).'.
      format(voigt.Dmin, voigt.Dmax, voigt.nDop), verb=2, indent=2)

  # Lorentz-width boundaries:
  if voigt.Lmin is None:
      voigt.Lmin = Lmin
  if voigt.Lmax is None:
      voigt.Lmax = Lmax
  pyrat.log.msg('Lorentz width limits: {:.1e} -- {:.1e} cm-1 ({:d} samples).'.
      format(voigt.Lmin, voigt.Lmax, voigt.nLor), verb=2, indent=2)


def calc_voigt(pyrat):
  """
  Wrapper to the Voigt-profile calculator.

  Determine the size of each voigt profile, find the ones that don't need
  to be recalculated (small Doppler/Lorentz width ratio) and get the profiles.
  """
  # Voigt object from pyrat:
  voigt = pyrat.voigt
  voigt.size  = np.zeros((voigt.nLor, voigt.nDop), np.int)
  voigt.index = np.zeros((voigt.nLor, voigt.nDop), np.int)
  # Calculate the half-size of the profiles:
  for i in np.arange(voigt.nLor):
      # Profile half-width in cm-1:
      pwidth = np.maximum(voigt.doppler, voigt.lorentz[i]) * voigt.extent
      # Width in number of spectral samples:
      psize = 2*np.asarray(pwidth/pyrat.spec.ownstep + 0.5, np.int) + 1
      # Clip to max and min values:
      psize = np.clip(psize, 3, 2*pyrat.spec.onwave+1)
      # Temporarily Set the size to 0 for not calculated profiles:
      # (sizes will be set in vp.grid())
      psize[np.where(voigt.doppler/voigt.lorentz[i] < voigt.DLratio)[0][1:]] = 0
      # Store half-size values for this Lorentz width:
      voigt.size[i] = psize//2
  pyrat.log.msg('Voigt half-sizes:\n{}'.format(voigt.size), verb=3, indent=2)

  pyrat.log.msg('Calculating Voigt profiles with Extent: {:.1f} widths.'.
      format(voigt.extent), verb=2, indent=2)
  # Allocate profile arrays (concatenated in a 1D array):
  voigt.profile = np.zeros(np.sum(2*voigt.size+1), np.double)

  # Calculate the Voigt profiles in C:
  vp.grid(voigt.profile, voigt.size, voigt.index,
          voigt.lorentz, voigt.doppler,
          pyrat.spec.ownstep, pyrat.verb)
  pyrat.log.msg('Voigt indices:\n{}'.format(voigt.index), verb=3, indent=2)
