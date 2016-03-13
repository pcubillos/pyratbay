import sys, os
import struct
import numpy as np
import scipy.interpolate as sip

from .. import tools     as pt
from .. import constants as pc

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../lib')
import vprofile as vp

def voigt(pyrat):
  """
  Driver to calculate a grid of Voigt profiles.
  """

  # Check if reading extinction-coefficient table or no TLI files:
  if (((pyrat.ex.extfile is not None) and os.path.isfile(pyrat.ex.extfile)) or
     pyrat.lt.nTLI == 0):
    pt.msg(pyrat.verb, "\nSkip Voigt-profile calculation.", pyrat.log)
    return

  pt.msg(pyrat.verb, "\nCalculate Voigt profiles:", pyrat.log)
  # Calculate Doppler and Lorentz-width boundaries:
  widthlimits(pyrat)

  # Make Voigt-width arrays:
  pyrat.voigt.doppler = np.logspace(np.log10(pyrat.voigt.Dmin),
                                   np.log10(pyrat.voigt.Dmax), pyrat.voigt.nDop)
  pyrat.voigt.lorentz = np.logspace(np.log10(pyrat.voigt.Lmin),
                                   np.log10(pyrat.voigt.Lmax), pyrat.voigt.nLor)

  # Calculate profiles:
  calcvoigt(pyrat)
  pt.msg(pyrat.verb, "Done.", pyrat.log)


def widthlimits(pyrat):
  """
  Calculate the boundaries for the Doppler and Lorentz widths.
  """

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
  # Minimum and maximum mass of molecules with line transitions:
  mmin = np.amin(pyrat.mol.mass[mols])
  mmax = np.amax(pyrat.mol.mass[mols])

  # Get wavenumber array boundaries:
  numin = np.amin(pyrat.spec.wn)
  numax = np.amax(pyrat.spec.wn)

  # Get max pressure:
  pmax = np.amax(pyrat.atm.press)
  # Get max collision diameter:
  cmax = (2.89/2.0 + np.amax(pyrat.mol.radius[mols])) * pc.A
  #cmax = 2.0*np.amax(pyrat.mol.radius[mols]) * pc.A

  # Calculate Doppler-width boundaries:
  if pyrat.voigt.Dmin is None:
    pyrat.voigt.Dmin = np.sqrt(2.0*pc.k*tmin/(mmax*pc.amu)) * numin / pc.c
  if pyrat.voigt.Dmax is None:
    pyrat.voigt.Dmax = np.sqrt(2.0*pc.k*tmax/(mmin*pc.amu)) * numax / pc.c
  pt.msg(pyrat.verb, "Doppler width limits: {:.3g} -- {:.3g}  cm-1".
                      format(pyrat.voigt.Dmin, pyrat.voigt.Dmax), pyrat.log, 2)

  # Calculate Lorentz-width boundaries:
  if pyrat.voigt.Lmin is None:
    pyrat.voigt.Lmin = pyrat.voigt.Dmin * pyrat.voigt.DLratio

  if pyrat.voigt.Lmax is None:
    pyrat.voigt.Lmax = (np.sqrt(2/(np.pi * pc.k * tmin *pc.amu)) * pmax / pc.c *
                        cmax**2.0 * np.sqrt(1.0/mmin + 1.0/2.01588))
  pt.msg(pyrat.verb, "Lorentz width limits: {:.3g} -- {:.3g}  cm-1".
                     format(pyrat.voigt.Lmin, pyrat.voigt.Lmax), pyrat.log, 2)


def calcvoigt(pyrat):
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
    psize = np.clip(psize, 3, 2*pyrat.spec.nwave+1)
    # Set the size to 0 for those that do not need to be calculated:
    psize[np.where(voigt.doppler/voigt.lorentz[i] < voigt.DLratio)[0][1:]] = 0
    # Store half-size values for this Lorentz width:
    voigt.size[i] = psize/2
  pt.msg(pyrat.verb, "Voigt half-sizes: \n{}".format(voigt.size), pyrat.log, 2)

  pt.msg(pyrat.verb, "Calculating Voigt profiles with Extent:  {:d} widths.".
                     format(voigt.extent), pyrat.log, 2)
  # Allocate profile arrays (concatenated in a 1D array):
  voigt.profile = np.zeros(np.sum(2*voigt.size+1), np.double)
  # Calculate the Voigt profiles in C:
  vp.grid(voigt.profile, voigt.size, voigt.index,
          voigt.lorentz, voigt.doppler,
          pyrat.spec.ownstep,  pyrat.verb)
  pt.msg(pyrat.verb, "Voigt indices:\n{}".format(voigt.index), pyrat.log, 2)
