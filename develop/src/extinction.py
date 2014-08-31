import sys, os
import numpy as np

sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/cfuncs/lib')
import ptools as pt
import pconstants as pc
import vprofile as vp

def voigt(pyrat):
  """
  Calculate a grid of Voigt profiles.

  Modification History:
  ---------------------
  2014-08-17  patricio  Initial version.
  """

  pt.msg(pyrat.verb, "\nCalculate Voigt profiles:")
  # Calculate Doppler and Lorentz-width boundaries:
  widthlimits(pyrat)

  # Make Voigt-width arrays:
  pyrat.voigt.doppler = np.logspace(np.log10(pyrat.voigt.Dmin),
                                   np.log10(pyrat.voigt.Dmax), pyrat.voigt.nDop)
  pyrat.voigt.lorentz = np.logspace(np.log10(pyrat.voigt.Lmin),
                                   np.log10(pyrat.voigt.Lmax), pyrat.voigt.nLor)

  # Calculate profiles:
  calcvoigt(pyrat)


def widthlimits(pyrat):
  """
  Calculate the boundaries for the Doppler and Lorentz widths.

  Modification History:
  ---------------------
  2014-08-17  patricio  Initial version.
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
  mols = np.unique(pyrat.iso.imol) # Moleciles with transitions
  mols = mols[np.where(mols>=0)]   # Remove -1's
  # Minimum and maximum mass of molecules with line transitions:
  mmin = np.amin(pyrat.mol.mass[mols])
  mmax = np.amax(pyrat.mol.mass[mols])

  # Get wavenumber array boundaries:
  numin = np.amin(pyrat.wn)
  numax = np.amax(pyrat.wn)

  # Get max pressure:
  pmax = np.amax(pyrat.atm.press)
  # Get max collision diameter:
  cmax = (2.89/2.0 + np.amax(pyrat.mol.radius[mols])) * pc.units['A']
  #cmax = 2.0*np.amax(pyrat.mol.radius[mols]) * pc.units['A']

  # Calculate Doppler-width boundaries:
  if pyrat.voigt.Dmin is None:
    pyrat.voigt.Dmin = np.sqrt(2.0*pc.k*tmin/(mmax*pc.u)) * numin / pc.c
  if pyrat.voigt.Dmax is None:
    pyrat.voigt.Dmax = np.sqrt(2.0*pc.k*tmax/(mmin*pc.u)) * numax / pc.c
  pt.msg(pyrat.verb, "Doppler width limits: {:.3g} -- {:.3g}  cm-1".format(
                                        pyrat.voigt.Dmin, pyrat.voigt.Dmax), 2)

  # Calculate Lorentz-width boundaries:
  if pyrat.voigt.Lmin is None:
    pyrat.voigt.Lmin = pyrat.voigt.Dmin * pyrat.voigt.DLratio

  if pyrat.voigt.Lmax is None:
    pyrat.voigt.Lmax = (np.sqrt(2/(np.pi * pc.k * tmin * pc.u)) * pmax / pc.c *
                        cmax**2.0 * np.sqrt(1.0/mmin + 1.0/2.01588))
  pt.msg(pyrat.verb, "Lorentz width limits: {:.3g} -- {:.3g}  cm-1".format(
                                        pyrat.voigt.Lmin, pyrat.voigt.Lmax), 2)


def calcvoigt(pyrat):
  """
  Wrapper to the Voigt-profile calculator.

  Determine the size of each Voigt profile, find the ones that don't need
  to be recalculated (small Doppler/Lorentz width ratio) and get the profiles.

  Modification History:
  ---------------------
  2014-08-24  patricio  Initial implementation.  Use modified functions from
                        the transit project (newprofile and voigtn).
  """
  # Voigt object from pyrat:
  Voigt = pyrat.voigt

  # Calculate the half-size of the profiles:
  Voigt.size = np.zeros((Voigt.nLor, Voigt.nDop), np.int)
  for i in np.arange(Voigt.nLor):
    # Profile half-width in cm-1:
    pwidth = np.maximum(Voigt.doppler, Voigt.lorentz[i]) * Voigt.width
    # Width in number of spectral samples:
    psize = 2*np.asarray(pwidth/pyrat.wnstep + 0.5, np.int) + 1
    # Clip to max and min values:
    psize = np.clip(psize, 3, 2*pyrat.nspec+1)
    # Set the size to -1 for those who are not being calculated:
    psize[np.where(Voigt.doppler/Voigt.lorentz[i] < Voigt.DLratio)[0][1:]] = 0
    # Store half-size values for this Lorentz width:
    Voigt.size[i] = psize/2

  pt.msg(pyrat.verb, "Calculating Voigt profiles with oversampling factor {:d} "
                 "and Nwidth factor {:d}.".format(Voigt.osamp, Voigt.width), 2)
  # Allocate profile array:
  Voigt.profile = np.zeros((Voigt.osamp, np.sum(2*Voigt.size+1)), np.double)
  # Calculate the Voigt profiles in C:
  vp.voigt(Voigt.profile, Voigt.lorentz, Voigt.doppler,
           Voigt.size,    Voigt.osamp,   pyrat.wnstep, pyrat.verb)



