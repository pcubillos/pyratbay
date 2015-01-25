import sys, os
import struct
import numpy as np
import scipy.interpolate as sip

sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/cfuncs/lib')
import ptools     as pt
import pconstants as pc
import vprofile   as vp
import extcoeff   as ec

def voigt(pyrat):
  """
  Calculate a grid of voigt profiles.

  Modification History:
  ---------------------
  2014-08-17  patricio  Initial version.
  """

  pt.msg(pyrat.verb, "\nCalculate voigt profiles:")
  # Calculate Doppler and Lorentz-width boundaries:
  widthlimits(pyrat)

  # Make voigt-width arrays:
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
  mols = np.unique(pyrat.iso.imol) # Molecules with transitions
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

  Determine the size of each voigt profile, find the ones that don't need
  to be recalculated (small Doppler/Lorentz width ratio) and get the profiles.

  Modification History:
  ---------------------
  2014-08-24  patricio  Initial implementation.  Use modified functions from
                        the transit project (newprofile and voigtn).
  """
  # voigt object from pyrat:
  voigt = pyrat.voigt

  # Calculate the half-size of the profiles:
  voigt.size = np.zeros((voigt.nLor, voigt.nDop), np.int)
  for i in np.arange(voigt.nLor):
    # Profile half-width in cm-1:
    pwidth = np.maximum(voigt.doppler, voigt.lorentz[i]) * voigt.width
    # Width in number of spectral samples:
    psize = 2*np.asarray(pwidth/pyrat.ownstep + 0.5, np.int) + 1
    # Clip to max and min values:
    psize = np.clip(psize, 3, 2*pyrat.nspec+1)
    # Set the size to 0 for those who do not need to be calculated:
    psize[np.where(voigt.doppler/voigt.lorentz[i] < voigt.DLratio)[0][1:]] = 0
    # Store half-size values for this Lorentz width:
    voigt.size[i] = psize/2

  pt.msg(pyrat.verb, "Voigt half-size: {}".format(voigt.size))
  print(pyrat.ownstep)
  pt.msg(pyrat.verb, "Calculating voigt profiles with Nwidth factor {:d}.".
                     format(voigt.width), 2)
  # Allocate profile arrays (concatenated in a 1D array):
  voigt.profile = np.zeros(np.sum(2*voigt.size+1), np.double)
  # Calculate the Voigt profiles in C:
  vp.voigt(voigt.profile, voigt.lorentz, voigt.doppler,
           voigt.size,    pyrat.ownstep,  pyrat.verb)


def opacity(pyrat):
  """
  Handle extinction-coefficient file (read/calculate/write).

  Modification History:
  ---------------------
  2015-01-19  patricio  pre-initial implementation.
  """
  if pyrat.ex.extfile is None:
    pt.msg(pyrat.verb-10, "No extinction coefficient table requested.")
    return

  # If file exists read:
  elif os.path.isfile(pyrat.ex.extfile):
    pt.msg(pyrat.verb, "Reading extinction-coefficient table file:"
                     "\n  '{:s}'".format(pyrat.ex.extfile))
    read_extinction(pyrat)

  # If it doesn't exist, calculate it:
  else:
    pt.msg(pyrat.verb, "Generating new extinction-coefficient table file:"
                     "\n  '{:s}'".format(pyrat.ex.extfile))
    calc_extinction(pyrat)


def read_extinction(pyrat):
  """
  Read an extinction-coefficient table from file.
  """
  ex = pyrat.ex
  # Open extinction coefficient file:
  f = open(pyrat.ex.extfile, "rb")

  # Read arrays lengths:
  ex.nmol    = struct.unpack('l', f.read(8))[0]
  ex.ntemp   = struct.unpack('l', f.read(8))[0]
  ex.nlayers = struct.unpack('l', f.read(8))[0]
  ex.nspec   = struct.unpack('l', f.read(8))[0]
  pt.msg(pyrat.verb-10, "File has {:d} molecules, {:d} temperature samples, "
                       "{:d} layers, and {:d} wavenumber samples.".
                       format(ex.nmol, ex.ntemp, ex.nlayers, ex.nspec), 2)

  # Read wavenumber, temperature, pressure, and isotope arrays:
  ex.molID = np.asarray(struct.unpack(str(ex.nmol )+'i',  f.read(4*ex.nmol )))
  pt.msg(pyrat.verb-15, "Molecules' IDs: {}".format(ex.mol), 2)

  ex.temp  = np.asarray(struct.unpack(str(ex.ntemp)+'d',  f.read(8*ex.ntemp)))
  pt.msg(pyrat.verb-15, "Temperatures (K): {}".
                         format(pt.pprint(ex.temp, fmt=np.int)), 2)

  ex.press = np.asarray(struct.unpack(str(ex.nlayers)+'d',f.read(8*ex.nlayers)))
  pt.msg(pyrat.verb-15, "Pressure layers (bar): {}".
                         format(pt.pprint(ex.press*1e-6,3)), 2)

  ex.wn    = np.asarray(struct.unpack(str(ex.nspec)+'d',  f.read(8*ex.nspec)))
  pt.msg(pyrat.verb-15, "Wavenumber array (cm-1): {}".
                         format(pt.pprint(ex.wn, 1)), 2)

  # Read extinction-coefficient data table:
  ndata = ex.nmol * ex.ntemp * ex.nlayers * ex.nspec
  data = np.asarray(struct.unpack('d'*ndata, f.read(8*ndata)))

  pyrat.ex.etable = np.reshape(data, (ex.nmol, ex.ntemp, ex.nlayers, ex.nspec))


def calc_extinction(pyrat):
  """
  """
  # Extinction-coefficient object:
  ex = pyrat.ex

  # Make temperature sample:
  if ex.tmin < 0.0 or ex.tmax < 0.0 or ex.tstep < 0.0:
    pt.error("All temperature variables: tmin ({:.1f}), tmax ({:.1f}), and "
       "tstep ({:.1f}) must be defined.".format(ex.tmin, ex.tmax, ex.tstep))
  ex.temp = np.arange(1.0*ex.tmin, ex.tmax, ex.tstep)
  ex.ntemp = len(ex.temp)
  pt.msg(pyrat.verb-15, "Temperature sample (K): {:s}".
                        format(pt.pprint(ex.temp)), 2)

  # Evaluate partition function:
  Ztable = np.zeros((pyrat.iso.niso, ex.ntemp), np.double)
  for i in np.arange(pyrat.lt.ndb):           # For each Database
    for j in np.arange(pyrat.lt.db[i].niso):  # For each isotope in DB
      zinterp = sip.interp1d(pyrat.lt.db[i].temp, pyrat.lt.db[i].z[j],
                             kind='cubic')
      Ztable[pyrat.lt.db[i].iiso+j] = zinterp(ex.temp)

  # Allocate wavenumber, pressure, and isotope arrays:
  ex.wn    = pyrat.wn
  ex.nspec = pyrat.nspec

  ex.molID = pyrat.mol.ID[np.unique(pyrat.iso.imol)]
  ex.nmol  = len(molID)

  ex.press = pyrat.atm.press
  ex.nlayers = pyrat.nlayers

  # Allocate extinction-coefficient array:
  ex.etable = np.zeros((ex.nmol, ex.ntemp, ex.nlayers, ex.nspec), np.double)

  # Compute extinction (do it in C):
  for nlayer
    for ntemp
      for nmol
      ex.etable[m, t, r] = extinction(pyrat, r, t, m)

  # Store values in file:


def extinction(pyrat, pressure, temp, molID=None):
  """
  Calculate the extinction coefficient over the wavenumber array for the
  given temperature and pressure.
  """
  # Unpack stuff:
  pressure = pyrat.press[r] # Layer pressure
  mm   = pyrat.atm.mm[r] # Mean molecular mass
  molq = pyrat.atm.q[r]  # Molecular abundance

  # Call coeff in C:
  ec.coeff(pyrat.voigt.profile, pyrat.voigt.size,
           pyrat.voigt.lorentz, pyrat.voigt.doppler,
           pyrat.wn, pyrat.own,
           molq, pyrat.mol.radius, pyrat.mol.mass,
           pyrat.iso.imol, pyrat.iso.mass,
           pressure, temp)

  # Return 





