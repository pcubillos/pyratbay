import sys, os
import struct
import numpy as np
import scipy.interpolate as sip

import ptools     as pt
import pconstants as pc
import extcoeff   as ec

def exttable(pyrat):
  """
  Handle extinction-coefficient table (read/calculate/write file).

  Modification History:
  ---------------------
  2015-01-19  patricio  pre-initial implementation.
  2015-01-25  patricio  intitial version.
  """

  # If the extinction file was not defined, skip this step:
  if pyrat.ex.extfile is None:
    pt.msg(pyrat.verb-10, "No extinction coefficient table requested.")
    return

  pt.msg(pyrat.verb, "\nBegin Extinction-coefficient table handling.")
  # If the extinction file exists, read it:
  if os.path.isfile(pyrat.ex.extfile):
    pt.msg(pyrat.verb, "Reading extinction-coefficient table file:"
                     "\n  '{:s}'".format(pyrat.ex.extfile), 2)
    read_extinction(pyrat)

  # If the extinction file doesn't exist, calculate it:
  else:
    pt.msg(pyrat.verb, "Generating new extinction-coefficient table file:"
                     "\n  '{:s}'".format(pyrat.ex.extfile), 2)
    calc_extinction(pyrat)


def read_extinction(pyrat):
  """
  Read an extinction-coefficient table from file.

  Modification History:
  ---------------------
  2015-01-25  patricio  intitial version.
  """
  ex = pyrat.ex               # Extinction-coefficient object
  f = open(ex.extfile, "rb")  # Open extinction coefficient file

  # Read arrays lengths:
  ex.nmol    = struct.unpack('l', f.read(8))[0]
  ex.ntemp   = struct.unpack('l', f.read(8))[0]
  ex.nlayers = struct.unpack('l', f.read(8))[0]
  ex.nspec   = struct.unpack('l', f.read(8))[0]
  pt.msg(pyrat.verb-10, "File has {:d} molecules, {:d} temperature samples, "
                       "{:d} layers, and {:d} wavenumber samples.".
                       format(ex.nmol, ex.ntemp, ex.nlayers, ex.nspec), 2)

  # Read wavenumber, temperature, pressure, and isotope arrays:
  ex.molID = np.asarray(struct.unpack(str(ex.nmol   )+'i', f.read(4*ex.nmol )))
  ex.temp  = np.asarray(struct.unpack(str(ex.ntemp  )+'d', f.read(8*ex.ntemp)))
  ex.press = np.asarray(struct.unpack(str(ex.nlayers)+'d',f.read(8*ex.nlayers)))
  ex.wn    = np.asarray(struct.unpack(str(ex.nspec  )+'d', f.read(8*ex.nspec)))

  pt.msg(pyrat.verb-15, "Molecules' IDs: {}".format(ex.molID), 2)
  pt.msg(pyrat.verb-15, "Temperatures (K): {}".
                         format(pt.pprint(ex.temp, fmt=np.int)), 2)
  pt.msg(pyrat.verb-15, "Pressure layers (bar): {}".
                         format(pt.pprint(ex.press*1e-6,3)), 2)
  pt.msg(pyrat.verb-15, "Wavenumber array (cm-1): {}".
                         format(pt.pprint(ex.wn, 1)), 2)

  # Read extinction-coefficient data table:
  ndata = ex.nmol * ex.ntemp * ex.nlayers * ex.nspec
  data = np.asarray(struct.unpack('d'*ndata, f.read(8*ndata)))

  pyrat.ex.etable = np.reshape(data, (ex.nmol, ex.ntemp, ex.nlayers, ex.nspec))


def calc_extinction(pyrat):
  """
  Compute the extinction-coefficient (per molecule) for a tabulated grid of
  temperatures and pressures, over a wavenumber array.

  Modification History:
  ---------------------
  2015-01-25  patricio  intitial version.
  """
  # Extinction-coefficient object:
  ex = pyrat.ex

  # Make temperature sample:
  if ex.tmin < 0.0 or ex.tmax < 0.0 or ex.tstep < 0.0:
    pt.error("All temperature variables: tmin ({:.1f}), tmax ({:.1f}), and "
       "tstep ({:.1f}) must be defined.".format(ex.tmin, ex.tmax, ex.tstep))

  ex.ntemp = int((ex.tmax-ex.tmin)/ex.tstep) + 1
  ex.temp  = np.linspace(ex.tmin, ex.tmin + (ex.ntemp-1)*ex.tstep, ex.ntemp)
  pt.msg(pyrat.verb-15, "Temperature sample (K): {:s}".
                         format(pt.pprint(ex.temp)), 2)

  # Evaluate the partition function at the given temperatures:
  pt.msg(pyrat.verb, "Interpolate partition function.", 2)
  ex.z = np.zeros((pyrat.iso.niso, ex.ntemp), np.double)
  for i in np.arange(pyrat.lt.ndb):           # For each Database
    for j in np.arange(pyrat.lt.db[i].niso):  # For each isotope in DB
      zinterp = sip.interp1d(pyrat.lt.db[i].temp, pyrat.lt.db[i].z[j],
                             kind='cubic')
      ex.z[pyrat.lt.db[i].iiso+j] = zinterp(ex.temp)

  # Allocate wavenumber, pressure, and isotope arrays:
  ex.wn    = pyrat.wn
  ex.nspec = pyrat.nspec

  ex.molID = pyrat.mol.ID[np.unique(pyrat.iso.imol)]
  ex.nmol  = len(ex.molID)

  ex.press = pyrat.atm.press
  ex.nlayers = pyrat.atm.nlayers

  # Allocate extinction-coefficient array:
  pt.msg(pyrat.verb, "Calculate extinction coefficient.", 2)
  ex.etable = np.zeros((ex.nmol, ex.ntemp, ex.nlayers, ex.nspec), np.double)

  # Compute extinction (in C):
  for r in np.arange(ex.nlayers):
    for t in np.arange(ex.ntemp):
      # Extinction coefficient for given temperature and pressure-layer:
      pt.msg(pyrat.verb, "\nR={}, T={}".format(r,t))
      extinction(pyrat, ex.etable[:,t,r], r, ex.temp[t], ex.z[:,t])

  # Store values in file:
  f = open(ex.extfile, "wb")

  # Write size of arrays:
  f.write(struct.pack("4l", ex.nmol, ex.ntemp, ex.nlayers, ex.nspec))
  # Write arrays:
  f.write(struct.pack(str(ex.nmol)   +"i", *list(ex.molID)))
  f.write(struct.pack(str(ex.ntemp)  +"d", *list(ex.temp) ))
  f.write(struct.pack(str(ex.nlayers)+"d", *list(ex.press)))
  f.write(struct.pack(str(ex.nspec)  +"d", *list(ex.wn)   ))
  # Write extinction-coefficient data:
  fmt = str(ex.nmol * ex.ntemp * ex.nlayers * ex.nspec) + "d"
  f.write(struct.pack(fmt, *list(pyrat.ex.etable.flatten())))
  f.close()
  pt.msg(pyrat.verb, "Extinction-coefficient table written to file:"
                     " '{:s}'.".format(ex.extfile), 2)


def extinction(pyrat, extcoeff, ilayer, temp, ziso, add=0):
  """
  Python wrapper for the extinction-coefficient calculation function.
  Compute the EC over the wavenumber array for the given temperature
  and pressure.

  Parameters:
  -----------
  pyrat: Pyrat Object
  extcoeff: 2D float ndarray
     Output extinction coefficient ().
  ilayer: Integer
     Index of the layer.
  temp: Float
     Layer's temperature.
  ziso: 1D float ndarray
     Isotopes' partition fuction at the given layer.
  add:  Boolean
     If True multiply the extinction coefficient by the density and co-add
     the contribution from all molecules.
  """
  # Unpack layer parameters:
  pressure = pyrat.atm.press[ilayer]  # Layer pressure
  molq     = pyrat.atm.q    [ilayer]  # Molecular abundance
  density  = pyrat.atm.d    [ilayer]  # Molecular density

  # Get mol index in extinction coefficient table for the isotopes:
  pyrat.iso.iext = np.zeros(pyrat.iso.niso, np.int)
  if pyrat.ex.extfile is not None:
    for i in np.arange(pyrat.iso.niso):
      pyrat.iso.iext[i] = np.where(pyrat.ex.molID ==
                                   pyrat.mol.ID[pyrat.iso.imol[i]])[0][0]

  # Calculate extinction-coefficient in C:
  ec.extinction(extcoeff,
                pyrat.voigt.profile, pyrat.voigt.size, pyrat.voigt.index,
                pyrat.voigt.lorentz, pyrat.voigt.doppler,
                pyrat.wn, pyrat.own, pyrat.odivisors,
                density, molq, pyrat.mol.radius, pyrat.mol.mass,
                pyrat.iso.imol, pyrat.iso.mass, pyrat.iso.ratio,
                ziso, pyrat.iso.iext,
                pyrat.lt.wn, pyrat.lt.elow, pyrat.lt.gf, pyrat.lt.isoid,
                pyrat.ex.ethresh, pressure, temp, add)

