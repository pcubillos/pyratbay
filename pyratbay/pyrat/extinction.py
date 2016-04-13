import sys, os
import struct
import numpy as np
import scipy.interpolate as sip

from .. import tools     as pt
from .. import constants as pc

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../lib')
import extcoeff   as ec

def exttable(pyrat):
  """
  Handle extinction-coefficient table (read/calculate/write file).
  """

  # If the extinction file was not defined, skip this step:
  if pyrat.ex.extfile is None:
    pt.msg(pyrat.verb-3, "No extinction-coefficient table requested.",
           pyrat.log)
    return

  pt.msg(pyrat.verb-3, "\nBegin extinction-coefficient table handling.",
         pyrat.log)
  # If the extinction file exists, read it:
  if os.path.isfile(pyrat.ex.extfile):
    pt.msg(pyrat.verb-4, "Reading extinction-coefficient table file:"
                     "\n  '{:s}'".format(pyrat.ex.extfile), pyrat.log, 2)
    read_extinction(pyrat)

  # If the extinction file doesn't exist, calculate it:
  else:
    pt.msg(pyrat.verb-4, "Generating new extinction-coefficient table file:"
                     "\n  '{:s}'".format(pyrat.ex.extfile), pyrat.log, 2)
    calc_extinction(pyrat)


def read_extinction(pyrat):
  """
  Read an extinction-coefficient table from file.
  """
  ex = pyrat.ex               # Extinction-coefficient object
  f = open(ex.extfile, "rb")  # Open extinction coefficient file

  # Read arrays lengths:
  ex.nmol    = struct.unpack('l', f.read(8))[0]
  ex.ntemp   = struct.unpack('l', f.read(8))[0]
  ex.nlayers = struct.unpack('l', f.read(8))[0]
  ex.nwave   = struct.unpack('l', f.read(8))[0]
  pt.msg(pyrat.verb-4, "File has {:d} molecules, {:d} temperature samples, "
           "{:d} layers, and {:d} wavenumber samples.".
           format(ex.nmol, ex.ntemp, ex.nlayers, ex.nwave), pyrat.log, 2)

  # Read wavenumber, temperature, pressure, and isotope arrays:
  # FINDME: pt.unpack
  ex.molID = np.asarray(struct.unpack(str(ex.nmol   )+'i', f.read(4*ex.nmol )))
  ex.temp  = np.asarray(struct.unpack(str(ex.ntemp  )+'d', f.read(8*ex.ntemp)))
  ex.press = np.asarray(struct.unpack(str(ex.nlayers)+'d',f.read(8*ex.nlayers)))
  ex.wn    = np.asarray(struct.unpack(str(ex.nwave  )+'d', f.read(8*ex.nwave)))

  pt.msg(pyrat.verb-4, "Molecules' IDs: {}".format(ex.molID),    pyrat.log, 2)
  pt.msg(pyrat.verb-4, "Temperatures (K): {}".
                         format(pt.pprint(ex.temp, fmt=np.int)), pyrat.log, 2)
  pt.msg(pyrat.verb-4, "Pressure layers (bar): {}".
                         format(pt.pprint(ex.press/pc.bar,3)),   pyrat.log, 2)
  pt.msg(pyrat.verb-4, "Wavenumber array (cm-1): {}".
                         format(pt.pprint(ex.wn, 1)), pyrat.log, 2)

  # Read extinction-coefficient data table:
  ndata = ex.nmol * ex.ntemp * ex.nlayers * ex.nwave
  data = np.asarray(struct.unpack('d'*ndata, f.read(8*ndata)))

  pyrat.ex.etable = np.reshape(data, (ex.nmol, ex.ntemp, ex.nlayers, ex.nwave))


def calc_extinction(pyrat):
  """
  Compute the extinction-coefficient (per species) for a tabulated grid of
  temperatures and pressures, over a wavenumber array.
  """
  # Extinction-coefficient object:
  ex = pyrat.ex

  # Make temperature sample:
  if (ex.tmin is None) or (ex.tmax is None):
    pt.error("Both, tmin ({}) and tmax ({}), must be defined to produce "
             "the extinction-coefficient table.".format(ex.tmin, ex.tmax),
             pyrat.log)
  # Temperature boundaries check:
  if (ex.tmin < pyrat.lt.tmin):
    pt.error("The extinction-coefficient table attempted to sample a "
             "temperature ({:.1f} K) below the lowest allowed TLI temperature "
             "({:.1f} K).".format(ex.tmin, pyrat.lt.tmin), pyrat.log)
  if (ex.tmax > pyrat.lt.tmax):
    pt.error("The extinction-coefficient table attempted to sample a "
             "temperature ({:.1f} K) above the highest allowed TLI temperature "
             "({:.1f} K).".format(ex.tmax, pyrat.lt.tmax), pyrat.log)
  # OK, now create the temperature array:
  ex.ntemp = int((ex.tmax-ex.tmin)/ex.tstep) + 1
  ex.temp  = np.linspace(ex.tmin, ex.tmin + (ex.ntemp-1)*ex.tstep, ex.ntemp)

  pt.msg(pyrat.verb-4, "Temperature sample (K): {:s}".
                         format(pt.pprint(ex.temp)), pyrat.log, 2)

  # Evaluate the partition function at the given temperatures:
  pt.msg(pyrat.verb-4, "Interpolate partition function.", pyrat.log, 2)
  ex.z = np.zeros((pyrat.iso.niso, ex.ntemp), np.double)
  for i in np.arange(pyrat.lt.ndb):           # For each Database
    for j in np.arange(pyrat.lt.db[i].niso):  # For each isotope in DB
      zinterp = sip.interp1d(pyrat.lt.db[i].temp, pyrat.lt.db[i].z[j],
                             kind='cubic')
      ex.z[pyrat.lt.db[i].iiso+j] = zinterp(ex.temp)

  # Allocate wavenumber, pressure, and isotope arrays:
  ex.wn    = pyrat.spec.wn
  ex.nwave = pyrat.spec.nwave

  ex.molID = pyrat.mol.ID[np.unique(pyrat.iso.imol)]
  ex.nmol  = len(ex.molID)

  ex.press = pyrat.atm.press
  ex.nlayers = pyrat.atm.nlayers

  # Allocate extinction-coefficient array:
  pt.msg(pyrat.verb-4, "Calculate extinction coefficient.", pyrat.log, 2)
  ex.etable = np.zeros((ex.nmol, ex.ntemp, ex.nlayers, ex.nwave), np.double)

  # Compute extinction (in C):
  for r in np.arange(ex.nlayers):
    for t in np.arange(ex.ntemp):
      # Extinction coefficient for given temperature and pressure-layer:
      pt.msg(pyrat.verb-4, "\nR={}, T={}".format(r,t), pyrat.log)
      extinction(pyrat, ex.etable[:,t,r], r, ex.temp[t], ex.z[:,t])

  # Store values in file:
  f = open(ex.extfile, "wb")

  # Write size of arrays:
  f.write(struct.pack("4l", ex.nmol, ex.ntemp, ex.nlayers, ex.nwave))
  # Write arrays:
  f.write(struct.pack(str(ex.nmol)   +"i", *list(ex.molID)))
  f.write(struct.pack(str(ex.ntemp)  +"d", *list(ex.temp) ))
  f.write(struct.pack(str(ex.nlayers)+"d", *list(ex.press)))
  f.write(struct.pack(str(ex.nwave)  +"d", *list(ex.wn)   ))
  # Write extinction-coefficient data:
  fmt = str(ex.nmol * ex.ntemp * ex.nlayers * ex.nwave) + "d"
  f.write(struct.pack(fmt, *list(pyrat.ex.etable.flatten())))
  f.close()
  pt.msg(pyrat.verb-3, "Extinction-coefficient table written to file:"
                       " '{:s}'.".format(ex.extfile), pyrat.log, 2)


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

  pyrat.iso.iext = np.zeros(pyrat.iso.niso, np.int)
  # Get species indices in extinction-coefficient table for the isotopes:
  if pyrat.ex.extfile is not None:
    for i in np.arange(pyrat.iso.niso):
      pyrat.iso.iext[i] = np.where(pyrat.ex.molID ==
                                   pyrat.mol.ID[pyrat.iso.imol[i]])[0][0]

  # Calculate extinction-coefficient in C:
  pt.msg(pyrat.verb-4, "Calculating extinction at layer {:3d}/{:d} "
         "(T={:6.1f} K, p={:.1e} bar).".format(ilayer+1, pyrat.atm.nlayers,
                                        temp, pressure/pc.bar), pyrat.log, 2)
  logtext = " "*800
  ec.extinction(extcoeff,
                pyrat.voigt.profile, pyrat.voigt.size, pyrat.voigt.index,
                pyrat.voigt.lorentz, pyrat.voigt.doppler,
                pyrat.spec.wn, pyrat.spec.own, pyrat.spec.odivisors,
                density, molq, pyrat.mol.radius, pyrat.mol.mass,
                pyrat.iso.imol, pyrat.iso.mass, pyrat.iso.ratio,
                ziso, pyrat.iso.iext,
                pyrat.lt.wn, pyrat.lt.elow, pyrat.lt.gf, pyrat.lt.isoid,
                pyrat.ex.ethresh, pressure, temp,
                logtext, pyrat.verb, add)
  pyrat.log.write(logtext.rstrip()[:-1])

