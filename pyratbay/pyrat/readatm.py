# Copyright (c) 2016-2017 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants   as sc
import scipy.integrate   as si
import scipy.interpolate as sip

from .. import tools      as pt
from .. import constants  as pc

def readatm(pyrat):
  """
  Read the atmospheric file, store variables in pyrat.
  """
  # Check atmfile:
  pt.msg(pyrat.verb-3, "\nReading atmospheric file: '{:s}'.".
                        format(pyrat.atmfile), pyrat.log)
  atmfile = open(pyrat.atmfile, "r")

  # Read keywords:
  iline = getkeywords(pyrat, atmfile)

  # Read molecular constant values:
  getconstants(pyrat)

  # Read atm data table (radius, pressure, temperature, abundances):
  atmfile.seek(iline, 0)
  getprofiles(pyrat, atmfile)

  atmfile.close()
  pt.msg(pyrat.verb-3, "Done.", pyrat.log)


def getkeywords(pyrat, atmfile):
  """
  Read keyword variables from the atmfile.
  """
  # User-input atmospheric-data object:
  atm = pyrat.inputs.atm

  atmfile.seek(0)
  while True:
    line = atmfile.readline().strip()

    # Stop when the per-layer data begins:
    if line == "@DATA":
      break

    # Skip empty and comment lines:
    elif line == '' or line.startswith('#'):
      pass

    # Radius, pressure, temperature, and abundance units of atm file:
    elif line == '@PRESSURE':
      atm.punits = atmfile.readline().strip()
    elif line == '@RADIUS':
      atm.runits = atmfile.readline().strip()
    elif line == '@TEMPERATURE':
      atm.tunits = atmfile.readline().strip()
    # Abundance by mass or number:
    elif line == '@ABUNDANCE':
      atm.qunits = atmfile.readline().strip()

    # Read in molecules:
    elif line == "@SPECIES":
      pyrat.mol.name = np.asarray(atmfile.readline().strip().split())
      pyrat.mol.nmol = len(pyrat.mol.name)

    else:
      pt.error("Atmosphere file has an unexpected line: \n'{:s}'".format(line),
               pyrat.log)

  iline = atmfile.tell()  # Current line position
  atm.qunits = pt.defaultp(atm.qunits, 'number',
       "Undefined abundance units in the input atmospheric file.  Assumed to "
       "be '{:s}'.", pyrat.wlog, pyrat.log)
  atm.punits = pt.defaultp(atm.punits, 'barye',
       "Undefined pressure units in the input atmospheric file.  Assumed to "
       "be '{:s}'.", pyrat.wlog, pyrat.log)

  pt.msg(pyrat.verb-4, "Species list: \n  {:s}".
                       format(str(pyrat.mol.name)),  pyrat.log, 2, si=4)
  if atm.qunits == "number":
    txt = "volume"
  else:
    txt = "mass"
  pt.msg(pyrat.verb-4, "Abundances are given by {:s} ({:s} mixing ratio).".
           format(atm.qunits, txt), pyrat.log, 2)
  pt.msg(pyrat.verb-4, "Unit factors: radius: {:s},  pressure: {:s},  "
    "temperature: {:s}".format(atm.runits, atm.punits, atm.tunits), pyrat.log,2)
  pt.msg(pyrat.verb-6, "Data starting position {:d}".format(iline), pyrat.log,2)
  return iline


def getconstants(pyrat):
  """
  Set molecules constant values (ID, mass, radius).
  """
  # Read file with molecular info:
  pt.msg(pyrat.verb-4, "Taking species constant parameters from: '{:s}'.".
         format(pyrat.molfile), pyrat.log, 2)
  molfile = open(pyrat.molfile, "r")
 
  # Skip comment and blank lines:
  line = molfile.readline().strip()
  while line == '' or line.startswith('#'):
    line = molfile.readline().strip()

  molID  = [] # Molecule ID
  symbol = [] # Molecule symbol
  mass   = [] # Molecule mass
  diam   = [] # Molecule diameter
  # Read Molecular values:
  while line != '' and not line.startswith('#'):  # Start reading species
    molinfo = line.split()
    # Extract info:
    molID .append(  int(molinfo[0]))
    symbol.append(      molinfo[1] )
    mass  .append(float(molinfo[2]))
    diam  .append(float(molinfo[3]))
    line = molfile.readline().strip()  # Read next line

  molID  = np.asarray(molID)
  symbol = np.asarray(symbol)
  mass   = np.asarray(mass)
  diam   = np.asarray(diam)

  # Check that all atmospheric species are listed in molfile:
  absent = np.setdiff1d(pyrat.mol.name, symbol)
  if len(absent) > 0:
    pt.error("These species: {:s} are not listed in the molecules info "
             "file: {:s}.\n".format(str(absent), pyrat.molfile), pyrat.log)

  # Set molecule's values:
  pyrat.mol.ID     = np.zeros(pyrat.mol.nmol, np.int)
  pyrat.mol.symbol = np.zeros(pyrat.mol.nmol, "|S15")
  pyrat.mol.mass   = np.zeros(pyrat.mol.nmol)
  pyrat.mol.radius = np.zeros(pyrat.mol.nmol)

  pt.msg(pyrat.verb-4, "Molecule   ID   Radius  Mass\n"
                       "                (A)     (gr/mol)", pyrat.log, 4)
  for i in np.arange(pyrat.mol.nmol):
    # Find the molecule in the list:
    imol = np.where(symbol == pyrat.mol.name[i])[0]
    # Set molecule ID:
    pyrat.mol.ID[i]     = molID [imol]
    # Set molecule symbol:
    pyrat.mol.symbol[i] = symbol[imol][0]
    # Set molecule mass:
    pyrat.mol.mass[i]   = mass  [imol]
    # Set molecule collision radius:
    pyrat.mol.radius[i] = diam[imol]/2.0 * pc.A
    pt.msg(pyrat.verb-4, "{:>10s}:  {:3d}  {:.3f}  {:8.4f}".
           format(pyrat.mol.name[i], pyrat.mol.ID[i],
                  pyrat.mol.radius[i]/pc.A, pyrat.mol.mass[i]), pyrat.log, 2)


def getprofiles(pyrat, atmfile):
  """
  Extract pressure, temperature, and mixing ratio data from the atmospheric
  file.
  """
  # User-input atmospheric-data object:
  atm = pyrat.inputs.atm

  # Read first line to count number of columns:
  datastart = atmfile.tell()
  line = atmfile.readline()
  ncolumns = len(line.split()) - pyrat.mol.nmol
  # Is the radius given in the atmosphere file?:
  rad = (ncolumns == 3)

  # Count number of layers:
  atm.nlayers = 1
  while True:
    line = atmfile.readline()
    if line == '' or line.startswith('#'):
      break
    atm.nlayers += 1
  pt.msg(pyrat.verb-4, "Number of layers in the input atmospheric file: {:d}".
                     format(atm.nlayers), pyrat.log, 2)

  # Initialize arrays:
  if rad:
    atm.runits = pt.defaultp(atm.runits, 'cm',
       "Undefined radius units in the input atmospheric file.  Assumed to "
       "be '{:s}'.", pyrat.wlog, pyrat.log)
    atm.radius = np.zeros( atm.nlayers)
  atm.press    = np.zeros( atm.nlayers)
  atm.temp     = np.zeros( atm.nlayers)
  atm.mm       = np.zeros( atm.nlayers)
  atm.q        = np.zeros((atm.nlayers, pyrat.mol.nmol))
  atm.d        = np.zeros((atm.nlayers, pyrat.mol.nmol))

  # Read table:
  nprofiles = pyrat.mol.nmol
  atmfile.seek(datastart, 0)
  for i in np.arange(atm.nlayers):
    data = atmfile.readline().split()
    if rad:
      atm.radius[i] = float(data[0])
    atm.press [i] = float(data[rad+0])
    atm.temp  [i] = float(data[rad+1])
    atm.q     [i] = np.asarray(data[rad+2:], float)

  # Sum of abundances per layer:
  sumq = np.sum(atm.q, axis=1)
  # FINDME: Add sumq != 1.0 warning.

  # Plot abundance profiles:
  if pyrat.verb >= 10:
    plt.figure(0)
    plt.clf()
    for i in np.arange(pyrat.mol.nmol):
      plt.loglog(atm.q[:,i], atm.press, label=pyrat.mol.name[i])
    plt.legend(loc="lower left")
    plt.ylim(atm.press[0], atm.press[-1])
    plt.ylabel("Pressure  ({:s})".format(atm.punits))
    plt.xlabel("Abundances")
    plt.draw()
    plt.savefig("atmprofiles.png")

  # Store values in CGS system of units:
  if rad:
    atm.radius *= pt.u(atm.runits)
  atm.press  *= pt.u(atm.punits)
  atm.temp   *= pt.u(atm.tunits)

  # Store the abundance as volume mixing ratio:
  if atm.qunits == "mass":
    atm.q = atm.q * atm.mm / pyrat.mol.mass

  # Calculate the mean molecular mass per layer:
  atm.mm = np.sum(atm.q*pyrat.mol.mass, axis=1)
  pt.msg(pyrat.verb-4, "Typical mean molecular mass: {:.3f} g mol-1.".
                        format(np.median(atm.mm)), pyrat.log, 2)

  # Calculate number density profiles for each molecule (in molecules cm-3):
  for i in np.arange(pyrat.mol.nmol):
    atm.d[:,i] = IGLdensity(atm.q[:,i], pyrat.mol.mass[i], atm.press, atm.temp)


def reloadatm(pyrat, temp, abund, radius=None):
  """
  Parameters:
  -----------
  pyrat: A Pyrat instance
  temp: 1D float ndarray
     Layer's temperature profile (in Kelvin) sorted from top to bottom.
  abund: 2D float ndarray
     Species mole mixing ratio profiles [nlayers, nmol].
  radius: 1D float ndarray
     Layer's altitude profile (in cm), same order as temp.
  """
  # Check that the dimensions match:
  if np.size(temp) != np.size(pyrat.atm.temp):
    pt.error("The temperature array size ({:d}) doesn't match the Pyrat's "
             "temperature size ({:d}).".
              format(np.size(temp), np.size(pyrat.atm.temp)), pyrat.log)
  if np.shape(abund) != np.shape(pyrat.atm.q):
    pt.error("The shape of the abundances array {:s} doesn't match the "
             "shape of the Pyrat's abundance size {:s}".
           format(str(np.shape(abund)), str(np.shape(pyrat.atm.q))), pyrat.log)

  # Check temperature boundaries:
  errorlog = ("One or more input temperature values lies out of the {:s} "
    "temperature boundaries (K): [{:6.1f}, {:6.1f}].\nHalted calculation.")
  if pyrat.ex.extfile is not None:
    if np.any(temp > pyrat.ex.tmax) or np.any(temp < pyrat.ex.tmin):
      pt.warning(pyrat.verb-1, errorlog.format("tabulated EC",
                                        pyrat.ex.tmin, pyrat.ex.tmax))
      return 0
  else:
    if (pyrat.lt.ntransitions > 0 and
        (np.any(temp > pyrat.lt.tmax) or np.any(temp < pyrat.lt.tmin))):
      pt.warning(pyrat.verb-1, errorlog.format("line-transition",
                                         pyrat.lt.tmin, pyrat.lt.tmax))
      return 0
    if (pyrat.cs.nfiles > 0 and
        (np.any(temp > pyrat.cs.tmax) or np.any(temp < pyrat.cs.tmin))):
      pt.warning(pyrat.verb-1, errorlog.format("cross-section",
                                         pyrat.cs.tmin, pyrat.cs.tmax))
      return 0

  # Put temperature and abundance data into the Pyrat object:
  pyrat.atm.temp = temp
  pyrat.atm.q    = abund

  # Mean molecular mass:
  pyrat.atm.mm = np.sum(pyrat.atm.q*pyrat.mol.mass, axis=1)

  # Number density (molecules cm-3):
  for i in np.arange(pyrat.mol.nmol):
    pyrat.atm.d[:,i] = IGLdensity(pyrat.atm.q[:,i], pyrat.mol.mass[i],
                                  pyrat.atm.press,  pyrat.atm.temp)

  # Take radius if provided, else use hydrostatic-equilibrium equation:
  if radius is not None:
    pyrat.atm.radius = radius
  else:
    # Compute gplanet from mass and radius if necessary/possible:
    if (pyrat.phy.gplanet is None and
      pyrat.phy.mplanet is not None and pyrat.phy.rplanet is not None):
      pyrat.phy.gplanet = pc.G * pyrat.phy.mplanet / pyrat.phy.rplanet**2

    # Check that the gravity variable is exists:
    if pyrat.phy.rplanet is None:
      pt.error("Undefined reference planetary radius (rplanet). Either provide "
        "the radius profile for the layers or the rplanet.", pyrat.log)
    if pyrat.phy.gplanet is None:
      pt.error("Undefined atmospheric gravity (gplanet).  Either provide "
        "the radius profile for the layers, the surface gravity, or the "
        "planetary mass (mplanet).", pyrat.log)
    if pyrat.refpressure is None:
      pt.error("Undefined reference pressure level (refpressure). Either "
         "provide the radius profile for the layers or refpressure.", pyrat.log)
    pyrat.atm.radius = pyrat.hydro(pyrat.atm.press, pyrat.atm.temp,
                          pyrat.atm.mm, pyrat.phy.gplanet, pyrat.phy.mplanet,
                          pyrat.refpressure, pyrat.phy.rplanet)
  # Check radii lie within Hill radius:
  rtop = np.where(pyrat.atm.radius > pyrat.phy.rhill)[0]
  if np.size(rtop) > 0:
    pyrat.atm.rtop = rtop[-1] + 1
    pt.warning(pyrat.verb-2, "The atmospheric pressure array extends "
       "beyond the Hill radius ({:.3e} km) at pressure {:.3e} bar (layer "
       "#{:d}).  Extinction beyond this layer will be neglected.".
        format(pyrat.phy.rhill/pc.km, pyrat.atm.press[pyrat.atm.rtop]/pc.bar,
               pyrat.atm.rtop), pyrat.log, pyrat.wlog)
  else:
    pyrat.atm.rtop = 0

  # Partition function:
  for i in np.arange(pyrat.lt.ndb):           # For each Database
    for j in np.arange(pyrat.lt.db[i].niso):  # For each isotope in DB
      zinterp = sip.interp1d(pyrat.lt.db[i].temp, pyrat.lt.db[i].z[j],
                             kind='slinear')
      pyrat.iso.z[pyrat.lt.db[i].iiso+j] = zinterp(pyrat.atm.temp)


def IGLdensity(abundance, mass, pressure, temperature):
  """
  Use the Ideal gas law to calculate the density in molecules cm-3.

  Parameters:
  -----------
  abundance: 1D ndarray
    Species volume mixing ratio profile.
  mass: Float
    Species mass (in gr/mol).
  pressure: 1D ndarray
    Atmospheric pressure profile (in barye units).
  temperature: 1D ndarray
    Atmospheric temperature (in kelvin).

  Returns:
  --------
  density: 1D float ndarray
     Atmospheric density in molecules per centimeter^3.
  """
  return abundance * pressure / (pc.k * temperature)
