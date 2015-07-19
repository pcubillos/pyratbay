import sys, os
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants   as sc
import scipy.integrate   as si
import scipy.interpolate as sip

import ptools     as pt
import pconstants as pc

def readatm(pyrat):
  """
  Read the atmospheric file, store variables in pyrat.
  """
  # Check atmfile:
  pt.msg(pyrat.verb, "\nReading atmosphere file: '{:s}'".format(pyrat.atmfile))
  atmfile = open(pyrat.atmfile, "r")

  # Read keywords:
  iline = getkeywords(pyrat, atmfile)

  # Read molecular constant values:
  getconstants(pyrat)

  # Read atm data table (radius, pressure, temperature, abundances):
  atmfile.seek(iline, 0)
  getprofiles(pyrat, atmfile)

  atmfile.close()
  pt.msg(pyrat.verb, "Done.", 0)


def getkeywords(pyrat, atmfile):
  """
  Read keyword variables from the atmfile.
  """
  # User-input atmospheric-data object:
  atm = pyrat.inputs.atm

  atmfile.seek(0)
  while True:
    line = atmfile.readline().strip()
    #print(line)

    # Stop when the per-layer data begins:
    if line == "@DATA":
      break

    # Skip empty and comment lines:
    elif line == '' or line.startswith('#'):
      pass

    # Atmospheric file info:
    elif line.startswith('n'):
      atm.info = line[1:].strip()

    # Radius, pressure, and temperature units of atm file:
    elif line == '@PRESSURE':
      atm.punits = atmfile.readline().strip()
    elif line == '@RADIUS':
      atm.runits = atmfile.readline().strip()
    elif line == '@TEMPERATURE':
      atm.tunits = atmfile.readline().strip()

    # Abundance by mass or number:
    elif line == '@ABUNDANCE':
      atm.abundance = atmfile.readline().strip() == "mass"

    # Read in molecules:
    elif line == "@SPECIES":
      pyrat.mol.name = np.asarray(atmfile.readline().strip().split())
      pyrat.mol.nmol = len(pyrat.mol.name)

    else:
      pt.error("Atmosphere file has an unexpected line: \n'{:s}'".format(line))

  iline = atmfile.tell()  # Current line position

  pt.msg(1, "Molecules list: \n{:s}".format(str(pyrat.mol.name)), 2)
  if atm.abundance:
    pt.msg(1, "Abundance is given by mass (mass mixing ratio)", 2)
  else:
    pt.msg(1, "Abundance is given by number (volume mixing ratio)", 2)
  pt.msg(1, "Unit factors: radius: {:s},  pressure: {:s},  temperature: {:s}".
          format(atm.runits, atm.punits, atm.tunits), 2)
  pt.msg(1, "Atm file info: {:s}".format(atm.info), 2)
  pt.msg(1, "Data starting position: {:d}".format(iline), 2)
  return iline


def getconstants(pyrat):
  """
  Set molecules constant values (ID, mass, radius).
  """
  # Read file with molecular info:
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
             "file: {:s}.\n".format(str(absent), pyrat.molfile))

  # Set molecule's values:
  pyrat.mol.ID     = np.zeros(pyrat.mol.nmol, np.int)
  pyrat.mol.symbol = np.zeros(pyrat.mol.nmol, "|S15")
  pyrat.mol.mass   = np.zeros(pyrat.mol.nmol)
  pyrat.mol.radius = np.zeros(pyrat.mol.nmol)

  pt.msg(pyrat.verb, "Molecule   ID   Radius  Mass\n"
                     "                (A)     (gr/mol)", 4)
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
    pyrat.mol.radius[i] = diam[imol]/2.0 * pc.units["A"]
    pt.msg(1, "{:>10s}:  {:3d}  {:.3f}  {:8.4f}".
           format(pyrat.mol.name[i], pyrat.mol.ID[i],
                  pyrat.mol.radius[i]/pc.units["A"], pyrat.mol.mass[i]), 2)


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
  pt.msg(pyrat.verb, "Number of layers in the input atmospheric file: {:d}".
                     format(atm.nlayers), 2)

  # Initialize arrays:
  if rad:
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
    atm.radius *= pc.units[atm.runits]
  atm.press  *= pc.units[atm.punits]
  atm.temp   *= pc.units[atm.tunits]

  # Store the abundance as volume mixing ratio:
  if atm.abundance:
    atm.q = atm.q * atm.mm / pyrat.mol.mass

  # Calculate the mean molecular mass per layer:
  atm.mm =     np.sum(atm.q*pyrat.mol.mass, axis=1)
  pt.msg(pyrat.verb-10, "Mean molecular mass array: {:s}".
                         format(str(atm.mm)), 2)

  # Calculate density profiles for each molecule:
  for i in np.arange(pyrat.mol.nmol):
    atm.d[:,i] = IGLdensity(atm.q[:,i], pyrat.mol.mass[i], atm.press, atm.temp)


def reloadatm(pyrat, temp, abund):
  """
  Parameters:
  -----------
  pyrat: A Pyrat instance
  temp: 1D float ndarray
     Array with a temperature profile in Kelvin (from top to bottom layer).
  abund: 2D float ndarray
     Array with the species mole mixing ratio profiles [nlayers, nmol].

  Notes:
  ------
  This code assumes that the input temperature and abundances correspond
  to the final sampling of the atmospheric layers (after makeradius).
  """
  # Check that the dimensions match:
  if np.size(temp) != np.size(pyrat.atm.temp):
    pt.error("The temperature array size ({:d}) doesn't match the Pyrat's "
             "temperature size ({:d}).".format(np.size(temp),
                                               np.size(pyrat.atm.temp)))
  if np.shape(abund) != np.shape(pyrat.atm.q):
    pt.error("The shape of the abundances array {:s} doesn't match the "
             "shape of the Pyrat's abundance size {:s}".format(
              str(np.shape(abund)), str(np.shape(pyrat.atm.q))))

  # Put temperature and abundance data into the Pyrat object:
  pyrat.atm.temp = temp
  pyrat.atm.q    = abund

  # Mean molecular mass:
  pyrat.atm.mm = np.sum(pyrat.atm.q*pyrat.mol.mass, axis=1)

  # Density:
  for i in np.arange(pyrat.mol.nmol):
    pyrat.atm.d[:,i] = IGLdensity(pyrat.atm.q[:,i], pyrat.mol.mass[i],
                                  pyrat.atm.press,  pyrat.atm.temp)

  # Radius:
  pyrat.atm.radius = hydro_equilibrium(pyrat.atm.press,    pyrat.atm.temp,
                                       pyrat.atm.mm,       pyrat.surfgravity,
                                       pyrat.pressurebase, pyrat.radiusbase)

  # Partition function:
  for i in np.arange(pyrat.lt.ndb):           # For each Database
    for j in np.arange(pyrat.lt.db[i].niso):  # For each isotope in DB
      zinterp = sip.interp1d(pyrat.lt.db[i].temp, pyrat.lt.db[i].z[j],
                             kind='slinear')
      pyrat.iso.z[pyrat.lt.db[i].iiso+j] = zinterp(pyrat.atm.temp)


def IGLdensity(abundance, mass, pressure, temperature):
  """
  Use the Ideal gas law to calculate the density.

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
  """
  return (mass * pc.u) * abundance * pressure / (pc.k * temperature)


def hydro_equilibrium(pressure, temperature, mu, g, p0=None, r0=None):
  """
  Calculate radii using the hydrostatic-equilibrium equation.

  Parameters:
  -----------
  pressure: 1D float ndarray
     Atmospheric pressure for each layer (in barye).
  temperature: 1D float ndarray
     Atmospheric temperature for each layer (in K).
  mu: 1D float ndarray
     Mean molecular mass for each layer (in g mol-1).
  g: Float
     Atmospheric gravity (in cm s-2).
  p0: Float
     Reference pressure level where radius(p0) = r0.
  r0: Flaot
     Reference radius level corresponding to p0.

  Returns:
  --------
  radius: 1D float ndarray
     Radius for each layer (in cm).

  Notes:
  ------
  If the reference values (p0 and r0) are not given, set radius = 0.0
  at the bottom of the matmosphere.
  """
  # Apply the HE equation:
  radius = si.cumtrapz((-pc.k*sc.N_A * temperature) / (mu*g), np.log(pressure))
  radius = np.concatenate(([0.0], radius))

  # Set absolute radii values if p0 and r0 are provided:
  if p0 is not None and r0 is not None:
    # Find current radius at p0:
    radinterp = sip.interp1d(pressure, radius, kind='slinear')
    r0_interp = radinterp(p0)
    # Set: radius(p0) = r0
    radius += r0 - r0_interp
  # Set radius = 0 at the bottom of the atmosphere:
  else:
    radius -= radius[-1]

  return radius
