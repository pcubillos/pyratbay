import numpy as np
import sys, os
import matplotlib.pyplot as plt

import ptools as pt
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

  Modification History:
  ---------------------
  2014-XX-XX  patricio  Initial implementation
  2015-01-19  patricio  Adapted to new file format.
  """
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
      pyrat.atmf.info = line[1:].strip()

    # Radius, pressure, and temperature units of atm file:
    elif line == '@PRESSURE':
      pyrat.atmf.punits = atmfile.readline().strip()
    elif line == '@RADIUS':
      pyrat.atmf.runits = atmfile.readline().strip()
    elif line == '@TEMPERATURE':
      pyrat.atmf.tunits = atmfile.readline().strip()

    # Abundance by mass or number:
    elif line == '@ABUNDANCE':
      pyrat.atmf.abundance = atmfile.readline().strip() == "mass"

    # Read in molecules:
    elif line == "@SPECIES":
      pyrat.mol.name = np.asarray(atmfile.readline().strip().split())
      pyrat.mol.nmol = len(pyrat.mol.name)

    else:
      pt.error("Atmosphere file has an unexpected line: \n'{:s}'".format(line))

  iline = atmfile.tell()  # Current line position

  pt.msg(1, "Molecules list: \n{:s}".format(str(pyrat.mol.name)), 2)
  if pyrat.atmf.abundance:
    pt.msg(1, "Abundance is given by mass (mass mixing ratio)", 2)
  else:
    pt.msg(1, "Abundance is given by number (volume mixing ratio)", 2)
  pt.msg(1, "Unit factors: radius: {:s},  pressure: {:s},  temperature: {:s}".
          format(pyrat.atmf.runits, pyrat.atmf.punits, pyrat.atmf.tunits), 2)
  pt.msg(1, "Atm file info: {:s}".format(pyrat.atmf.info), 2)
  pt.msg(1, "Data starting position: {:d}".format(iline), 2)
  return iline


def getconstants(pyrat):
  """
  Set molecules constant values (ID, mass, radius).

  Modification History:
  ---------------------
  2014-??-??  patricio  Initial implemetation
  2014-08-17  patricio  Added Documentation. Adapted to new molecules.dat file.
                        Store molecule symbol and ID.
  2015-01-19  patricio  Updated molecules.dat reading.
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
  Extract data from atmospheric file into pyrat object.

  Modification History:
  ---------------------
  2014-XX-XX  patricio  Initial implementation.
  2015-01-19  patricio  Adapted to new atmosphere file format.
  """
  # Read first line to count number of columns:
  datastart = atmfile.tell()
  line = atmfile.readline()
  ncolumns = len(line.split()) - pyrat.mol.nmol
  # Is the radius given in the atmosphere file?:
  rad = (ncolumns == 3)

  # Count number of layers:
  pyrat.atmf.nlayers = 1
  while True:
    line = atmfile.readline()
    if line == '' or line.startswith('#'):
      break
    pyrat.atmf.nlayers += 1
  pt.msg(pyrat.verb, "Number of layers in the atmospheric file: {:d}".
                     format(pyrat.atmf.nlayers), 2)

  # Initialize arrays:
  if rad:
    pyrat.atmf.radius = np.zeros( pyrat.atmf.nlayers)
  pyrat.atmf.press  = np.zeros( pyrat.atmf.nlayers)
  pyrat.atmf.temp   = np.zeros( pyrat.atmf.nlayers)
  pyrat.atmf.mm     = np.zeros( pyrat.atmf.nlayers)
  pyrat.atmf.q      = np.zeros((pyrat.atmf.nlayers, pyrat.mol.nmol))
  pyrat.atmf.d      = np.zeros((pyrat.atmf.nlayers, pyrat.mol.nmol))

  # Read table:
  nprofiles = pyrat.mol.nmol
  atmfile.seek(datastart, 0)
  for i in np.arange(pyrat.atmf.nlayers):
    data = atmfile.readline().split()
    if rad:
      pyrat.atmf.radius[i] = float(data[0])
    pyrat.atmf.press [i] = float(data[rad+0])
    pyrat.atmf.temp  [i] = float(data[rad+1])
    pyrat.atmf.q     [i] = np.asarray(data[rad+2:], float)

  # Sum of abundances per layer:
  sumq = np.sum(pyrat.atmf.q, axis=1)
  # FINDME: Add sumq != 1.0 warning.

  # Plot abundance profiles:
  if pyrat.verb >= 10:
    plt.figure(0)
    plt.clf()
    for i in np.arange(pyrat.mol.nmol):
      plt.loglog(pyrat.atmf.q[:,i], pyrat.atmf.press, label=pyrat.mol.name[i])
    plt.legend(loc="lower left")
    plt.ylim(pyrat.atmf.press[0], pyrat.atmf.press[-1])
    plt.ylabel("Pressure  ({:s})".format(pyrat.atmf.punits))
    plt.xlabel("Abundances")
    plt.draw()
    plt.savefig("atmprofiles.png")

  # Store values in CGS system of units:
  if rad:
    pyrat.atmf.radius *= pc.units[pyrat.atmf.runits]
  pyrat.atmf.press  *= pc.units[pyrat.atmf.punits]
  pyrat.atmf.temp   *= pc.units[pyrat.atmf.tunits]

  # Store the abundance as volume mixing ratio:
  if pyrat.atmf.abundance:
    pyrat.atmf.q = pyrat.atmf.q * pyrat.atmf.mm / pyrat.mol.mass

  # Calculate the mean molecular mass per layer:
  pyrat.atmf.mm =     np.sum(pyrat.atmf.q*pyrat.mol.mass, axis=1)
  pt.msg(pyrat.verb-10, "Mean molecular mass array: {:s}".
                        format(str(pyrat.atmf.mm)), 2)

  # Calculate density profiles for each molecule:
  for i in np.arange(pyrat.mol.nmol):
    pyrat.atmf.d[:,i] = IGLdensity(pyrat.atmf.q[:,i], pyrat.mol.mass[i],
                                  pyrat.atmf.press, pyrat.atmf.temp)


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

  Modification History:
  ---------------------
  2014-06-08  patricio  Initial implementation.
  """
  return (mass * pc.u) * abundance * pressure / (pc.k * temperature)
