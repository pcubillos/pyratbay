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
  pt.msg(pyrat.verb, "\nReading atmosphere file: '%s'"%pyrat.atmfile, 0)
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
  # Initialize molecules object:

  atmfile.seek(0)
  while True:
    iline = atmfile.tell()  # Current line position
    line = atmfile.readline().strip()
    #print(line)
    # Skip empty and comment lines:
    if line == '' or line.startswith('#'):
      pass

    # Abundance by mass or number:
    elif line.startswith('q'):
      pyrat.atmf.abundance = (line[1:].strip() == 'mass')

    # Atmospheric file info:
    elif line.startswith('n'):
      pyrat.atmf.info = line[1:].strip()

    # Radius offset level:
    elif line.startswith('z'):
      pyrat.atmf.roffset = float(line[1:].strip())

    # Radius, pressure, and temperature units of atm file:
    elif line.startswith('u'):
      if line[1] == 'r':
        pyrat.atmf.runits = line[2:].strip()
      if line[1] == 'p':
        pyrat.atmf.punits = line[2:].strip()
      if line[1] == 't':
        pyrat.atmf.tunits = line[2:].strip()

    # Read in molecules with abundance profiles:
    elif line.startswith('i'):
      pyrat.mol.name = np.asarray(line[1:].strip().split())
      pyrat.mol.nmol = len(pyrat.mol.name)

    # Read in molecules without abundance profiles:
    elif line.startswith('f'):
      remnmol, remnfactor = line[1:].strip().split()
      pyrat.mol.nmol += 1
      pyrat.mol.name = np.concatenate((pyrat.mol.name, [remnmol])) 
      pyrat.atmf.remainder = np.concatenate((pyrat.atmf.remainder,
                                            [float(remnfactor)])) 
    else:
      break

  pt.msg(1, "Molecules list: \n%s"%(str(pyrat.mol.name)), 2)
  if pyrat.atmf.abundance:
    pt.msg(1, "Abundance is by mass (mass mixing ratio)", 2)
  else:
    pt.msg(1, "Abundance is by number (volume mixing ratio)", 2)
  pt.msg(1, "Unit factors: radius: %s,  pressure: %s,  temperature: %s"%(
         pyrat.atmf.runits, pyrat.atmf.punits, pyrat.atmf.tunits), 2)
  pt.msg(1, "Atm file info: %s"%pyrat.atmf.info, 2)
  pt.msg(1, "Radius offset: %.2f %s"%(pyrat.atmf.roffset, pyrat.atmf.runits), 2)
  pt.msg(1, "Data starting position: %s"%iline, 2)
  return iline


def getconstants(pyrat):
  """
  Set molecules constant values (mass, radius).

  Modification History:
  ---------------------
  2014-??-??  patricio  Initial implemetation
  2014-08-17  patricio  Added Documentation. Adapted to new molecules.dat file.
                        Store molecule symbol and ID.
  """
  molfile = open("../inputs/molecules.dat", "r")
 
  # Read Molecular name aliases:
  alias = [] # Molecule alias
  amol  = [] # Molecule name
  line = molfile.readline().strip()
  while line == '' or line.startswith('#'):
    line = molfile.readline().strip()
  while line != '' and not line.startswith('#'):
    alias.append(line.split()[0]) # Alias name
    amol.append(line.split()[1])  # Molecule name
    line = molfile.readline().strip()
  alias = np.asarray(alias)

  # Read Molecular values:
  molID  = [] # Molecule ID
  symbol = [] # Molecule symbol
  mass   = [] # Molecule mass
  diam   = [] # Molecule diameter
  while line == '' or line.startswith('#'):
    line = molfile.readline().strip()
  while line != '' and not line.startswith('#'):
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

  # Set molecule's values:
  pyrat.mol.ID     = np.zeros(pyrat.mol.nmol, np.int)
  pyrat.mol.symbol = np.zeros(pyrat.mol.nmol, "|S15")
  pyrat.mol.mass   = np.zeros(pyrat.mol.nmol)
  pyrat.mol.radius = np.zeros(pyrat.mol.nmol)

  pt.msg(pyrat.verb, "Molecule   ID   Radius  Mass\n"
                     "                (A)     (gr/mol)", 4)
  for i in np.arange(pyrat.mol.nmol):
    # Check if name is an alias:
    if pyrat.mol.name[i] in alias:
      molname = amol[np.where(alias == pyrat.mol.name[i])[0]]
    else:
      molname = pyrat.mol.name[i]

    # Find the molecule in the list:
    imol = np.where(symbol == molname)[0]
    # Set molecule ID:
    pyrat.mol.ID[i]     = molID [imol]
    # Set molecule symbol:
    pyrat.mol.symbol[i] = symbol[imol][0]
    # Set molecule mass:
    pyrat.mol.mass[i]   = mass  [imol]
    # Set molecule radius:
    pyrat.mol.radius[i] = diam  [imol]/2.0
    pt.msg(1, "{:>10s}:  {:3d}  {:.3f}  {:8.4f}".format(pyrat.mol.name[i],
                  pyrat.mol.ID[i], pyrat.mol.radius[i], pyrat.mol.mass[i]), 2)


def getprofiles(pyrat, atmfile):
  """
  Extract data from atmospheric file into pyrat object.
  """
  datastart = atmfile.tell()
  # Count number of layers:
  while (True):
    line = atmfile.readline()
    if line == '' or line.startswith('#'):
      break
    pyrat.atmf.layers += 1
  pt.msg(pyrat.verb, "Number of layers in Atm. file: %d"%pyrat.atmf.layers, 2)

  # Initialize arrays:
  pyrat.atmf.rad   = np.zeros(pyrat.atmf.layers)
  pyrat.atmf.press = np.zeros(pyrat.atmf.layers)
  pyrat.atmf.temp  = np.zeros(pyrat.atmf.layers)
  pyrat.atmf.mm    = np.zeros(pyrat.atmf.layers)
  pyrat.atmf.q     = np.zeros((pyrat.atmf.layers, pyrat.mol.nmol))
  pyrat.atmf.d     = np.zeros((pyrat.atmf.layers, pyrat.mol.nmol))

  # Read table:
  nprofiles = pyrat.mol.nmol - len(pyrat.atmf.remainder)
  atmfile.seek(datastart, 0)
  for i in np.arange(pyrat.atmf.layers):
    data = atmfile.readline().split()
    pyrat.atmf.rad  [i] = float(data[0])
    pyrat.atmf.press[i] = float(data[1])
    pyrat.atmf.temp [i] = float(data[2])
    pyrat.atmf.q    [i,0:nprofiles] = np.asarray(data[3:], float)

  # abundance of molecules with the remainder:
  sumq = np.sum(pyrat.atmf.q, axis=1)
  for i in np.arange(len(pyrat.atmf.remainder)):
    pyrat.atmf.q[:,nprofiles+i] = pyrat.atmf.remainder[i]*(1.0-sumq)

  # Plot abundance profiles:
  if pyrat.verb >= 10:
    plt.figure(0)
    plt.clf()
    for i in np.arange(pyrat.mol.nmol):
      plt.loglog(pyrat.atmf.q[:,i], pyrat.atmf.press, label=pyrat.mol.name[i])
    plt.legend(loc="lower left")
    plt.ylim(pyrat.atmf.press[0], pyrat.atmf.press[-1])
    plt.ylabel("Pressure  (%s)"%pyrat.atmf.punits)
    plt.xlabel("Abundances")
    plt.draw()
    plt.savefig("atmprofiles.png")

  # Store values in CGS system of units:
  pyrat.atmf.rad    = ((pyrat.atmf.rad + pyrat.atmf.roffset) * 
                                                   pc.units[pyrat.atmf.runits])
  pyrat.atmf.press *= pc.units[pyrat.atmf.punits]
  pyrat.atmf.temp  *= pc.units[pyrat.atmf.tunits]

  # Calculate the mean molecular mass per layer:
  if pyrat.atmf.abundance:  # Abundance by mass:
    pyrat.atmf.mm = 1.0/np.sum(pyrat.atmf.q/pyrat.mol.mass, axis=1)
  else:                    # Abundance by number:
    pyrat.atmf.mm =     np.sum(pyrat.atmf.q*pyrat.mol.mass, axis=1)

  pt.msg(pyrat.verb-10, "Mean molecular mass: %s"%str(pyrat.atmf.mm), 2)
  # Store the abundance as volume mixing ratio:
  if pyrat.atmf.abundance:
    pyrat.atmf.q = pyrat.atmf.q * pyrat.atmf.mm / pyrat.mol.mass

  # Calculate density profiles for each molecule:
  for i in np.arange(pyrat.mol.nmol):
    pyrat.atmf.d[:,i] = IGLdensity(pyrat.atmf.q[:,i], pyrat.mol.mass[i],
                                  pyrat.atmf.press, pyrat.atmf.temp)


def IGLdensity(abundance, mass, pressure, temperature):
  """
  Use the Ideal gas law to calculate pressure

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
  return (mass * pc.u) * abundance  * pressure / (pc.k * temperature) 
