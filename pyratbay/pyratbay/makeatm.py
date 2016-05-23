__all__ = ["writeatm", "readatm", "uniform", "makeatomic", "readatomic",
           "makepreatm", "TEA2pyrat"]

import os
import numpy as np

from .. import tools as pt

# Get Pyrat-Bay inputs dir:
thisdir = os.path.dirname(os.path.realpath(__file__))
indir   = thisdir + "/../../inputs/"


def writeatm(atmfile, pressure, temperature, species, abundances,
             punits, header):
  """
  Write an atmospheric file following the Pyrat format.

  Parameters
  ----------
  atmfile: String
     Name of output atmospheric file.
  pressure: 1D float ndarray
     Monotonously decreasing pressure profile (in barye).
  temperature: 1D float ndarray
     Temperature profile for pressure layers (in Kelvin).
  species: 1D string ndarray
     List of atmospheric species.
  abundances: 1D float ndarray
     The species mole mixing ratio.
  punits:  String
     Pressure units of output.
  header:  String
     Header message (comment) to include at the top of the file.
  """
  # Open file for writing:
  f = open(atmfile, "w")

  f.write(header)

  # Set the values units:
  f.write("# Abundance units (by number or mass):\n@ABUNDANCE\nnumber\n")
  f.write("# Pressure units:\n@PRESSURE\n{:s}\n".format(punits))
  f.write("# Temperatures are always Kelvin.\n\n")

  # Write the species names:
  f.write("# Atmospheric composition:\n"
          "@SPECIES\n" +
          "  ".join(["{:<s}".format(mol) for mol in species]) + '\n\n')

  # Write the per-layer data:
  f.write("# Pressure  Temperature  " +
          "".join(["{:<14s}".format(mol) for mol in species]) + "\n")
  f.write("@DATA\n")

  pressure = pressure/pt.u(punits)
  # Write data for each layer:
  nlayers = len(pressure)
  for i in np.arange(nlayers):
    # Pressure, and Temperature:
    f.write("{:10.4e}  {:11.3f}  ".format(pressure[i], temperature[i]))
    # Species mole mixing ratios:
    f.write("  ".join(["{:12.6e}".format(ab) for ab in abundances[i]]) + "\n")
  f.close()


def readatm(atmfile):
  """
  Read a Pyrat atmospheric file.
  """

  atmfile = open(atmfile, "r")
  while True:
    line = atmfile.readline().strip()

    # Stop when the per-layer data begins:
    if line == "@DATA":
      break

    # Skip empty and comment lines:
    elif line == '' or line.startswith('#'):
      pass

    # Radius, pressure, and temperature units of atm file:
    elif line == '@PRESSURE':
      print("Pressure units: {:s}".   format(atmfile.readline().strip()))
    elif line == '@RADIUS':
      print("Radius units: {:s}".     format(atmfile.readline().strip()))
    elif line == '@TEMPERATURE':
      print("Temperature units: {:s}".format(atmfile.readline().strip()))
    # Abundance by mass or number:
    elif line == '@ABUNDANCE':
      print("Abundance units: {:s}".  format(atmfile.readline().strip()))
    # Read in molecules:
    elif line == "@SPECIES":
      species = np.asarray(atmfile.readline().strip().split())
      nspecies = len(species)
    else:
      print("Atmosphere file has an unexpected line: \n'{:s}'".format(line))

  # Read first line to count number of columns:
  datastart = atmfile.tell()
  line = atmfile.readline()
  # Is there a column for the radius:
  rad = (len(line.split()) - nspecies == 3)

  # Count number of layers:
  nlayers = 1
  while True:
    line = atmfile.readline()
    if line == '' or line.startswith('#'):
      break
    nlayers += 1

  # Initialize arrays:
  if rad:
    radius = np.zeros(nlayers)
  press    = np.zeros(nlayers)
  temp     = np.zeros(nlayers)
  mm       = np.zeros(nlayers)
  q        = np.zeros((nlayers, nspecies))

  # Read table:
  nprofiles = nspecies
  atmfile.seek(datastart, 0)
  for i in np.arange(nlayers):
    data = atmfile.readline().split()
    if rad:
      radius[i] = float(data[0])
    press [i] = float(data[rad+0])
    temp  [i] = float(data[rad+1])
    q     [i] = np.asarray(data[rad+2:], float)

  if rad:
    return species, press, temp, q, radius
  return species, press, temp, q


def uniform(atmfile, pressure, temperature, species, abundances, punits="bar"):
  """
  Generate an atmospheric file with uniform abundances.

  Parameters
  ----------
  atmfile: String
     Name of output atmospheric file.
  pressure: 1D float ndarray
     Monotonously decreasing pressure profile (in punits).
  temperature: 1D float ndarray
     Temperature profile for pressure layers (in Kelvin).
  species: 1D string ndarray
     List of atmospheric species.
  abundances: 1D float ndarray
     The species mole mixing ratio.
  punits:  String
     Pressure units.
  """

  # Safety checks:
  nlayers = len(pressure)

  if len(temperature) != nlayers:
    pt.error("Pressure array length ({:d}) and Temperature array length "
             "({:d}) don't match.".format(nlayers, len(temperature)))
  if len(species) != len(abundances):
    pt.error("Species array length ({:d}) and Abundances array length ({:d}) "
             "don't match.".format(len(species), len(abundances)))
  # FINDME: Check pressure array is monotonously decreasing.

  # Expand abundances to 2D array:
  abund = np.tile(abundances, (nlayers,1))
  # File header:
  header = ("# This is an atmospheric file with pressure, temperature, \n"
            "# and uniform mole mixing ratio profiles.\n\n")

  writeatm(atmfile, pressure, temperature, species, abund, punits, header)


def makeatomic(solar, afile, xsolar=1.0, swap=None):
  """
  Generate an (atomic) elemental-abundances file by scaling the
  solar abundances from Asplund et al. (2009).
  http://adsabs.harvard.edu/abs/2009ARA%26A..47..481A

  Parameters
  ----------
  solar: String
     Input solar abundances filename.
  afile: String
     Output filename of modified atomic abundances.
  xsolar: Integer
     Multiplication factor for metal abundances (everything
     except H and He).
  swap: 2-element string tuple
      Swap the abundances of the given elements.

  Uncredited developers
  ---------------------
  Jasmina Blecic
  """
  # Read the Asplund et al. (2009) solar elementa abundances:
  index, symbol, dex, name, mass = readatomic(solar)
  # Count the number of elements:
  nelements = len(symbol)

  # Scale the metals aundances:
  imetals = np.where((symbol != "H") & (symbol != "He"))
  dex[imetals] += np.log10(xsolar)

  # Swap abundances if requested:
  if swap is not None:
    dex0 = dex[np.where(symbol == swap[0])]
    dex[np.where(symbol == swap[0])] = dex[np.where(symbol == swap[1])]
    dex[np.where(symbol == swap[1])] = dex0

  # Save data to file
  with open(afile, "w") as f:
    # Write header
    f.write("# Elemental abundances:\n"
            "# atomic number, symbol, dex abundance, name, molar mass (u).\n")
    # Write data
    for i in np.arange(nelements):
      f.write("{:3d}  {:2s}  {:5.2f}  {:10s}  {:12.8f}\n".
              format(index[i], symbol[i], dex[i], name[i], mass[i]))


def readatomic(afile):
  """
  Read an elemental (atomic) composition file.

  Parameters
  ----------
  afile: String
    File with atomic composition.

  Returns
  -------
  anum: 1D integer ndarray
     Atomic number (except for Deuterium, which has anum=0).
  symbol: 1D string ndarray
     Elemental chemical symbol.
  dex: 1D float ndarray
     Logarithmic number-abundance, scaled to log(H) = 12.
  name: 1D string ndarray
     Element names.
  mass: 1D float ndarray
     Elemental mass in amu.

  Uncredited developers
  ---------------------
  Jasmina Blecic
  """
  # Allocate arrays:
  nelements = 84  # Fixed number
  anum   = np.zeros(nelements, np.int)
  symbol = np.zeros(nelements, '|S2')
  dex    = np.zeros(nelements, np.double)
  name   = np.zeros(nelements, '|S20')
  mass   = np.zeros(nelements, np.double)

  # Open-read file:
  with open(afile, 'r') as f:
    # Read-discard first two lines (header):
    f.readline()
    f.readline()
    # Store data into the arrays:
    for i in np.arange(nelements):
      anum[i], symbol[i], dex[i], name[i], mass[i] = f.readline().split()

  return anum, symbol, dex, name, mass


def makepreatm(pressure, temp, afile, elements, species, patm):
  """
  Generate a pre-atm file for TEA containing the elemental abundances
  at each atmospheric layer.

  Parameters
  ----------
  pressure: String
     Pressure atmospheric profile (bar).
  abun_file: String
     Name of the abundances file.
  in_elem: String
     String containing input elemental species.
  out_spec: String
     String containing output molecular species.
  pre_atm: String
     Pre-atmospheric filename.
  temp: 1D float array
     Temperature atmospheric profile (in K).

  Uncredited developers
  ---------------------
  Jasmina Blecic
  """
  # Number of layers:
  nlayers = len(pressure)

  # Get the elemental abundace data:
  index, symbol, dex, name, mass = readatomic(afile)
  # Take only the elements we need:
  iatoms = np.in1d(symbol, elements)

  # Lisf of fractional number of molecules relative to hydrogen:
  nfrac = (10.0**(dex-dex[np.where(symbol=="H")]))[iatoms].tolist()

  # Write pre-atm file:
  f = open(patm, 'w')

  # Pre-atm header with basic instructions
  f.write("# TEA pre-atmosphere file.\n"
          "# Units: pressure (bar), temperature (K), abundance (unitless).\n\n")

  # Load defaults:
  sdefaults = {}
  with open(indir+"TEA_gdata_defaults.txt", "r") as d:
    for line in d:
      sdefaults[line.split()[0]] = line.split()[1]

  # Check species names:
  for i in np.arange(len(species)):
    if species[i].find("_") < 0:
      species[i] = sdefaults[species[i]]

  # Write species names
  f.write('#SPECIES\n{:s}\n\n'.format(" ".join(species)))

  # Write TEA data:
  f.write("#TEADATA\n#Pressure          Temp  " +
          "  ".join(["{:>12s}".format(atom) for atom in symbol[iatoms]]) + "\n")
  for i in np.arange(nlayers):
    # Pressure, and Temperature:
    f.write("{:10.4e}     {:>8.2f}  ".format(pressure[i], temp[i]))
    # Elemental abundance list:
    f.write("  ".join(["{:12.6e}".format(abun) for abun in nfrac]) + "\n")
  f.close()


def TEA2pyrat(teafile, atmfile):
  """
  Format a TEA atmospheric file into a Pyrat atmospheric file.

  Paramters
  ---------
  teafile:  String
     Input TEA atmospheric file.
  atmfile:  String
     Output Pyrat atmospheric file.
  """
  # Open read the TEA file:
  f = open(teafile, "r")
  tea = np.asarray(f.readlines())
  f.close()

  # Line-indices where the species and data are:
  ispec = np.where(tea == "#SPECIES\n")[0][0] + 1  # Species list
  idata = np.where(tea == "#TEADATA\n")[0][0] + 2  # data starting line

  # Read and clean species names:
  species = tea[ispec].split()
  nspecies = len(species)
  for i in np.arange(nspecies):
    species[i] = species[i].rpartition("_")[0]

  # Extract per-layer data:
  nlayers = len(tea) - idata
  temperature = np.zeros(nlayers)
  pressure    = np.zeros(nlayers)
  abundance   = np.zeros((nlayers, nspecies))
  for i in np.arange(nlayers):
    data = np.asarray(tea[idata+i].split(), np.double)
    pressure[i], temperature[i] = data[0:2]
    abundance[i,:] = data[2:]

  # TEA pressure units are always bars:
  punits = "bar"
  pressure = pressure * pt.u(punits)  # Set in barye units
  # File header:
  header = "# TEA atmospheric file formatted for Pyrat.\n\n"

  writeatm(atmfile, pressure, temperature, species, abundance, punits, header)
