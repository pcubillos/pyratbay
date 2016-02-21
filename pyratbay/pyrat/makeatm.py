import numpy as np

from .. import tools     as pt

def uniform(atmfile, pressure, temperature, species, abundances):
  """
  Generate an atmospheric file with uniform abundances.

  Parameters
  ----------
  atmfile: String
     Name of output atmospheric file.
  pressure: 1D float ndarray
     Monotonously decreasing pressure profile (in bar).
  temperature: 1D float ndarray
     Temperature profile for pressure layers (in K).
  species: 1D string ndarray
     List of atmospheric species.
  abundances: 1D float ndarray
     The species mole mixing ratio.
  """

  # Safety checks:
  if len(pressure) != len(temperature):
    pt.error("Pressure array length ({:d}) and Temperature array length "
             "({:d}) don't match.".format(len(pressure), len(temperature)))
  if len(species) != len(abundances):
    pt.error("Species array length ({:d}) and Abundances array length ({:d}) "
             "don't match.".format(len(species), len(abundances)))
  # FINDME: Check pressure array is monotonously decreasing.

  # Write to file:
  f = open(atmfile, "w")

  f.write("# This is an atmospheric file with pressure, temperature, \n"
          "# and uniform mole mixing ratio profiles.\n\n")

  # Set the values units:
  f.write("# Abundance units (by number or mass):\n@ABUNDANCE\nnumber\n")
  f.write("# Pressure units:\n@PRESSURE\nbar\n")
  f.write("# Temperatures are always Kelvin.\n\n")

  # Write the species names:
  f.write("# Atmospheric composition:\n"
          "@SPECIES\n" +
          "  ".join(["{:<s}".format(mol) for mol in species]) + '\n\n')

  # Write the per-layer data:
  # Column's header:
  f.write("# Pressure  Temperature  " +
        "".join(["{:<14s}".format(mol) for mol in species]) + "\n")
  f.write("@DATA\n")

  # For each layer write TEA data
  for i in np.arange(len(pressure)):
      # Pressure, and Temperature:
      f.write("{:10.4e}  {:11.3f}  ".format(pressure[i], temperature[i]))
      # Species mole mixing ratios:
      f.write("  ".join(["{:12.6e}".format(ab) for ab in abundances]) + "\n")
  f.close()
