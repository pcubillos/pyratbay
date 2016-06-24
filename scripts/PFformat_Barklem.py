#!/usr/bin/env python

# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys
import numpy as np

def main():
  """
  Format  Partridge and Schwenke H2O partition function file.

  Usage
  -----
  Execute from the Shell:
  ./PFformat_Barklem.py fileIn Species

  Parameters
  ----------
  fileIn: String
     Input partition-function filename.
  Species: String
     Species name.

  Notes
  -----
  Download the partition-function file to the working directory:
    TBD 
  """
  # Parse arguments:
  fileIn  = sys.argv[1]
  species = sys.argv[2]
  fileOut = "PF_Barklem_{:s}.dat".format(species)

  # Read and extract data from files:
  f = open(fileIn, "r")
  lines = f.readlines()
  f.close()

  # extract temperature array:
  temp = np.asarray(lines[2].split()[2:], np.double)

  iso  = []
  data = []
  # Find the species and extract dthe PF data:
  for i in np.arange(4, len(lines)):
    line = lines[i].strip()
    if line.startswith("{:s}_".format(species)):
      iso .append(line.split()[0])
      data.append(line.split()[1:])

  data = np.asarray(data, np.double).T
  # Number of isotopes and temperature samples:
  niso  = len(iso)
  ntemp = len(temp)

  # Output file name:
  fout = open(fileOut, "w")

  fout.write(
  "# This file incorporates the tabulated {:s} Partition-function data from\n"
  "#   Barklem et al.\n\n")

  fout.write("@ISOTOPES\n            ")
  for j in np.arange(niso):
      fout.write("  {:10s}".format(iso[j]))
  fout.write("\n\n")

  fout.write("# Temperature (K), partition function for each isotope:\n")

  fout.write("@DATA\n")
  for i in np.arange(ntemp):
    fout.write(" {:7.2e} ".format(temp[i]))
    for j in np.arange(niso):
      fout.write("  {:10.6f}".format(data[i,j]))
    fout.write("\n")

  fout.close()


if __name__ == "__main__":
  main()
