#!/usr/bin/env python

# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys, os
import re 
import numpy as np

def main():
  """
  Format the Partridge and Schwenke H2O partition function file.

  Usage:
  ------
  Execute from the Shell:
  ./PFformat_PandS_H2O.py [fileIn] [fileOut]

  Parameters:
  -----------
  fileIn: String
     Input P&S H2O partition-function filename.
  fileOut: String
     Output PF filename.

  Notes:
  ------
  Download the partition-function file to the working directory:
    http://kurucz.harvard.edu/molecules/h2o/h2opartfn.dat
  """
  # Parse arguments:
  if len(sys.argv) > 1:
    fileIn = sys.argv[1:]

  #if len(sys.argv) > 2:
  #  fileOut = sys.argv[2]
  #else:
  fileOut = "PF_Exomol_{:s}.dat"

  # Read and extract data from files:
  lines    = []
  isotopes = []
  data     = []
  for j in np.arange(len(fileIn)):
    with open(fileIn[j], "r") as f:
      lines = f.readlines()

    # Get info from file name:
    s = os.path.split(fileIn[j])[1].split("_")[0].split("-")
    molecule = ""
    iso = ""
    for i in np.arange(len(s)):
      match = re.match(r"([0-9]+(?:.[0-9]+)?)([a-z]+)", s[i], re.I)
      molecule += match.group(2)
      iso      += match.group(1)[-1:]

    print(iso)
    # Extract the isotopes array:
    isotopes.append(iso)

    # Number of temperature samples:
    ntemp = len(lines)
    # Allocate arrays:
    temp = np.zeros(ntemp, np.double)
    z    = np.zeros(ntemp, np.double)
    for i in np.arange(ntemp):
      temp[i], z[i] = lines[i].split()
    data.append(z)

  print(isotopes)
  # Number of isotopes:
  niso  = len(isotopes)
  data = np.asarray(data).T

  # Output file name:
  fout = open(fileOut.format(molecule), "w")

  fout.write(
  "# This file incorporates the tabulated {:s} Partition-function data from\n"
  "# Exomol\n\n".format(molecule))

  fout.write("@ISOTOPES\n            ")
  for j in np.arange(niso):
      fout.write("  {:10s}".format(isotopes[j]))
  fout.write("\n\n")

  fout.write("# Temperature (K), partition function for each isotope:\n")

  fout.write("@DATA\n")
  for i in np.arange(ntemp):
    fout.write(" {:7.1f} ".format(temp[i]))
    for j in np.arange(niso):
      fout.write("  {:10.3f}".format(data[i,j]))
    fout.write("\n")

  fout.close()


if __name__ == "__main__":
  main()
