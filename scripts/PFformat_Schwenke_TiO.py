#!/usr/bin/env python

import sys
import numpy as np

def main():
  """
  Format the Schwenke TiO partition function file:

  Usage:
  ------
  Execute from the Shell:
  ./PFformat_Schwenke_TiO.py [fileIn] [fileOut]

  Parameters:
  -----------
  fileIn: String
     Input Schwenke TiO partition-function filename.
  fileOut: String
     Output PF filename.

  Download the partition-function file to the working directory:
    http://kurucz.harvard.edu/molecules/tio/tiopart.dat
  """

  # Parse arguments:
  if len(sys.argv) > 1:
    fileIn = sys.argv[1]
  else:
    fileIn = "tiopart.dat"

  if len(sys.argv) > 2:
    fileOut = sys.argv[2]
  else:
    fileOut = "PF_Schwenke_TiO.dat"

  # Read and extract data from files:
  f = open(fileIn, "r")
  lines = f.readlines()
  f.close()

  # Number of header lines in file:
  offset = 1

  # Extract the isotopes array:
  iso = lines[0].split()[1:]
  # Number of isotopes:
  niso  = len(iso)

  # Number of temperature samples:
  ntemp = len(lines) - offset

  # Allocate arrays:
  temp = np.zeros( ntemp,        np.double)
  data = np.zeros([ntemp, niso], np.double)

  for i in np.arange(ntemp):
    info = lines[i+offset].split()
    temp[i] = info[0]
    data[i] = info[1:]

  # Output file name:
  fout = open(fileOut, "w")

  fout.write(
  "# This file incorporates the tabulated TiO Partition-function data from\n"
  "#   http://kurucz.harvard.edu/molecules/tio/tiopart.dat\n\n")

  fout.write("@ISOTOPES\n              ")
  for j in np.arange(niso):
      fout.write("  {:10s}".format(iso[j]))
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
