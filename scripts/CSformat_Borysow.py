#!/usr/bin/env python

# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys, os
import numpy as np

def main():
  """
  Format Borysow's H2-He CIA data file.

  Usage
  -----
  Execute from the Shell:
  ./CSformat_Borysow.py fileIn fileOut

  Examples
  --------
  $ ./CSformat_Borysow.py ciah2he_dh_quantmech CIA_Borysow_H2He_1000-7000K_0.5-400um.dat  H2 He
  $ ./CSformat_Borysow.py final_CIA_LT.dat     CIA_Borysow_H2H2_0060-0350K_0.6-1000um.dat H2 H2
  $ ./CSformat_Borysow.py final_CIA_HT.dat     CIA_Borysow_H2H2_0400-1000K_0.6-1000um.dat H2 H2
  $ ./CSformat_Borysow.py CIA.H2H2.Yi          CIA_Borysow_H2H2_1000-7000K_0.6-0500um.dat H2 H2

  Download the CIA tabulated data from A. Borysow's  webpage:
    http://www.astro.ku.dk/~aborysow/programs/ciah2he_dh_quantmech
  """

  # Parse arguments:
  fileIn, fileOut, mol1, mol2 = sys.argv[1:]

  # Read and extract data from files:
  f = open(fileIn, "r")
  lines = f.readlines()
  f.close()

  # Extract temperature arrays:
  temp = lines[1].split()[1:]
  ntemp = len(temp)
  # Remove the trailing 'K':
  for i in np.arange(ntemp):
    temp[i] = temp[i][:-1]
  temp = np.asarray(temp, np.double)

  # Number of wavenumber samples:
  nwave = len(lines) - 3

  # Allocate arrays:
  wn   = np.zeros( nwave,        np.double)
  data = np.zeros((nwave,ntemp), np.double)

  # Extract wavenumber and extinction arrays:
  for i in np.arange(nwave):
    info = lines[i+3].split()
    wn  [i] = info[0]
    data[i] = info[1:]


  # Output file name:
  fout = open(fileOut, "w")

  # Write comments:
  fout.write("# This file contains the reformated {:s}-{:s} CIA data from:\n"
             "#  http://www.astro.ku.dk/~aborysow/programs/{:s}\n\n".
              format(mol1, mol2, os.path.basename(fileIn)))
  # Write header:
  fout.write("@SPECIES\n{:s} {:s}\n\n".format(mol1, mol2))
  fout.write("@TEMPERATURES\n         ")
  for j in np.arange(len(temp)):
      fout.write("      {:4.0f}".format(temp[j]))
  fout.write("\n\n")

  # Write the data:
  fout.write("# Wavenumber in cm-1, CIA coefficients in cm-1 amagat-2:\n")
  fout.write("@DATA\n")
  for i in np.arange(nwave):
    fout.write(" {:7.1f} ".format(wn[i]))
    for j in np.arange(ntemp):
      fout.write(" {:.3e}".format(data[i,j]))
    fout.write("\n")

  fout.close()


if __name__ == "__main__":
  main()
