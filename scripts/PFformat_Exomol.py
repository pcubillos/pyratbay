#!/usr/bin/env python

# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys, os
import re 
import numpy as np

def main():
  """
  Format Exomol partition function files.

  Usage:
  ------
  Execute from the Shell:
  ./PFformat_Exomol.py file1 [file2] ... [fileN]

  Parameters:
  -----------
  file1--fileN: String
     Input Exomol partition-function filenames.

  Notes:
  ------
  As far as I've seen, each input file contains a single isotope.
  More exotic formats will break the code.

  Exomol database: http://www.exomol.com/
  """
  # Parse arguments:
  if len(sys.argv) > 1:
    fileIn = sys.argv[1:]
  # Output file:
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
      match = re.match(r"([0-9]+)([a-z]+)([0-9]*)", s[i], re.I)
      N = 1 if match.group(3) == "" else int(match.group(3))
      molecule += match.group(2) + match.group(3)
      iso      += match.group(1)[-1:] * N

    # Extract the isotopes array:
    isotopes.append(iso)

    # Number of temperature samples:
    ntemp = len(lines)
    # Allocate arrays:
    temp = np.zeros(ntemp, np.double)
    z    = np.zeros(ntemp, np.double)
    for i in np.arange(ntemp):
      temp[i], z[i] = lines[i].split()
    if data != [] and len(z) != len(data[-1]):
      print("Warning! Lengths of PF files do not match!\n")
      # FINDME: I should stop here, but for the moment, I just roll with it
    data.append(z)

  # Number of isotopes:
  niso  = len(isotopes)
  # Temporary patch (zero pad):
  maxlen = 0
  for i in np.arange(niso):
    maxlen = np.amax((maxlen, len(data[i])))
  for i in np.arange(niso):
    d = np.zeros(maxlen)
    d[0:len(data[i])] = data[i]
    data[i] = d

  data = np.asarray(data).T
  # Output file name:
  fout = open(fileOut.format(molecule), "w")

  fout.write(
  "# This file incorporates the tabulated {:s} partition-function data\n"
  "# from Exomol\n\n".format(molecule))

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
  print("\nWriten partition-function file:\n  '{:s}'\nfor {:s} molecule, "
   "with isotopes {:s},\nand temperature range {:.0f}K--{:.0f}K.".
      format(fileOut.format(molecule), molecule, isotopes, temp[0], temp[-1]))


if __name__ == "__main__":
  main()
