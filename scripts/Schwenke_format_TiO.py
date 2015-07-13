#!/usr/bin/env python
# Format the Schwenke TiO partition function file:

# Download the partition-function file to the working directory:
#   http://kurucz.harvard.edu/molecules/tio/tiopart.dat

import numpy as np

# Read and extract data from files:
f = open("tiopart.dat", "r")
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
fout = open("TiO_PF_Schwenke.dat", "w")

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

