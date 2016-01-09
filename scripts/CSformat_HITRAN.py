#!/usr/bin/env python

import sys, os
import numpy as np

def main():
  """
  Format HITRAN CIA data (Richard et al. 2012) for use in Transit.

  Usage
  -----
  To run this code, execute from the Shell:
   ./HITRAN_CIA_format.py fileIN [fileOUT] [temp_step] [wnumber_step]

  Parameters
  ----------
  fileIN: String
     File name of the input CIA file as given by Richard et al. (2012)
  fileOUT: String
     File name of the ouput CIA file for use with Transit.
  temp_step: Float
     Optional temperature sampling step to thin the output array.
  wnumber_step: Float
     Optional wavenumber sampling step to thin the output array.

  Examples
  --------
  For H2-H2 run from the Shell prompt:
    ./HITRAN_CIA_format.py  H2-H2_2011.cia  CIA_HITRAN_H2H2_0200-3000K_1.0-500um.dat  50  10

  Previous (uncredited) developers
  --------------------------------
  - Dylan Bruce (UCF)
  """

  # Loschmidt number (cm-3):
  N0 = 2.6867774e19

  # Input/output files:
  filein  = sys.argv[1]
  fileout = None
  if len(sys.argv) > 2:
    fileout = sys.argv[2]

  # User-defined sampling rates:
  tstep = None           # For temperature
  if len(sys.argv) > 3:
    tstep = float(sys.argv[3])

  wstep = None           # For wavenumber
  if len(sys.argv) > 4:
    wstep = float(sys.argv[4])

  # Read and extract data from files:
  f = open(filein, "r")
  lines = f.readlines()
  f.close()


  # Find all the header lines:
  headers = []
  for line in lines:
    if line.startswith(lines[0][0:20]):  # First 20 chars contain the species
      headers.append(line)
  # Save the species names:
  species = lines[0][0:20].strip().split("-")


  # Some files have multiple sets of opacities,
  # Count number of sets:
  wnmin, wnmax = [], []
  temps  = []
  nwaves = []

  wmin, wmax = -1.0, -1.0
  for header in headers:
    # New set:
    if float(header[20:30]) != wmin or float(header[30:40]) != wmax:
      wnmin. append( float(header[20:30]))
      wnmax. append( float(header[30:40]))
      nwaves.append(   int(header[40:47]))
      temps. append([float(header[47:54])])
      wmin = wnmin[-1]
      wmax = wnmax[-1]
    # Same current set:
    else:
      temps[-1].append(float(header[47:54]))

  # Number of sets:
  nsets = len(nwaves)

  # Extract the opacity data:
  wave = []
  data = []
  ntemps = []
  istart = 1
  for i in np.arange(nsets):
    # Initialize arrays:
    ntemps.append(len(temps[i]))
    data.append( np.zeros((nwaves[i], ntemps[i]), np.double) )
    wave.append( np.zeros(nwaves[i], np.double))

    # Write down wavenumbers and coefficients into arrays:
    for j in np.arange(ntemps[i]):
      for k in np.arange(nwaves[i]):
        line = lines[istart + k].split()
        if j == 0:
          wave[i][k] = line[0]
        data[i][k, j] = line[1]
      istart += nwaves[i] + 1


  # Thin the arrays if requested:
  # Note: So far used for H2-H2 and H2-He files.
  if tstep is not None  and  wstep is not None  and nsets == 1:
    # Thinned arrays:
    Tthin  = np.arange(np.amin(temps[0]), np.amax(temps[0])+1, tstep)
    wnthin = np.arange(np.amin(wave[0]),  np.amax(wave[0]) +1, wstep)
    # Indices corresponding to the thinned arrays:
    itemp = np.where(np.in1d(temps[0], Tthin ))[0]
    iwn   = np.where(np.in1d(wave[0],  wnthin))[0]
    # Slice the data arrays:
    temps[0] = temps[0][itemp]
    wave [0] = wave [0][iwn]
    data [0] = data [0][iwn,:]
    data [0] = data [0][:, itemp]
    # Update array sizes:
    ntemps = [len(temps[0])]
    nwaves = [len(wave [0])]

  # Scale the opacity to the correct units (cm5 molecule-1 --> cm-1 amagat-2):
  for i in np.arange(nsets):
    data[i] *= N0**2


  # Write to the output file:
  for i in np.arange(nsets):
    if fileout is None or nsets > 1:
      fileout = "CIA_HITRAN_{:s}-{:s}_{:06.1f}-{:06.1f}K_{:.1f}-{:.1f}um.dat".\
                 format(species[0],   species[1],
                        temps[i][0],  temps[i][-1],
                        1e4/wnmax[i], 1e4/wnmin[i])

    fout = open(fileout, "w")

    # Write comments:
    fout.write("# This file contains the reformated {:s}-{:s} CIA data from:\n"
               "#  Richard et al. (2012), HITRAN file: {:s}\n\n".
               format(species[0], species[1], filein))
    # Write header:
    fout.write("@SPECIES\n{:s} {:s}\n\n".format(species[0], species[1]))
    fout.write("@TEMPERATURES\n          ")
    for j in np.arange(ntemps[i]):
        fout.write("      {:4.0f}".format(temps[i][j]))
    fout.write("\n\n")

    # Write down the data:
    fout.write("# Wavenumber in cm-1, CIA coefficients in cm-1 amagat-2:\n")
    fout.write("@DATA\n")
    for j in np.arange(nwaves[i]):
      fout.write("  {:7.1f} ".format(wave[i][j]))
      for k in np.arange(ntemps[i]):
        fout.write(" {:.3e}".format(data[i][j,k]))
      fout.write("\n")

    fout.close()

if __name__ == "__main__":
  main()
