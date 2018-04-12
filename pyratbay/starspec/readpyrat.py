# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["readpyrat"]


import numpy as np


def readpyrat(specfile, wn=True):
  """
  Read a pyrat spectrum file.

  Parameters
  ----------
  specfile: String
     Path to output Transit spectrum file to read.
  wn: Boolean
     If True convert wavelength to wavenumber.

  Return
  ------
  wave: 1D float ndarray
     The spectrum's wavenumber or wavelength array.
  spectrum: 1D float ndarray
     The pyrat spectrum.
  """
  f = open(specfile, "r")

  # Count number of lines in file:
  f.seek(0)
  # Ignore first two lines of comments:
  l = f.readline()
  l = f.readline()
  ndata = 0
  for line in f:
    ndata += 1

  # Initialize arrays:
  wave     = np.zeros(ndata, np.double)
  spectrum = np.zeros(ndata, np.double)

  # Return to begining of file:
  f.seek(0)
  f.readline()
  f.readline()
  for i in np.arange(ndata):
    l = f.readline().strip().split()
    wave    [i] = np.double(l[ 0])
    spectrum[i] = np.double(l[-1])

  # Convert wavelength (micron) to wavenumber (cm-1):
  if wn:
    wave = 1e4/wave

  return wave, spectrum
