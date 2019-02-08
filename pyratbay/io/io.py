# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["read_pyrat", "write_opacity", "read_opacity"]

import struct

import numpy as np


def read_pyrat(specfile, wn=True):
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
  with open(specfile, "r") as f:
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


def write_opacity(ofile, molID, temp, press, wn, etable):
  """
  Write an opacity table as a binary file.

  Parameters
  ----------
  ofile: String
      Path to a Pyrat Bay opacity file.
  molID: 1D integer ndarray
      molecule ID.
  temp: 1D float ndarray
      Temperature (Kelvin degree).
  press: 1D float ndarray
      Pressure (barye).
  wn: 1D float ndarray
      Wavenumber (cm-1).
  etable: 4D float ndarray
      Tabulated opacities (cm-1).
  """
  # Get array shapes:
  nmol    = len(molID)
  ntemp   = len(temp)
  nlayers = len(press)
  nwave   = len(wn)

  with open(ofile, "wb") as f:
      # Size of arrays:
      f.write(struct.pack("4l", nmol, ntemp, nlayers, nwave))
      # Arrays:
      f.write(struct.pack(str(nmol)   +"i", *list(molID)))
      f.write(struct.pack(str(ntemp)  +"d", *list(temp) ))
      f.write(struct.pack(str(nlayers)+"d", *list(press)))
      f.write(struct.pack(str(nwave)  +"d", *list(wn)   ))
      # Write opacity data in chunks to avoid memory crashes:
      fmt = str(ntemp * nlayers * nwave) + "d"
      for i in np.arange(nmol):
          f.write(struct.pack(fmt, *list(etable[i].flatten())))


def read_opacity(ofile):
  """
  Read an opacity table from file.

  Parameters
  ----------
  ofile: String
      Path to a Pyrat Bay opacity file.

  Returns
  -------
  sizes: 4-element integer tuple
      Sizes of the dimensions of the opacity table:
      (nmol, ntemp, nlayers, nwave)
  arrays: 4-element 1D ndarray tuple
      The dimensions of the opacity table:
      - molecule ID (integer, unitless, see inputs/molecules.dat)
      - temperature (float, Kelvin)
      - pressure    (float, barye)
      - wavenumber  (float, cm-1)
  etable: 4D float ndarray tuple
      The tabulated opacities (cm-1), of shape [nmol, ntemp, nlayers, nwave].
  """
  with open(ofile, "rb") as f:
      # Read arrays lengths:
      nmol    = struct.unpack('l', f.read(8))[0]
      ntemp   = struct.unpack('l', f.read(8))[0]
      nlayers = struct.unpack('l', f.read(8))[0]
      nwave   = struct.unpack('l', f.read(8))[0]

      # Read wavenumber, temperature, pressure, and isotope arrays:
      molID = np.asarray(struct.unpack(str(nmol   )+'i', f.read(4*nmol   )))
      temp  = np.asarray(struct.unpack(str(ntemp  )+'d', f.read(8*ntemp  )))
      press = np.asarray(struct.unpack(str(nlayers)+'d', f.read(8*nlayers)))
      wn    = np.asarray(struct.unpack(str(nwave  )+'d', f.read(8*nwave  )))

      # Read extinction-coefficient data table:
      ndata = nmol * ntemp * nlayers * nwave
      etable = np.asarray(struct.unpack('d'*ndata, f.read(8*ndata))).reshape(
                          (nmol, ntemp, nlayers, nwave))

  return ((nmol, ntemp, nlayers, nwave),
          (molID, temp, press, wn),
          etable)
