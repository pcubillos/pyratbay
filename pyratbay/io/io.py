# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["write_spectrum", "read_spectrum",
           "write_opacity", "read_opacity"]

import struct

import numpy as np



def write_spectrum(wl, spectrum, filename, path, wlunits='um'):
  """
  Write a Pyrat spectrum to file.


  Parameters
  ----------
  wl: 1D float iterable
      Wavelength array  in cm units.
  spectrum: 1D float iterable
      Spectrum array. (rp/rs)**2 for transmission (unitless),
      planetary flux for emission (erg s-1 cm-2 cm units).
  filename: String
      Output file name.
  path: String
      Observing mode: 'transit' for transmission, 'eclipse' for
      emission.
  wlunits: String
      Output units for wavelength.

  Examples
  --------
  >>> # See read_spectrum() examples.
  """
  # Need to import here to avoid circular imports:
  from .. import tools as pt

  if filename is None:
      return

  # Type of spectrum and units:
  if path == "transit":
      spectype  = "(Rp/Rs)**2"
      specunits = "unitless"
  elif path == "eclipse":
      spectype  = "Flux"
      specunits = "erg s-1 cm-2 cm"
  else:
      raise ValueError("Input 'path' argument must be either 'transit' "
                       "or 'eclipse'.")

  # Wavelength units in brackets:
  wl = wl/pt.u(wlunits)
  # Precision of 5 decimal places (or better if needed):
  precision = -np.floor(np.log10(np.amin(np.abs(np.ediff1d(wl)))))
  precision = int(np.clip(precision+1, 5, np.inf))
  buff = precision + 5

  # Open-write file:
  with open(filename, "w") as f:
      # Write header:
      f.write("# {:>{:d}s}   {:>15s}\n".format("Wavelength", buff, spectype))
      f.write("# {:>{:d}s}   {:>15s}\n".format(wlunits, buff, specunits))
      # Write the spectrum values:
      for wave, flux in zip(wl, spectrum):
          f.write("{:>{:d}.{:d}f}   {:.9e}\n".
                  format(wave, buff+2, precision, flux))


def read_spectrum(filename, wn=True):
  """
  Read a Pyrat spectrum file.

  Parameters
  ----------
  filename: String
     Path to output Transit spectrum file to read.
  wn: Boolean
     If True convert wavelength to wavenumber.

  Return
  ------
  wave: 1D float ndarray
     The spectrum's wavenumber (in cm units) or wavelength array (in
     the input file's units).
  spectrum: 1D float ndarray
     The spectrum in the input file.

  Examples
  --------
  >>> import pyratbay.io as io
  >>> # Write a spectrum to file:
  >>> nwave = 7
  >>> wl = np.linspace(1.1, 1.7, nwave) * 1e-4
  >>> spectrum = np.ones(nwave)
  >>> io.write_spectrum(wl, spectrum,
  >>>     filename='sample_spectrum.dat', path='transit', wlunits='um')
  >>> # Take a look at the output file:
  >>> with open('sample_spectrum.dat', 'r') as f:
  >>>     print("".join(f.readlines()))
  # Wavelength        (Rp/Rs)**2
  #         um          unitless
       1.10000   1.000000000e+00
       1.20000   1.000000000e+00
       1.30000   1.000000000e+00
       1.40000   1.000000000e+00
       1.50000   1.000000000e+00
       1.60000   1.000000000e+00
       1.70000   1.000000000e+00
  >>> # Now, read from file (getting wavenumber array):
  >>> wn, flux = io.read_spectrum('sample_spectrum.dat')
  >>> print(wn)
  [9090.90909091 8333.33333333 7692.30769231 7142.85714286 6666.66666667
   6250.         5882.35294118]
  >>> print(flux)
  [1. 1. 1. 1. 1. 1. 1.]
  >>> # Read from file (getting wavelength array):
  >>> wl, flux = io.read_spectrum('sample_spectrum.dat', wn=False)
  >>> print(wl)
  [1.1 1.2 1.3 1.4 1.5 1.6 1.7]
  >>> print(flux)
  [1. 1. 1. 1. 1. 1. 1.]
  """
  # Need to import here to avoid circular imports:
  from .. import tools as pt
  with open(filename, "r") as f:
    # Count number of lines in file:
    f.seek(0)
    # Get wavelength units from header:
    l = f.readline()
    l = f.readline()
    wlunits = l.split()[1]

    wave, spectrum = np.array([line.strip().split() for line in f], np.double).T

    # Convert wavelength to wavenumber (always in cm-1):
    if wn:
        wave = 1.0/(wave*pt.u(wlunits))

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
