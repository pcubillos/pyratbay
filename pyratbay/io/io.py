# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["write_spectrum", "read_spectrum",
           "write_opacity", "read_opacity",
           "write_pf", "read_pf",
           "write_cs", "read_cs",
          ]

import os
import struct

import numpy as np

import pyratbay.constants as pc


def write_spectrum(wl, spectrum, filename, type, wlunits='um'):
  """
  Write a Pyrat spectrum to file.

  Parameters
  ----------
  wl: 1D float iterable
      Wavelength array in cm units.
  spectrum: 1D float iterable
      Spectrum array. (rp/rs)**2 for transmission (unitless),
      planetary flux for emission (erg s-1 cm-2 cm units).
  filename: String
      Output file name.
  type: String
      Data type:
      'transit' for transmission,
      'eclipse' for emission,
      'filter' for a instrumental filter transmission.
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
  if type == "transit":
      spectype  = "(Rp/Rs)**2"
      specunits = "unitless"
  elif type == "eclipse":
      spectype  = "Flux"
      specunits = "erg s-1 cm-2 cm"
  elif type == "filter":
      spectype  = "transmission"
      specunits = "unitless"
  else:
      raise ValueError("Input 'type' argument must be 'transit', 'eclipse',"
                       " or 'filter'.")

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
  Read a Pyrat spectrum file, a plain text file with two-columns: the
  wavelength and signal.  If wn is true, this function converts
  wavelength to wavenumber in cm-1.  The very last comment line sets
  the wavelength units (the first string following a blank, e.g., the
  string '# um' sets the wavelength units as microns).
  If the units are not defined, assume wavelength units are microns.

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
  >>>     filename='sample_spectrum.dat', type='transit', wlunits='um')
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
  # Extract data:
  data = np.loadtxt(filename, unpack=True)
  wave, spectrum = data[0], data[1]

  if not wn:
      return wave, spectrum

  # Check 'header' (last comment line) for wavelength units:
  with open(filename, "r") as f:
      for line in f:
          info = line
          if not line.strip().startswith('#') and line.strip() != '':
              break

  # Get wavelength units from last line of comments:
  if len(info.split()) > 1:
      wlunits = info.split()[1]
  else:
      wlunits = 'um'
  if wlunits not in pc.validunits:
      wlunits = 'um'

  # Convert wavelength to wavenumber in cm-1:
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


def write_pf(pffile, pf, isotopes, temp, header=None):
    """
    Write a partition-function file in Pyrat Bay format.

    Parameters
    ----------
    pffile: String
        Output partition-function file.
    pf: 2D float iterable
        Partition-function data (of shape [niso, ntemp]).
    isotopes: 1D string iterable
        Isotope names.
    temp: 1D float iterable
        Temperature array.
    header: String
        A header for the partition-function file (must be as comments).

    Examples
    --------
    >>> # See read_pf() examples.
    """
    if len(isotopes) != np.shape(pf)[0]:
        raise ValueError('Shape of the partition-function array does not '
                         'match with the number of isotopes.')
    if len(temp) != np.shape(pf)[1]:
        raise ValueError('Shape of the partition-function array does not '
                         'match with the number of temperature samples.')

    # Write output file:
    with open(pffile, "w") as f:
        if header is not None:
            f.write(header)
        f.write("@ISOTOPES\n           "
              + "  ".join(["{:>11s}".format(iso) for iso in isotopes])
              + "\n\n")

        f.write("# Temperature (K), partition function for each isotope:\n")
        f.write("@DATA\n")
        for t, z in zip(temp, pf.T):
            f.write("  {:7.1f}  ".format(t)
                  + "  ".join("{:.5e}".format(d) for d in z) + "\n")


def read_pf(pffile):
    r"""
    Read a partition-function file.

    Parameters
    ----------
    pffile: String
        Partition function file to read.

    Returns
    -------
    pf: 2D float ndarray
        The partition function data (of shape [niso, ntemp]).
    isotopes: List of strings
         The names of the tabulated isotopes.
    temp: 1D float ndarray
        Array with temperature sample.

    Examples
    --------
    >>> import pyratbay.io as io
    >>> # Generate some mock PF data and write to file:
    >>> pffile = 'PF_Exomol_NH3.dat'
    >>> isotopes = ['4111', '5111']
    >>> temp   = np.linspace(10,100,4)
    >>> pf     = np.array([np.logspace(0,3,4),
    >>>                    np.logspace(1,4,4)])
    >>> header = '# Mock partition function for NH3.\n'
    >>> io.write_pf(pffile, pf, isotopes, temp, header)

    >>> # Now, read it back:
    >>> pf, iso, temp = io.read_pf(pffile)
    >>> for item in [iso, temp, pf]:
    >>>     print(item)
    ['4111' '5111']
    [ 10.  40.  70. 100.]
    [[1.e+00 1.e+01 1.e+02 1.e+03]
     [1.e+01 1.e+02 1.e+03 1.e+04]]
    """
    if not os.path.isfile(pffile):
        raise ValueError("Partition-function file '{:s}' does not exist.".
                         format(pffile))

    with open(pffile, "r") as f:
        lines = f.readlines()

    nlines = len(lines)
    lines  = iter(lines)
    for i,line in enumerate(lines):
        line = line.strip()
        # Stop when the tabulated data begins:
        if line == "@DATA":
            break
        # Read isotopes:
        if line == "@ISOTOPES":
            isotopes = np.asarray(next(lines).split())

    # Number of samples:
    niso = len(isotopes)
    ntemp = nlines - i - 2
    # Allocate arrays:
    temp = np.zeros(ntemp, np.double)
    pf   = np.zeros((niso, ntemp), np.double)

    # Read the data:
    for i,line in enumerate(lines):
        info = line.split()
        temp[i] = info[0]
        pf[:,i] = info[1:]

    return pf, isotopes, temp


def write_cs(csfile, cs, species, temp, wn, header=None):
    """
    Write a cross-section file in Pyrat Bay format.

    Parameters
    ----------
    csfile: String
        Output cross-section file.
    cs: 2D float iterable
        Cross-section opacity in units of cm-1 amagat^-N, with N the
        number of species, of shape [ntemp, nwave].
    species: 1D string iterable
        Species names.
    temp: 1D float iterable
        Temperature array in Kelvin degree.
    wn: 1D float iterable
        Wavenumber array in cm-1.
    header: String
        A header for the cross-section file (must be as comments).

    Examples
    --------
    >>> # See read_cs() examples.
    """
    if len(temp) != np.shape(cs)[0]:
        raise ValueError('Shape of the cross-section array does not '
                         'match the number of temperature samples.')
    if len(wn) != np.shape(cs)[1]:
        raise ValueError('Shape of the cross-section array does not '
                         'match the number of wavenumber samples.')

    with open(csfile, "w") as f:
        if header is not None:
            f.write(header)
        f.write("@SPECIES\n"
              + "  ".join(["{:s}".format(spec) for spec in species])
              + "\n\n")
        f.write("@TEMPERATURES\n                "
              + "      ".join(["{:4.0f}".format(t) for t in temp])
              + "\n\n")

        f.write("# Wavenumber in cm-1, opacity in cm-1 amagat-{:d}:\n".
                format(len(species)))
        f.write("@DATA\n")
        for wave, data in zip(wn, cs.T):
            f.write("  {:7.1f}  ".format(wave)
                  + " ".join("{:.3e}".format(d) for d in data) + "\n")


def read_cs(csfile):
    r"""
    Read a cross-section file.

    Parameters
    ----------
    csfile: String
        Partition function file to read.

    Returns
    -------
    cs: 2D float ndarray
        Cross-section opacity in units of cm-1 amagat^-N, with N the
        number of species, of shape [ntemp, nwave].
    species: 1D string list
        Species names.
    temp: 1D float ndarray
        Temperature array in Kelvin degree.
    wn: 1D float ndarray
        Wavenumber array in cm-1.

    Examples
    --------
    >>> import pyratbay.io as io
    >>> # Generate some mock PF data and write to file:
    >>> csfile = 'CS_Mock-HITRAN_H2-H2.dat'
    >>> species = ['H2', 'H2']
    >>> temp = np.linspace(100, 1000, 3)
    >>> wn   = np.arange(10, 15, 1.0)
    >>> cs   = np.array([np.logspace( 0,-4,5),
    >>>                  np.logspace(-1,-5,5),
    >>>                  np.logspace(-2,-6,5)])
    >>> header = '# Mock cross-section for H2-H2.\n'
    >>> io.write_cs(csfile, cs, species, temp, wn, header)
    >>> # Now, read it back:
    >>> cs, species, temp, wn = io.read_cs(csfile)
    >>> for item in [species, temp, wn, cs]:
    >>>     print(item)
    ['H2', 'H2']
    [ 100.  550. 1000.]
    [10. 11. 12. 13. 14.]
    [[1.e+00 1.e-01 1.e-02 1.e-03 1.e-04]
     [1.e-01 1.e-02 1.e-03 1.e-04 1.e-05]
     [1.e-02 1.e-03 1.e-04 1.e-05 1.e-06]]
    """
    if not os.path.isfile(csfile):
        raise ValueError("Cross-section file '{:s}' does not exist.".
                         format(csfile))

    with open(csfile, "r") as f:
        lines = f.readlines()

    # Number of header lines (to skip when reading the tabulated data):
    nlines = len(lines)
    lines  = iter(lines)
    for i,line in enumerate(lines):
        line = line.strip()
        # Stop when the tabulated data begins:
        if line == "@DATA":
            break
        # Get species:
        elif line == "@SPECIES":
            species = next(lines).split()
        # Get the sampled temperatures:
        elif line == "@TEMPERATURES":
            temp = np.array(next(lines).split(), np.double)

    # Number of samples:
    nwave = nlines - i - 3
    ntemp = len(temp)
    # Allocate arrays:
    wn = np.zeros(nwave, np.double)
    cs = np.zeros((ntemp, nwave), np.double)

    # Read the data:
    for i,line in enumerate(lines):
        info = line.split()
        wn[i]   = info[0]
        cs[:,i] = info[1:]

    return cs, species, temp, wn
