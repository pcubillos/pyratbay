# Copyright (c) 2016-2020 Patricio Cubillos.
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE).

__all__ = [
    'save_pyrat', 'load_pyrat',
    'write_atm', 'read_atm',
    'write_spectrum', 'read_spectrum',
    'write_opacity',  'read_opacity',
    'write_pf', 'read_pf',
    'write_cs', 'read_cs',
    'read_pt',
    ]

import os
import struct
import pickle

import numpy as np

from .. import constants as pc
from .. import tools as pt


def save_pyrat(pyrat, pfile=None):
    """
    Save a pyrat instance into a pickle file.

    Parameters
    ----------
    pyrat: A Pyrat instance
        Object to save.
    pfile: String
        Name of output file.  Default to the pyrat logname (changing
        the extension to '.pickle').
    """
    if pfile is None:
        pfile = os.path.splitext(pyrat.log.logname)[0] + '.pickle'
        print('Saving pyrat instance to: {}'.format(pfile))
    # Reset values to reduce pickle size:
    with pt.tmp_reset(pyrat, 'spec.own', 'voigt.profile', 'log.file',
            'ex.ec', 'ex.etable', 'ret.posterior',
            lt=pyrat.lt.clone_new(pyrat)):
        with open(pfile, 'wb') as f:
            pickle.dump(pyrat, f, pickle.HIGHEST_PROTOCOL)


def load_pyrat(pfile):
    """
    Load a pyrat instance from a pickle file.

    Parameters
    ----------
    pfile: String
        Name of input pickle file.

    Returns
    -------
    pyrat: A Pyrat instance
        Loaded object.
    """
    with open(pfile, 'rb') as f:
        pyrat = pickle.load(f)
    pyrat.log.verb = -1
    pyrat.setup_spectrum()
    pyrat.log.verb = pyrat.verb
    # Recover MCMC posterior:
    if pt.isfile(pyrat.ret.mcmcfile) == 1:
        with np.load(pyrat.ret.mcmcfile) as d:
            Z = d['Z']
            Zchain = d['Zchain']
        burnin = pyrat.ret.burnin
        ipost = np.ones(len(Z), bool)
        for c in np.unique(Zchain):
            ipost[np.where(Zchain == c)[0][0:burnin]] = False
        ipost[np.where(Zchain == -1)] = False
        pyrat.ret.posterior = Z[ipost]

    return pyrat


def write_atm(atmfile, pressure, temperature, species, abundances,
              punits, header, radius=None, runits=None):
    r"""
    Write an atmospheric file following the Pyrat format.

    Parameters
    ----------
    atmfile: String
        Name of output atmospheric file.
    pressure: 1D float ndarray
        Monotonously decreasing pressure profile (in barye).
    temperature: 1D float ndarray
        Temperature profile for pressure layers (in Kelvin).
    species: 1D string ndarray
        List of atmospheric species.
    abundances: 2D float ndarray
        The species mole mixing ratio (of shape [nlayers,nspecies]).
    punits:  String
        Pressure units of output.
    header:  String
        Header message (comment) to include at the top of the file.
    radius: 1D float ndarray
        Monotonously increasing radius profile (in cm).
    runits:  String
        Radius units of output.

    Examples
    --------
    >>> import numpy as np
    >>> import pyratbay.io as io
    >>> import pyratbay.atmosphere as pa

    >>> atmfile = 'WASP-00b.atm'
    >>> nlayers = 5
    >>> pressure    = pa.pressure('1e-8 bar', '1e2 bar', nlayers)
    >>> temperature = pa.tmodels.Isothermal(nlayers)(1500.0)
    >>> species     = "H2 He H2O".split()
    >>> abundances  = [0.8499, 0.15, 1e-4]
    >>> qprofiles = pa.uniform(pressure, temperature, species, abundances)
    >>> io.write_atm(atmfile, pressure, temperature, species, qprofiles,
    >>>     punits='bar', header='# Example atmospheric file:\n')
    >>> # Print output file:
    >>> with open(atmfile, 'r') as f:
    >>>     print(f.read())
    # Example atmospheric file:
    # Abundance units (by number or mass):
    @ABUNDANCE
    number
    # Pressure units:
    @PRESSURE
    bar
    # Temperatures units:
    @TEMPERATURE
    kelvin
    # Atmospheric composition:
    @SPECIES
    H2  He  H2O

    # Pressure  Temperature  H2            He            H2O
    @DATA
    1.0000e-08     1500.000  8.499000e-01  1.500000e-01  1.000000e-04
    3.1623e-06     1500.000  8.499000e-01  1.500000e-01  1.000000e-04
    1.0000e-03     1500.000  8.499000e-01  1.500000e-01  1.000000e-04
    3.1623e-01     1500.000  8.499000e-01  1.500000e-01  1.000000e-04
    1.0000e+02     1500.000  8.499000e-01  1.500000e-01  1.000000e-04
    """
    f = open(atmfile, "w")
    f.write(header)

    # Set the values units:
    f.write("# Abundance units (by number or mass):\n@ABUNDANCE\nnumber\n")
    f.write("# Pressure units:\n@PRESSURE\n{:s}\n".format(punits))
    if radius is not None:
        f.write("# Radius units:\n@RADIUS\n{:s}\n".format(runits))
    f.write("# Temperatures units:\n@TEMPERATURE\nkelvin\n")

    # Write the species names:
    f.write("# Atmospheric composition:\n@SPECIES\n" +
            "  ".join(["{:<s}".format(mol) for mol in species]) + '\n\n')
    # Write the per-layer data:
    if radius is not None:
        f.write("# Radius    Pressure    Temperature  ")
    else:
        f.write("# Pressure  Temperature  ")
    f.write("".join(["{:<14s}".format(mol) for mol in species]) + "\n")
    f.write("@DATA\n")

    pressure = pressure/pt.u(punits)
    if radius is not None:
        radius = radius/pt.u(runits)

    # Write data for each layer:
    nlayers = len(pressure)
    for i in np.arange(nlayers):
        # (radius,) pressure, and temperature:
        if radius is not None:
            f.write("{:10.4e}  ".format(radius[i]))
        f.write("{:10.4e}  {:11.3f}  ".format(pressure[i], temperature[i]))
        # Species mole mixing ratios:
        f.write("  ".join(["{:12.6e}".format(ab) for ab in abundances[i]])
                + "\n")
    f.close()


def read_atm(atmfile):
    r"""
    Read a Pyrat atmospheric file.

    Parameters
    ----------
    atmfile: String
       File path to a Pyrat Bay's atmospheric file.

    Returns
    -------
    units: 4-element string tuple
        Units for pressure, temperature, abundance, and radius as given
        in the atmospheric file.
    species: 1D string ndarray
        The list of species names read from the atmospheric file (of
        size nspec).
    press: 1D float ndarray
        The atmospheric pressure profile (of size nlayers). The
        file's @PRESSURE keyword indicates the ouptput units.
    temp: 1D float ndarray
        The atmospheric temperature profile (of size nlayers). The
        file's @TEMPERATURE keyword indicates the ouptput units.
    q: 2D float ndarray
        The mixing ratio profiles of the atmospheric species (of size
        [nlayers,nspec]).  The file's @ABUNDANCE indicates the output
        units.
    radius: 1D float ndarray
        The atmospheric altiture profile (of size nlayers).  None if the
        atmospheric file does not contain a radius profile.
        The file's @RADIUS keyword indicates the output units.

    Examples
    --------
    >>> # Continuing example from io.write_atm():
    >>> import pyratbay.io as io

    >>> atmfile = 'WASP-00b.atm'
    >>> units, specs, pressure, temp, q, rad = io.read_atm(atmfile)
    >>> print(units, specs, pressure, temp, q, rad, sep='\n')
    ('bar', 'kelvin', 'number', None)
    ['H2' 'He' 'H2O']
    [1.0000e-08 3.1623e-06 1.0000e-03 3.1623e-01 1.0000e+02]
    [1500. 1500. 1500. 1500. 1500.]
    [[8.499e-01 1.500e-01 1.000e-04]
     [8.499e-01 1.500e-01 1.000e-04]
     [8.499e-01 1.500e-01 1.000e-04]
     [8.499e-01 1.500e-01 1.000e-04]
     [8.499e-01 1.500e-01 1.000e-04]]
    None
    """
    atmfile = open(atmfile, "r")
    punits, runits, tunits, qunits, species = None, None, None, None, None

    while True:
        line = atmfile.readline().strip()
        # Stop when the per-layer data begins:
        if line == "@DATA":
            break
        # Skip empty and comment lines:
        elif line == '' or line.startswith('#'):
            pass
        # Extract units, and species from header:
        elif line == '@PRESSURE':
            punits = atmfile.readline().strip()
        elif line == '@RADIUS':
            runits = atmfile.readline().strip()
        elif line == '@TEMPERATURE':
            tunits = atmfile.readline().strip()
        elif line == '@ABUNDANCE':
            qunits = atmfile.readline().strip()
        elif line == '@SPECIES':
            species = np.asarray(atmfile.readline().strip().split())
        else:
            raise ValueError("Atmosphere file has unexpected line: \n'{:s}'".
                             format(line))

    if punits is None:
        raise ValueError("Atmospheric file does not have '@PRESSURE' header")
    if tunits is None:
        raise ValueError("Atmospheric file does not have '@TEMPERATURE' header")
    if qunits is None:
        raise ValueError("Atmospheric file does not have '@ABUNDANCE' header")
    if species is None:
        raise ValueError("Atmospheric file does not have '@SPECIES' header")

    nspecies = len(species)

    # Read first line to count number of columns:
    datastart = atmfile.tell()
    line = atmfile.readline()
    # Is there a column for the radius:
    rad = len(line.split()) - nspecies == 3

    if rad and runits is None:
        raise ValueError("Atmospheric file does not have '@RADIUS' header")

    # Count number of layers:
    nlayers = 1
    while True:
        line = atmfile.readline()
        if line == '' or line.startswith('#'):
            break
        nlayers += 1

    # Initialize arrays:
    if rad:
        radius = np.zeros(nlayers, np.double)
    else:
        radius = None
    press = np.zeros(nlayers, np.double)
    temp  = np.zeros(nlayers, np.double)
    q     = np.zeros((nlayers, nspecies), np.double)

    # Read table:
    atmfile.seek(datastart, 0)
    for i in np.arange(nlayers):
        data = atmfile.readline().split()
        if rad:
            radius[i] = data[0]
        press[i] = data[rad+0]
        temp [i] = data[rad+1]
        q    [i] = data[rad+2:]

    return (punits, tunits, qunits, runits), \
           species, press, temp, q, radius


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
  if not hasattr(pc, wlunits):
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


def read_pt(ptfile):
    r"""
    Read a pressure and temperature profile from a file.

    Parameters
    ----------
    ptfile: String
        Input file with pressure (in bars, first column) and temperature
        profiles (in Kelvin degree, second column).

    Returns
    -------
    pressure: 1D float ndarray
        Pressure profile in barye.
    temperature: 1D float ndarray
        Temperature profile in Kelvin.

    Examples
    --------
    >>> import pyratbay.io as io
    >>> ptfile = 'pt_profile.dat'
    >>> temp  = np.array([100.0, 150.0, 200.0, 175.0, 150.0])
    >>> press = np.array([1e-6,  1e-4,  1e-2,  1e0,   1e2])
    >>> with open(ptfile, 'w') as f:
    >>>     for p,t in zip(press, temp):
    >>>         f.write('{:.3e}  {:5.1f}\n'.format(p, t))
    >>> pressure, temperature = io.read_pt(ptfile)
    >>> for p,t in zip(pressure, temperature):
    >>>     print('{:.1e} barye  {:5.1f} K'.format(p, t))
    1.0e+00 barye  100.0 K
    1.0e+02 barye  150.0 K
    1.0e+04 barye  200.0 K
    1.0e+06 barye  175.0 K
    1.0e+08 barye  150.0 K
    """
    pressure, temperature = np.loadtxt(ptfile, usecols=(0,1), unpack=True)
    pressure *= pc.bar
    return pressure, temperature

