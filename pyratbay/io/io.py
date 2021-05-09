# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'save_pyrat',
    'load_pyrat',
    'write_atm',
    'read_atm',
    'write_spectrum',
    'read_spectrum',
    'write_opacity',
    'read_opacity',
    'write_pf',
    'read_pf',
    'write_cs',
    'read_cs',
    'read_pt',
    'read_atomic',
    'read_molecs',
    'read_isotopes',
    'import_xs',
    'import_tea',
    'export_pandexo',
    ]

import os
import pickle

import numpy as np
import mc3

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
        print(f'Saving pyrat instance to: {pfile}')
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
        with np.load(pyrat.ret.mcmcfile) as mcmc:
            posterior, zchain, zmask = mc3.utils.burn(mcmc)
        pyrat.ret.posterior = posterior
    return pyrat


def write_atm(atmfile, pressure, temperature, species=None, abundances=None,
        radius=None, punits='bar', runits=None, header=None):
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
    radius: 1D float ndarray
        Monotonously increasing radius profile (in cm).
    punits:  String
        Pressure units of output.
    runits:  String
        Radius units of output.
    header:  String
        Header message (comment) to include at the top of the file.

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
    # Pressure units:
    @PRESSURE
    bar
    # Temperatures units:
    @TEMPERATURE
    kelvin
    # Abundance units (mixing ratio):
    @ABUNDANCE
    volume

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
    if (species is None) != (abundances is None):
        raise ValueError('Both species and abundances must be defined')

    f = open(atmfile, "w")
    if header is not None:
        f.write(header)

    # Set the values units:
    f.write(f"# Pressure units:\n@PRESSURE\n{punits}\n")
    f.write("# Temperatures units:\n@TEMPERATURE\nkelvin\n")
    if radius is not None:
        f.write(f"# Radius units:\n@RADIUS\n{runits}\n")
    # At the moment, only take volume MR (mass if theres's popular demand)
    if species is not None:
        f.write("# Abundance units (mixing fraction):\n@ABUNDANCE\nvolume\n")
        # Species names:
        f.write("\n# Atmospheric composition:\n@SPECIES\n" +
                "  ".join([f"{mol:<s}" for mol in species]) + '\n')

    # Write the per-layer data:
    if radius is not None:
        f.write("\n# Radius    Pressure    Temperature  ")
    else:
        f.write("\n# Pressure  Temperature  ")
    if species is not None:
        f.write("".join([f"{mol:<14s}" for mol in species]))
    f.write("\n@DATA\n")

    pressure = pressure/pt.u(punits)
    if radius is not None:
        radius = radius/pt.u(runits)

    # Write data for each layer:
    nlayers = len(pressure)
    for i in range(nlayers):
        if radius is not None:
            f.write(f"{radius[i]:10.4e}  ")
        f.write(f"{pressure[i]:10.4e}  {temperature[i]:11.3f}  ")
        if species is not None:
            f.write(f"  ".join([f"{q:12.6e}" for q in abundances[i]]))
        f.write('\n')
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
            continue
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
            raise ValueError(f"Atmosphere file has unexpected line: \n'{line}'")

    if punits is None:
        raise ValueError("Atmospheric file does not have '@PRESSURE' header")
    if tunits is None:
        raise ValueError("Atmospheric file does not have '@TEMPERATURE' header")

    if qunits is None and species is not None:
        raise ValueError("Atmospheric file does not have '@ABUNDANCE' header")
    if species is None and qunits is not None:
        raise ValueError("Atmospheric file does not have '@SPECIES' header")

    has_radius = runits is not None
    has_q = species is not None

    nspecies = len(species) if has_q else 0
    nrad = int(has_radius)

    # Read first line to count number of columns:
    datastart = atmfile.tell()
    line = atmfile.readline()
    ncolumns = len(line.split())

    if ncolumns == 3 + nspecies and runits is None:
        raise ValueError("Atmospheric file does not have '@RADIUS' header")

    if ncolumns != 2 + nrad + nspecies:
        rad_txt = ", 1 column for radius" if has_radius else ""
        q_txt = f", {nspecies} columns for abundances" if has_q else ""

        raise ValueError(
            f"Inconsistent number of columns ({ncolumns}) in '@DATA', "
             "expected 2 columns for temperature and pressure"
            f"{rad_txt}{q_txt}")

    # Count number of layers:
    nlayers = 1
    while True:
        line = atmfile.readline()
        if line == '' or line.startswith('#'):
            break
        nlayers += 1

    # Initialize arrays:
    radius = np.zeros(nlayers, np.double) if has_radius else None
    press = np.zeros(nlayers, np.double)
    temp  = np.zeros(nlayers, np.double)
    q = np.zeros((nlayers, nspecies), np.double) if has_q else None

    # Read table:
    atmfile.seek(datastart, 0)
    for i in np.arange(nlayers):
        data = atmfile.readline().split()
        if has_radius:
            radius[i] = data[0]
        press[i] = data[nrad+0]
        temp [i] = data[nrad+1]
        if has_q:
            q[i] = data[nrad+2:]

    return (punits, tunits, qunits, runits), \
           species, press, temp, q, radius


def write_spectrum(wl, spectrum, filename, type, wlunits='um'):
  """
  Write a spectrum to file.

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
      - 'transit' for transmission
      - 'emission' for emission
      - 'filter' for a instrumental filter transmission
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
  elif type == "emission":
      spectype  = "Flux"
      specunits = "erg s-1 cm-2 cm"
  elif type == "filter":
      spectype  = "transmission"
      specunits = "unitless"
  else:
      raise ValueError(
          "Input 'type' argument must be 'transit', 'emission', or 'filter'.")

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


def write_opacity(ofile, species, temp, press, wn, opacity):
    """
    Write an opacity table as a binary file.

    Parameters
    ----------
    ofile: String
        Output filename where to save the opacity data.
        File extension must be .npz
    species: 1D string iterable
        Species names.
    temp: 1D float ndarray
        Temperature array (Kelvin degree).
    press: 1D float ndarray
        Pressure array (barye).
    wn: 1D float ndarray
        Wavenumber array (cm-1).
    opacity: 4D float ndarray
        Tabulated opacities (cm2 molecule-1) of shape
        [nspec, ntemp, nlayers, nwave].
    """
    np.savez(ofile,
        species=species, temperature=temp, pressure=press, wavenumber=wn,
        opacity=opacity)


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
        (nspec, ntemp, nlayers, nwave)
    arrays: 4-element 1D ndarray tuple
        The dimensions of the opacity table:
        - species     (string, the species names)
        - temperature (float, Kelvin)
        - pressure    (float, barye)
        - wavenumber  (float, cm-1)
    opacity: 4D float ndarray tuple
        The tabulated opacities (cm2 molecule-1), of shape
        [nspec, ntemp, nlayers, nwave].
    """
    with np.load(ofile) as f:
        species = f['species']
        temp    = f['temperature']
        press   = f['pressure']
        wn      = f['wavenumber']
        opacity = f['opacity']

    # Arrays lengths:
    nspec = len(species)
    ntemp = len(temp)
    nlayers = len(press)
    nwave = len(wn)

    return ((nspec, ntemp, nlayers, nwave),
            (species, temp, press, wn),
            opacity)


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


def read_atomic(afile):
    """
    Read an elemental (atomic) composition file.

    Parameters
    ----------
    afile: String
        File with atomic composition.

    Returns
    -------
    atomic_num: 1D integer ndarray
        Atomic number (except for Deuterium, which has anum=0).
    symbol: 1D string ndarray
        Elemental chemical symbol.
    dex: 1D float ndarray
        Logarithmic number-abundance, scaled to log(H) = 12.
    name: 1D string ndarray
        Element names.
    mass: 1D float ndarray
        Elemental mass in amu.

    Uncredited developers
    ---------------------
    Jasmina Blecic
    """
    # Allocate arrays:
    nelements = 84  # Fixed number
    atomic_num = np.zeros(nelements, np.int)
    symbol = np.zeros(nelements, '|U2')
    dex    = np.zeros(nelements, np.double)
    name   = np.zeros(nelements, '|U20')
    mass   = np.zeros(nelements, np.double)

    # Open-read file:
    with open(afile, 'r') as f:
        # Read-discard first two lines (header):
        f.readline()
        f.readline()
        # Store data into the arrays:
        for i in range(nelements):
            atomic_num[i], symbol[i], dex[i], name[i], mass[i] = \
                f.readline().split()

    return atomic_num, symbol, dex, name, mass


def read_molecs(file):
    r"""
    Read a molecules file to extract their symbol, mass, and diameter.

    Parameters
    ----------
    file: String
        The molecule file path.

    Returns
    -------
    symbol: 1D string ndarray
        The molecule's name.
    mass: 1D float ndarray
        The mass of the molecules (in g mol-1).
    diam: 1D float ndarray
        The collisional diameter of the molecules (in Angstrom).

    Notes
    -----
    In all truthfulness, these are species, not only molecules, as the
    file also contain elemental particles.

    Examples
    --------
    >>> import pyratbay.io as io
    >>> import pyratbay.constants as pc
    >>> names, mass, diam = io.read_molecs(pc.ROOT+'pyratbay/data/molecules.dat')
    >>> names = list(names)
    >>> print(f"H2O: mass = {mass[names.index('H2O')]} g mol-1, "
    >>>       f"diameter = {diam[names.index('H2O')]} Angstrom.")
    H2O: mass = 18.01528 g mol-1, diameter = 3.2 Angstrom.
    """
    symbol = [] # Molecule symbol
    mass   = [] # Molecule mass
    diam   = [] # Molecule diameter

    for line in open(file, 'r'):
        # Skip comment and blank lines:
        if line.strip() == '' or line.strip().startswith('#'):
            continue
        info = line.split()
        symbol.append(info[0])
        mass.append(info[1])
        diam.append(info[2])

    symbol = np.asarray(symbol)
    mass = np.asarray(mass, np.double)
    diam = np.asarray(diam, np.double)

    return symbol, mass, diam


def read_isotopes(file):
    r"""
    Read an isotopes file to extract their molecule, hitran name,
    exomol name, isotopic ratio, and mass.

    Parameters
    ----------
    file: String
        The isotope file path.

    Returns
    -------
    mol_ID: 1D integer ndarray
        HITRAN molecule ID.
    mol: 1D string ndarray
        Molecule names.
    hitran_iso: 1D string ndarray
        Isotope name as in HITRAN database.
    exomol_iso: 1D string ndarray
        Isotope name based on exomol database.
    iso_ratio: 1D float ndarray
        Isotopic ratios.
    iso_mass: 1D float ndarray
        The mass of the molecules (in g mol-1).

    Examples
    --------
    >>> import pyratbay.io as io
    >>> import pyratbay.constants as pc
    >>> ID, mol, hit_iso, exo_iso, ratio, mass = \
    >>>     io.read_isotopes(pc.ROOT+'pyratbay/data/isotopes.dat')
    >>> print("H2O isotopes:\n iso    iso    isotopic  mass"
    >>>                    "\n hitran exomol ratio     g/mol")
    >>> for i in range(len(mol)):
    >>>     if mol[i] == 'H2O':
    >>>         print(f" {hit_iso[i]:6} {exo_iso[i]:6} "
    >>>               f"{ratio[i]:.3e} {mass[i]:.4f}")
    H2O isotopes:
    iso    iso    isotopic  mass
    hitran exomol ratio     g/mol
    161    116    9.973e-01 18.0106
    181    118    1.999e-03 20.0148
    171    117    3.719e-04 19.0148
    162    126    3.107e-04 19.0168
    182    000    6.230e-07 21.0211
    172    000    1.158e-07 20.0211
    262    226    2.420e-08 20.0210
    282    000    0.000e+00 22.0000
    272    000    0.000e+00 21.0000
    """
    mol_ID, mol, hitran_iso, exomol_iso, ratio, mass = [], [], [], [], [], []
    for line in open(file, 'r'):
        # Skip comment and blank lines:
        if line.strip() == '' or line.strip().startswith('#'):
            continue
        info = line.split()
        mol_ID.append(info[0])
        mol.append(info[1])
        hitran_iso.append(info[2])
        exomol_iso.append(info[3])
        ratio.append(info[4])
        mass.append(info[5])

    mol_ID = np.asarray(mol_ID, int)
    mol = np.asarray(mol)
    hitran_iso = np.asarray(hitran_iso)
    exomol_iso = np.asarray(exomol_iso)
    ratio = np.asarray(ratio, np.double)
    mass = np.asarray(mass, np.double)

    return mol_ID, mol, hitran_iso, exomol_iso, ratio, mass


def import_xs(filename, source, read_all=True, ofile=None):
    """
    Read a cross-section opacity file from an external source.

    Parameters
    ----------
    filename: String
        The opacity pickle file to read.
    source: String
        The cross-section source: exomol or taurex (see note below).
    read_all: Bool
        If True, extract all contents in the file: cross-section,
        pressure, temperature, and wavenumber.
        If False, extract only the cross-section data.
    ofile: String
        If not None, store Exomol XS data into a Pyratbay opacity
        format.

    Returns
    -------
    xs: 3D float ndarray
        Opacity cross-section in cm2 molecule-1.
        with shape [npress, ntemp, nwave].
    pressure: 1D float ndarray
        Pressure sample of the opacity file (in barye units)
    temperature: 1D float ndarray
        Temperature sample of the opacity file (in Kelvin degrees units).
    wavenumber: 1D float ndarray
        Wavenumber sample of the opacity file (in cm-1 units).
    species: String
        The species name.

    Notes
    -----
    - exomol cross sections (Chubb et al. 2020, AA) can be found here:
    http://www.exomol.com/data/data-types/opacity/
    - taurex cross sections (Al-Refaie et al. 2019) can be found here:
    https://taurex3-public.readthedocs.io/en/latest/user/taurex/quickstart.html

    Examples
    --------
    >>> # For this example, you'll need to have/download the following
    >>> # file into the current folder:
    >>> # http://www.exomol.com/db/H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL__R15000_0.3-50mu.xsec.TauREx.h5
    >>> import pyratbay.io as io
    >>> filename = '1H2-16O__POKAZATEL__R15000_0.3-50mu.xsec.TauREx.h5'
    >>> xs_H2O, press, temp, wn, species = io.import_xs(filename, 'exomol')
    """
    try:
        import h5py
    except ModuleNotFoundError as e:
        if source == 'exomol':
            raise e

    if source == 'exomol':
        with h5py.File(filename, 'r') as xs_data:
            xs = np.array(xs_data['xsecarr'])
            if read_all or ofile is not None:
                pressure    = np.array(xs_data['p']) * pc.bar
                temperature = np.array(xs_data['t'])
                wavenumber  = np.array(xs_data['bin_edges'])
                species     = [xs_data['mol_name'][0]]

    elif source == 'taurex':
        with open(filename, 'rb') as f:
            xs_data = pickle.load(f)
            xs = xs_data['xsecarr']
            if read_all or ofile is not None:
                pressure    = xs_data['p'] * pc.bar
                temperature = xs_data['t']
                wavenumber  = xs_data['wno']
                species     = [xs_data['name']]

    else:
        raise ValueError("Invalid cross-section source type.")

    if ofile is not None:
        nlayers, ntemp, nwave = np.shape(xs)
        xs_pb = np.swapaxes(xs,0,1).reshape(1,ntemp,nlayers,nwave)
        write_opacity(ofile, species, temperature, pressure, wavenumber, xs_pb)

    if read_all:
        return xs, pressure, temperature, wavenumber, species[0]
    return xs


def import_tea(teafile, atmfile, req_species=None):
    """
    Format a TEA atmospheric file into a Pyrat atmospheric file.

    Paramters
    ---------
    teafile:  String
        Input TEA atmospheric file.
    atmfile:  String
        Output Pyrat atmospheric file.
    req_species: List of strings
        The requested species for output.  If None, request all species
        in teafile.
    """
    # Open read the TEA file:
    with open(teafile, "r") as f:
        tea = np.asarray(f.readlines())

    # Line-indices where the species and data are:
    ispec = np.where(tea == "#SPECIES\n")[0][0] + 1  # Species list
    idata = np.where(tea == "#TEADATA\n")[0][0] + 2  # data starting line

    # TEA--Pyrat molecules names dictionary:
    tea_to_pyrat = {}
    for line in open(pc.ROOT+"pyratbay/data/TEA_gdata_defaults.txt", "r"):
        pyrat_name, tea_name = line.split()
        tea_to_pyrat[tea_name] = pyrat_name

    # Read and clean species names:
    species = tea[ispec].split()
    nspecies = len(species)
    for i in range(nspecies):
        species[i] = tea_to_pyrat[species[i]]

    if req_species is None:
        req_species = species

    # Species indices corresponding to req_species:
    nreqspecies = len(req_species)
    sindex = np.zeros(len(req_species), int)
    for i in range(len(req_species)):
        sindex[i] = species.index(req_species[i])

    # Extract per-layer data:
    nlayers = len(tea) - idata
    temperature = np.zeros(nlayers)
    pressure    = np.zeros(nlayers)
    abundance   = np.zeros((nlayers, nreqspecies))
    for i in range(nlayers):
        data = np.asarray(tea[idata+i].split(), np.double)
        pressure[i], temperature[i] = data[0:2]
        abundance[i,:] = data[2:][sindex]

    # TEA pressure units are always bars:
    punits = "bar"
    pressure = pressure * pt.u(punits)  # Set in barye units
    # File header:
    header = "# TEA atmospheric file formatted for Pyrat.\n\n"

    write_atm(atmfile, pressure, temperature, req_species, abundance,
        punits=punits, header=header)


def export_pandexo(pyrat, baseline, transit_duration,
    Vmag=None, Jmag=None, Hmag=None, Kmag=None, metal=0.0,
    instrument=None, n_transits=1, resolution=None,
    noise_floor=0.0, sat_level=80.0,
    save_file=True):
    """
    Parameters
    ----------
    pyrat: A Pyrat instance
        Pyrat object from which to extract the system physical properties.
    baseline: Float or string
        Total observing time in sec (float) or with given units (string).
    transit_duration: Float or string
        Transit/eclipse duration in sec (float) or with given units (string).
    metal: Float
        Stellar metallicity as log10(Fe/H).
    Vmag: Float
        Stellar magnitude in the Johnson V band.
        Only one of Vmag, Jmag, Hmag, or Kmag should be defined.
    Jmag: Float
        Stellar magnitude in the Johnson J band.
        Only one of Vmag, Jmag, Hmag, or Kmag should be defined.
    Hmag: Float
        Stellar magnitude in the Johnson H band.
        Only one of Vmag, Jmag, Hmag, or Kmag should be defined.
    Kmag: Float
        Stellar magnitude in the Johnson Kband.
        Only one of Vmag, Jmag, Hmag, or Kmag should be defined.
    instrument: String or list of strings or dict
        Observing instrument to simulate.
        If None, this function returns the input dictionary.
    n_transits: Integer
        Number of transits/eclipses.
    resolution: Float
        Approximate output spectral sampling R = 0.5*lambda/delta-lambda.
    sat_level: Float
        Saturation level in percent of full well.
    noise_floor: Float or string
        Noise-floor level in ppm at all wavelengths (if float) or
        wavelength dependent (if string, filepath).
    save_file: Bool or string
        If string, store pandexo output pickle file with this filename.
        If True, store pandexo output with default name based on
        the pyrat object's output filename.

    Returns
    -------
    pandexo_sim: dict
        Output from pandexo.engine.justdoit.run_pandexo().
        Note this dict has R=None, noccultations=1 (as suggested in pandexo).
    wavelengths: List of 1D float arrays
        Wavelengths of simulated observed spectra for each instrument.
        Returned only if instrument is not None.
    spectra: List of 1D float arrays
        Simulated observed spectra for each instrument.
        Returned only if instrument is not None.
    uncertainties: List of 1D float arrays
        Uncertainties of simulated observed spectra for each instrument.
        Returned only if instrument is not None.

    Examples
    --------
    >>> import pyratbay as pb
    >>> import pyratbay.io as io

    >>> pyrat = pb.run('demo_spectrum-transmission.cfg')
    >>> instrument = 'NIRCam F322W2'
    >>> #instrument = jdi.load_mode_dict(instrument)
    >>> baseline = '4.0 hour'
    >>> transit_duration = '2.0 hour'
    >>> resolution = 100.0
    >>> n_transits = 2
    >>> Jmag = 8.0
    >>> metal = 0.0

    >>> pandexo_sim, wls, spectra, uncerts = io.export_pandexo(
    >>>     pyrat, baseline, transit_duration,
    >>>     n_transits=n_transits,
    >>>     resolution=resolution,
    >>>     instrument=instrument,
    >>>     Jmag=Jmag,
    >>>     metal=metal)
    """
    import pandexo.engine.justdoit as jdi
    import pandexo.engine.justplotit as jpi

    if isinstance(baseline, str):
        baseline = pt.get_param(baseline)
    if isinstance(transit_duration, str):
        transit_duration = pt.get_param(transit_duration)

    ref_wave = {
        'Vmag':0.55,
        'Jmag':1.25,
        'Hmag':1.6,
        'Kmag':2.22,
        }
    mags = {
        'Vmag':Vmag,
        'Jmag':Jmag,
        'Hmag':Hmag,
        'Kmag':Kmag,
        }
    mag = {key:val for key,val in mags.items() if val is not None}
    if len(mag) != 1:
        raise ValueError(
            f'Exactly one of {list(mags.keys())} should be defined')
    band_mag, mag = mag.popitem()

    exo_dict = jdi.load_exo_dict()

    exo_dict['observation']['sat_level'] = sat_level
    exo_dict['observation']['sat_unit'] = '%'
    exo_dict['observation']['noccultations'] = 1
    exo_dict['observation']['R'] = None

    exo_dict['planet']['transit_duration'] = transit_duration
    exo_dict['planet']['td_unit'] = 's'
    exo_dict['observation']['baseline'] = baseline
    exo_dict['observation']['baseline_unit'] = 'total'
    exo_dict['observation']['noise_floor'] = noise_floor

    # Stellar flux from erg s-1 cm-2 cm to erg s-1 cm-2 Hz-1:
    starflux = {'f':pyrat.spec.starflux/pc.c, 'w':1.0/pyrat.spec.wn}
    exo_dict['star']['type'] = 'user'
    exo_dict['star']['starpath'] = starflux
    exo_dict['star']['w_unit'] = 'cm'
    exo_dict['star']['f_unit'] = 'erg/cm2/s/Hz'
    exo_dict['star']['mag'] = mag
    exo_dict['star']['ref_wave'] = ref_wave[band_mag]

    exo_dict['star']['temp'] = pyrat.phy.tstar
    exo_dict['star']['metal'] = metal
    exo_dict['star']['logg'] = np.log10(pyrat.phy.gstar)
    exo_dict['star']['radius'] = pyrat.phy.rstar/pc.rsun
    exo_dict['star']['r_unit'] = 'R_sun'

    if pyrat.od.path == 'transit':
        exo_dict['planet']['f_unit'] = 'rp^2/r*^2'
        spectrum = pyrat.spec.spectrum
    elif pyrat.od.path == 'eclipse':
        exo_dict['planet']['f_unit'] = 'fp/f*'
        rprs = pyrat.phy.rplanet/pyrat.phy.rstar
        spectrum = pyrat.spec.spectrum/pyrat.spec.starflux * rprs**2

    exo_dict['planet']['type'] ='user'
    exo_dict['planet']['exopath'] = {'f':spectrum, 'w':1.0/pyrat.spec.wn}
    exo_dict['planet']['w_unit'] = 'cm'
    exo_dict['planet']['radius'] = pyrat.phy.rplanet
    exo_dict['planet']['r_unit'] = 'cm'

    if instrument is None:
        return exo_dict

    if isinstance(instrument, str):
        instrument = [instrument]

    if save_file is True:
        output_path = os.path.dirname(pyrat.log.logname)
        output_file = os.path.basename(pyrat.log.logname).replace(
            '.log', '_pandexo.p')
    elif isinstance(save_file, str):
        output_path = os.path.dirname(save_file)
        output_file = os.path.basename(save_file)
        save_file = True
    else:
        save_file = False

    pandexo_sim = jdi.run_pandexo(
        exo_dict,
        instrument,
        save_file=save_file,
        output_path=output_path,
        output_file=output_file,
        num_cores=pyrat.ncpu)

    if isinstance(pandexo_sim, list):
        pandexo_sim = [sim[list(sim.keys())[0]] for sim in pandexo_sim]
    else:
        pandexo_sim = [pandexo_sim]

    wavelengths, spectra, uncerts = [], [], []
    for sim in pandexo_sim:
        wl, spec, unc = jpi.jwst_1d_spec(
            sim, R=resolution, num_tran=n_transits, plot=False)
        wavelengths += wl
        spectra += spec
        uncerts += unc

    return pandexo_sim, wavelengths, spectra, uncerts

