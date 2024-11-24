# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'save_pyrat',
    'load_pyrat',
    'write_atm',
    'read_atm',
    'write_spectrum',
    'read_spectrum',
    'write_spectra',
    'read_spectra',
    'write_opacity',
    'read_opacity',
    'write_pf',
    'read_pf',
    'write_cs',
    'read_cs',
    'read_pt',
    'write_observations',
    'read_observations',
    'read_atomic',
    'read_molecs',
    'read_isotopes',
    'import_xs',
    'import_tea',
    'export_pandexo',
]

from decimal import Decimal
import inspect
import os
import pickle

import numpy as np
import mc3
import h5py

from .. import constants as pc
from .. import tools as pt
from .. import spectrum as ps


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
    with pt.tmp_reset(
            pyrat,
            'voigt.profile', 'log.file',
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
    pyrat.set_spectrum()
    pyrat.log.verb = pyrat.verb
    # Recover MCMC posterior:
    if pt.isfile(pyrat.ret.mcmcfile) == 1:
        with np.load(pyrat.ret.mcmcfile) as mcmc:
            posterior, zchain, zmask = mc3.utils.burn(mcmc)
        pyrat.ret.posterior = posterior
    return pyrat


def write_atm(
        atmfile, pressure, temperature, species=None, abundances=None,
        radius=None, punits='bar', runits=None, header=None,
    ):
    r"""
    Write an atmospheric file following the Pyrat format.

    Parameters
    ----------
    atmfile: String
        Name of output atmospheric file.
    pressure: 1D float ndarray
        Monotonously decreasing pressure profile (in bar).
    temperature: 1D float ndarray
        Temperature profile for pressure layers (in Kelvin).
    species: 1D string ndarray
        List of atmospheric species.
    abundances: 2D float ndarray
        The species mole mixing ratio (of shape [nlayers,nspecies]).
    radius: 1D float ndarray
        Monotonously increasing radius profile (in cm).
    punits: String
        Pressure units of output.
    runits: String
        Radius units of output.
    header: String
        Header message (comment) to include at the top of the file.

    Examples
    --------
    >>> import numpy as np
    >>> import pyratbay.io as io
    >>> import pyratbay.atmosphere as pa

    >>> atmfile = 'WASP-00b.atm'
    >>> nlayers = 5
    >>> pressure = pa.pressure('1e-8 bar', '1e2 bar', nlayers)
    >>> temperature = pa.tmodels.Isothermal(pressure)(1500.0)
    >>> species = "H2 He H2O".split()
    >>> abundances = [0.8499, 0.15, 1e-4]
    >>> vmr = pa.uniform(abundances, nlayers)
    >>> io.write_atm(atmfile, pressure, temperature, species, vmr,
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
        f.write("\n# Radius      Pressure    Temperature  ")
    else:
        f.write("\n# Pressure  Temperature  ")
    if species is not None:
        f.write("".join([f"{mol:<14s}" for mol in species]))
    f.write("\n@DATA\n")

    pressure = pressure*pc.bar/pt.u(punits)
    if radius is not None:
        radius = radius/pt.u(runits)

    # Write data for each layer:
    nlayers = len(pressure)
    for i in range(nlayers):
        if radius is not None:
            f.write(f"{radius[i]:12.6e}  ")
        f.write(f"{pressure[i]:12e}  {temperature[i]:9.3f}  ")
        if species is not None:
            f.write("  ".join([f"{q:12.6e}" for q in abundances[i]]))
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
    vmr: 2D float ndarray
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
    ('bar', 'kelvin', 'volume', None)
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
    punits, runits, tunits, vmr_units, species = None, None, None, None, None

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
            vmr_units = atmfile.readline().strip()
        elif line == '@SPECIES':
            species = np.asarray(atmfile.readline().strip().split())
        else:
            raise ValueError(f"Atmosphere file has unexpected line: \n'{line}'")

    if punits is None:
        raise ValueError("Atmospheric file does not have '@PRESSURE' header")
    if tunits is None:
        raise ValueError("Atmospheric file does not have '@TEMPERATURE' header")

    if vmr_units is None and species is not None:
        raise ValueError("Atmospheric file does not have '@ABUNDANCE' header")
    if species is None and vmr_units is not None:
        raise ValueError("Atmospheric file does not have '@SPECIES' header")

    has_radius = runits is not None
    has_vmr = species is not None

    nspecies = len(species) if has_vmr else 0
    nrad = int(has_radius)

    # Read first line to count number of columns:
    datastart = atmfile.tell()
    line = atmfile.readline()
    ncolumns = len(line.split())

    if ncolumns == 3 + nspecies and runits is None:
        raise ValueError("Atmospheric file does not have '@RADIUS' header")

    if ncolumns != 2 + nrad + nspecies:
        rad_txt = ", 1 column for radius" if has_radius else ""
        q_txt = f", {nspecies} columns for abundances" if has_vmr else ""

        raise ValueError(
            f"Inconsistent number of columns ({ncolumns}) in '@DATA', "
             "expected 2 columns for temperature and pressure"
            f"{rad_txt}{q_txt}"
        )

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
    temp = np.zeros(nlayers, np.double)
    vmr = np.zeros((nlayers, nspecies), np.double) if has_vmr else None

    # Read table:
    atmfile.seek(datastart, 0)
    for i in range(nlayers):
        data = atmfile.readline().split()
        if has_radius:
            radius[i] = data[0]
        press[i] = data[nrad+0]
        temp [i] = data[nrad+1]
        if has_vmr:
            vmr[i] = data[nrad+2:]

    units = (punits, tunits, vmr_units, runits)
    return units, species, press, temp, vmr, radius


def write_spectra(spectra, wl, temperatures, filename):
    """
    Write flux spectra as function of wavelength and temperature to file.

    Parameters
    ----------
    spectra: 1D float iterable
        Flux spectra arrays (erg s-1 cm-2 cm units).
    wl: 1D float iterable
        Wavelength array in microns.
    temperatures: 1D float iterable
        Effective temperature array in Kelvin degrees.
    filename: String
        Output file name.
    """
    # Precision of 5 decimal places (or better if needed):
    min_diff = np.amin(np.abs(np.ediff1d(wl)))
    dec_place = Decimal(min_diff).adjusted()
    dec = np.clip(1 - dec_place, 5, 20)

    ntemps, nwave = np.shape(spectra)
    with open(filename, 'w') as f:
        # Write header:
        f.write(
            "# Temperature units: K\n"
            "# Wavelength units: um\n"
            "# Flux units: erg s-1 cm-2 cm\n"
        )
        f.write(f"\n@TEMPERATURES\n{'':{dec+4}}")
        for temp in temperatures:
            f.write(f"{temp:>15.1f}")
        f.write('\n\n@SPECTRA')
        # Write the spectrum values:
        for i in range(nwave):
            f.write(f"\n{wl[i]:>{dec+4}.{dec}f}")
            for j in range(ntemps):
                f.write(f"{spectra[j,i]:15.7e}")


def read_spectra(filename):
    """
    Write flux spectra as function of wavelength and temperature to file.

    Parameters
    ----------
    filename: String
        Input file name.

    Returns
    -------
    spectra: 1D float iterable
        Flux spectra arrays (erg s-1 cm-2 cm units).
    wn: 1D float iterable
        Wavenumber array in cm-1.
    temperatures: 1D float iterable
        Effective temperature array in Kelvin degrees.
    """
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f.readlines()]

    # Old style:
    if '@SPECTRA' not in lines:
        wn, spectra = read_spectrum(filename)
        temperatures = np.array([0])
        return spectra, wn, temperatures

    # Remove comments and blanks:
    lines = [
        line for line in lines
        if not line.startswith('#') and line != ''
    ]

    itemp = lines.index('@TEMPERATURES')
    temperatures = np.array(lines[itemp+1].split(), np.double)
    ntemps = len(temperatures)

    iflux = lines.index('@SPECTRA') + 1
    nwave = len(lines) - iflux
    data = np.zeros((nwave, ntemps+1))
    for i in range(nwave):
        data[i] = lines[i+iflux].split()
    spectra = data[:, 1:].T
    wn = 1.0/data[:,0]/pc.um

    return spectra, wn, temperatures


def write_spectrum(wl, spectrum, filename, type):
    """
    Write a spectrum to file.

    Parameters
    ----------
    wl: 1D float iterable
        Wavelength array in micron units.
    spectrum: 1D float iterable
        Spectrum array. (rp/rs)**2 for transmission (unitless),
        planetary flux for emission (erg s-1 cm-2 cm units).
    filename: String
        Output file name.
    type: String
        Data type (only used for header comments):
        - 'transit' for transit spectra
        - 'eclipse' for secondary eclipse spectra
        - 'emission' for emission flux
        - 'filter' for a instrumental filter transmission

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
        spectype  = "Fp/Fs"
        specunits = "unitless"
    elif type == "emission":
        spectype  = "Flux"
        specunits = "erg s-1 cm-2 cm"
    elif type == "filter":
        spectype  = "transmission"
        specunits = "unitless"
    else:
        raise ValueError(
            "Input 'type' argument must be 'transit', 'eclipse', "
            "'emission', or 'filter'"
        )

    # Precision of 5 decimal places (or better if needed):
    precision = -np.floor(np.log10(np.amin(np.abs(np.ediff1d(wl)))))
    precision = int(np.clip(precision+1, 5, np.inf))
    buff = precision + 5

    # Open-write file:
    with open(filename, 'w') as f:
        # Write header:
        f.write(f'# {"Wavelength":>{buff:d}s}   {spectype:>15s}\n')
        f.write(f"# {'um':>{buff:d}s}   {specunits:>15s}\n")
        # Write the spectrum values:
        for wave, flux in zip(wl, spectrum):
            f.write(f"{wave:>{buff+2:d}.{precision:d}f}   {flux:.9e}\n")


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
    >>> wl = np.linspace(1.1, 1.7, nwave)
    >>> spectrum = np.ones(nwave)
    >>> io.write_spectrum(wl, spectrum, 'sample_spectrum.dat', type='transit')
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
    wave, spectrum = np.loadtxt(filename, unpack=True)

    # Convert wavelength (um) to wavenumber (cm-1)
    if wn:
        wave = 1.0 / (wave*pc.um)

    return wave, spectrum



def write_opacity(ofile, species, temp, press, wn, opacity):
    """
    Write an opacity table as a binary npz file.

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
        Pressure array (bar).
    wn: 1D float ndarray
        Wavenumber array (cm-1).
    opacity: 4D float ndarray
        Tabulated opacities (cm2 molecule-1) of shape
        [nspec, ntemp, nlayers, nwave].
    """
    units = {
        'temperature': 'K',
        'pressure': 'bar',
        'wavenumber': 'cm-1',
        'cross section': 'cm2 molecule-1',
    }
    np.savez(
        ofile,
        species=species,
        temperature=temp,
        pressure=press,
        wavenumber=wn,
        opacity=opacity,
        units=units,
    )


def read_opacity(ofile, extract='all'):
    """
    Read an opacity table from file.
    Compatible with petitRADTRANS3 opacity files as well.

    Parameters
    ----------
    ofile: String
        Path to a Pyrat Bay opacity file.
    extract: String
        Information to extract, select between:
        - 'arrays' for the species, temperature, pressure, and wavenumber
        - 'opacity' for the opacity grid
        - 'all' to get both and the array dimensions.

    Returns
    -------
    sizes: 4-element integer tuple
        Sizes of the dimensions of the opacity table:
        (nspec, ntemp, nlayers, nwave)
    units: dict
        The physical units for the different quantities
    arrays: 4-element 1D ndarray tuple
        The dimensions of the opacity table:
        - species (string, the species names)
        - temperature (float, K)
        - pressure (float, bar)
        - wavenumber (float, cm-1)
    opacity: 4D float ndarray tuple
        The tabulated opacities (cm2 molecule-1), of shape
        [nspec, ntemp, nlayers, nwave].
    """
    if ofile.endswith('petitRADTRANS.h5'):
        with h5py.File(ofile, 'r') as f:
            species = list(f['mol_name'])
            species = np.array([species[0].decode('utf-8')])
            temp = np.array(f['t'])
            press = np.array(f['p'])
            wn = np.array(f['bin_edges'])
            if extract in ['opacity', 'all']:
                opacity = np.array(f['xsecarr'])
                # Same format as pyratbay files: (nmol, npress, ntemp, nwave)
                opacity = np.swapaxes(opacity, 0, 1)
                opacity = np.expand_dims(opacity, axis=0)
            units = {
                'temperature': 'K',
                'pressure': 'bar',
                'wavenumber': 'cm-1',
                'cross section': 'cm2 molecule-1',
            }
    else:
        with np.load(ofile, allow_pickle=True) as f:
            species = f['species']
            temp = f['temperature']
            press = f['pressure']
            wn = f['wavenumber']
            if extract in ['opacity', 'all']:
                opacity = f['opacity']
            units = np.ndarray.item(f['units']) if 'units' in f else None

    # If it does not have units, must be pyratbay<2.0, where pressures
    # were stored in barye units, thus, need to convert to bars
    if units is None:
        press /= pc.bar
        units = {
            'temperature': 'K',
            'pressure': 'bar',
            'wavenumber': 'cm-1',
            'cross section': 'cm2 molecule-1',
        }

    if extract == 'opacity':
        return opacity
    if extract == 'arrays':
        return (species, temp, press, wn)
    if extract == 'all':
        shape = np.shape(opacity)
        return (
            shape,
            units,
            (species, temp, press, wn),
            opacity,
        )


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
        f.write("@ISOTOPES\n            "
              + "  ".join(["{:13s}".format(iso) for iso in isotopes])
              + "\n\n")

        f.write("# Temperature (K), partition function for each isotope:\n")
        f.write("@DATA\n")
        for t, z in zip(temp, pf.T):
            f.write("  {:7.1f}   ".format(t)
                  + "  ".join("{:.7e}".format(d) for d in z) + "\n")


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
        Pressure profile in bar.
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
    >>>     print('{:.1e} bar  {:5.1f} K'.format(p, t))
    1.0e-06 bar  100.0 K
    1.0e-04 bar  150.0 K
    1.0e-02 bar  200.0 K
    1.0e+00 bar  175.0 K
    1.0e+02 bar  150.0 K
    """
    pressure, temperature = np.loadtxt(ptfile, usecols=(0,1), unpack=True)
    return pressure, temperature


def write_observations(
        obs_file, inst_names, wl, wl_half_width,
        depth=None, depth_err=None, depth_units='none',
    ):
    r"""
    Write an observation file for use in pyrat bay.
    These can be a combination of tophat pass bands (non-zero wl) or
    paths to files with tabulated passbands (zero wl value).

    Parameters
    ----------
    obs_file: String
        Name of observation file.
    inst_names: 1D string iterable
        Name for the bandpasses.  If this is of type str, use
        the same name for all bands.
    wl: 1D string ndarray
        Central wavelength of the bands (in microns).
        If wl is zero for a data point, assume that the inst_name
        is a file path to a tabulated pass band.
    wl_half_width: 2D float ndarray
        Bandpass half-width (in microns).
    data: 1D float ndarray
        Transit of eclipse depth corresponding to  each band.
        If not None, include this data into the observation file.
    depth_err: 1D float ndarray
        Depth uncertainties.
    depth_units:  String
        Units of input depth data.

    Examples
    --------
    >>> import pyratbay.io as io

    >>> # Observation file with only the pass bands:
    >>> wl = [2.144, 2.333, 2.523]
    >>> half_widths = [0.095, 0.095, 0.095]
    >>> io.write_observations('obs_file.txt', 'HST', wl, half_widths)

    >>> # Observation file with pass bands and data:
    >>> data = np.array([329.6, 344.5, 301.4])
    >>> uncert = np.array([20.4, 21.9, 23.5])
    >>> io.write_observations(
    >>>     'obs_file.txt', 'HST', wl, half_widths, data, uncert, 'ppm',
    >>> )
    """
    default_header = inspect.cleandoc(
        """
        # Passband info could be (1) a path to a filter file or
        # (2) a tophat filter defined by a central wavelength, half-width,
        # and optionally a name

        # Comment lines (like this one) and blank lines are ignored,
        # central-wavelength and half-width units are always microns

        # @DEPTH_UNITS sets the depth and uncert units (none, percent, ppt, ppm)
        # and also indicates that there's data and uncertainties to read
        # as two columns before the passband info
        """
    )
    # Consistency checks TBD:
    ndata = len(wl)
    if isinstance(inst_names, str):
         inst_names = np.tile(inst_names, ndata)

    has_data = depth is not None and depth_err is not None

    # Position of least significant digit
    wl_dec = [Decimal(width).adjusted() for width in wl_half_width]
    wl_dec += [
        Decimal(wl_diff).adjusted()
        for wl_diff in np.ediff1d(sorted(wl))
    ]
    wl_dec = np.clip(-np.amin(wl_dec) + 2, 1, 10)
    wl_len = wl_dec + 4

    depth_header1 = ''
    depth_header2 = ''
    if has_data:
        depth_dec = np.amin([Decimal(err).adjusted() for err in depth_err])
        depth_dec = np.clip(-depth_dec + 2, 1, 10)
        depth_len = depth_dec + 6
        depth_header1 = f'\n\n@DEPTH_UNITS\n{depth_units}'
        depth_header2 = 'depth    depth_err'

    with open(obs_file, 'w') as f:
        f.write(default_header)

        f.write(
            f'{depth_header1}'
            f'\n\n# {depth_header2}     wl  half_width    instrument\n@DATA\n'
        )
        for i in range(ndata):
            if has_data:
                f.write(
                    f'{depth[i]:{depth_len}.{depth_dec}f}  '
                    f'{depth_err[i]:{depth_len}.{depth_dec}f}   '
                )
            # Passband file
            if wl[i]==0.0:
                f.write(f'{inst_names[i]}\n')
            # Tophat
            else:
                f.write(
                    f'{wl[i]:{wl_len}.{wl_dec}f}  '
                    f'{wl_half_width[i]:{wl_len}.{wl_dec}f}    '
                    f'{inst_names[i]}\n'
                )


def read_observations(obs_file):
    r"""
    Read an observations file.

    Parameters
    ----------
    obs_file: String
        Path to file containing observations info, see Notes below.

    Returns
    -------
    filters: List
        Filter passband objects.
    depth: 1D string list
        The transit or eclipse depths for each filter.
    depth_err: 1D float ndarray
        The depth uncertainties.

    Notes
    -----
    An obs_file contains passband info (indicated by a '@DATA' flag),
    one line per passband. The passband info could be:
    (1) a path to a file containing the spectral response, or
    (2) a tophat filter defined by a central wavelength, half-width,
    and optionally a name.

    A @DEPTH_UNITS flag sets the depth and uncert units (which can be
    set to: none, percent, ppt, ppm).
    This flag also indicates that there's data and uncerts to read
    as two columns before the passband info.

    Comment and blank lines are ignored,
    central-wavelength and half-width units are always microns.

    Examples
    --------
    >>> import pyratbay.io as io
    >>> # File including depths and uncertainties:
    >>> obs_file = 'observations.dat'
    >>> bands, depth, depth_err = io.read_observations(obs_file)

    >>> # File including only the passband info:
    >>> obs_file = 'filters.dat'
    >>> bands = io.read_observations(obs_file)
    """
    if not os.path.isfile(obs_file):
        raise ValueError(f"Observation file '{obs_file}' does not exist")

    # Skip comment and blank lines:
    lines = []
    for line in open(obs_file, 'r'):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        lines.append(line)

    if '@DATA' not in lines:
        raise ValueError("Observation file does not have a '@DATA' header")

    # Number of header lines (to skip when reading the tabulated data):
    i = 0
    depth_units = None
    nlines = len(lines)
    for i in range(nlines):
        line = lines[i]
        # Stop when the tabulated data begins:
        if line == "@DATA":
            i += 1
            break
        # Get the sampled temperatures:
        elif line == "@DEPTH_UNITS":
            depth_units = lines[i+1]

    nobs = nlines - i
    has_data = depth_units is not None

    filters = []
    depth = np.zeros(nobs)
    depth_err = np.zeros(nobs)

    # Read the data:
    for j in range(nobs):
        info = lines[i+j].split()
        ndata = 0
        if has_data:
            depth[j] = info[0]
            depth_err[j] = info[1]
            ndata = 2

        if len(info) - ndata == 1:
            filter_file = info[ndata].replace('{ROOT}', pc.ROOT)
            filter_file = filter_file.replace('{FILTERS}', pc.FILTERS)
            filter_file = os.path.realpath(filter_file)
            filters.append(ps.PassBand(filter_file))
        elif len(info) - ndata == 2:
            wl0, half_width = np.array(info[ndata:ndata+2], np.double)
            filters.append(ps.Tophat(wl0, half_width))
        elif len(info) - ndata == 3:
            wl0, half_width = np.array(info[ndata:ndata+2], np.double)
            name = info[ndata+2]
            filters.append(ps.Tophat(wl0, half_width, name=name))
        else:
            error_msg = 'Invalid number of values in obs_file'
            if has_data and len(info) in [1,2]:
                error_msg += (
                    ", perhaps remove the '@DEPTH_UNITS' flag if there's no "
                    "depth/uncert data"
                )
            elif not has_data and len(info) in [4,5]:
                error_msg += ", perhaps the '@DEPTH_UNITS' flag is missing"
            raise ValueError(error_msg)

    if has_data:
        depth *= pt.u(depth_units)
        depth_err *= pt.u(depth_units)
        return filters, depth, depth_err
    return filters


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
    atomic_num = np.zeros(nelements, int)
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
    Read a molecules file to extract their names, masses, and radii.
    The output also includes the ions denoted by a '-' and '+'
    character appended at the end of the species names.

    Parameters
    ----------
    file: String
        The molecule file path.

    Returns
    -------
    names: 1D string ndarray
        The molecules' names.
    masses: 1D float ndarray
        The mass of the molecules (in g mol-1).
    radii: 1D float ndarray
        The collisional radius of the molecules (in angstrom).

    Notes
    -----
    In all truthfulness, these are species, not only molecules, as the
    file also contain elemental particles.

    Examples
    --------
    >>> import pyratbay.io as io
    >>> import pyratbay.constants as pc
    >>> names, masses, radii = io.read_molecs(
    >>>     pc.ROOT+'pyratbay/data/molecules.dat')
    >>> names = list(names)
    >>> print(f"H2O: mass = {masses[names.index('H2O')]} g mol-1, "
    >>>       f"radius = {radii[names.index('H2O')]} angstrom.")
    H2O: mass = 18.015 g mol-1, radius = 1.6 Angstrom.
    """
    names = [] # Molecule names
    masses = [] # Molecule masses
    radii = [] # Molecule radii

    for line in open(file, 'r'):
        # Skip comment and blank lines:
        if line.strip() == '' or line.strip().startswith('#'):
            continue
        info = line.split()
        names.append(info[0])
        masses.append(info[1])
        radii.append(info[2])

    electron_index = names.index('e-')
    names.pop(electron_index)
    e_mass = masses.pop(electron_index)
    e_radius = radii.pop(electron_index)

    names = np.array(
        names
        + ['e-']
        + [f'{name}-' for name in names]
        + [f'{name}+' for name in names])
    masses = np.array(
        masses + [e_mass] + masses + masses,
        np.double)
    radii = np.array(
        radii + [e_radius] + radii + radii,
        np.double)

    return names, masses, radii


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
        Pressure sample of the opacity file (in bars)
    temperature: 1D float ndarray
        Temperature sample of the opacity file (in Kelvin degrees).
    wavenumber: 1D float ndarray
        Wavenumber sample of the opacity file (in cm-1).
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
                pressure = np.array(xs_data['p'])
                temperature = np.array(xs_data['t'])
                wavenumber = np.array(xs_data['bin_edges'])
                species = [xs_data['mol_name'][0]]

    elif source == 'taurex':
        with open(filename, 'rb') as f:
            xs_data = pickle.load(f)
            xs = xs_data['xsecarr']
            if read_all or ofile is not None:
                pressure = xs_data['p']
                temperature = xs_data['t']
                wavenumber = xs_data['wno']
                species = [xs_data['name']]

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
    pressure = np.zeros(nlayers)
    abundance = np.zeros((nlayers, nreqspecies))
    for i in range(nlayers):
        data = np.asarray(tea[idata+i].split(), np.double)
        pressure[i], temperature[i] = data[0:2]
        abundance[i,:] = data[2:][sindex]

    # TEA pressure units are always bars:
    punits = "bar"
    # File header:
    header = "# TEA atmospheric file formatted for Pyrat.\n\n"

    write_atm(
        atmfile, pressure, temperature, req_species, abundance,
        punits=punits, header=header,
    )


def export_pandexo(
        pyrat, baseline, transit_duration,
        Vmag=None, Jmag=None, Hmag=None, Kmag=None, metal=0.0,
        instrument=None, n_transits=1, resolution=None,
        noise_floor=0.0, sat_level=80.0,
        save_file=True,
    ):
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
    exo_dict['star']['logg'] = pyrat.phy.log_gstar
    exo_dict['star']['radius'] = pyrat.phy.rstar/pc.rsun
    exo_dict['star']['r_unit'] = 'R_sun'

    if pyrat.od.rt_path == 'transit':
        exo_dict['planet']['f_unit'] = 'rp^2/r*^2'
        spectrum = pyrat.spec.spectrum
    elif pyrat.od.rt_path == 'emission':
        exo_dict['planet']['f_unit'] = 'fp/f*'
        rprs = pyrat.atm.rplanet/pyrat.phy.rstar
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
        num_cores=pyrat.ncpu,
    )

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

