# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'get_tips_molname',
    'check_exomol_files',
    'tips',
    'exomol_pf',
    'exomol_states',
    'kurucz',
]

import pickle
import bz2

import numpy as np
from scipy.interpolate import CubicSpline
import mc3.utils as mu

from ... import io as io
from ... import constants as pc
from ... import tools as pt


def get_tips_molname(molID):
    """
    Get the TIPS molecule name for given molecule ID.

    Parameters
    ----------
    molID: Integer
        HITRAN molecule ID. See for example: https://hitran.org/lbl/

    Returns
    -------
    molname: String
        Name of molecule.

    Examples
    --------
    >>> import pyratbay.opacity.partitions as pf
    >>> print(pf.get_tips_molname(1), pf.get_tips_molname(6))
    H2O CH4
    """
    with open(pc.ROOT+'pyratbay/data/tips_2021.pkl', 'rb') as p:
        data = pickle.load(p)
    if molID not in data['mol_ID']:
        raise ValueError(
            f'TIPS 2021 database does not contain molecule ID: {molID}')
    return data['mol_ID'][molID]


def tips(molecule, isotopes=None, outfile=None, db_type='as_tips'):
    """
    Extract TIPS 2021 partition-function values for given molecule.
    If requested, write the partition-function into a file for use
    with Pyrat Bay.
    References:
        Gamache et al. (2017), JQSRT, 203, 70.
        Gamache et al. (2021), JQSRT, 271, 107713.

    Parameters
    ----------
    molecule: String
        Name of the molecule.
    isotopes: String or list of strings
        If not None, only extract the requested isotopes.
    outfile: String
        If not None, save output to file.
        If outfile == 'default', save output to file named as
        PF_tips_molecule.dat
    db_type: String
        If db_type == 'as_exomol', return isotopic names following
        the exomol notation.

    Returns
    -------
    pf: 2D float ndarray
        TIPS partition function for input molecule.
    isotopes: 1D string list
        List of isotopes.
    temp: 1D float ndarray
        Partition-function temperature samples (K).

    Examples
    --------
    >>> import pyratbay.opacity.partitions as pf
    >>> pf_data, isotopes, temp = pf.tips('H2O', outfile='default')

    Written partition-function file:
      'PF_tips_H2O.dat'
    for molecule H2O, with isotopes ['161', '181', '171', '162', '182', '172', '262', '282', '272'],
    and temperature range 1--6000 K.
    """
    with open(pc.ROOT+'pyratbay/data/tips_2021.pkl', 'rb') as p:
        data = pickle.load(p)
    if molecule not in data:
        raise ValueError(f"Molecule '{molecule}' is not in TIPS database.")

    if isotopes is None:
        isotopes = list(data[molecule])
    isotopes = [isotopes] if isinstance(isotopes, str) else isotopes

    for iso in isotopes:
        if iso not in data[molecule]:
            raise ValueError(
                f"Molecule '{molecule}' does not have isotope '{iso}'"
            )

    # Extend to maximum tmax
    tips_temp = data['temp']
    ntemps = [
        len(data[molecule][iso])
        for iso in isotopes
    ]
    ntemp_max = np.amax(ntemps)
    temp = tips_temp[0:ntemp_max]

    # Compute partition function:
    niso = len(isotopes)
    pf = np.zeros((niso, ntemp_max), np.double)
    for i,iso in enumerate(isotopes):
        iso = isotopes[i]
        ntemp = ntemps[i]
        part = data[molecule][iso]
        iso_temp = tips_temp[0:ntemp]
        # Extrapolate with cubic spline in log-pf:
        if ntemp == np.amax(ntemps):
            pf[i] = part
        else:
            pf[i,0:ntemp] = part
            thin = 10
            spline = CubicSpline(
                iso_temp[::thin],
                np.log(part[::thin]),
                bc_type='not-a-knot',
            )
            pf[i,ntemp:] = np.exp(spline(tips_temp[ntemp:ntemp_max]))


    # Get exomol isotope names if requested:
    if db_type == 'as_exomol':
        ID, molecs, hitran, exomol, iso_ratio, iso_mass = \
            io.read_isotopes(pc.ROOT+'pyratbay/data/isotopes.dat')
        iso_map = {
            exo.item(): hit.item()
            for mol, hit, exo in zip(molecs, exomol, hitran)
            if mol==molecule
        }
        isotopes = [iso_map[iso] for iso in isotopes]

    # Write output file:
    if outfile == 'default':
        outfile = f'PF_tips_{molecule}.dat'

    if outfile is not None:
        header = f'# Tabulated {molecule} partition-function from TIPS.\n\n'
        io.write_pf(outfile, pf, isotopes, temp, header)

        print(
            f"\nWritten partition-function file:\n  '{outfile}'"
            f"\nfor molecule {molecule}, with isotopes {isotopes},"
            f"\nand temperature range {temp[0]:.0f}--{temp[-1]:.0f} K."
        )
    return pf, isotopes, temp


def check_exomol_files(files):
    """
    Check that all input exomol files are of the same type.
    Check that all refer to a same molecule.
    Collect molecule and isotopes names.

    Parameters
    ----------
    files: List of strings
        A list of Exomol files.

    Returns
    -------
    file_type: String
        Whether all input files are .pf files (return 'pf'),
        all input files are .states or .states.bz2 files (return 'states'),
        or else return ''.
    molecule: String
        Molecule's name.
    isotopes: List of strings
        List of isotope names.
    """
    are_pf = np.all([
        file.strip().endswith('.pf')
        for file in files
    ])
    are_states = np.all([
        file.strip().endswith('.states') or file.strip().endswith('.states.bz2')
        for file in files
    ])

    if are_pf:
        file_type = 'pf'
    elif are_states:
        file_type = 'states'
    else:
        file_type = ''

    # Get molecule and isotopes, check all files of same molecule
    isotopes = []
    molecule = ''
    for file in files:
        mol, iso = pt.get_exomol_mol(file)
        if molecule != '' and molecule != mol:
            raise ValueError('All files must correspond to the same molecule')
        molecule = mol
        isotopes.append(iso)

    return file_type, molecule, isotopes



def exomol_pf(files, outfile=None):
    """
    Extract ExoMol partition-function values from input files.
    If requested, write the partition-function into a file for use
    with Pyrat Bay.

    Parameters
    ----------
    files: String or List of strings
        Input Exomol ilenames.  Files must either all correspond to .pf
        files or all correspond to .states files.
        For multiple isotopes, all files must correspond to a same molecule.
    outfile: String
         If not None, save output to file.  If outfile == 'default',
         save output to file named as PF_exomol_molecule.dat

    Returns
    -------
    pf: 2D float ndarray
        TIPS partition function for input molecule.
    isotopes: 1D string list
        List of isotopes.
    temps: 1D float ndarray
        Partition-function temperature samples (K).

    Examples
    --------
    >>> import pyratbay.opacity.partitions as pf
    >>>
    >>> # Extract data from Exomol .pf files
    >>> # wget https://www.exomol.com/db/HCN/1H-12C-14N/Harris/1H-12C-14N__Harris.pf
    >>> # wget https://www.exomol.com/db/HCN/1H-13C-14N/Larner/1H-13C-14N__Larner.pf
    >>> files = ['1H-12C-14N__Harris.pf', '1H-13C-14N__Larner.pf']
    >>> pf_data, isotopes, temps = pf.exomol_pf(files)
    """
    # Put into list if necessary:
    if isinstance(files, str):
        files = [files]

    # Make sure input files are consistent of the same type
    file_type, molecule, isotopes = check_exomol_files(files)
    if file_type != 'pf':
        error = "All input files must be exomol '.pf' files"
        raise ValueError(error)

    # Read PF data
    data, temps = [], []
    for file in files:
        temp, z = np.loadtxt(file).T
        data.append(z)
        temps.append(temp)

    niso = len(files)
    minlen = min(len(temp) for temp in temps)
    maxlen = max(len(temp) for temp in temps)
    ntemp = minlen
    if minlen != maxlen:
        for temp in temps:
            if np.any(temp[0:minlen] - temps[0][0:minlen] != 0):
                error = 'Temperature sampling in PF files are not compatible'
                raise ValueError(error)
        warning = 'Length of PF files do not match. Trimming to shorter size'
        with mu.Log() as log:
            log.warning(warning)
    pf = np.zeros((niso, ntemp), np.double)
    for i,z in enumerate(data):
        pf[i] = z[0:minlen]
    temps = temps[0][0:minlen]

    # Write output file:
    if outfile == 'default':
        outfile = f'PF_exomol_{molecule}.dat'

    if outfile is not None:
        header = (
            f'# This file incorporates the tabulated {molecule} '
            'partition-function data\n# from Exomol\n\n'
        )
        io.write_pf(outfile, pf, isotopes, temps, header)
        print(
            f"\nWritten partition-function file:\n  '{outfile}'\n"
            f"for molecule {molecule}, with isotopes {isotopes},\n"
            f"and temperature range {temps[0]:.0f}--{temps[-1]:.0f} K."
        )
    return pf, isotopes, temps


def exomol_states(files, tmin, tmax, tstep, outfile=None):
    """
    Extract ExoMol partition-function values from input files.
    If requested, write the partition-function into a file for use
    with Pyrat Bay.

    Parameters
    ----------
    files: String or List of strings
        Input Exomol ilenames.  Files must either all correspond to .pf
        files or all correspond to .states files.
        For multiple isotopes, all files must correspond to a same molecule.
    tmin: Float
        Mimimum temperature to sample the partitions.
        Required to sample from .state files only.
    tmax: Float
        Maximum temperature to sample the partitions.
        Required to sample from .state files only.
    tstep: Float
        Temperature step at which to sample the temperature array
        Required to sample from .state files only.
    outfile: String
         If not None, save output to file.  If outfile == 'default',
         save output to file named as PF_exomol_molecule.dat

    Returns
    -------
    pf: 2D float ndarray
        TIPS partition function for input molecule.
    isotopes: 1D string list
        List of isotopes.
    temps: 1D float ndarray
        Partition-function temperature samples (K).

    Examples
    --------
    >>> import pyratbay.opacity.partitions as pf
    >>>
    >>> # Extract data from Exomol .states files
    >>> # wget https://www.exomol.com/db/HCN/1H-12C-14N/Harris/1H-12C-14N__Harris.states.bz2
    >>> # wget https://www.exomol.com/db/HCN/1H-13C-14N/Larner/1H-13C-14N__Larner.states.bz2
    >>> files = [
    >>>     '1H-12C-14N__Harris.states.bz2',
    >>>     '1H-13C-14N__Larner.states.bz2',
    >>> ]
    >>> pf, isotopes, temps = pf.exomol_states(
    >>>     files, tmin=5.0, tmax=5000.0, tstep=5.0,
    >>> )
    """
    # Put into list if necessary:
    if isinstance(files, str):
        files = [files]

    # Make sure input files are consistent of the same type
    file_type, molecule, isotopes = check_exomol_files(files)
    if file_type != 'states':
        error = "All input files must be exomol '.states' files"
        raise ValueError(error)

    # Read PF data
    C2 = pc.h * pc.c / pc.k
    ntemps = int((tmax-tmin)/tstep) + 1
    temps = np.linspace(tmin, tmin + (ntemps-1)*tstep, ntemps)

    nfiles = len(files)
    pf = np.zeros((nfiles, ntemps))
    for i,state in enumerate(files):
        if state.endswith('.bz2'):
            with bz2.open(state, 'rt') as file:
                lines = file.readlines()
        else:
            with open(state, 'r') as file:
                lines = file.readlines()
        nlines = len(lines)

        energy = np.zeros(nlines)
        degeneracy = np.zeros(nlines, int)
        for j in range(nlines):
            energy[j], degeneracy[j] = lines[j][12:32].split()

        for j in range(ntemps):
            pf[i,j] = np.sum(degeneracy*np.exp(-C2*energy/temps[j]))

    # Write output file:
    if outfile == 'default':
        outfile = f'PF_exomol_{molecule}.dat'

    if outfile is not None:
        header = (
            f'# This file incorporates the tabulated {molecule} '
            'partition-function data\n# from Exomol\n\n'
        )
        io.write_pf(outfile, pf, isotopes, temps, header)
        print(
            f"\nWritten partition-function file:\n  '{outfile}'\n"
            f"for molecule {molecule}, with isotopes {isotopes},\n"
            f"and temperature range {temps[0]:.0f}--{temps[-1]:.0f} K."
        )
    return pf, isotopes, temps


def kurucz(pf_file, outfile=None, type_flag='as_exomol'):
    """
    Extract Kurucz partition-function values from input file.
    If requested, write the partition-function into a file for use
    with Pyrat Bay.

    Parameters
    ----------
    pf_file: String
        Input partition-function from Kurucz webpage.  Currently only H2O
        and TiO are available (probably there's no need for any other support).
        Files can be downloaded from these links:
          http://kurucz.harvard.edu/molecules/h2o/h2opartfn.dat
          http://kurucz.harvard.edu/molecules/tio/tiopart.dat
    outfile: String
        If not None, save output to file.
        If outfile == 'default', save output to file named as
        PF_kurucz_molecule.dat

    Returns
    -------
    pf: 2D float ndarray
        TIPS partition function for input molecule.
    isotopes: 1D string list
        List of isotopes.
    temp: 1D float ndarray
        Partition-function temperature samples (K).

    Examples
    --------
    >>> # First, download kurucz data to current dictory, e.g.:
    >>> # wget http://kurucz.harvard.edu/molecules/h2o/h2opartfn.dat
    >>> # wget http://kurucz.harvard.edu/molecules/tio/tiopart.dat

    >>> import pyratbay.opacity.partitions as pf
    >>> pf_data, isotopes, temp = pf.kurucz('h2opartfn.dat', outfile='default')

    Written partition-function file:
      'PF_kurucz_H2O.dat'
    for molecule H2O, with isotopes ['1H1H16O', '1H1H17O', '1H1H18O', '1H2H16O'],
    and temperature range 10--6000 K.

    >>> pf_data, isotopes, temp = pf.kurucz('tiopart.dat', outfile='default')

    Written partition-function file:
      'PF_kurucz_TiO.dat'
    for molecule TiO, with isotopes ['66', '76', '86', '96', '06'],
    and temperature range 10--6000 K.
    """
    if 'h2o' in pf_file:
        molecule = 'H2O'
        url = 'http://kurucz.harvard.edu/molecules/h2o/h2opartfn.dat'
        isotopes = ['1H1H16O', '1H1H17O', '1H1H18O', '1H2H16O']
        skiprows = 6
    elif 'tio' in pf_file:
        molecule = 'TiO'
        url = 'http://kurucz.harvard.edu/molecules/tio/tiopart.dat'
        isotopes = ['66', '76', '86', '96', '06']
        skiprows = 1
    else:
        print('Invalid Kurucz partition-function file.')

    # Read and extract data from files:
    data = np.loadtxt(pf_file, skiprows=skiprows, unpack=True)

    # Allocate arrays:
    temp = data[0]
    pf   = data[1:]

    # Write output file:
    if outfile == 'default':
        outfile = 'PF_kurucz_{:s}.dat'.format(molecule)

    if outfile is not None:
        header = ('# This file incorporates the tabulated {:s} '
                  'partition-function data\n# from {:s}\n\n'.
                  format(molecule, url))
        io.write_pf(outfile, pf, isotopes, temp, header)
        print("\nWritten partition-function file:\n  '{:s}'\nfor molecule "
              "{:s}, with isotopes {},\nand temperature range {:.0f}--{:.0f} "
              "K.".format(outfile, molecule, isotopes, temp[0], temp[-1]))
    return pf, isotopes, temp
