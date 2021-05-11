# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'get_tips_molname',
    'tips',
    'exomol',
    'kurucz',
    ]

import pickle

import numpy as np
from numpy.core.numeric import isscalar

from ... import io as io
from ... import constants as pc
from ... import tools as pt
import mc3.utils as mu


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
    with open(pc.ROOT+'pyratbay/data/tips_2017.pkl', 'rb') as p:
        data = pickle.load(p)
    if molID not in data['mol_ID']:
        raise ValueError('TIPS 2017 database does not contain molecule ID: {}'
            .format(molID))
    return data['mol_ID'][molID]


def tips(molecule, isotopes=None, outfile=None, db_type='as_tips'):
    """
    Extract TIPS 2017 partition-function values for given molecule.
    If requested, write the partition-function into a file for use
    with Pyrat Bay.
    Reference: Gamache et al. (2017), JQSRT, 203, 70.

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
    and temperature range 1--5000 K.
    """
    with open(pc.ROOT+'pyratbay/data/tips_2017.pkl', 'rb') as p:
        data = pickle.load(p)
    if molecule not in data:
        raise ValueError(f"Molecule '{molecule}' is not in TIPS database.")

    if isotopes is None:
        isotopes = list(data[molecule].keys())
    isotopes = [isotopes] if isscalar(isotopes) else isotopes

    for iso in isotopes:
        if iso not in data[molecule]:
            raise ValueError(
                f"Molecule '{molecule}' does not have isotope '{iso}'")

    ntemp = np.amin([data['ntemp'][molecule][iso] for iso in data[molecule]])
    temp = data['temp'][0:ntemp]

    # Compute partition function:
    niso  = len(isotopes)
    pf = np.zeros((niso, ntemp), np.double)
    for i,iso in enumerate(isotopes):
        pf[i] = data[molecule][iso][0:ntemp]

    # Get exomol isotope names if requested:
    if db_type == 'as_exomol':
        ID, molecs, hitran, exomol, iso_ratio, iso_mass = \
            io.read_isotopes(pc.ROOT+'pyratbay/data/isotopes.dat')
        iso_map = {exo:hit for mol, hit, exo in zip(molecs, exomol, hitran)
                   if mol==molecule}
        isotopes = [iso_map[iso] for iso in isotopes]

    # Write output file:
    if outfile == 'default':
        outfile = f'PF_tips_{molecule}.dat'

    if outfile is not None:
        header = f'# Tabulated {molecule} partition-function from TIPS.\n\n'
        io.write_pf(outfile, pf, isotopes, temp, header)

        print(f"\nWritten partition-function file:\n  '{outfile}'"
              f"\nfor molecule {molecule}, with isotopes {isotopes},"
              f"\nand temperature range {temp[0]:.0f}--{temp[-1]:.0f} K.")
    return pf, isotopes, temp


def exomol(pf_files, outfile=None):
    """
    Extract ExoMol partition-function values from input files.
    If requested, write the partition-function into a file for use
    with Pyrat Bay.

    Parameters
    ----------
    pf_files: String or List of strings
        Input Exomol partition-function filenames.  If there are
        multiple isotopes, all of them must correspond to the same
        molecule.
    outfile: String
        If not None, save output to file.
        If outfile == 'default', save output to file named as
        PF_exomol_molecule.dat

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
    >>> # First, download ExoMol data to current dictory, e.g.:
    >>> # wget http://www.exomol.com/db/NH3/14N-1H3/BYTe/14N-1H3__BYTe.pf
    >>> # wget http://www.exomol.com/db/NH3/15N-1H3/BYTe-15/15N-1H3__BYTe-15.pf
    >>> import pyratbay.opacity.partitions as pf
    >>> # A single file:
    >>> pf_data, isotopes, temp = pf.exomol('14N-1H3__BYTe.pf',
    >>>     outfile='default')
    Written partition-function file:
      'PF_exomol_NH3.dat'
    for molecule NH3, with isotopes ['4111'],
    and temperature range 1--1600 K.

    >>> # Multiple files (isotopes) for a molecule:
    >>> pf_data, isotopes, temp = pf.exomol(
    >>>     ['14N-1H3__BYTe.pf', '15N-1H3__BYTe-15.pf'], outfile='default')

    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Warning:
        Length of PF files do not match.  Trimming to shorter size.
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    Written partition-function file:
      'PF_exomol_NH3.dat'
    for molecule NH3, with isotopes ['4111', '5111'],
    and temperature range 1--1600 K.
    """
    # Put into list if necessary:
    if isinstance(pf_files, str):
        pf_files = [pf_files]

    # Read and extract data from files:
    isotopes = []
    data, temps = [], []
    molecule = ''
    for pf_file in pf_files:
        # Get info from file name:
        mol, iso = pt.get_exomol_mol(pf_file)

        # Check all files correspond to the same molecule.
        if molecule == '':
            molecule = mol
        elif molecule != mol:
            raise ValueError('All files must correspond to the same molecule.')

        isotopes.append(iso)
        # Read data:
        temp, z = np.loadtxt(pf_file).T
        data.append(z)
        temps.append(temp)

    # Number of isotopes:
    niso = len(isotopes)
    # Check temp sampling:
    minlen = min(len(temp) for temp in temps)
    maxlen = max(len(temp) for temp in temps)
    ntemp = minlen
    if minlen != maxlen:
        for temp in temps:
            if np.any(temp[0:minlen] - temps[0][0:minlen] != 0):
                raise ValueError(
                    'Temperature sampling in PF files are not compatible.')
        with mu.Log() as log:
            log.warning('Length of PF files do not match.  Trimming to '
                        'shorter size.')

    pf = np.zeros((niso, ntemp), np.double)
    for i,z in enumerate(data):
        pf[i] = z[0:minlen]
    temp = temps[i][0:minlen]

    # Write output file:
    if outfile == 'default':
        outfile = f'PF_exomol_{molecule}.dat'

    if outfile is not None:
        header = (f'# This file incorporates the tabulated {molecule} '
                  'partition-function data\n# from Exomol\n\n')
        io.write_pf(outfile, pf, isotopes, temp, header)
        print("\nWritten partition-function file:\n  '{:s}'\nfor molecule "
              "{:s}, with isotopes {},\nand temperature range {:.0f}--{:.0f} "
              "K.".format(outfile, molecule, isotopes, temp[0], temp[-1]))
    return pf, isotopes, temp


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
