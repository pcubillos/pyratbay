# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'make_tli',
]

import os
import sys
import time
import struct

import numpy as np
import mc3.utils as mu

from . import linelist
from .. import constants as pc
from .. import tools as pt
from .. import version as ver


def make_tli(
        dblist, pflist, dbtype, tlifile, wl_low,  wl_high, wl_units,
        log=None,
    ):
    """
    Create a transition-line-information (TLI) file.

    Parameters
    ----------
    dblist: List of strings
        Opacity databases to read.
    pflist: List of strings
        Partition function for each of the databases.
    dbtype: List of strings
        Database type of each database.
    tlifile: String
        Output TLI file name.
    wl_low: Float
        Lower wavelength boundary to consider, units given by wl_units.
    wl_high: Float
        High wavelength boundary to consider, units given by wl_units.
    wl_units: String
        Wavelength units (when not specified in wl_low nor wl_high).
    log: Log object
        An mc3.utils.Log instance to log screen outputs to file.
    """
    if log is None:
        log = mu.Log(verb=2)

    # Input-not-found error messages:
    if tlifile is None:
        log.error('Undefined TLI file (tlifile).')

    if wl_low is None:
        log.error('Undefined low wavelength boundary (wl_low)')
    if wl_high is None:
        log.error('Undefined high wavelength boundary (wl_high)')

    if dblist is None:
        log.error('There are no input database files (dblist)')
    if dbtype is None:
        log.error('There are no input database types (dbtype)')
    if pflist is None:
        log.error('There are no partition-function inputs (pflist)')

    # Check number of files match:
    if isinstance(dblist, str):
        dblist = [dblist]
    nfiles = len(dblist)

    if isinstance(pflist, str):
        pflist = [pflist]
    if len(pflist) == 1:
        pflist = [pflist[0] for _ in range(nfiles)]
    if isinstance(dbtype, str):
        dbtype = [dbtype]
    if len(dbtype) == 1:
        dbtype = [dbtype[0] for _ in range(nfiles)]

    if nfiles != len(pflist) or nfiles != len(dbtype):
        log.error(
            f'The number of Line-transition files ({nfiles}) does not match '
            f'the number of partition-function files ({len(pflist)}) or '
            f'database-type files ({len(dbtype)})'
        )

    # Driver routine to read the databases:
    db_readers = {
        dbname.lower(): getattr(linelist,dbname)
        for dbname in pc.dbases
    }
    dblist = [
        os.path.realpath(dbase.replace('{ROOT}', pc.ROOT))
        for dbase in dblist
    ]

    databases = []
    db_names = []
    log.head('\nReading input database files:')
    for (dbase, pf, dtype) in zip(dblist, pflist, dbtype):
        if dtype not in db_readers:
            log.error(
                f"Unknown type '{dtype}' for database '{dbase}'.  "
                f"Select from: {str(pc.dbases)}"
            )
        log.head(f'- {dbase}')
        databases.append(db_readers[dtype](dbase, pf, log))
        db_names.append(databases[-1].name)
    log.msg(f'There are {nfiles} input database file(s).')

    # Open output file:
    tli = open(tlifile, 'wb')

    # Get the machine endian type (big/little):
    if sys.byteorder == 'big':
        endian = 'b'
    if sys.byteorder == 'little':
        endian = 'l'

    # Start storing TLI header values:
    header = struct.pack('s', endian.encode())
    header += struct.pack('3h', ver.LR_VER, ver.LR_MIN, ver.LR_REV)

    # Boundaries in wavenumber space (in cm-1):
    wn_low = 1.0 / wl_high / pt.u(wl_units)
    wn_high = 1.0 / wl_low / pt.u(wl_units)

    # Add initial and final wavenumber boundaries (in cm-1):
    header += struct.pack('2d', wn_low, wn_high)

    Ndb = len(np.unique(db_names))
    header += struct.pack('h', Ndb)
    tli.write(header)

    log.msg(
        f'\nOS endianness:  {sys.byteorder}\n'
        f'Initial TLI wavelength ({wl_units}): {wl_low:7.3f} ({wn_high:9.3f} cm-1)\n'
        f'Final TLI wavelength ({wl_units}):   {wl_high:7.3f} ({wn_low:9.3f} cm-1)\n'
        f'There are {Ndb} different database(s).'
    )


    log.msg('\nReading and writting partition function info.')
    idb = 1         # Database correlative number
    niso_total = 0  # Cumulative number of isotopes
    accum = [0]     # Cumulative number of isotopes per database
    db_names = []
    # Loop through the partition files (if more than one) and write the
    # data to a processed TLI file:
    for db in databases:
        # Skip if we already stored the pf info of this DB:
        if db.name in db_names:
            continue
        db_names.append(db.name)

        # Get partition function values:
        temp, partition, pf_iso = db.getpf(log.verb)
        iso_names = db.isotopes
        iso_mass = db.mass
        iso_ratio = db.isoratio

        # Number of temperature samples and isotopes:
        ntemp = len(temp)
        niso = len(iso_names)

        # Extract partition-function info sorted by iso_names:
        pf = np.zeros((niso, ntemp), np.double)
        for part,iso in zip(partition, pf_iso):
            # Ignore PF isotopes that don't exist in isotopes.dat:
            if iso not in iso_names:
                continue
            idx = iso_names.index(iso)
            pf[idx] = part

        # Store length of and database name:
        name = db.name
        pack = struct.pack(f'h{len(name)}s', len(name), name.encode('utf-8'))
        tli.write(pack)
        # Store the molecule name:
        mol = db.molecule
        pack = struct.pack(f'h{len(mol)}s', len(mol), mol.encode('utf-8'))
        tli.write(pack)
        # Store the number of temperature samples and isotopes:
        tli.write(struct.pack('hh', ntemp, niso))
        log.msg(
            f"Database ({idb}/{Ndb}): '{db.name}' ({db.molecule} molecule)",
            indent=2,
        )
        log.msg(
            f'Number of temperatures: {ntemp}\n'
            f'Number of isotopes: {niso}',
            indent=4,
        )

        # Write the temperature array:
        tli.write(struct.pack(f'{ntemp}d', *temp))
        log.msg(
            'Temperatures (K): '
            f'[{temp[0]:6.1f}, {temp[1]:6.1f}, ..., {temp[-1]:6.1f}]',
            indent=4,
        )

        # For each isotope, write partition function information.
        for j in range(niso):
            iname = iso_names[j]
            log.msg(f"Isotope ({j+1}/{niso}): '{iname}'", indent=4)

            # Store length of isotope name, isotope name, and isotope mass:
            pack = struct.pack(
                f'h{len(iname)}s', len(iname), str(iname).encode('utf-8'),
            )
            tli.write(pack)
            tli.write(struct.pack('d', iso_mass[j]))
            tli.write(struct.pack('d', iso_ratio[j]))

            # Write the partition function per isotope:
            tli.write(struct.pack(f'{ntemp}d', *pf[j]))
            log.msg(
                f'Mass (u):        {iso_mass[j]:8.4f}\n'
                f'Isotopic ratio:  {iso_ratio[j]:8.4g}\n'
                f'Part. Function:  '
                f'[{pf[j,0]:.2e}, {pf[j,1]:.2e}, ..., {pf[j,-1]:.2e}]',
                indent=6,
            )

        # Calculate cumulative number of isotopes per database:
        niso_total += niso
        idb += 1
        accum.append(niso_total)

    log.msg(f'Cumulative number of isotopes per database: {accum}')

    log.head('\nExtracting line transition info.')
    wnumber = np.array([], np.double)
    gf = np.array([], np.double)
    elow = np.array([], np.double)
    isoID = np.array([], int)
    # Read from file and write the transition info:
    for db in databases:
        # Get database index:
        idb = db_names.index(db.name)

        ti = time.time()
        transitions = db.dbread(wn_low, wn_high, log.verb)
        tf = time.time()

        if transitions is None:
            continue

        wnumber = np.concatenate((wnumber, transitions[0]))
        gf = np.concatenate((gf, transitions[1]))
        elow = np.concatenate((elow, transitions[2]))
        isoID = np.concatenate((isoID, transitions[3]+accum[idb]))

        unique_iso = np.unique(transitions[3])
        log.msg(
            f'Isotope in-database indices: {unique_iso}\n'
            f'Isotope correlative indices: {unique_iso+accum[idb]}',
            indent=2,
        )
        log.debug('Reading time: {tf-ti:8.3f} seconds', indent=2)


    # Total number of transitions:
    ntransitions = np.size(wnumber)
    # Number of transitions per isotope:
    ntrans_iso = np.bincount(isoID)
    ntrans_iso = ntrans_iso[np.where(ntrans_iso>0)]  # Remove zeroes

    # Sort by isotope ID:
    ti = time.time()
    isort = np.argsort(isoID)
    # Sort each isotope by wavenumber:
    ihi = 0
    for ntrans in ntrans_iso:
        ilo = ihi
        ihi += ntrans
        wn_sort = np.argsort(wnumber[isort][ilo:ihi])
        isort[ilo:ihi] = isort[ilo:ihi][wn_sort]
    tf = time.time()

    # Actual sorting:
    wnumber = wnumber[isort]
    gf = gf[isort]
    elow = elow[isort]
    isoID = isoID[isort]

    log.debug(f'Sort time:    {tf-ti:8.3f} seconds', indent=2)
    ntrans_str = '  '.join([f'{val:,d}' for val in ntrans_iso])
    log.msg(f'\nTransitions per isotope:\n[{ntrans_str}]')

    # Pack:
    tli.write(struct.pack('i', ntransitions))
    wn_min = np.amin(wnumber)
    wn_max = np.amax(wnumber)
    wl_min = 1.0 / wn_max / pc.um
    wl_max = 1.0 / wn_min / pc.um
    log.msg(
        f'\nWriting {ntransitions:,d} transition lines '
        f'between wavenumbers {wn_min:.2f} and {wn_max:.2f} cm-1 '
        f'({wl_min:.3f} -- {wl_max:.3f} um).'
    )

    # Write the number of transitions for each isotope:
    niso = len(ntrans_iso)
    tli.write(struct.pack('i', niso))
    tli.write(struct.pack(str(niso)+'i', *list(ntrans_iso)))

    # Write the Line-transition data:
    ti = time.time()
    transinfo  = struct.pack(str(ntransitions)+'d', *list(wnumber))
    transinfo += struct.pack(str(ntransitions)+'h', *list(isoID))
    transinfo += struct.pack(str(ntransitions)+'d', *list(elow))
    transinfo += struct.pack(str(ntransitions)+'d', *list(gf))
    tf = time.time()
    log.debug(f'Packing time: {tf-ti:8.3f} seconds')

    ti = time.time()
    tli.write(transinfo)
    tf = time.time()
    log.debug(f'Writing time: {tf-ti:8.3f} seconds')

    log.head(f"Generated TLI file: '{tlifile}'.")
    tli.close()
    log.close()
