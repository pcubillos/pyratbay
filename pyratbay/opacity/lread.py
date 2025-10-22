# Copyright (c) 2021-2025 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'make_tli',
]

import os
import time
import struct
import sys

import numpy as np
import mc3.utils as mu

from . import linelist
from .. import constants as pc
from .. import tools as pt
from .. import version as ver


def pack_str(file, string):
    """Convenience function to pack length of and then a string"""
    size = len(string)
    file.write(struct.pack(f'h{size}s', size, string.encode('utf-8')))


def pack_array(file, array, format, size=None):
    """Convenience function to pack an array of given data type"""
    if size is None:
        size = len(array)
    file.write(struct.pack(f'{size}{format}', *list(array)))


def make_tli(
        dblist, pflist, dbtype, tlifile, wl_low, wl_high, wl_units,
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
        dbname.lower(): getattr(linelist, dbname)
        for dbname in pc.dbases
    }
    dblist = [
        os.path.realpath(dbase.replace('{ROOT}', pc.ROOT))
        for dbase in dblist
    ]

    databases = []
    db_names = []
    unique_dbs = []
    log.head('\nReading input database files:')
    for (dbase, pf, dtype) in zip(dblist, pflist, dbtype):
        if dtype not in db_readers:
            log.error(
                f"Unknown type '{dtype}' for database '{dbase}'.  "
                f"Select from: {str(pc.dbases)}"
            )
        log.head(dbase, indent=2)
        db = db_readers[dtype](dbase, pf, log)
        databases.append(db)
        db_names.append(db.name)
        if db.name not in unique_dbs:
            unique_dbs.append(db.name)
    log.msg(f'There are {nfiles} input database file(s).\n\n')

    # Boundaries in wavenumber space (in cm-1):
    wn_low = 1.0 / wl_high / pt.u(wl_units)
    wn_high = 1.0 / wl_low / pt.u(wl_units)

    # Output file:
    tli = {}
    tli['version'] = f'{ver.LR_VER}.{ver.LR_MIN}.{ver.LR_REV}'
    tli['wn_units'] = 'cm-1'
    tli['wn_min'] = wn_low
    tli['wn_max'] = wn_high

    n_databases = len(unique_dbs)
    tli['n_databases'] = n_databases

    log.msg(
        f'Initial wavelength: {wl_low:7.3f} {wl_units} ({wn_high:9.3f} cm-1)\n'
        f'Final wavelength:   {wl_high:7.3f} {wl_units} ({wn_low:9.3f} cm-1)\n'
        f'There are {n_databases} different database(s).'
    )

    log.head('\nExtracting line transition info.')
    tli['databases'] = []
    for i,db_name in enumerate(unique_dbs):
        dbase = {}
        tli['databases'].append(dbase)
        ti = time.time()
        wn = []
        gf = []
        elow = []
        iso_id = []
        for db in databases:
            if db.name != db_name:
                continue
            this_db = db
            transitions = db.dbread(wn_low, wn_high, log.verb)
            if transitions is None:
                continue

            wn.append(transitions[0])
            gf.append(transitions[1])
            elow.append(transitions[2])
            iso_id.append(transitions[3])

        db = this_db
        wn = np.concatenate(wn)
        gf = np.concatenate(gf)
        elow = np.concatenate(elow)
        iso_id = np.concatenate(iso_id)
        tf = time.time()
        log.debug(f'Reading time: {tf-ti:8.3f} seconds', indent=2)

        ntransitions = np.size(wn)
        # iso_id are indices as in db.isotopes
        # iso_idx are indices 0-N, after filtering out isotopes with no lines
        unique_iso, iso_idx, ntrans_iso = np.unique(
            iso_id, return_inverse=True, return_counts=True,
        )

        # Sort by isotope ID, then each isotope by wavenumber:
        ti = time.time()
        isort = np.argsort(iso_id)
        ihi = 0
        for ntrans in ntrans_iso:
            ilo = ihi
            ihi += ntrans
            wn_sort = np.argsort(wn[isort][ilo:ihi])
            isort[ilo:ihi] = isort[ilo:ihi][wn_sort]
        tf = time.time()

        wn = wn[isort]
        gf = gf[isort]
        elow = elow[isort]
        iso_id = iso_id[isort]
        iso_idx = iso_idx[isort]
        log.debug(f'Sort time:    {tf-ti:8.3f} seconds', indent=2)

        dbase['name'] = db.name
        dbase['molecule'] = db.molecule
        dbase['n_lines'] = ntransitions
        dbase['n_lines_iso'] = ntrans_iso
        dbase['iso_id'] = iso_idx
        dbase['wn'] = wn
        dbase['elow'] = elow
        dbase['gf'] = gf

        # Filter out isotopes with no line transitions:
        iso_names = np.array(db.isotopes)[unique_iso]
        iso_mass = np.array(db.mass)[unique_iso]
        iso_ratio = np.array(db.isoratio)[unique_iso]

        temp, partition, pf_iso = db.getpf(log.verb)
        iso_match = np.isin(iso_names, pf_iso)
        if np.any(~iso_match):
            log.error(
                'No partition functions found for these isotopes of the '
                f'{db.molecule} line list: {iso_names[~iso_match]}'
            )

        # Filter and sort PF by isotopes in iso_names:
        pf_idx = [pf_iso.index(iso) for iso in iso_names]
        pf = partition[pf_idx]

        # Store the number of temperature samples and isotopes:
        dbase['temperatures'] = temp
        dbase['isotopes'] = iso_names
        dbase['iso_mass'] = iso_mass
        dbase['iso_ratio'] = iso_ratio
        dbase['partition'] = pf

        # Report info for each isotope
        wl_min = 1.0 / np.amax(wn) / pc.um
        wl_max = 1.0 / np.amin(wn) / pc.um
        n_iso = len(iso_names)
        log.msg(
            f"Database ({i+1}/{n_databases}): {repr(db.name)} "
            f"({db.molecule} molecule)\n"
            f'Number of isotopes with line transitions: {n_iso}',
            indent=2,
        )
        log.msg("idx  isotope    mass (u)    fraction       n_lines", indent=2)
        for j in range(n_iso):
            name = f'{repr(str(iso_names[j])):10s}'
            mass = iso_mass[j]
            ratio = iso_ratio[j]
            ratio = f'{ratio:9.7f}' if ratio >= 1e-5 else f'{ratio:.3e}'
            ntrans = ntrans_iso[j]
            log.msg(
                f"{j+1:3d}  {name}  {mass:7.3f}    {ratio}  {ntrans:11,d}",
                indent=2,
            )

        log.msg(
            f'Total: {ntransitions:,d} line transitions '
            f'between {wl_min:.3f} -- {wl_max:.3f} um\n\n'
            f'Number of temperatures: {len(temp)}\n'
            '  Temperatures (K): '
            f'[{temp[0]:6.1f}, {temp[1]:6.1f}, ..., {temp[-1]:6.1f}]',
            indent=2,
        )
        for j in range(n_iso):
            log.msg(
                f'Partition Function ({str(iso_names[j])}):  '
                f'[{pf[j,0]:.2e}, {pf[j,1]:.2e}, ..., {pf[j,-1]:.2e}]',
                indent=4,
            )

    # Store to file
    ti = time.time()
    tli_file = open(tlifile, 'wb')
    endian = sys.byteorder[0]
    tli_file.write(struct.pack('s', endian.encode('utf-8')))
    tli_file.write(struct.pack('3h', ver.LR_VER, ver.LR_MIN, ver.LR_REV))
    tli_file.write(struct.pack('2d', tli['wn_min'], tli['wn_max']))
    tli_file.write(struct.pack('h', tli['n_databases']))

    dbases = tli['databases']
    for dbase in dbases:
        pack_str(tli_file, dbase['name'])
        pack_str(tli_file, dbase['molecule'])
        ntemp = len(dbase['temperatures'])
        niso = len(dbase['isotopes'])
        tli_file.write(struct.pack('hh', ntemp, niso))
        pack_array(tli_file, dbase['temperatures'], 'd')
        for j,iso in enumerate(dbase['isotopes']):
            pack_str(tli_file, str(iso))
            tli_file.write(struct.pack('d', dbase['iso_mass'][j]))
            tli_file.write(struct.pack('d', dbase['iso_ratio'][j]))
            pack_array(tli_file, dbase['partition'][j], 'd')

    n_lines = np.sum([dbase['n_lines'] for dbase in dbases])
    tli_file.write(struct.pack('i', n_lines))

    n_lines_iso = np.concatenate([dbase['n_lines_iso'] for dbase in dbases])
    tli_file.write(struct.pack('i', len(n_lines_iso)))
    for dbase in dbases:
        pack_array(tli_file, dbase['n_lines_iso'], 'i')

    for dbase in dbases:
        pack_array(tli_file, dbase['wn'], 'd')
    for dbase in dbases:
        pack_array(tli_file, dbase['iso_id'], 'h')
    for dbase in dbases:
        pack_array(tli_file, dbase['elow'], 'd')
    for dbase in dbases:
        pack_array(tli_file, dbase['gf'], 'd')
    tli_file.close()
    tf = time.time()
    log.debug(f'Writing time: {tf-ti:8.3f} seconds')
    log.head(f"Generated TLI file: '{tlifile}'.")
    log.close()
