# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import sys
import numpy as np

from .. import constants as pc
from .. import tools as pt
from . import objects as o


def read_tli(pyrat):
    """
    Main driver to read the line transition data from TLI files.
    """
    pyrat.lt  = pyrat.lt.clone_new(pyrat)
    pyrat.iso = o.Isotopes()
    # Count number of TLI files:
    if pyrat.lt.tlifile is None:
        pyrat.log.head("\nNo line transition file to read.")
        return

    pyrat.log.head("\nReading line transition info.")
    # TLI file object:
    tlis = []
    # Index of first database in TLI file:
    dbindex = [0]

    # Read data bases header info:
    for tlifile in pyrat.lt.tlifile:
        pyrat.log.head(f"Read TLI file: '{tlifile}'.", indent=2)
        tli = open(tlifile, "rb")
        tlis.append(tli)
        dbindex.append(dbindex[-1] + read_header(pyrat, tli))

    # Set link to molecules' indices:
    setimol(pyrat)

    # Read line-transition data (if there's no extinction-coefficient table):
    if pt.isfile(pyrat.ex.extfile) != 1:
        for tli, dbi in zip(tlis, dbindex):
            read_linetransition(pyrat, tli, dbi)
            tli.close()

    pyrat.log.msg(f"Read a total of {pyrat.lt.ntransitions:,d} line "
        "transitions.", indent=2)
    pyrat.log.head("Line-transition done.\n")


def read_header(pyrat, linefile):
    """
    Open TLI file and read database header info.

    Returns:
    --------
    Ndb: Integer
        Number of databases in this TLI file.
    """
    log = pyrat.log
    # Read header:
    endian = linefile.read(1).decode()
    if endian == 'b':
        log.msg("TLI data storage: Big endian", indent=2)
    elif endian == 'l':
        log.msg("TLI data storage: Little endian", indent=2)

    # Compare endianness:
    endianness = sys.byteorder
    if endianness[0:1] != endian:
        log.error(
            f"Incompatible endianness between TLI file ({endian}) and "
            f"Pyrat ({endianness[0:1]}).")

    # Read TLI version:
    TLI_ver, TLI_min, TLI_rev = pt.unpack(linefile, 3, "h")
    log.msg("TLI version: {:d}.{:d}.{:d}.".
                  format(TLI_ver, TLI_min, TLI_rev), indent=2)
    if TLI_ver != 6 or TLI_min not in [1,2,3,4,5]:
        log.error(
            "Incompatible TLI version.  The TLI file must be created "
            "with Lineread version 6.1-6.5.")

    # Read initial and final wavenumber from TLI:
    lt_wni, lt_wnf = pt.unpack(linefile, 2, "d")
    log.msg(f"TLI wavenumber range (cm-1): [{lt_wni:.1f}, {lt_wnf:.1f}]",
        indent=2)
    # Check TLI and pyrat wavelength ranges:
    checkrange(pyrat, lt_wni, lt_wnf)

    # Read number of data bases:
    Ndb = pt.unpack(linefile, 1, "h")
    log.msg(f"Number of data bases: {Ndb}", indent=2)

    # Cumulative isotope index:
    acumiso = pyrat.iso.niso

    for i in range(Ndb):
        db = o.Database()
        # Read Database name:
        lenDBname = pt.unpack(linefile, 1, "h")
        db.name   = pt.unpack(linefile, lenDBname, "s")
        log.msg(f"Data base name: '{db.name}'", indent=2)
        # Read Molecule name:
        lenMolec   = pt.unpack(linefile, 1,        "h")
        db.molname = pt.unpack(linefile, lenMolec, "s")
        log.msg(f"Molecule name: '{db.molname}'", indent=2)
        # Read temperature array:
        db.ntemp, db.niso =  pt.unpack(linefile, 2,        "h")
        db.temp = np.asarray(pt.unpack(linefile, db.ntemp, "d"))
        # Update temperature boundaries:
        pyrat.lt.tmin = np.amax((pyrat.lt.tmin, db.temp[ 0]))
        pyrat.lt.tmax = np.amin((pyrat.lt.tmax, db.temp[-1]))
        log.msg("Temperature range: "
               f"{db.temp[0]:4.1f} -- {db.temp[-1]:4.1f} K.", indent=2)

        # Allocate arrays for isotopic info:
        name    = np.zeros(db.niso, 'U20')
        mass    = np.zeros(db.niso)
        ratio   = np.zeros(db.niso)
        dbindex = np.zeros(db.niso, np.int)
        db.z    = np.zeros((db.niso, db.ntemp))

        # Store per-isotope info:
        for j in range(db.niso):
            dbindex[j] = i + pyrat.lt.ndb
            lenIsoName = pt.unpack(linefile, 1,          "h")
            name[j]    = pt.unpack(linefile, lenIsoName, "s")
            mass[j]    = pt.unpack(linefile, 1,          "d")
            ratio[j]   = pt.unpack(linefile, 1,          "d")
            db.z[j]    = np.asarray(pt.unpack(linefile, db.ntemp, "d"))

            # Print info to screen:
            log.msg(f'Isotope: {name[j]},  mass: {mass[j]:.4f} u,  '
                f'isotopic ratio: {ratio[j]:.5g}', indent=3)
            log.debug(f'Z = [{db.z[j,0]:.2e}, {db.z[j,1]:.2e}, ..., '
                f'{db.z[j,-1]:.2e}]', indent=4)

        # Add the number of isotopes read:
        pyrat.iso.niso += db.niso
        log.debug(f"Number of isotopes: {pyrat.iso.niso}", indent=2)

        # Store name, mass in isotopes structure:
        pyrat.iso.name    = np.concatenate((pyrat.iso.name,    name))
        pyrat.iso.mass    = np.concatenate((pyrat.iso.mass,    mass))
        pyrat.iso.ratio   = np.concatenate((pyrat.iso.ratio,   ratio))
        pyrat.iso.dbindex = np.concatenate((pyrat.iso.dbindex, dbindex))

        # Set isotope correlative index for DB:
        db.iiso  = acumiso
        acumiso += db.niso
        pyrat.lt.db.append(db)
        log.msg(f"DB index: {i+pyrat.lt.ndb}, Cumulative Iso: {acumiso}\n\n",
            indent=2)

    # Keep count of number of databases:
    pyrat.lt.ndb += Ndb

    # Return the pointer position in lineinfo file:
    return Ndb


def read_linetransition(pyrat, linefile, dbindex):
    """
    Read the databases line transition info.
    """
    # Read the number of line transitions:
    nTransitions = pt.unpack(linefile, 1, "i")
    # Read the number of isotopes in line-transition array:
    nIso = pt.unpack(linefile, 1, "i")
    # Read the number of transitions per isotope:
    NisoTran = np.atleast_1d(pt.unpack(linefile, nIso, "i"))

    # Position where the line-transition data begins:
    init_wl  = linefile.tell()
    init_iso = init_wl  + nTransitions*pc.dreclen  # Init pos of isoID data
    init_el  = init_iso + nTransitions*pc.sreclen  # Init pos of Elow data
    init_gf  = init_el  + nTransitions*pc.dreclen  # Init pos of gf data

    # Count the number of transitions:
    linefile.seek(0, 2)
    endrec = linefile.tell()
    nrec = (endrec - init_wl)*1.0 / pc.tlireclen
    if nrec != nTransitions:
        pyrat.log.error(
            f'The remaining data file size ({nrec:.1f}) does not '
            f'correspond to the number of transitions ({nTransitions}).')
    pyrat.log.msg(f'There are {nTransitions:,d} line transitions in TLI file.',
        indent=2)

    # Allocate arrays:
    wn    = np.zeros(nTransitions)
    isoid = np.zeros(nTransitions, np.short)
    elow  = np.zeros(nTransitions)
    gf    = np.zeros(nTransitions)

    # Number of line-transitions offset for a given isotope:
    offset = 0
    start  = init_wl
    nlt    = 0  # Total number of line-transitions read
    for i in range(nIso):
        # Search lower and higher line-transition indices to read:
        ifirst = pt.binsearch(
            linefile, pyrat.spec.wnlow,  start, NisoTran[i], upper=False)
        ilast  = pt.binsearch(
            linefile, pyrat.spec.wnhigh, start, NisoTran[i], upper=True)
        # Add offset for this isotope:
        ifirst += offset
        ilast  += offset

        if ifirst < 0 or ilast < 0:
            start  += NisoTran[i]*pc.dreclen
            offset += NisoTran[i]
            continue

        # Print boundaries:
        linefile.seek(ifirst*pc.dreclen + init_wl, 0)
        first = pt.unpack(linefile, 1, 'd')
        linefile.seek(ilast*pc.dreclen  + init_wl, 0)
        last = pt.unpack(linefile, 1, 'd')
        pyrat.log.debug(
            f"Found initial transition ({ifirst:8d}):  {first:13.4f} cm-1",
            indent=2)
        pyrat.log.debug(
            f"Found final   transition ({ilast:8d}):  {last:13.4f} cm-1",
            indent=2)

        # Number of transitions to read:
        nread = ilast - ifirst + 1
        # Get isotope ID:
        i0 = pt.binsearch(linefile, 0, start, NisoTran[i], upper=False)
        linefile.seek((i0+offset)*pc.sreclen + init_iso, 0)
        isoID = pt.unpack(linefile, 1, 'h')

        # Read data into arrays:
        linefile.seek(ifirst*pc.dreclen + init_wl,  0)
        wn   [nlt:nlt+nread] = pt.unpack(linefile, nread, "d")

        linefile.seek(ifirst*pc.sreclen + init_iso, 0)
        isoid[nlt:nlt+nread] = pt.unpack(linefile, nread, 'h')

        linefile.seek(ifirst*pc.dreclen + init_el,  0)
        elow [nlt:nlt+nread] = pt.unpack(linefile, nread, 'd')

        linefile.seek(ifirst*pc.dreclen + init_gf,  0)
        gf   [nlt:nlt+nread] = pt.unpack(linefile, nread, 'd')
        pyrat.log.msg(
            f"Read {nread:11,d} transitions for isotope "
            f"{pyrat.lt.db[dbindex].iiso+isoID:2d}.", indent=4)

        start  += NisoTran[i]*pc.dreclen
        offset += NisoTran[i]
        nlt    += nread

    # Add the pre-existing number of isotopes:
    isoid += pyrat.lt.db[dbindex].iiso

    pyrat.lt.wn    = np.concatenate((pyrat.lt.wn,    wn   [0:nlt]))
    pyrat.lt.gf    = np.concatenate((pyrat.lt.gf,    gf   [0:nlt]))
    pyrat.lt.elow  = np.concatenate((pyrat.lt.elow,  elow [0:nlt]))
    pyrat.lt.isoid = np.concatenate((pyrat.lt.isoid, isoid[0:nlt]))
    pyrat.lt.ntransitions += nlt


def checkrange(pyrat, wn_low, wn_high):
    """
    Display a warning if line database spectral range does not completely
    include the pyrat spectral range.

    Parameters:
    -----------
    pyrat: Object
    wn_low: Float
        Database's wavenumber lower boundary (in cm^-1).
    wn_high: Float
        Database's wavenumber higher boundary (in cm^-1).
    """
    # Print out warning if ranges dont overlap:
    if wn_low > pyrat.spec.wnhigh or wn_high < pyrat.spec.wnlow:
        pyrat.log.warning(
            f"TLI wavenumber range ({wn_low:.2f} - {wn_high:.2f} cm^-1) does "
            f"not overlap with Pyrat wavenumber range "
            f"({pyrat.spec.wnlow:.2f} - {pyrat.spec.wnhigh:.2f} cm^-1).")
    # Print out warning if TLI range is smaller than the pyrat required range:
    elif wn_low > pyrat.spec.wnlow or wn_high < pyrat.spec.wnhigh:
        pyrat.log.warning(
            f"TLI wavenumber range ({wn_low:.2f} - {wn_high:.2f} cm^-1) does "
            f"not cover the full Pyrat wavenumber range "
            f"({pyrat.spec.wnlow:.2f} - {pyrat.spec.wnhigh:.2f} cm^-1).")


def setimol(pyrat):
    """
    Set the molecule index for the list of isotopes.
    """
    # Allocate imol array:
    pyrat.iso.imol = np.zeros(pyrat.iso.niso, np.int)
    # For each isotope:
    for i in range(pyrat.iso.niso):
        # Get molecule name from database object:
        molname = pyrat.lt.db[pyrat.iso.dbindex[i]].molname
        # Set index:
        if molname not in pyrat.mol.symbol:  # Isotope's not in molecule list
            pyrat.iso.imol[i] = -1
        else:
            pyrat.iso.imol[i] = np.where(pyrat.mol.symbol == molname)[0]
    pyrat.log.msg(f"Isotope's molecule indices:\n  {pyrat.iso.imol}", indent=2)
    # Report missing species:
    imiss = np.unique(pyrat.iso.dbindex[np.where(pyrat.iso.imol < 0)])
    for i in imiss:
        iso_miss = pyrat.iso.name[np.where(pyrat.iso.dbindex==i)]
        pyrat.log.warning(
            f"The species '{pyrat.lt.db[i].molname}' for isotopes "
            f"{iso_miss} is not present in the atmosphere.")

