# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import numpy as np

from .  import objects   as o
from .. import tools     as pt
from .. import constants as pc


def readlinedb(pyrat):
  """
  Main driver to read the line transition data from TLI files.
  """
  # Count number of TLI files:
  if pyrat.linedb is None:
    pyrat.lt.nTLI = 0
    pyrat.log.msg("\nNo line transition file to read.")
  else:
    pyrat.lt.nTLI = len(pyrat.linedb)
    pyrat.log.msg("\nReading line transition info.")

  # TLI file object:
  TLI = []
  # Index of first database in TLI file:
  dbindex = [0]

  # Read data bases header info:
  for n in np.arange(pyrat.lt.nTLI):
    # Open-read TLI data base:
    TLI.append(open(pyrat.linedb[n], "r"))
    pyrat.log.msg("Read TLI file: '{:s}'.".format(pyrat.linedb[n]), indent=2)
    # Read headers info:
    dbindex.append(dbindex[-1] + readheader(pyrat, TLI[n]))

  # Set link to molecules' indices:
  setimol(pyrat)

  # Read line-transition data (if there's no extinction-coefficient table):
  if (pyrat.ex.extfile is None) or (not os.path.isfile(pyrat.ex.extfile)):
    for n in np.arange(pyrat.lt.nTLI):
      readlinetransition(pyrat, TLI[n], dbindex[n])

  pyrat.log.msg("Read a total of {:,d} line transitions.".
                format(pyrat.lt.ntransitions), verb=2, indent=2)
  pyrat.log.msg("Line-transition done.\n")


def readheader(pyrat, linefile):
  """
  Open TLI file and read database header info.

  Returns:
  --------
  Ndb: Integer
     Number of databases in this TLI file.
  """
  # Read header:
  endian = linefile.read(1)
  if endian == 'b':
    pyrat.log.msg("TLI data storage: Big endian", verb=2, indent=2)
  elif endian == 'l':
    pyrat.log.msg("TLI data storage: Little endian", verb=2, indent=2)

  # Compare endianness:
  endianness = sys.byteorder
  if endianness[0:1] != endian:
    pyrat.log.error("Incompatible endianness between TLI file ({:s}) and "
                    "Pyrat ({:s}).".format(endian, endianness[0:1]))

  # Read TLI version:
  TLI_ver, TLI_min, TLI_rev = pt.unpack(linefile, 3, "h")
  pyrat.log.msg("TLI version: {:d}.{:d}.{:d}.".
                format(TLI_ver, TLI_min, TLI_rev), verb=2, indent=2)
  if (TLI_ver != 6) or (TLI_min not in [1,2,3,4]):
    pyrat.log.error("Incompatible TLI version.  The TLI file must be created "
                    "with Lineread version 6.1-6.4.")

  # Read initial and final wavenumber from TLI:
  lt_wni, lt_wnf = pt.unpack(linefile, 2, "d")
  pyrat.log.msg("TLI wavenumber range (cm-1): [{:.1f}, {:.1f}]".
                format(lt_wni, lt_wnf), verb=2, indent=2)
  # Check TLI and pyrat wavelength ranges:
  checkrange(pyrat, lt_wni, lt_wnf)

  # Read number of data bases:
  Ndb = pt.unpack(linefile, 1, "h")
  pyrat.log.msg("Number of data bases: {:d}".format(Ndb), verb=2, indent=2)

  # Cumulative isotope index:
  acumiso = pyrat.iso.niso

  for i in np.arange(Ndb):
    db = o.Database()
    # Read Database name:
    lenDBname = pt.unpack(linefile, 1,         "h")
    db.name   = pt.unpack(linefile, lenDBname, "s")
    pyrat.log.msg("Data base name: '{:s}'".format(db.name), verb=2, indent=2)
    # Read Molecule name:
    lenMolec   = pt.unpack(linefile, 1,        "h")
    db.molname = pt.unpack(linefile, lenMolec, "s")
    pyrat.log.msg("Molecule name: '{:s}'".format(db.molname), verb=2, indent=2)
    # Read temperature array:
    db.ntemp, db.niso =  pt.unpack(linefile, 2,        "h")
    db.temp = np.asarray(pt.unpack(linefile, db.ntemp, "d"))
    # Update temperature boundaries:
    pyrat.lt.tmin = np.amax((pyrat.lt.tmin, db.temp[ 0]))
    pyrat.lt.tmax = np.amin((pyrat.lt.tmax, db.temp[-1]))
    pyrat.log.msg("Temperature range: {:4.1f} -- {:4.1f} K.".
                  format(db.temp[0], db.temp[-1]), verb=2, indent=2)

    # Allocate arrays for isotopic info:
    name    = np.zeros(db.niso, pc.strfmt)
    mass    = np.zeros(db.niso)
    ratio   = np.zeros(db.niso)
    dbindex = np.zeros(db.niso, np.int)
    db.z    = np.zeros((db.niso, db.ntemp))

    # Store per-isotope info:
    for j in np.arange(db.niso):
      dbindex[j] = i + pyrat.lt.ndb
      lenIsoName = pt.unpack(linefile, 1,          "h")
      name[j]    = pt.unpack(linefile, lenIsoName, "s")
      mass[j]    = pt.unpack(linefile, 1,          "d")
      ratio[j]   = pt.unpack(linefile, 1,          "d")
      db.z[j]    = np.asarray(pt.unpack(linefile, db.ntemp, "d"))

      # Print info to screen:
      pyrat.log.msg("Isotope: {:s},  mass: {:.4f} u,  isotopic ratio: {:.5g}".
                    format(name[j], mass[j], ratio[j]), verb=2, indent=3)
      pyrat.log.msg("Z = [{:.2e}, {:.2e}, ..., {:.2e}]".
                    format(db.z[j,0], db.z[j,1], db.z[j,-1]), verb=3, indent=4)

    # Add the number of isotopes read:
    pyrat.iso.niso += db.niso
    pyrat.log.msg("Number of isotopes: {:d}".format(pyrat.iso.niso),
                  verb=3, indent=2)

    # Store name, mass in isotopes structure:
    pyrat.iso.name    = np.concatenate((pyrat.iso.name,    name))
    pyrat.iso.mass    = np.concatenate((pyrat.iso.mass,    mass))
    pyrat.iso.ratio   = np.concatenate((pyrat.iso.ratio,   ratio))
    pyrat.iso.dbindex = np.concatenate((pyrat.iso.dbindex, dbindex))

    # Set isotope correlative index for DB:
    db.iiso  = acumiso
    acumiso += db.niso
    pyrat.lt.db.append(db)
    pyrat.log.msg("DB index: {:d}, Cumulative Iso: {:d}\n\n".
                  format(i+pyrat.lt.ndb, acumiso), verb=2, indent=2)

  # Keep count of number of databases:
  pyrat.lt.ndb  += Ndb

  # Return the pointer position in lineinfo file:
  return Ndb


def readlinetransition(pyrat, linefile, dbindex):
  """
  Read the databases line transition info.
  """
  # Read the number of line transitions:
  nTransitions = pt.unpack(linefile, 1,    "i")
  # Read the number of isotopes in line-transition array:
  nIso         = pt.unpack(linefile, 1,    "i")
  # Read the number of transitions per isotope:
  NisoTran     = np.atleast_1d(pt.unpack(linefile, nIso, "i"))

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
    pyrat.log.error("The remaining data file size ({:.1f}) does not correspond "
             "to the number of transitions ({:d}).".format(nrec, nTransitions))
  pyrat.log.msg("There are {:,d} line transitions in TLI file.".
                format(nTransitions), verb=2, indent=2)

  # Allocate arrays:
  wn    = np.zeros(nTransitions)
  isoid = np.zeros(nTransitions, np.short)
  elow  = np.zeros(nTransitions)
  gf    = np.zeros(nTransitions)

  # Number of line-transitions offset for a given isotope:
  offset = 0
  start  = init_wl
  nlt     = 0  # Total number of line-transitions read
  for i in np.arange(nIso):
    # Search lower and higher line-transition indices to read:
    ifirst = pt.binsearch(linefile, pyrat.spec.wnlow,  start, NisoTran[i]-1,
                          False)
    ilast  = pt.binsearch(linefile, pyrat.spec.wnhigh, start, NisoTran[i]-1,
                          True)
    # Add offset for this isotope:
    ifirst += offset
    ilast  += offset

    # Print boundaries:
    linefile.seek(ifirst*pc.dreclen + init_wl, 0)
    pyrat.log.msg("Found initial transition ({:8d}):  {:13.4f} cm-1".
                  format(ifirst, pt.unpack(linefile, 1, 'd')), verb=3, indent=2)
    # Pyrat low-wavelength boundary in microns:
    linefile.seek(ilast*pc.dreclen  + init_wl, 0)
    pyrat.log.msg("Found final   transition ({:8d}):  {:13.4f} cm-1".
                  format(ilast, pt.unpack(linefile, 1, 'd')), verb=3, indent=2)

    # Number of transitions to read:
    nread = ilast - ifirst + 1
    # Get isotope ID:
    i0 = pt.binsearch(linefile, 0, start, NisoTran[i]-1, False)
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
    pyrat.log.msg("Read {:11,d} transitions for isotope {:2d}.".
        format(nread, pyrat.lt.db[dbindex].iiso+isoID), verb=2, indent=4)

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
  if (wn_low > pyrat.spec.wnhigh) or (wn_high < pyrat.spec.wnlow):
    pyrat.log.warning("TLI wavenumber range ({:.2f} - {:.2f} cm^-1) "
       "does not overlap with Pyrat wavenumber range ({:.2f} - {:.2f} cm^-1).".
        format(wn_low, wn_high, pyrat.spec.wnlow, pyrat.spec.wnhigh))
  # Print out warning if TLI range is smaller than the pyrat required range:
  elif (wn_low > pyrat.spec.wnlow) or  (wn_high < pyrat.spec.wnhigh):
    pyrat.log.warning("TLI wavenumber range ({:.2f} - {:.2f} cm^-1) "
       "does not cover the full Pyrat wavenumber range ({:.2f} - {:.2f} "
       "cm^-1).".format(wn_low, wn_high, pyrat.spec.wnlow, pyrat.spec.wnhigh))


def setimol(pyrat):
  """
  Set the molecule index for the list of isotopes.
  """
  # Allocate imol array:
  pyrat.iso.imol = np.zeros(pyrat.iso.niso, np.int)
  # For each isotope:
  for i in np.arange(pyrat.iso.niso):
    # Get molecule name from database object:
    molname = pyrat.lt.db[pyrat.iso.dbindex[i]].molname
    # Set index:
    if molname not in pyrat.mol.symbol:  # Isotope's not in molecule list
      pyrat.iso.imol[i] = -1
    else:
      pyrat.iso.imol[i] = np.where(pyrat.mol.symbol == molname)[0]
  pyrat.log.msg("Isotope's molecule indices:\n  {:s}".
                format(str(pyrat.iso.imol)), verb=2, indent=2)
  # Report missing species:
  imiss = np.unique(pyrat.iso.dbindex[np.where(pyrat.iso.imol < 0)])
  for i in imiss:
    pyrat.log.warning("The species '{:s}' for isotopes {:s} is not "
       "present in the atmosphere.".format(pyrat.lt.db[i].molname,
       pyrat.iso.name[np.where(pyrat.iso.dbindex==i)]))

