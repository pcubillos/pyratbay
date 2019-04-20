# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["makeTLI"]

import sys
import time
import struct

import numpy as np

from .  import database
from .. import constants as pc
from .. import tools     as pt
from .. import VERSION as ver


def makeTLI(dblist, pflist, dbtype, tlifile,
            wllow,  wlhigh, wlunits, log):
  """
  Driver function to create a TLI file.

  Parameters
  ----------
  dblist: List of strings
      Opacity databases to read.
  pflist: List of strings
      Partition function for each of the databases.
  dbtype: List of strings
      Type of each database.
  tlifile: String
      Output TLI file name.
  wllow: String or float
      Lower wavelength boundary to consider. If float, assume units
      from wlunits input.  Otherwise, wllow sets the value and units
      (for example: '1.0 um').
  wlhigh: String or float
      High wavelength boundary to consider. If float, assume units
      from wlunits input.  Otherwise, wlhigh sets the value and units.
  wlunits: String
      Wavelength units (when not specified in wllow nor wlhigh).
  log: Log object
      An MCcubed.utils.Log instance to log screen outputs to file.
  """
  # Input-not-found error messages:
  if tlifile is None:
      log.error('No output TLI file specified.')

  if wllow is None:
      log.error("Unspecified low wavelength boundary (wllow).")
  if wlhigh is None:
      log.error("Unspecified high wavelength boundary (wlhigh).")

  wllow  = pt.getparam('wllow',  wllow,  wlunits, log)
  wlhigh = pt.getparam('wlhigh', wlhigh, wlunits, log)

  if dblist is None:
      log.error("There are no input database files ('dblist').")
  if dbtype is None:
      log.error("There are no input database types ('dbtype').")
  if pflist is None:
      log.error("There are no partition-function inputs ('pflist').")

  # Check number of files match:
  nfiles = len(dblist)
  if (nfiles != len(pflist)) or (nfiles != len(dbtype)):
      log.error("The number of Line-transition files ({:d}) does not match the "
                "number of partition-function files ({:d}) or database-type "
                "files ({:d}).".format(nfiles, len(pflist), len(dbtype)))

  # Driver routine to read the databases:
  db_readers = {dbname:getattr(database,dbname) for dbname in pc.dbases}

  databases = []
  db_names = []
  log.msg("\nReading input database files:")
  for (dbase, pf, dtype) in zip(dblist, pflist, dbtype):
      if dtype not in db_readers:
          log.error("Unknown type '{:s}' for database '{:s}'.  Select "
                    "from: {:s}".format(dtype, dbase, str(pc.dbases)))
      log.msg("- {:s}".format(dbase))
      databases.append(db_readers[dtype](dbase, pf, log))
      db_names.append(databases[-1].name)
  log.msg("There are {:d} input database file(s).".format(nfiles), verb=2)

  # Open output file:
  tli = open(tlifile, "wb")

  # Get the machine endian type (big/little):
  if sys.byteorder == 'big':
      endian = 'b'
  if sys.byteorder == 'little':
      endian = 'l'

  # Start storing TLI header values:
  header  = struct.pack("s", endian.encode())
  header += struct.pack("3h", ver.LR_VER, ver.LR_MIN, ver.LR_REV)

  # Boundaries in wavenumber space (in cm-1):
  wnlow  = 1.0/wlhigh
  wnhigh = 1.0/wllow

  # Add initial and final wavenumber boundaries (in cm-1):
  header += struct.pack("2d", wnlow, wnhigh)

  Ndb = len(np.unique(db_names))
  header += struct.pack("h", Ndb)
  tli.write(header)

  units = pt.u(wlunits)
  log.msg("\nOS endianness:  {:s}\n"
          "Initial TLI wavelength ({:s}): {:7.3f}  ({:9.3f} cm-1)\n"
          "Final   TLI wavelength ({:s}): {:7.3f}  ({:9.3f} cm-1)".
          format(sys.byteorder,
                 wlunits, wllow/units,  wnhigh,
                 wlunits, wlhigh/units, wnlow), verb=2)
  log.msg("There are {:d} different database(s).".format(Ndb), verb=2)


  log.msg("\nReading and writting partition function info.", verb=2)
  idb = 1         # Database correlative number:
  niso_total = 0  # Cumulative number of isotopes
  acum = [0]      # Cumulative number of isotopes per database
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
      iso_mass  = db.mass
      iso_ratio = db.isoratio

      # Number of temperature samples and isotopes:
      ntemp = len(temp)
      niso  = len(iso_names)

      # Partition-function sorted according to iso_names:
      pf = np.zeros((niso, ntemp), np.double)
      for p,iso in zip(partition, pf_iso):
          idx = iso_names.index(iso)
          pf[idx] = p

      # Store length of and database name:
      tli.write(struct.pack("h{:d}s".format(len(db.name)),
                            len(db.name), db.name.encode("utf-8")))
      # Store the molecule name:
      tli.write(struct.pack("h{:d}s".format(len(db.molecule)),
                            len(db.molecule), db.molecule.encode("utf-8")))
      # Store the number of temperature samples and isotopes:
      tli.write(struct.pack("hh", ntemp, niso))
      log.msg("Database ({:d}/{:d}): '{:s}' ({:s} molecule)".
              format(idb, Ndb, db.name, db.molecule), verb=2, indent=2)
      log.msg("Number of temperatures: {:d}\n"
              "Number of isotopes: {:d}".format(ntemp, niso), verb=2, indent=4)

      # Write the temperature array:
      tli.write(struct.pack("{:d}d".format(ntemp), *temp))
      log.msg("Temperatures (K): [{:6.1f}, {:6.1f}, ..., {:6.1f}]".
              format(temp[0], temp[1], temp[-1]), verb=2, indent=4)

      # For each isotope, write partition function information.
      for j in range(niso):
          iname = iso_names[j]
          log.msg("Isotope ({:d}/{:d}): '{:s}'".format(j+1, niso, iname),
                  verb=2, indent=4)

          # Store length of isotope name, isotope name, and isotope mass:
          tli.write(struct.pack("h{:d}s".format(len(iname)),
                                len(iname), str(iname).encode("utf-8")))
          tli.write(struct.pack("d", iso_mass[j]))
          tli.write(struct.pack("d", iso_ratio[j]))

          # Write the partition function per isotope:
          tli.write(struct.pack("{:d}d".format(ntemp), *pf[j]))
          log.msg("Mass (u):        {:8.4f}\n"
                  "Isotopic ratio:  {:8.4g}\n"
                  "Part. Function:  [{:.2e}, {:.2e}, ..., {:.2e}]".
                  format(iso_mass[j], iso_ratio[j],
                         pf[j,0], pf[j,1], pf[j,-1]), verb=2, indent=6)

      # Calculate cumulative number of isotopes per database:
      niso_total += niso
      idb += 1
      acum.append(niso_total)

  log.msg("Cumulative number of isotopes per database: {}".format(acum), verb=2)


  log.msg("\nExtracting line transition info.")
  wnumber = np.array([], np.double)
  gf      = np.array([], np.double)
  elow    = np.array([], np.double)
  isoID   = np.array([], np.int)
  # Read from file and write the transition info:
  for db in databases:
      # Get database index:
      idb = db_names.index(db.name)

      ti = time.time()
      transitions = db.dbread(wnlow, wnhigh, log.verb)
      tf = time.time()

      if transitions is None:
          continue

      wnumber = np.concatenate((wnumber, transitions[0]))
      gf      = np.concatenate((gf,      transitions[1]))
      elow    = np.concatenate((elow,    transitions[2]))
      isoID   = np.concatenate((isoID,   transitions[3]+acum[idb]))

      unique_iso = np.unique(transitions[3])
      log.msg("Isotope in-database indices: {}".format(unique_iso),
              verb=2, indent=2)
      log.msg("Isotope correlative indices: {}".format(unique_iso+acum[idb]), 
              verb=2, indent=2)
      log.msg("Reading time: {:8.3f} seconds".format(tf-ti), verb=3, indent=2)


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
      ilo  = ihi
      ihi += ntrans
      wnsort = np.argsort(wnumber[isort][ilo:ihi])
      isort[ilo:ihi] = isort[ilo:ihi][wnsort]
  tf = time.time()

  # Actual sorting:
  wnumber = wnumber[isort]
  gf      = gf     [isort]
  elow    = elow   [isort]
  isoID   = isoID  [isort]

  log.msg("Sort time:    {:8.3f} seconds".format(tf-ti), verb=3, indent=2)
  log.msg("\nTransitions per isotope:\n{}".format(ntrans_iso), verb=2)

  # Pack:
  tli.write(struct.pack("i", ntransitions))
  log.msg("\nWriting {:,d} transition lines.".format(ntransitions), verb=2)
  # Write the number of transitions for each isotope:
  niso = len(ntrans_iso)
  tli.write(struct.pack("i", niso))
  tli.write(struct.pack(str(niso)+"i", *list(ntrans_iso)))

  # Write the Line-transition data:
  ti = time.time()
  transinfo  = struct.pack(str(ntransitions)+"d", *list(wnumber))
  transinfo += struct.pack(str(ntransitions)+"h", *list(isoID))
  transinfo += struct.pack(str(ntransitions)+"d", *list(elow))
  transinfo += struct.pack(str(ntransitions)+"d", *list(gf))
  tf = time.time()
  log.msg("Packing time: {:8.3f} seconds".format(tf-ti), verb=3)

  ti = time.time()
  tli.write(transinfo)
  tf = time.time()
  log.msg("Writing time: {:8.3f} seconds".format(tf-ti), verb=3)

  log.msg("Generated TLI file: '{:s}'.".format(tlifile))
  tli.close()
  log.close()
