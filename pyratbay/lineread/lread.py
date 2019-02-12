# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["makeTLI", "parser"]

import sys
import os
import time
import argparse
import struct
from datetime import date
if sys.version_info.major == 3:
    import configparser
else:
    import ConfigParser as configparser

import numpy as np
import matplotlib.pyplot as plt

from .  import db
from .. import constants as pc
from .. import VERSION as ver

sys.path.append(pc.ROOT + "/modules/MCcubed")
import MCcubed.utils as mu


def parser(cfile=None):
  """
  Read and process the command line arguments.
  """
  # Parser to process a configuration file:
  cparser = argparse.ArgumentParser(description=__doc__, add_help=False,
                           formatter_class=argparse.RawDescriptionHelpFormatter)
  # Add config file option:
  cparser.add_argument("-c", "--cfile",
                       help="Configuration filename (string).", metavar="FILE")
  # remaining_argv contains all other command-line-arguments:
  args, remaining_argv = cparser.parse_known_args()

  # Get parameters from configuration file (if exists):
  if args.cfile:
    cfile = args.cfile

  # Parse the configuration-file arguments:
  if cfile is not None:
    if not os.path.isfile(cfile):
      print("Configuration file '{:s}' does not exist.".format(cfile))
      sys.exit(0)
    config = configparser.SafeConfigParser()
    config.read([cfile])
    if "pyrat" not in config.sections():
      print("Invalid configuration file: '{:s}'. The configuration-file "
            "section must be 'pyrat'.".format(cfile))
      sys.exit(0)
    defaults = dict(config.items("pyrat"))
    # Store these arguments as lists:
    if "dblist" in defaults:
      defaults["dblist"] = defaults["dblist"].split()
    if "pflist" in defaults:
      defaults["pflist"] = defaults["pflist"].split()
    if "dbtype" in defaults:
      defaults["dbtype"] = defaults["dbtype"].split()
  else:
    defaults = {}

  # Inherit options from cparser:
  parser = argparse.ArgumentParser(parents=[cparser])

  # General Options:
  parser.add_argument("-v", "--verbose-level", action="store",
                       help="Verbosity level (integer) [default: %(default)s].",
                       dest="verb", type=int, default=2)
  parser.add_argument("-q", "--quiet",         action="store_false",
                       help="Set verbosity level to 0.",
                       dest="verb")
  # Database Options:
  group = parser.add_argument_group("Database Options")
  group.add_argument("-o", "--outfile",        action  = "store",
                     help="Output filename (string) [default: '%(default)s'].",
                     dest= "outfile",          default = "output.tli")
  group.add_argument("-d", "--database",       action="append",
                     help="Path (string) to the input line-transition "
                          "database file(s).",
                     dest="dblist")
  group.add_argument("-p", "--partition",      action="append",
                     help="Path (string) to auxiliary partition-function "
                          "file(s).",
                     dest="pflist")
  group.add_argument("-t", "--dbtype",         action="append",
                     help="Database type (string).  'ps' for Partridge & "
                          "Schwenke's H2O; 'hit' for HITRAN and HITEMP; or "
                          "'ts' for Schwenke's TiO.",
                     dest="dbtype")
  # Wavelength Options:
  group = parser.add_argument_group("Wavelength Options")
  group.add_argument("-i", "--wl-init",       action="store",
                     help="Initial wavelength (microns) [default: "
                          "%(default)s].",
                     dest="iwl", type=float, default=0.01)
  group.add_argument("-f", "--wl-final",      action="store",
                     help="Final wavelength (microns) [default: %(default)s].",
                     dest="fwl", type=float, default=999.9)
  parser.set_defaults(**defaults)
  args = parser.parse_args(remaining_argv)

  return args


def makeTLI(dblist=None, pflist=None, dbtype=None, outfile=None,
            iwl=None,    fwl=None,    verb=None,   cfile=None):
  """
  Driver function to create a TLI file.

  Uncredited developers
  ---------------------
  Patricio Rojo  (Cornell U.)
  Madison Stemm  (UCF)
  """
  # Open log file (removing .tli extension):
  logname = outfile.replace(".tli", ".log")
  log = mu.Log(logname, verb=verb, width=80)

  # Welcome message:
  log.msg("{:s}\n"
          "  Lineread.\n"
          "  Version {:d}.{:d}.{:d}.\n"
          "  Copyright (c) 2016-{:d} Patricio Cubillos and collaborators.\n"
          "  Lineread is (temporarily) proprietaty software (see LICENSE).\n"
          "{:s}\n\n".format(log.sep, ver.LR_VER, ver.LR_MIN, ver.LR_REV,
                            date.today().year, log.sep))

  # Input-not-found error messages:
  if dblist is None:
    log.error("There are no input database files ('dblist').")
  if dbtype is None:
    log.error("There are no input database types ('dbtype').")
  if pflist is None:
    log.error("There are no partition-function inputs ('pflist').")
  # Defaulted-inputs warning messages:
  if iwl == 0.01:
    log.warning("Using default initial wavelength boundary: iwl = 0.01 um.")
  if fwl == 999.9:
    log.warning("Using default final wavelength boundary: fwl = 999.9 um.")

  # Number of files:
  Nfiles = len(dblist)

  # Double-check the number of files:
  Npf   = len(pflist)
  Ntype = len(dbtype)
  if (Nfiles != Npf) or (Nfiles != Ntype):
    log.error("The number of Line-transition files ({:d}) does not match the "
        "number of partition-function files ({:d}) or database-type "
        "files ({:d}).".format(Nfiles, Npf, Ntype))

  # Driver routine to read the databases:
  driver = []
  for i in np.arange(Nfiles):
    if   dbtype[i] == "ps":
      driver.append(db.pands(      dblist[i], pflist[i], log))
    elif dbtype[i] == "hit":
      driver.append(db.hitran(     dblist[i], pflist[i], log))
    elif dbtype[i] == "ts":
      driver.append(db.tioschwenke(dblist[i], pflist[i], log))
    elif dbtype[i] == "vo":
      driver.append(db.voplez(     dblist[i], pflist[i], log))
    elif dbtype[i] == "vald":
      driver.append(db.vald(       dblist[i], pflist[i], log))
    elif dbtype[i] == "emol":
      driver.append(db.exomol(     dblist[i], pflist[i], log))
    elif dbtype[i] == "repack":
      driver.append(db.repack(     dblist[i], pflist[i], log))
    else:
      log.error("Unknown Database type ({:d}): '{:s}'.  Select from: {:s}".
                format(i+1, dbtype[i], pc.dbases))
    log.msg("Reading input database file '{:s}'.".format(dblist[i]))
  log.msg("There are {:d} input database file(s).".format(Nfiles), verb=2)

  # Open output file:
  TLIout  = open(outfile, "wb")

  # Get the machine endian type (big/little):
  if sys.byteorder == 'big':
    endian = b'b'
  if sys.byteorder == 'little':
    endian = b'l'

  # Start storing TLI header values:
  header = endian
  header += struct.pack("3h", ver.LR_VER, ver.LR_MIN, ver.LR_REV)

  # Boundaries in wavenumber space (in cm-1):
  iwn = 1.0/(fwl*pc.um)
  fwn = 1.0/(iwl*pc.um)

  # Add initial and final wavenumber boundaries (in cm-1):
  header += struct.pack("2d", iwn, fwn)

  # Get number of databases:
  DBnames = []  # Database names
  DBskip  = []  # Index of repeated databases
  for i in np.arange(Nfiles):
    dbname = driver[i].name
    if dbname in DBnames:
      DBskip.append(i)  # Ommit repeated databases
    else:
      DBnames.append(dbname)
  Ndb = len(DBnames)
  header += struct.pack("h", Ndb)
  TLIout.write(header)

  log.msg("\nOS endianness:  {:s}\n"
          "Initial TLI wavelength (um): {:7.3f}  ({:9.3f} cm-1)\n"
          "Final   TLI wavelength (um): {:7.3f}  ({:9.3f} cm-1)".
          format(sys.byteorder, iwl, fwn,  fwl, iwn), verb=2)
  log.msg("There are {:d} different database(s).".format(Ndb), verb=2)
  log.msg("List of databases:\n{}".format(DBnames), verb=2)

  # Partition info:
  totIso = 0                   # Cumulative number of isotopes
  acum = np.zeros(Ndb+1, int)  # Cumul. number of isotopes per database


  log.msg("\nReading and writting partition function info.", verb=2)
  # Database correlative number:
  idb = 0
  # Loop through the partition files (if more than one) and write the
  # data to a processed TLI file:
  for i in np.arange(Nfiles):
    # Skip if we already stored the PF info of this DB:
    if i in DBskip:
      continue

    # Get partition function values:
    Temp, partDB, PFiso = driver[i].getpf(verb)
    isoNames     = driver[i].isotopes
    iso_mass     = driver[i].mass
    iso_ratio    = driver[i].isoratio

    # Number of temperature samples:
    Ntemp = len(Temp)
    # Number of isotopes:
    Niso  = len(isoNames)

    # Partition-function sorted according to isoNames:
    PF = np.zeros((Niso, Ntemp), np.double)
    for j in np.arange(np.shape(partDB)[0]):
      idx = isoNames.index(PFiso[j])
      PF[idx] = partDB[j]

    # DB and molecule name lengths:
    lenDBname = len(DBnames[idb])
    lenMolec  = len(driver[i].molecule)

    # Store length of database name, database name, number of temperatures,
    #  and number of isotopes in TLI file:
    TLIout.write(struct.pack("h{:d}s".format(lenDBname),
                             lenDBname, DBnames[idb].encode("utf-8")))
    # Store the molecule name:
    TLIout.write(struct.pack("h{:d}s".format(lenMolec),
                             lenMolec, driver[i].molecule.encode("utf-8")))
    # Store the number of temperature samples and isotopes:
    TLIout.write(struct.pack("hh", Ntemp, Niso))
    log.msg("Database ({:d}/{:d}): '{:s}' ({:s} molecule)".format(
            idb+1, Ndb, DBnames[idb], driver[i].molecule), verb=2, indent=2)
    log.msg("Number of temperatures: {:d}\n"
            "Number of isotopes: {:d}".format(Ntemp, Niso), verb=2, indent=4)

    # Write the temperature array:
    TLIout.write(struct.pack("{:d}d".format(Ntemp), *Temp))
    log.msg("Temperatures (K): [{:6.1f}, {:6.1f}, ..., {:6.1f}]".
            format(Temp[0], Temp[1], Temp[Ntemp-1]), verb=2, indent=4)

    # For each isotope, write partition function information.
    # Keep a tally of isotopes for multiple databse support:
    for j in np.arange(Niso):
      log.msg("Isotope ({:d}/{:d}): '{:s}'".format(j+1, Niso, isoNames[j]),
              verb=2, indent=4)

      # Store length of isotope name, isotope name, and isotope mass:
      lenIsoname = len(isoNames[j])
      Iname = str(isoNames[j])
      TLIout.write(struct.pack("h{:d}s".format(lenIsoname), lenIsoname,
                   Iname.encode()))
      TLIout.write(struct.pack("d", iso_mass[j]))
      TLIout.write(struct.pack("d", iso_ratio[j]))

      # Write the partition function per isotope:
      TLIout.write(struct.pack("{:d}d".format(Ntemp), *PF[j]))
      log.msg("Mass (u):        {:8.4f}\n"
              "Isotopic ratio:  {:8.4g}\n"
              "Part. Function:  [{:.2e}, {:.2e}, ..., {:.2e}]".
              format(iso_mass[j], iso_ratio[j],
                  PF[j,0], PF[j,1], PF[j,Ntemp-1]), verb=2, indent=6)

    # Calculate cumulative number of isotopes per database:
    totIso += Niso
    idb += 1
    acum[idb] = totIso

  # Cumulative number of isotopes:
  log.msg("Cumulative number of isotopes per DB: {}".format(acum), verb=2)
  log.msg("Lineread done.")


  log.msg("\nExtracting line transition info.")
  wnumber = np.array([], np.double)
  gf      = np.array([], np.double)
  elow    = np.array([], np.double)
  isoID   = np.array([], np.int)
  # Read from file and write the transition info:
  for dbase in np.arange(Nfiles):
    # Get database index:
    dbname = driver[dbase].name
    idb = DBnames.index(dbname)

    # Read databases:
    ti = time.time()
    transDB = driver[dbase].dbread(iwn, fwn, verb, pflist[dbase])
    tf = time.time()

    if transDB is None:
      continue

    wnumber = np.concatenate((wnumber, transDB[0]))
    gf      = np.concatenate((gf,      transDB[1]))
    elow    = np.concatenate((elow,    transDB[2]))
    isoID   = np.concatenate((isoID,   transDB[3]+acum[idb]))

    log.msg("Isotope in-database indices: {}".format(np.unique(transDB[3])),
            verb=2, indent=2)
    log.msg("Isotope correlative indices: {}".
                    format(np.unique(transDB[3]+acum[idb])), verb=2, indent=2)
    log.msg("Reading time: {:8.3f} seconds".format(tf-ti), verb=3, indent=2)


  # Total number of transitions:
  nTransitions = np.size(wnumber)
  # Number of transitions per isotope:
  Nisotran = np.bincount(isoID)
  Nisotran = Nisotran[np.where(Nisotran>0)]  # Remove zeroes

  # Sort by isotope ID:
  ti = time.time()
  isort   = np.argsort(isoID)
  # Sort each isotope by wavenumber:
  ihi = 0
  for j in np.arange(len(Nisotran)):
    ilo  = ihi
    ihi += Nisotran[j]
    wnsort = np.argsort(wnumber[isort][ilo:ihi])
    isort[ilo:ihi] = isort[ilo:ihi][wnsort]
  tf = time.time()
  log.msg("Sort time:    {:8.3f} seconds".format(tf-ti), verb=3, indent=2)
  log.msg("Done.")

  # Actual sorting:
  wnumber = wnumber[isort]
  gf      = gf     [isort]
  elow    = elow   [isort]
  isoID   = isoID  [isort]

  log.msg("\nTransitions per isotope:\n{}".format(Nisotran), verb=2)

  # FINDME: Implement well this:
  if False:
    plt.figure(0)
    plt.clf()
    plt.plot(isoID)
    plt.xlabel("Line index")
    plt.ylabel("Isotope ID")
    plt.ylim(-0.1, np.amax(isoID)+0.1)
    plt.savefig("ID.png")

    plt.clf()
    plt.plot(wnumber, ",")
    plt.xlabel("Line index")
    plt.ylabel("Wavelength  (um)")
    plt.savefig("wavelength.png")

  # Pack:
  # Write the number of transitions:
  TLIout.write(struct.pack("i", nTransitions))
  log.msg("\nWriting {:,d} transition lines.".format(nTransitions), verb=2)
  # Write the number of transitions for each isotope:
  Niso = len(Nisotran)
  # Note that nIso may differ from accumiso, since accum iso accounts for
  # all the existing isotopes for an species, whereas nIso accounts only
  # for the isotopes that do have line transitions in the given range.
  TLIout.write(struct.pack("i", Niso))
  TLIout.write(struct.pack(str(Niso)+"i", *list(Nisotran)))

  # Write the Line-transition data:
  ti = time.time()
  transinfo  = struct.pack(str(nTransitions)+"d", *list(wnumber))
  transinfo += struct.pack(str(nTransitions)+"h", *list(isoID))
  transinfo += struct.pack(str(nTransitions)+"d", *list(elow))
  transinfo += struct.pack(str(nTransitions)+"d", *list(gf))
  tf = time.time()
  log.msg("Packing time: {:8.3f} seconds".format(tf-ti), verb=3)

  ti = time.time()
  TLIout.write(transinfo)
  tf = time.time()
  log.msg("Writing time: {:8.3f} seconds".format(tf-ti), verb=3)

  TLIout.close()
  log.msg("Generated TLI file: '{:s}'.".format(outfile))

  if len(log.warnings) > 0:
    # Write all warnings to file:
    wpath, wfile = os.path.split(os.path.realpath(logname))
    wfile = "{:s}/warnings_{:s}".format(wpath, wfile)
    with open(wfile, "w") as warnings:
      warnings.write("Warnings log:\n\n{:s}\n".format(log.sep))
      warnings.write("\n\n{:s}\n".format(log.sep).join(log.warnings))
    # Report it:
    flag = len(log.warnings) > 1
    log.msg("\n{:s}\nThere {:s} {:d} warning{:s} raised.  See '{:s}'.".
            format(log.sep, ['was','were'][flag], len(log.warnings),
                   ["","s"][flag], wfile), verb=0)

  log.close()
