#!/usr/bin/env python

# FINDME a LICENSE

__all__ = ["makeTLI", "parser"]

import sys, os
import time
import ConfigParser, argparse
import struct
import numpy as np
import matplotlib.pyplot as plt

from .. import tools as pt
from .. import constants as pc
from .  import db
from .. import VERSION as ver

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
      pt.error("Configuration file '{:s}' does not exist.".format(cfile))
    config = ConfigParser.SafeConfigParser()
    config.read([cfile])
    if "lineread" not in config.sections():
      pt.error("Invalid configuration file: '{:s}'.".format(cfile))
    defaults = dict(config.items("lineread"))
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

  # FINDME: For some reason, I need to import 'db' again.
  #         Why has db dissapeared from the namespace?
  from . import db
  args = locals()
  sys.argv = ["ipython"]

  # Process configuration-file arguments:
  cargs = parser(cfile=cfile)

  # Unpack parameters that have not been defined already:
  for key in args.keys():
    if args[key] is None:
      exec("{:s} = cargs.{:s}".format(key, key))

  # Open log file (removing .tli extension):
  logname = outfile.replace(".tli", ".log")
  log = open(logname, "w")
  # Warning log:
  wlog = []

  # Welcome message:
  pt.msg(1, "{:s}\n  Lineread.\n"
            "  Version {:d}.{:d}.{:d}.\n"
            "  Copyright (c) 2016 Patricio Cubillos and collaborators.\n"
            "  Lineread is open-source software under the FINDME license.\n"
            "{:s}\n\n".format(pt.sep, ver.LR_VER, ver.LR_MIN,
                                      ver.LR_REV, pt.sep), log)

  # Input-not-found error messages:
  if dblist is None:
    pt.error("There are no input database files ('dbfile').", log)
  if dbtype is None:
    pt.error("There are no input database types ('dbtype').", log)
  if pflist is None:
    pt.error("There are no partition-function inputs ('pflist').", log)
  # Defaulted-inputs warning messages:
  if iwl == 0.01:
    pt.warning("Using default initial wavelength boundary: iwl = 0.01 um.",
               wlog, log)
  if fwl == 999.9:
    pt.warning("Using default final wavelength boundary: fwl = 999.9 um.",
               wlog, log)

  # Number of files:
  Nfiles = len(dblist)

  # Double-check the number of files:
  Npf   = len(pflist)
  Ntype = len(dbtype)
  if (Nfiles != Npf) or (Nfiles != Ntype):
    pt.error("The number of Line-transition files ({:d}) does not match the "
        "number of partition-function files ({:d}) or database-type "
        "files ({:d}).".format(Nfiles, Npf, Ntype), log)

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
    else:
      pt.error("Unknown Database type ({:d}): '{:s}'".format(i+1, dbtype[i]),
               log)
    pt.msg(verb-10, "File {:d}, database name: '{:s}'".
                                             format(i+1, driver[i].name), log)

  pt.msg(verb, "Beginning to write the TLI file: '{:s}'".format(outfile), log)
  # Open output file:
  TLIout  = open(outfile, "wb")

  # Get the machine endian type (big/little):
  if sys.byteorder == 'big':
    endian = 'b'
  if sys.byteorder == 'little':
    endian = 'l'

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

  pt.msg(verb-8, "Endianness:   {:s}\n"
                 "TLI Version:  {:d}.{:d}.{:d}\n"
                 "Initial wavelength (um): {:7.3f}  ({:9.3f} cm-1)\n"
                 "Final wavelength (um):   {:7.3f}  ({:9.3f} cm-1)".
                 format(endian, ver.LR_VER, ver.LR_MIN, ver.LR_REV,
                        iwl, fwn,  fwl, iwn), log)
  pt.msg(verb-8, "There are {:d} databases in {:d} files.".
                    format(Ndb, Nfiles), log)
  pt.msg(verb-9, "List of databases:\n{}".format(DBnames), log)

  # Partition info:
  totIso = 0                   # Cumulative number of isotopes
  acum = np.zeros(Ndb+1, int)  # Cumul. number of isotopes per database

  pt.msg(verb-2, "Reading and writting partition function info:", log)
  # Database correlative number:
  idb = 0
  # Loop through the partition files (if more than one) and write the
  # data to a processed TLI file:
  for i in np.arange(Nfiles):
    # Skip if we already stored the PF info of this DB:
    if i in DBskip:
      continue

    # Get partition function values:
    Temp, partDB = driver[i].getpf(verb)
    isoNames     = driver[i].isotopes
    iso_mass     = driver[i].mass
    iso_ratio    = driver[i].isoratio

    # Number of temperature samples:
    Ntemp = len(Temp)
    # Number of isotopes:
    Niso  = len(isoNames)

    # DB and molecule name lengths:
    lenDBname = len(DBnames[idb])
    lenMolec  = len(driver[i].molecule)

    # Store length of database name, database name, number of temperatures,
    #  and number of isotopes in TLI file:
    TLIout.write(struct.pack("h{:d}c".format(lenDBname),
                             lenDBname, *DBnames[idb]))
    # Store the molecule name:
    TLIout.write(struct.pack("h{:d}c".format(lenMolec),
                             lenMolec, *driver[i].molecule))
    # Store the number of temperature samples and isotopes:
    TLIout.write(struct.pack("hh", Ntemp, Niso))
    pt.msg(verb-8, "\nDatabase ({:d}/{:d}): '{:s}' ({:s} molecule)\n"
                   "  Number of temperatures: {:d}\n"
                   "  Number of isotopes: {:d}".
                    format(idb+1, Ndb, DBnames[idb], driver[i].molecule,
                           Ntemp, Niso), log)

    # Write the temperature array:
    TLIout.write(struct.pack("{:d}d".format(Ntemp), *Temp))
    pt.msg(verb-8, "Temperatures: [{:6.1f}, {:6.1f}, ..., {:6.1f}]".
                    format(Temp[0], Temp[1], Temp[Ntemp-1]), log, 2)

    # For each isotope, write partition function information.
    # Keep a tally of isotopes for multiple databse support:
    for j in np.arange(Niso):
      pt.msg(verb-9, "Isotope ({:d}/{:d}): '{:s}'".
                      format(j+1, Niso, isoNames[j]), log, 2)

      # Store length of isotope name, isotope name, and isotope mass:
      lenIsoname = len(isoNames[j])
      Iname = str(isoNames[j])
      TLIout.write(struct.pack("h{:d}c".format(lenIsoname), lenIsoname, *Iname))
      TLIout.write(struct.pack("d", iso_mass[j]))
      TLIout.write(struct.pack("d", iso_ratio[j]))

      # Write the partition function per isotope:
      TLIout.write(struct.pack("{:d}d".format(Ntemp), *partDB[j]))
      pt.msg(verb-9, "Mass (u):        {:8.4f}\n"
                     "Isotopic ratio:  {:8.4g}\n"
                     "Part. Function:  [{:.2e}, {:.2e}, ..., {:.2e}]".
                     format(iso_mass[j], iso_ratio[j],
                         partDB[j,0], partDB[j,1], partDB[j,Ntemp-1]), log, 4)

    # Calculate cumulative number of isotopes per database:
    totIso += Niso
    idb += 1
    acum[idb] = totIso

  # Cumulative number of isotopes:
  pt.msg(verb-5, "Cumulative number of isotopes per DB: {}".format(acum), log)
  pt.msg(verb, "Done.", log)

  pt.msg(verb, "\nWriting transition info to TLI file:", log)
  wnumber = np.array([], np.double)
  gf      = np.array([], np.double)
  elow    = np.array([], np.double)
  isoID   = np.array([], np.int)
  # Read from file and write the transition info:
  for db in np.arange(Nfiles):
    # Get database index:
    dbname = driver[db].name
    idb = DBnames.index(dbname)

    # Read databases:
    ti = time.time()
    transDB = driver[db].dbread(iwn, fwn, verb, pflist[db])
    tf = time.time()
    pt.msg(verb-3, "Reading time: {:8.3f} seconds".format(tf-ti), log)

    wnumber = np.concatenate((wnumber, transDB[0]))
    gf      = np.concatenate((gf,      transDB[1]))
    elow    = np.concatenate((elow,    transDB[2]))
    isoID   = np.concatenate((isoID,   transDB[3]+acum[idb]))

    pt.msg(verb-8, "Isotope in-database indices: {}".
                    format(np.unique(transDB[3])), log)
    pt.msg(verb-8, "Isotope correlative indices: {}\n\n".
                    format(np.unique(transDB[3]+acum[idb])), log)

  # Total number of transitions:
  nTransitions = np.size(wnumber)

  # Sort the line transitions (by isotope, then by wavenumber):
  ti = time.time()
  isort = sorted(zip(np.arange(nTransitions), isoID, wnumber),
                 key=lambda x:(x[1], x[2]))
  isort = list(zip(*isort)[0])
  tf = time.time()
  pt.msg(verb-3, "Sort time: {:9.7f} seconds".format(tf-ti), log)

  wnumber = wnumber[isort]
  gf      = gf     [isort]
  elow    = elow   [isort]
  isoID   = isoID  [isort]

  # Calculate the total number of transitions per isotope:
  Nisotran = np.bincount(isoID)
  Nisotran = Nisotran[np.where(Nisotran>0)]  # Remove zeroes
  pt.msg(verb-5, "Transitions per isotope:\n{}".format(Nisotran), log)

  # FINDME: Implement well this:
  if False:
    plt.figure(0)
    plt.clf()
    plt.plot(isoID)
    plt.xlabel("Line index")
    plt.ylabel("Isotope ID")
    plt.savefig("ID.png")

    plt.clf()
    plt.plot(wnumber, ",")
    plt.xlabel("Line index")
    plt.ylabel("Wavelength  (um)")
    plt.savefig("wavelength.png")

  # Pack:
  transinfo = ""
  # Write the number of transitions:
  TLIout.write(struct.pack("i", nTransitions))
  pt.msg(verb-3, "\nWriting {:d} transition lines.".format(nTransitions), log)
  # Write the number of transitions for each isotope:
  Niso = len(Nisotran)
  # Note that nIso may differ from accumiso, since accum iso accounts for
  # all the existing isotopes for an species, whereas nIso accounts only
  # for the isotopes that do have line transitions in the given range.
  TLIout.write(struct.pack("i", Niso))
  TLIout.write(struct.pack(str(Niso)+"i", *list(Nisotran)))

  # Write the Line-transition data:
  ti = time.time()
  transinfo += struct.pack(str(nTransitions)+"d", *list(wnumber))
  transinfo += struct.pack(str(nTransitions)+"h", *list(isoID))
  transinfo += struct.pack(str(nTransitions)+"d", *list(elow))
  transinfo += struct.pack(str(nTransitions)+"d", *list(gf))
  tf = time.time()
  pt.msg(verb-3, "Packing time: {:8.3f} seconds".format(tf-ti), log)

  ti = time.time()
  TLIout.write(transinfo)
  tf = time.time()
  pt.msg(verb-3, "Writing time: {:8.3f} seconds".format(tf-ti), log)

  TLIout.close()
  pt.msg(verb, "\nGenerated output TLI file: '{:s}'".format(outfile), log)

  if len(wlog) > 0:
    # Write all warnings to file:
    wpath, wfile = os.path.split(os.path.realpath(logname))
    wfile = "{:s}/warnings_{:s}".format(wpath, wfile)
    warns = open(wfile, "w")
    warns.write("Warnings log:\n\n{:s}\n".format(pt.sep))
    warns.write("\n\n{:s}\n".format(pt.sep).join(wlog))
    warns.close()
    # Report it:
    pt.warning("There was(were) {:d} warning(s) raised.\nSee '{:s}'.".
                format(len(wlog), wfile), [], log)

  pt.msg(verb, "Done.\n", log)
  log.close()

  sys.exit(0)
