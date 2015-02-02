#!/usr/bin/env python

# ****************************** START LICENSE ******************************
#
# FINDME a LICENSE
#
# ******************************* END LICENSE ******************************

import matplotlib.pyplot as plt

import numpy as np
import scipy.constants as sc
import sys, os, time

import ConfigParser
import argparse
import struct
import heapq as hq

# Main dir (where lineread.py is located):
maindir = os.path.dirname(os.path.realpath(__file__))
# Add paths to Python a C folders:
sys.path.append(maindir + '/pysrc/')
sys.path.append(maindir + '/csrc/lib/')

import pconstants as pc
import ptools     as pt
import db_pands   as ps
import db_hitran  as hit
import db_tioschwenke as ts


def parseargs():
  """
  Read and process the command line arguments.

  Returns:
  --------
  args: Namespace object
     Object with the command line arguments.

  Modification History:
  ---------------------
  2013        madison   Initial implementation.        
  2014-03-06  patricio  Updated from optparse to argparse. Added documentation
                        and config_file option.    pcubillos@fulbrightmail.org
  2015-02-01  patricio  
  """
  # Parser to process a configuration file:
  cparser = argparse.ArgumentParser(description=__doc__, add_help=False,
                           formatter_class=argparse.RawDescriptionHelpFormatter)
  # Add config file option:
  cparser.add_argument("-c", "--config_file",
                       help="Configuration filename (string).", metavar="FILE")
  # remaining_argv contains all other command-line-arguments:
  args, remaining_argv = cparser.parse_known_args()

  # Get parameters from configuration file (if exists):
  if args.config_file:
    if not os.path.isfile(args.config_file):
      pt.error("Configuration file '{:s}' does not exist.".
                format(args.config_file))
    config = ConfigParser.SafeConfigParser()
    config.read([args.config_file])
    if "Parameters" not in config.sections():
      pt.error("Invalid configuration file: '{:s}'.".format(args.config_file))
    defaults = dict(config.items("Parameters"))
    # Store these arguments as lists:
    if "db_list" in defaults:
      defaults["db_list"]   = defaults["db_list"].split()
    if "part_list" in defaults:
      defaults["part_list"] = defaults["part_list"].split()
    if "dbtype" in defaults:
      defaults["dbtype"]    = defaults["dbtype"].split()
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
  group.add_argument("-o", "--output",         action  = "store",
                     help="Output filename (string) [default: '%(default)s'].",
                     dest= "output",           default = "output.tli")
  group.add_argument("-d", "--database",       action="append",
                     help="Path (string) to the input line-transition "
                          "database file(s).",
                     dest="db_list")
  group.add_argument("-p", "--partition",      action="append",
                     help="Path (string) to auxiliary partition-function "
                          "file(s).",
                     dest="part_list")
  group.add_argument("-t", "--dbtype",         action="append",
                     help="Database type (string).  'ps' for Partridge & "
                          "Schwenke's H2O; 'hit' for HITRAN and HITEMP; or "
                          "'ts' for Schwenke's TiO.",
                     choices=('ps', 'hit', 'ts'),
                     dest="dbtype")
  # Wavelength Options:
  group = parser.add_argument_group("Wavelength Options")
  group.add_argument("-i", "--wl-init",       action="store",
                     help="Initial wavelength (microns) [default: "
                          "%(default)s].",
                     dest="iwav", type=float, default=1.0)
  group.add_argument("-f", "--wl-final",      action="store",
                     help="Final wavelength (microns) [default: %(default)s].",
                     dest="fwav", type=float, default=2.0)
  parser.set_defaults(**defaults)
  args = parser.parse_args(remaining_argv)

  return args


if __name__ == "__main__":
  """
  Main function to create a TLI file.

  Usage:
  ------
  Execute from the shell:
  ./lineread.py [--option <args>]

  To dysplay the list of command-line arguments execute:
  ./lineread.py --help

  Modification History:
  ---------------------
  2013-10-21  madison   Initial version based on P. Rojo's lineread C code.
                                                 madison.stemm@ucf.edu
  2014-03-05  patricio  Added documentation and updated Madison's code.
                                                 pcubillos@fulbrightmail.org
  2014-07-27  patricio  Updated to version 5.0
  2015-02-01  patricio  
  """

  # Process command-line-arguments:
  cla = parseargs()
  # Unpack parameters:
  verb       = cla.verb
  dblist     = cla.db_list
  pflist     = cla.part_list
  dbtype     = cla.dbtype 
  outputfile = cla.output

  # Number of files:
  Nfiles = len(dblist)
  # Driver routine to read the databases:
  driver = []
  for i in np.arange(Nfiles):
    if   dbtype[i] == "ps":
      driver.append(ps.pands(dblist[i], pflist[i]))
    elif dbtype[i] == "hit":
      driver.append(hit.hitran(dblist[i], pflist[i]))
    elif dbtype[i] == "ts":
      driver.append(ts.tioschwenke(dblist[i], pflist[i]))
    else:
      pt.error("Unknown Database type ({:d}): '{:s}'".format(i+1, dbtype[i]))
    pt.msg(verb-10, "File {:d}, database name: '{:s}'".
                       format(i+1, driver[i].name))

  pt.msg(verb, "Beginning to write the TLI file: '{:s}'".format(outputfile))
  # Open output file:
  TLIout  = open(outputfile, "wb")

  # Get the machine endian type (big/little):
  endianness = sys.byteorder

  # Hardcoded implementation of lineread's "magic number" check for endianness
  # derived from binary encoding the letters TLI into binary location
  # and checking the order
  if endianness == 'big':
    endian = 'b'
  if endianness == 'little':
    endian = 'l'

  # TLI header: Tells endianness of binary, TLI version, and number of
  # databases used.
  header = endian
  header += struct.pack("3h", pc.TLI_VER, pc.LR_VER, pc.LR_REV)

  # Boundaries in wavenumber space (in cm-1):
  iwn = 1.0/(cla.fwav*pc.MTC)
  fwn = 1.0/(cla.iwav*pc.MTC)

  # Add initial and final wavenumber boundaries:
  header += struct.pack("2d", iwn, fwn)

  # Get number of databases:
  DBnames = [] # Database names
  DBskip  = [] # Index of repeated databases
  for i in np.arange(Nfiles):
    dbname = driver[i].name
    if dbname in DBnames:
      DBskip.append(i) # Ommit repeated databases
    else:
      DBnames.append(dbname) 
  Ndb = len(DBnames)
  header += struct.pack("h", Ndb)
  TLIout.write(header)

  pt.msg(verb-8, "Endianness:   {:s}\n"
                 "TLI Version:  {:d}\n"
                 "LR Version:   {:d}\n"
                 "LR Revision:  {:d}\n"
                 "Initial wavelength (um): {:7.3f}  ({:9.3f} cm-1)\n"
                 "Final wavelength (um):   {:7.3f}  ({:9.3f} cm-1)".
                 format(endian, pc.TLI_VER, pc.LR_VER, pc.LR_REV,
                        cla.iwav, fwn,  cla.fwav, iwn))
  pt.msg(verb-8, "There are {:d} databases in {:d} files.".
                    format(Ndb, Nfiles))
  pt.msg(verb-9, "List of databases: {}".format(DBnames))

  # Partition info:
  totIso = 0                  # Cumulative number of isotopes
  acum = np.zeros(Ndb+1, int) # Cumul. number of isotopes per database

  pt.msg(verb-2, "Reading and writting partition function info:")
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
    isoNames  = driver[i].isotopes
    iso_mass  = driver[i].mass
    iso_ratio = driver[i].isoratio

    #pt.msg(verb, "Isotopes mass: {}".format(iso_mass))

    # Number of temperature samples:
    Ntemp = len(Temp)
    # Number of isotopes:
    Niso  = len(isoNames)

    # DB and molecule name lengths:
    lenDBname = len(DBnames[idb])
    lenMolec  = len(driver[i].molecule)

    # Store length of database name, database name, number of temperatures,
    #  and number of isotopes in TLI file:
    TLIout.write(struct.pack("h%dc"%lenDBname,  lenDBname, *DBnames[idb]))
    # Store the molecule name:
    TLIout.write(struct.pack("h%dc"%lenMolec, lenMolec, *driver[i].molecule))
    # Store the number of temperature samples and isotopes:
    TLIout.write(struct.pack("hh", Ntemp, Niso))
    pt.msg(verb-8, "Database ({:d}/{:d}): '{:s}' ({:s} molecule)\n"
                   "  Number of temperatures: {:d}\n"
                   "  Number of isotopes: {:d}".
                    format(idb+1, Ndb, DBnames[idb], driver[i].molecule,
                           Ntemp, Niso))

    # Write the temperature array:
    TLIout.write(struct.pack("%dd"%Ntemp, *Temp))
    pt.msg(verb-8, "Temperatures: [{:6.1f}, {:6.1f}, ..., {:6.1f}]".
                    format(Temp[0], Temp[1], Temp[Ntemp-1]), 2)

    # For each isotope, write partition function information.
    # Keep a tally of isotopes for multiple databse support:
    for j in np.arange(Niso):
      pt.msg(verb-9, "Isotope ({:d}/{:d}): '{:s}'".
                      format(j+1, Niso, isoNames[j]), 2)

      # Store length of isotope name, isotope name, and isotope mass:
      lenIsoname = len(isoNames[j])
      Iname = str(isoNames[j])
      TLIout.write(struct.pack("h%dc"%lenIsoname, lenIsoname, *Iname))
      TLIout.write(struct.pack("d", iso_mass[j]))
      TLIout.write(struct.pack("d", iso_ratio[j]))

      # Write the partition function per isotope:
      TLIout.write(struct.pack("%dd"%Ntemp, *partDB[j]))
      pt.msg(verb-9, "Mass (u):        {:8.4f}\n"
                     "Isotopic ratio:  {:8.4g}\n"
                     "Part. Function:  [{:.2e}, {:.2e}, ..., {:.2e}]".
                      format(iso_mass[j], iso_ratio[j],
                             partDB[j,0], partDB[j,1], partDB[j,Ntemp-1]), 4)

    # Calculate cumulative number of isotopes per database:
    totIso += Niso
    idb += 1
    acum[idb] = totIso

  # Cumulative number of isotopes:
  pt.msg(verb-5, "Cumulative number of isotopes per DB: {}".format(acum))
  pt.msg(verb, "Done.")

  pt.msg(verb, "\nWriting transition info to TLI file:")
  wnumber = []
  gf      = []
  elow    = []
  isoID   = []
  # Read from file and write the transition info:
  for db in np.arange(Nfiles):
    # Get database index:
    dbname = driver[db].name
    idb = DBnames.index(dbname)

    # Read databases:
    ti = time.time()
    transDB = driver[db].dbread(iwn, fwn, cla.verb, pflist[db])
    tf = time.time()
    pt.msg(verb-3, "Reading time: {:8.3f} seconds".format(tf-ti))
    
    wnumber = np.concatenate((wnumber, transDB[0]))
    gf      = np.concatenate((gf,      transDB[1]))
    elow    = np.concatenate((elow,    transDB[2]))
    isoID   = np.concatenate((isoID,   transDB[3]+acum[idb]))

    pt.msg(verb-8, "Isotpe in-database indices: {}".
                    format(np.unique(transDB[3])))
    pt.msg(verb-8, "Isotpe correlative indices: {}".
                    format(np.unique(transDB[3]+acum[idb])))

  # Total number of transitions:
  nTransitions = np.size(wnumber)

  # Sort the line transitions (by isotope, then by wavenumber):
  ti = time.time()
  isort = sorted(zip(np.arange(nTransitions), isoID, wnumber),
                 key=lambda x:(x[1], x[2]))
  isort = list(zip(*isort)[0])
  tf = time.time()
  pt.msg(verb-3, "Sort time: {:9.7f} seconds".format(tf-ti))
  wnumber = wnumber[isort]
  gf      = gf     [isort]
  elow    = elow   [isort]
  isoID   = isoID  [isort]

  # FINDME: Implement well this:
  if True:
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

  pt.msg(verb-3, "Writing {:d} transition lines.".format(nTransitions))
  ti = time.time()
  transinfo += struct.pack(str(nTransitions)+"d", *list(wnumber))
  transinfo += struct.pack(str(nTransitions)+"h", *list(isoID))
  transinfo += struct.pack(str(nTransitions)+"d", *list(elow))
  transinfo += struct.pack(str(nTransitions)+"d", *list(gf))
  tf = time.time()
  pt.msg(verb-3, "Packing time: {:8.3f} seconds".format(tf-ti))

  ti = time.time()
  TLIout.write(transinfo)
  tf = time.time()
  pt.msg(verb-3, "Writing time: {:8.3f} seconds".format(tf-ti))

  TLIout.close()
  pt.msg(verb, "Done.\n")
  sys.exit(0)

