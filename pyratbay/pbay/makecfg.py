# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["makeTEA"]

import os, sys
import argparse, ConfigParser
import numpy as np
import scipy.constants as sc

#import reader as rd
#import constants as c

pbdir  = os.path.dirname(os.path.realpath(__file__))
TEAdir = os.path.realpath(pbdir + "/../../modules/TEA/")

def makeTransit(cfile, tepfile, shareOpacity):
  """
  Make the transit configuration file.

  Parameters:
  -----------
  cfile: String
     BART configuration file.
  tepfile: String
     A TEP file.
  """

  # Known transit arguments:
  known_args = ["radlow",  "radhigh",  "raddelt",  "radfct",
                "wllow",   "wlhigh",   "wlfct",
                "wnlow",   "wnhigh",   "wndelt",  "wnfct", "wnosamp",
                "tlow",    "thigh",    "tempdelt",
                "allowq",  "nwidth",
                "opacityfile",
                "molfile",
                "toomuch", "tauiso", "outtau", "taulevel", "modlevel",
                "starrad", "transparent",
                "solution", "raygrid",
                "gsurf", "refpress", "refradius",
                "cloudrad", "cloudfct", "cloudext",
                "verb",
                "savefiles",
                "outtoomuch", "outsample", "outmod", "outflux", "outintens",
                "linedb", "csfile", "orbpars", "orbparsfct"]

  # Name of the configuration-file section:
  section = "MCMC"
  # Read BART configuration file:
  Bconfig = ConfigParser.SafeConfigParser()
  Bconfig.read([cfile])

  # Keyword names of the arguments in the BART configuration file:
  args = Bconfig.options(section)

  # transit configuration filename:
  tcfile = open(Bconfig.get(section, "tconfig"), "w")

  # Keyword for the atmospheric file is different in transit:
  tcfile.write("atm {:s}\n".format(Bconfig.get(section, "atmfile")))

  # Default molfile path:
  if "molfile" not in args:
    tcfile.write("molfile {:s}\n".format(
      os.path.realpath(filedir + "/../modules/transit/inputs/molecules.dat")))

  # Calculate gsurf and refradius from the tepfile:
  tep = rd.File(tepfile)
  mplanet = float(tep.getvalue('Mp')[0]) * c.Mjup
  rplanet = float(tep.getvalue('Rp')[0]) * c.Rjup

  # Planetary radius reference level in km:
  Bconfig.set(section, "refradius", "{:.2f}".format(rplanet * 1e-3))
  # Planetary surface gravity in cm/s2:
  Bconfig.set(section, "gsurf", "{:.1f}".format(100*sc.G*mplanet/rplanet**2))
  # Add these keywords:
  args = np.union1d(args, ["refradius", "gsurf"])

  # Reformat the cross-section inputs for Transit:
  cs = Bconfig.get(section, "csfile")
  Bconfig.set(section, "csfile", ",".join(cs.split()))

  # Print the known arguments to file:
  for key in np.intersect1d(args, known_args):
    # FINDME: Why am I splitting?
    values = Bconfig.get(section, key).split("\n")
    for val in values:
      tcfile.write("{:s} {:s}\n".format(key, val))

  if shareOpacity:
    tcfile.write("shareOpacity \n")
  tcfile.close()


def makeMCMC(cfile, MCMC_cfile, logfile):
  """
  Reformat configuration file to remove relative paths.  This output 
  configuration file is used by the BART's MCMC program.

  Parameters:
  -----------
  cfile: String
     BART configuration file.
  MCMC_cfile: String
     Reformated configuration file.
  logfile: String
     Default logfile argument if not specified in cfile.
  """

  # Name of the configuration-file section:
  section = "MCMC"
  # Open input BART configuration file:
  Bconfig = ConfigParser.SafeConfigParser()
  Bconfig.optionxform = str
  Bconfig.read([cfile])
  # Keyword names of the arguments in the BART configuration file:
  args = Bconfig.options(section)

  # Known arguments that may have a path:
  input_args = ["tep_name", "kurucz", "molfile", "filter", "linedb",
                "csfile", "loc_dir"]
  output_args = ["tconfig",    "atmfile",   "opacityfile", "press_file",
                 "abun_basic", "abun_file", "preatm_file", "outflux",
                 "outmod",     "savemodel", "logfile"]

  # Set default logfile:
  if "logfile" not in args:
    Bconfig.set(section, "logfile", logfile)
    args.append("logfile")

  # Inputs should always exist:
  for arg in np.intersect1d(args, input_args):
    # Split multiple values (if more than one), get full path, join back:
    values = Bconfig.get(section, arg).split("\n")
    for v in np.arange(len(values)):
      values[v] = os.path.realpath(values[v])
    Bconfig.set(section, arg, "\n".join(values))

  # Outputs are stored/copied into loc_dir:
  for arg in np.intersect1d(args, output_args):
    Bconfig.set(section, arg,
                Bconfig.get(section, "loc_dir") + "/" +
                os.path.basename(Bconfig.get(section, arg)))

  # Add mpi:
  Bconfig.set(section, "mpi", "True")

  # Add func:
  Bconfig.set(section, "func", "hack BARTfunc {:s}".format(filedir))

  # Params is a special case:
  params = Bconfig.get(section, "params")
  # It may or not be a file path:
  if os.path.isfile(params):
    Bconfig.set(section, "params", os.path.realpath(params))

  # Write the configuration file for use by MC3:
  with open(MCMC_cfile, 'w') as configfile:
    Bconfig.write(configfile)


def makeTEA(cfile=None, maxiter=100, save_headers=False, save_outputs=False,
            doprint=False, times=False, location_TEA=TEAdir, abun_file=None,
            location_out="./TEA"):
  """
  Make a TEA configuration file.

  Parameters
  ----------
  cfile: String
     BART configuration file
  TEAdir: String
     Default TEA directory.
  """
  # Open New ConfigParser:
  config = ConfigParser.SafeConfigParser()
  config.add_section('TEA')
  config.set("TEA", "maxiter",      str(maxiter))
  config.set("TEA", "save_headers", str(save_headers))
  config.set("TEA", "save_outputs", str(save_outputs))
  config.set("TEA", "doprint",      str(doprint))
  config.set("TEA", "times",        str(times))
  config.set("TEA", "location_TEA", str(location_TEA))
  config.set("TEA", "location_out", str(location_out))
  config.set("TEA", "abun_file",    str(abun_file))

  # Override with input ConfigParser values:
  if cfile is not None:
    Bconfig = ConfigParser.SafeConfigParser()
    Bconfig.read([cfile])

    keys = ["maxiter", "save_headers", "save_outputs", "doprint",
            "times", "location_TEA", "abun_file", "location_out"]
    # Set TEA default arguments:
    for i in np.arange(len(keys)):
      if Bconfig.has_option("PBAY", keys[i]):
        config.set("TEA", keys[i], Bconfig.get("PBAY", keys[i]))

  # For completion:
  config.add_section('PRE-ATM')
  config.set("PRE-ATM", "PT_file",        "None")
  config.set("PRE-ATM", "pre_atm_name",   "None")
  config.set("PRE-ATM", "input_elem",     "None")
  config.set("PRE-ATM", "output_species", "None")

  # Write TEA configuration file:
  with open("TEA.cfg", 'w') as configfile:
    config.write(configfile)
