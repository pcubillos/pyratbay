# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os, sys
import ConfigParser as cp
import numpy as np

__all__ = ["makeTEA"]


pbdir  = os.path.dirname(os.path.realpath(__file__))
TEAdir = os.path.realpath(pbdir + "/../../modules/TEA/")


def makeTEA(cfile=None, maxiter=300, savefiles=False,
            verb=1, times=False, abun_file=None,
            location_out="./TEA", guess=None, ncpu=1):
  """
  Make a TEA configuration file.

  Parameters
  ----------
  cfile: String
     BART configuration file
  TEAdir: String
     Default TEA directory.
  """
  # Open New Config parser:
  config = cp.SafeConfigParser()
  config.add_section('TEA')
  config.set("TEA", "maxiter",      str(maxiter))
  config.set("TEA", "savefiles",    str(savefiles))
  config.set("TEA", "verb",         str(verb))
  config.set("TEA", "times",        str(times))
  config.set("TEA", "guess",        str(guess))
  config.set("TEA", "ncpu",         str(ncpu))
  config.set("TEA", "location_out", str(location_out))
  config.set("TEA", "abun_file",    str(abun_file))

  # Override with input Config parser values:
  if cfile is not None:
    Bconfig = cp.SafeConfigParser()
    Bconfig.read([cfile])

    keys = ["maxiter", "savefiles", "verb",
            "times", "location_TEA", "abun_file", "ncpu", "location_out"]
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
