# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
if sys.version_info.major == 3:
  import configparser
else:
  import ConfigParser as configparser
import numpy as np

__all__ = ["makeTEA"]


pbdir  = os.path.dirname(os.path.realpath(__file__))
TEAdir = os.path.realpath(pbdir + "/../../modules/TEA/")


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
  # Open New Config parser:
  config = configparser.SafeConfigParser()
  config.add_section('TEA')
  config.set("TEA", "maxiter",      str(maxiter))
  config.set("TEA", "save_headers", str(save_headers))
  config.set("TEA", "save_outputs", str(save_outputs))
  config.set("TEA", "doprint",      str(doprint))
  config.set("TEA", "times",        str(times))
  config.set("TEA", "location_TEA", str(location_TEA))
  config.set("TEA", "location_out", str(location_out))
  config.set("TEA", "abun_file",    str(abun_file))

  # Override with input Config parser values:
  if cfile is not None:
    Bconfig = configparser.SafeConfigParser()
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
