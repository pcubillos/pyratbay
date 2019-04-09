# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ['parse']

import os
import sys
import argparse
if sys.version_info.major == 3:
    import configparser
else:
    import ConfigParser as configparser
import numpy as np
from datetime import date

from .  import tools as pt
from .. import VERSION as ver
from .. import constants as pc

sys.path.append(pc.ROOT + "modules/MCcubed")
import MCcubed.utils as mu


def parse(cfile):
  """
  Read the command line arguments.

  Parameters
  ----------
  cfile: String
      A Pyrat Bay configuration file.

  Returns
  -------
  args: Namespace
      Object storing the attributes defined in this function, with
      the values given in cfile.
  log: Log object
      An MCcubed.utils.Log instance to log screen outputs to file.
  """
  if cfile is None:
      print("No configuration file specified.")
      sys.exit(0)
  elif not os.path.isfile(cfile):
      print("Configuration file '{:s}' does not exist.".format(cfile))
      sys.exit(0)

  config = configparser.ConfigParser()
  config.optionxform = str  # Enable case-sensitive variable names
  config.read([cfile])
  if "pyrat" not in config.sections():
      print("\nInvalid configuration file: '{:s}', no [pyrat] section.".
            format(cfile))
      sys.exit(0)
  defaults = dict(config.items("pyrat"))

  # Use argparse only to define the data types and defaults:
  parser = argparse.ArgumentParser(
               formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("--verb", type=int, default=None)

  pt.addarg("dblist",      parser, pt.parray, None)
  pt.addarg("pflist",      parser, pt.parray, None)
  pt.addarg("dbtype",      parser, pt.parray, None)

  pt.addarg("linedb",      parser, pt.parray, None)
  pt.addarg("tlifile",     parser, pt.parray, None)
  pt.addarg("csfile",      parser, pt.parray, None)
  pt.addarg("molfile",     parser, str,       None)
  pt.addarg("extfile",     parser, str,       None)
  # Spectrum sampling options:
  pt.addarg("wlunits",     parser, str,       None)
  pt.addarg("wllow",       parser, str,       None)
  pt.addarg("wlhigh",      parser, str,       None)
  pt.addarg("wnunits",     parser, str,       None)
  pt.addarg("wnlow",       parser, str,       None)
  pt.addarg("wnhigh",      parser, str,       None)
  pt.addarg("wnstep",      parser, str,       None)
  pt.addarg("wnosamp",     parser, int,       None)
  pt.addarg("resolution",  parser, np.double, None)
  # Atmospheric sampling options:
  pt.addarg("tmodel",      parser, str,       None)
  pt.addarg("tpars",       parser, pt.parray, None)
  pt.addarg("radlow",      parser, str,       None)
  pt.addarg("radhigh",     parser, str,       None)
  pt.addarg("radstep",     parser, str,       None)
  pt.addarg("runits",      parser, str,       None)
  pt.addarg("punits",      parser, str,       None)
  pt.addarg("nlayers",     parser, int,       None)
  pt.addarg('ptop',        parser, str,       None)
  pt.addarg('pbottom',     parser, str,       None)
  # Variables for TEA calculations:
  pt.addarg('atmfile',     parser, str,       None)
  pt.addarg('species',     parser, pt.parray, None)
  pt.addarg('uniform',     parser, pt.parray, None)
  pt.addarg('ptfile',      parser, str,       None)
  pt.addarg('solar',       parser, str,       None)
  pt.addarg('xsolar',      parser, np.double, None)
  pt.addarg('atomicfile',  parser, str,       None)
  pt.addarg('patm',        parser, str,       None)
  pt.addarg('elements',    parser, pt.parray, None)
  # Extinction options:
  pt.addarg("tmin",        parser, np.double, None)
  pt.addarg("tmax",        parser, np.double, None)
  pt.addarg("tstep",       parser, np.double, None)
  pt.addarg("ethresh",     parser, np.double, 1e-15)
  pt.addarg('ncpu',        parser, int,       None)
  # Voigt-profile options:
  pt.addarg("vextent",     parser, np.double, None)
  pt.addarg("Dmin",        parser, np.double, None)
  pt.addarg("Dmax",        parser, np.double, None)
  pt.addarg("nDop",        parser, np.int,    None)
  pt.addarg("Lmin",        parser, np.double, None)
  pt.addarg("Lmax",        parser, np.double, None)
  pt.addarg("nLor",        parser, np.int,    None)
  pt.addarg("DLratio",     parser, np.double, None)
  # Hazes and clouds options:
  pt.addarg("hazes",       parser, pt.parray, None)
  pt.addarg("hpars",       parser, pt.parray, None)
  pt.addarg("rayleigh",    parser, pt.parray, None)
  pt.addarg("rpars",       parser, pt.parray, None)
  pt.addarg("fpatchy",     parser, np.double, None)
  # Alkali opacity options:
  pt.addarg("alkali",      parser, pt.parray, None)
  # Optical depth options:
  pt.addarg("path",        parser, str,       None)
  pt.addarg("maxdepth",    parser, np.double, None)
  pt.addarg("raygrid",     parser, pt.parray, None)
  pt.addarg("quadrature",  parser, int,       None)
  # Data options:
  pt.addarg("runmode",     parser, str,       None)
  pt.addarg("data",        parser, pt.parray, None)
  pt.addarg("uncert",      parser, pt.parray, None)
  pt.addarg("filter",      parser, pt.parray, None)
  # Retrieval options:
  pt.addarg("retflag",     parser, pt.parray, None)
  pt.addarg("bulk",        parser, pt.parray, None)
  pt.addarg("molmodel",    parser, pt.parray, None)
  pt.addarg("molfree",     parser, pt.parray, None)
  pt.addarg("molpars",     parser, pt.parray, None)
  pt.addarg("qcap",        parser, np.double, 0.99)
  pt.addarg("params",      parser, pt.parray, None)
  pt.addarg("stepsize",    parser, pt.parray, None)
  pt.addarg('tlow',        parser, np.double, 0.0)
  pt.addarg("thigh",       parser, np.double, np.inf)
  # Retrieval variables:
  pt.addarg('pmin',        parser, pt.parray, None)
  pt.addarg('pmax',        parser, pt.parray, None)
  pt.addarg('prior',       parser, pt.parray, None)
  pt.addarg('priorlow',    parser, pt.parray, None)
  pt.addarg('priorup',     parser, pt.parray, None)
  pt.addarg('walk',        parser, str,       'snooker')
  pt.addarg('nsamples',    parser, eval,      int(1e5))
  pt.addarg('nchains',     parser, int,       7)
  pt.addarg('burnin',      parser, eval,      0)
  pt.addarg('thinning',    parser, int,       1)
  pt.addarg('grbreak',     parser, np.double, 0.0)
  pt.addarg('grnmin',      parser, eval,      0.5)
  pt.addarg('resume',      parser, str,       False, action='store_true')
  # System physical parameters:
  pt.addarg("starspec",    parser, str,       None)
  pt.addarg("kurucz",      parser, str,       None)
  pt.addarg("marcs",       parser, str,       None)
  pt.addarg("phoenix",     parser, str,       None)
  pt.addarg("rstar",       parser, str,       None)
  pt.addarg("gstar",       parser, np.double, None)
  pt.addarg("tstar",       parser, np.double, None)
  pt.addarg("mstar",       parser, str,       None)
  pt.addarg("rplanet",     parser, str,       None)
  pt.addarg("refpressure", parser, str,       None)
  pt.addarg("mplanet",     parser, str,       None)
  pt.addarg("gplanet",     parser, np.double, None)
  pt.addarg("smaxis",      parser, str,       None)
  pt.addarg("tint",        parser, np.double, None)
  # Output file options:
  pt.addarg("outspec",     parser, str,       None)
  pt.addarg("logfile",     parser, str,       None)
  # Output plotting:
  pt.addarg('logxticks',   parser, pt.parray, None)
  pt.addarg('yran',        parser, pt.parray, None)


  # Set the defaults from the configuration file:
  parser.set_defaults(**defaults)
  # Set values from command line:
  args, unknowns = parser.parse_known_args()
  args.configfile = cfile

  args.logfile = pt.path(args.logfile)
  log = mu.Log(logname=args.logfile, verb=args.verb, width=80,
               append=args.resume)

  # Welcome message:
  log.msg("{:s}\n"
          "  Python Radiative Transfer in a Bayesian framework (Pyrat Bay).\n"
          "  Version {:d}.{:d}.{:d}.\n"
          "  Copyright (c) 2016-{:d} Patricio Cubillos and collaborators.\n"
          "  Pyrat Bay is (temporarily) proprietaty software (see LICENSE).\n"
          "{:s}\n\n".format(log.sep, ver.PBAY_VER, ver.PBAY_MIN,
                            ver.PBAY_REV, date.today().year, log.sep), verb=0)

  log.msg("Read command-line arguments from configuration file: '{:s}'".
          format(cfile))

  return args, log
