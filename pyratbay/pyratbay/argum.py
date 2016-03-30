import sys, os
import ConfigParser, argparse
import numpy as np

from .. import tools as pt
from .. import VERSION   as ver


def parse(wlog):
  """
  Read the command line arguments.
  """

  # Parser to process a configuration file:
  cparser = argparse.ArgumentParser(description=__doc__, add_help=False,
                           formatter_class=argparse.RawDescriptionHelpFormatter)
  # Add config file option:
  cparser.add_argument("-c", "--cfile",
                       help="Configuration filename (string).", metavar="FILE")
  # remaining_argv contains all other command-line-arguments:
  args, remaining_argv = cparser.parse_known_args()

  if args.cfile is None:
    pt.error("No configuration file specified.")
  elif not os.path.isfile(args.cfile):
    pt.error("Configuration file '{:s}' does not exist.".format(args.cfile))

  config = ConfigParser.SafeConfigParser()
  config.read([args.cfile])
  if "pbay" not in config.sections():
    pt.error("Invalid configuration file: '{:s}'.".format(args.cfile))

  defaults = dict(config.items("pbay"))

  # Inherit options from cparser:
  parser = argparse.ArgumentParser(parents=[cparser])

  # Pressure layers options:
  group = parser.add_argument_group("Pressure layer sampling")
  group.add_argument("--nlayers", dest="nlayers",
           help="Number of atmospheric layers [default: %(default)s]",
           action="store", type=int, default=100)
  group.add_argument("--punits",  dest="punits",
           help="Pressure units [default: %(default)s]",
           action="store", type=str, default="---")
  group.add_argument("--ptop",    dest="ptop",
           help="Pressure at the top of the atmosphere (punits) "
                "[default: %(default)s]",
           action="store", type=str, default="1e-8 bar")
  group.add_argument("--pbottom", dest="pbottom",
           help="Pressure at the botom of the atmosphere (punits) "
                "[default: %(default)s]",
           action="store", type=str, default="100 bar")

  # Physical variables of the system:
  group = parser.add_argument_group("Physical parameters")
  group.add_argument("--rstar",   dest="rstar",
           help="Stellar radius [default: %(default)s]",
           action="store", type=str,       default="1 rsun")
  group.add_argument("--tstar",   dest="tstar",
           help="Stellar effective temperature (kelvin) "
                "[default: %(default)s]",
           action="store", type=np.double, default=5600.0)
  group.add_argument("--smaxis",  dest="smaxis",
           help="Planetary semi-major axis (runits) [default: %(default)s]",
           action="store", type=str,       default=None)
  group.add_argument("--tint",    dest="tint",
           help="Planetary internal temperature (kelvin) "
                "[default: %(default)s]",
           action="store", type=np.double, default=100.0)
  group.add_argument("--rplanet", dest="rplanet",
           help="Planetary radius (runits) [default: %(default)s]",
           action="store", type=str,       default="1.0 rearth")

  # Temperature profile options:
  group = parser.add_argument_group("Temperature-profile options")
  group.add_argument("--tmodel", dest="tmodel",
           help="Temperature-profile model [default: %(default)s]",
           action="store", type=str,       default="TCEA")
  group.add_argument("--tparams", dest="tparams",
           help="Temperature-profile parameters [default: %(default)s]",
           action="store", type=pt.parray, default=None)

  # Atmospheric file:
  group = parser.add_argument_group("Atmospheric-file options")
  group.add_argument("--atmfile",    dest="atmfile",
           help="Atmospheric file [default: %(default)s]",
           action="store", type=str,       default=None)
  group.add_argument("--species",    dest="species",
           help="Atmospheric species for the atmospheric file "
                "[default: %(default)s]",
           action="store", type=pt.parray, default=None)
  # Uniform decides uniform-abundances vs TEA-abundances:
  group.add_argument("--uniform",    dest="uniform",
           help="Species mole mixing ratios for uniform-abundance "
                "profiles [default: %(default)s]",
           action="store", type=pt.parray, default=None)
  # VAriables for TEA calculations:
  group.add_argument("--solar",      dest="solar",
           help="Solar composition file (solar atomic composition) "
                "[default: %(default)s]",
           action="store", type=str,       default=None)
  group.add_argument("--xsolar",     dest="xsolar",
           help="Metallicity scaling factor [default: %(default)s]",
           action="store", type=np.double, default=1.0)
  group.add_argument("--atomicfile", dest="atomicfile",
           help="Atomic-composition file (atomic composition) "
                "[default: %(default)s]",
           action="store", type=str,       default=None)
  group.add_argument("--patm",       dest="patm",
           help="Pre-atmospheric file (atomic abundances per layer) "
                "[default: %(default)s]",
           action="store", type=str,       default=None)
  group.add_argument("--elements",   dest="elements",
           help="Atmospheric elements for the pre-atmospheric file "
                "[default: %(default)s]",
           action="store", type=pt.parray, default=None)

  # Output files options:
  group = parser.add_argument_group("Output files")
  group.add_argument("--logfile", dest="logfile",
           help="Log file name.",
           action="store", type=str,       default=None)

  parser.set_defaults(**defaults)
  args = parser.parse_args(remaining_argv)

  # Get logfile:
  if args.logfile is not None:
    log = open(args.logfile, "w")

  # Welcome message:
  pt.msg(1, "{:s}\n"
            "  Python Radiative Transfer in a Bayesian framework (Pyrat Bay).\n"
            "  Version {:d}.{:d}.{:d}.\n"
            "  Copyright (c) 2016 Patricio Cubillos and collaborators.\n"
            "  Pyrat Bay is open-source software under the RR license.\n"
            "{:s}\n\n".format(pt.sep, ver.PBAY_VER, ver.PBAY_MIN,
                                      ver.PBAY_REV, pt.sep), log)

  # Throw warnings for the default input variables:
  if args.punits == "---":
    args.punits = "bar"
    pt.warning("punits input variable defaulted to: '{:s}'.".
               format(args.punits), wlog, log)

  return args, log
