# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys, os
import ConfigParser, argparse
import numpy as np

from .. import tools as pt
from .. import VERSION as ver

rootdir = os.path.realpath(os.path.dirname(__file__) + "/../../")

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
  if "pyrat" not in config.sections():
    pt.error("Invalid configuration file: '{:s}', no [pyrat] section.".
             format(args.cfile))

  cfile = args.cfile
  defaults = dict(config.items("pyrat"))

  # Inherit options from cparser:
  parser = argparse.ArgumentParser(parents=[cparser])

  # Pressure-layer options:
  group = parser.add_argument_group("Pressure layer sampling")
  group.add_argument("--nlayers", dest="nlayers",
           help="Number of atmospheric layers [default: 100]",
           action="store", type=int, default=None)
  group.add_argument("--punits",  dest="punits",
           help="Pressure units [default: bar]",
           action="store", type=str, default=None)
  group.add_argument("--ptop",    dest="ptop",
           help="Pressure at the top of the atmosphere (punits) "
                "[default: 1e-8 bar]",
           action="store", type=str, default=None)
  group.add_argument("--pbottom", dest="pbottom",
           help="Pressure at the botom of the atmosphere (punits) "
                "[default: 100 bar]",
           action="store", type=str, default=None)
  group.add_argument("--runmode", dest="runmode",
           help="Run mode flag.  Select from: 'tli', 'pt', 'atmosphere', "
                "'opacity', 'spectrum', 'mcmc'.",
           action="store", type=str, default=None)
  # Physical variables of the system:
  group = parser.add_argument_group("Physical parameters")
  group.add_argument("--radunits",  dest="radunits",
           help="Radius units [default: cm]",
           action="store", type=str,       default=None)
  group.add_argument("--rstar",   dest="rstar",
           help="Stellar radius (default units: radunits).",
           action="store", type=str,       default=None)
  group.add_argument("--tstar",   dest="tstar",
           help="Stellar effective temperature (kelvin).",
           action="store", type=np.double, default=None)
  group.add_argument("--gstar",   dest="gstar",
           help="Stellar surface gravity (cm s-2).",
           action="store", type=np.double, default=None)
  group.add_argument("--smaxis",  dest="smaxis",
           help="Orbital semi-major axis (default units: radunits).",
           action="store", type=str,       default=None)
  group.add_argument("--tint",    dest="tint",
           help="Planetary internal temperature (kelvin) [default: 100].",
           action="store", type=np.double, default=None)
  group.add_argument("--mplanet",    dest="mplanet",
           help="Planetary mass (default units: g).",
           action="store", type=str, default=None)
  group.add_argument("--gplanet",    dest="gplanet",
           help="Planetary surface gravity (cm s-2).",
           action="store", type=np.double, default=None)
  group.add_argument("--rplanet", dest="rplanet",
           help="Planetary radius (radunits) [default: %(default)s]",
           action="store", type=str,       default=None)
  group.add_argument("--starspec",   dest="starspec",
           help="Stellar spectrum model file.",
           action="store", type=str, default=None)
  group.add_argument("--kurucz",   dest="kurucz",
           help="Kurucz stellar model file.",
           action="store", type=str, default=None)

  # Temperature profile options:
  group = parser.add_argument_group("Temperature-profile options")
  group.add_argument("--tmodel", dest="tmodel",
           help="Temperature-profile model [default: %(default)s]",
           action="store", type=str,       default=None)
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
  group.add_argument("--uniform",    dest="uniform",
           help="Species mole mixing ratios for uniform-abundance "
                "profiles [default: %(default)s]",
           action="store", type=pt.parray, default=None)
  # Variables for TEA calculations:
  group.add_argument("--solar",      dest="solar",
           help="Solar composition file (solar atomic composition) "
                "[default: %(default)s]",
           action="store", type=str,       default=None)
  group.add_argument("--xsolar",     dest="xsolar",
           help="Metallicity scaling factor [default: %(default)s]",
           action="store", type=np.double, default=None)
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
  # Extinction-coefficient variables:
  group.add_argument("--extfile",    dest="extfile",
           help="Extinction-coefficient table file.",
           action="store", type=str,       default=None)
  # Retrieval variables:
  group = parser.add_argument_group("Retrieval options")
  group.add_argument("--data",   dest="data",
           help="Transit or eclipse depths.",
           action="store", type=pt.parray, default=None)
  group.add_argument("--uncert",   dest="uncert",
           help="Transit or eclipse depth uncertainties.",
           action="store", type=pt.parray, default=None)
  group.add_argument("--params", dest="params",
           help="Filename or list of initial-guess model-fitting "
                "parameter [required]",
           action="store", type=pt.parray, default=None)
  group.add_argument("--pmin", dest="pmin",
           help="Filename or list of parameter lower boundaries.",
           action="store", type=pt.parray, default=None)
  group.add_argument("--pmax", dest="pmax",
           help="Filename or list of parameter upper boundaries.",
           action="store", type=pt.parray, default=None)
  group.add_argument("--stepsize", dest="stepsize",
           help="Filename or list with proposal jump scale.",
           action="store", type=pt.parray, default=None)
  group.add_argument("--prior", dest="prior",
           help="Filename or list of parameter priors.",
           action="store", type=pt.parray, default=None)
  group.add_argument("--priorlow", dest="priorlow",
           help="Filename or list of parameter prior lower uncertainties.",
           action="store", type=pt.parray, default=None)
  group.add_argument("--priorup", dest="priorup",
           help="Filename or list of parameter prior upper uncertainties.",
           action="store", type=pt.parray, default=None)
  group.add_argument("--thigh",   dest="thigh",
           help="Upper boundary for temperature sampling (kelvin).",
           action="store", type=np.double, default=np.inf)
  group.add_argument("--tlow",   dest="tlow",
           help="Lower boundary for temperature sampling (kelvin).",
           action="store", type=np.double, default=0.0)
  group.add_argument("--walk",     dest="walk",
           help="Random walk algorithm, select from: ['mrw', 'demc', "
                "'snooker']. [default: %(default)s]",
           action="store", type=str,  default="snooker")
  group.add_argument("--nsamples", dest="nsamples",
           help="Number of MCMC samples [default: %(default)s]",
           action="store", type=eval,  default=int(1e5))
  group.add_argument("--nchains", dest="nchains",
           help="Number of chains [default: %(default)s]",
           action="store", type=int,  default=7)
  group.add_argument("--burnin", dest="burnin",
           help="Number of burn-in iterations per chain [default: %(default)s]",
           action="store", type=eval,  default=0)
  group.add_argument("--thinning", dest="thinning",
           help="Chains thinning factor (use every thinning-th iteration) "
                "[default: %(default)s]",
           action="store", type=int,  default=1)
  group.add_argument("--bulk",   dest="bulk",
           help="Bulk-abundance atmospheric species",
           action="store", type=pt.parray, default=None)
  group.add_argument("--molscale",   dest="molscale",
           help="Variable-abundance atmospheric species",
           action="store", type=pt.parray, default=None)
  group.add_argument("--filter",     dest="filter",
           help="Waveband filter filenames.",
           action="store", type=pt.parray, default=None)
  # Output files options:
  group = parser.add_argument_group("Output files")
  group.add_argument("--logfile", dest="logfile",
           help="Log file name.",
           action="store", type=str,       default=None)

  parser.set_defaults(**defaults)
  args = parser.parse_args(remaining_argv)
  args.cfile = cfile

  # Get logfile:
  if args.logfile is not None:
    log = open(args.logfile, "w")
  else:
    log = None

  # Welcome message:
  pt.msg(1, "{:s}\n"
            "  Python Radiative Transfer in a Bayesian framework (Pyrat Bay).\n"
            "  Version {:d}.{:d}.{:d}.\n"
            "  Copyright (c) 2016 Patricio Cubillos and collaborators.\n"
            "  Pyrat Bay is (temporarily) proprietaty software (see LICENSE).\n"
#            "  Pyrat Bay is open-source software under the RR license.\n"
            "{:s}\n\n".format(pt.sep, ver.PBAY_VER, ver.PBAY_MIN,
                                      ver.PBAY_REV, pt.sep), log)

  return args, log


def checkpressure(args, log, wlog):
  """
  Check the input arguments to calculate the pressure profile.
  """
  args.punits = pt.defaultp(args.punits, "bar",
    "punits input variable defaulted to '{:s}'.",   wlog, log)
  args.nlayers = pt.defaultp(args.nlayers, 100,
    "Number of atmospheric-model layers defaulted to {:d}.",         wlog, log)
  args.ptop = pt.defaultp(args.ptop, "1e-8 bar",
    "Atmospheric-model top-pressure boundary defaulted to {:s}.",    wlog, log)
  args.pbottom = pt.defaultp(args.pbottom, "100 bar",
    "Atmospheric-model bottom-pressure boundary defaulted to {:s}.", wlog, log)


def checktemp(args, log, wlog):
  """
  Check the input arguments to calculate the temperature profile.
  """
  if args.tmodel is None:
    pt.error("Undefined temperature model (tmodel).", log)
  if args.tparams is None:
    pt.error("Undefined temperature-model parameters (tparams).", log)

  if args.tmodel == "TCEA":
    if len(args.tparams) != 5:
      pt.error("Wrong number of parameters ({:d}) for the TCEA temperature "
               "model (5).".format(len(args.tparams)), log)
    if args.rstar is None:
      pt.error("Undefined stellar radius (rstar).", log)
    if args.tstar is None:
      pt.error("Undefined stellar temperature (tstar).", log)
    if args.smaxis is None:
      pt.error("Undefined orbital semi-major axis (smaxis).", log)
    if (args.gplanet is None and
        (args.mplanet is None or args.rplanet is None)):
      pt.error("Undefined planetary surface gravity (gplanet).", log)
    args.tint = pt.defaultp(args.tint, "100",
      "Planetary internal temperature defaulted to {:s} K.", wlog, log)
    args.radunits = pt.defaultp(args.radunits, "cm",
                   "radunits input variable defaulted to '{:s}'.", wlog, log)

  elif args.tmodel == "isothermal":
    if len(args.tparams) != 1:
      pt.error("Wrong number of parameters ({:d}) for the isothermal "
               "temperature model (1).".format(len(args.tparams)), log)


def checkatm(args, log, wlog):
  """
  Check the input arguments to calculate the atmospheric model.
  """
  if args.runmode != "atmosphere" and args.atmfile is None:
    pt.error("Undefined atmospheric file (atmfile).", log)
  if args.species is None:
    pt.error("Undefined atmospheric species list (species).", log)
  args.punits = pt.defaultp(args.punits, "bar",
    "punits input variable defaulted to '{:s}'.",   wlog, log)

  # Uniform-abundances profile:
  if args.uniform is not None:
    if len(args.uniform) != len(args.species):
      pt.error("Number of uniform abundances ({:d}) does not match the "
               "number of species ({:d}).".
                format(len(args.uniform), len(args.species)), log)
    return
  else:  # TEA abundances:
    if args.elements is None:
      pt.error("Undefined atmospheric atomic-composition list (elements).", log)
    args.solar = pt.defaultp(args.solar, rootdir+"/inputs/AsplundEtal2009.txt",
      "Solar-abundances file defaulted to '{:s}'.", wlog, log)
    args.atomicfile = pt.defaultp(args.atomicfile, "./atomic.tea",
      "Atomic-composition file defaulted to '{:s}'.", wlog, log)
    args.patm = pt.defaultp(args.patm, "./preatm.tea",
      "Pre-atmospheric file defaulted to '{:s}'.", wlog, log)
    args.xsolar = pt.defaultp(args.xsolar, 1.0,
      "Solar-metallicity scaling factor defaulted to {:.2f}.", wlog, log)


def checkinputs(args, log, wlog):
  """
  Check that the input values (args) make sense.
  """

  if args.runmode is None:
    pt.error("Undefined run mode (runmode), select from: 'tli', 'pt', "
             "'atmosphere', 'opacity', 'spectrum', 'mcmc'.", log)
  if args.runmode not in ['tli', 'pt', 'atmosphere', 'opacity',
                          'spectrum', 'mcmc']:
    pt.error("Invalid runmode ({:s}).  Must select one from: 'tli', 'pt', "
        "'atmosphere', 'opacity', 'spectrum', 'mcmc'.".format(args.runmode))


  # Stellar model:
  if args.starspec is not None:
    # Check file exists
    pass
  elif args.kurucz is not None:
    # Check file exists
    # Check gstar exists
    pass
  else:
    #pt.error("Stellar spectrum model was not specified.")
    pass


