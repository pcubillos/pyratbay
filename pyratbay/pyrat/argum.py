# Copyright (c) 2016-2017 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys, os
import argparse, ConfigParser
import numpy as np
import scipy.constants   as sc
import scipy.interpolate as si
import scipy.special     as ss
import multiprocessing   as mpr

from .. import tools      as pt
from .. import constants  as pc
from .. import wine       as w
from .. import starspec   as sspec
from .. import atmosphere as atm
from .. import VERSION    as ver

from .  import haze      as hz
from .  import alkali    as al

rootdir = os.path.realpath(os.path.dirname(__file__) + "/../../")
sys.path.append(rootdir + "/pyratbay/lib/")
import pt as PT

def parse(pyrat, log=None):
  """
  Parse the command-line arguments into the pyrat object

  Parameters
  ----------
  pyrat: Object
     A Pyrat instance where to store the CLA.
  """

  # Parse configuration file:
  cparser = argparse.ArgumentParser(description=__doc__, add_help=False,
                           formatter_class=argparse.RawDescriptionHelpFormatter)
  # Add config file option:
  cparser.add_argument("-c", "--configfile",
                       help="Specify config file", metavar="FILE")
  cparser.add_argument("-v",  "--verb",  dest="verb",
                       help="Verbosity level [default: %(default)s]",
                       action="store", type=int, default=2)
  # remaining_argv contains all other command-line-arguments:
  args, remaining_argv = cparser.parse_known_args()

  # Get parameters from configuration file (if exists):
  cfile = args.configfile # The configuration file
  #if cfile is None:
  #  pt.exit(message="Undefined configuration file.")
  if cfile is not None and not os.path.isfile(cfile):
    pt.error("Configuration file: '{:s}' not found.".format(cfile))
  if cfile:
    config = ConfigParser.SafeConfigParser()
    config.optionxform = str  # Enable case-sensitive variable names
    config.read([cfile])
    defaults = dict(config.items("pyrat"))
  else:
    defaults = {}

  # Inherit options from cparser:
  parser = argparse.ArgumentParser(parents=[cparser])  #, add_help=False) ??
  # Process pyrat Options:
  group = parser.add_argument_group("Input Files Options")
  pt.addarg("atmfile",     group, str,       None,
      "Atmospheric file [default: %(default)s]")
  pt.addarg("linedb",      group, pt.parray, None,
      "Line database files [default: %(default)s]")
  pt.addarg("csfile",      group, pt.parray, None,
      "Cross-section files [default: %(default)s]")
  pt.addarg("molfile",     group, str,       None,
      "Molecular info file [default: 'pyrat/inputs/molecules.dat']")
  pt.addarg("extfile",     group, str,       None,
      "Extinction-coefficient table file [default: %(default)s]")
  # Spectrum sampling options:
  group = parser.add_argument_group("Spectrum Sampling Options")
  pt.addarg("wlunits",     group, str,       None,
      "Wavelength (input) units [default: um]")
  pt.addarg("wllow",       group, str,       None,
      "Wavelength low boundary [default: %(default)s]")
  pt.addarg("wlhigh",      group, str,       None,
      "Wavelength high boundary [default: %(default)s]")
  pt.addarg("wnunits",     group, str,       None,
      "Wavenumber (input) inverse units [default: cm]")
  pt.addarg("wnlow",       group, str,       None,
      "Wavenumber low boundary [default: %(default)s]")
  pt.addarg("wnhigh",      group, str,       None,
      "Wavenumber high boundary [default: %(default)s]")
  pt.addarg("wnstep",      group, str,       None,
      "Wavenumber sampling step [default: 1.0 cm]")
  pt.addarg("wnosamp",     group, int,       None,
      "Wavenumber oversampling factor [default: 2160]")
  # Atmospheric sampling options:
  group = parser.add_argument_group("Atmosphere Sampling Options")
  pt.addarg("radlow",      group, str,       None,
      "Atmospheric radius low boundary [default: Atmospheric file value]")
  pt.addarg("radhigh",     group, str,       None,
      "Atmospheric radius high boundary [default: Atmospheric file value]")
  pt.addarg("radstep",     group, str,       None,
      "Atmospheric radius sampling step [default: Atmospheric file value]")
  pt.addarg("radunits",    group, str,       None,
      "Radius (user) units [default: km]")
  pt.addarg("plow",        group, str,       None,
      "Atmospheric pressure low boundary (overrides radius high boundary) "
      "[default: %(default)s]")
  pt.addarg("phigh",       group, str,       None,
      "Atmospheric pressure high boundary (overrides radius low boundary) "
      "[default: %(default)s]")
  pt.addarg("nlayers",     group, np.int,    None,
      "Number of atmospheric layers [default: %(default)s]")
  pt.addarg("punits",      group, str,       None,
      "Pressure (user) units [default: bar]")
  # Extinction options:
  group = parser.add_argument_group("Extinction Calculations Options")
  pt.addarg("tmin",        group, np.double, None,
      "Minimum extinction-coefficient grid temperature.")
  pt.addarg("tmax",        group, np.double, None,
      "Maximum extinction-coefficient grid temperature.")
  pt.addarg("tstep",       group, np.double, None,
      "Temperature sample step interval in Kelvin [default: 100]")
  pt.addarg("ethresh",     group, np.double, 1e-9,
      "Extinction-coefficient threshold [default: %(default)s]")
  pt.addarg("nproc",       group, int,       1,
      "Number of processors [default: %(default)s]")
  # Voigt-profile options:
  group = parser.add_argument_group("Voigt-profile  Options")
  pt.addarg("vextent",     group, np.double, None,
      "Extent of Voigt profile in number of Voigt widths [default: 20]")
  pt.addarg("Dmin",        group, np.double, None,
      "Minimum Doppler-width to sample in cm-1 [default: 1.0e-03]")
  pt.addarg("Dmax",        group, np.double, None,
      "Maximum Doppler-width to sample in cm-1 [default: 0.25]")
  pt.addarg("nDop",        group, np.int,    None,
      "Number of Doppler-width samples [default: 40]")
  pt.addarg("Lmin",        group, np.double, None,
      "Minimum Lorentz width to sample in cm-1 [default: 1.0e-04]")
  pt.addarg("Lmax",        group, np.double, None,
      "Maximum Lorentz width to sample in cm-1 [default: 10.0]")
  pt.addarg("nLor",        group, np.int,    None,
      "Number of Lorentz-width samples [default: 40]")
  pt.addarg("DLratio",     group, np.double, None,
      "Minimum Doppler/Lorentz-width ratio to re-calculate a Voigt profile "
      "[default: 0.1]")
  # Hazes and clouds options:
  group = parser.add_argument_group("Hazes and Clouds Options")
  pt.addarg("hazes",       group, pt.parray, None,
      "Haze models [default: %(default)s].")
  pt.addarg("hpars",       group, pt.parray, None,
      "Haze model fitting parameters.")
  # Alkali opacity options:
  group = parser.add_argument_group("Alkali Options")
  pt.addarg("alkali",      group, pt.parray, None,
      "Alkali absorption models [default: %(default)s].")
  # Optical depth options:
  group = parser.add_argument_group("Optical Depth Options")
  pt.addarg("path",        group, str,       None,
      "Observing geometry. Select between: 'transit' or 'eclipse'.")
  pt.addarg("maxdepth",    group, np.double, None,
      "Maximum optical depth to calculate [default: 10]")
  pt.addarg("raygrid",     group, pt.parray, None,
      "Incident angles over day-side hemisphere for intensity integration."
      "Values in degrees between 0 and 90 [default: 0, 20, 40, 60, 80]")
  pt.addarg("quadrature",  group, int,       None,
      "Polynomial degree for quadrature-integration over day-side hemisphere.")
  # Data options:
  group = parser.add_argument_group("Data options")
  pt.addarg("runmode",     group, str,       None,
      "Run mode flag.  Select from: tli, pt, atmosphere, spectrum, "
      "opacity, or mcmc.")
  pt.addarg("data",        group, pt.parray, None,
      "Transit or eclipse depth uncertainties.")
  pt.addarg("uncert",      group, pt.parray, None,
      "Transit or eclipse depth uncertainties.")
  pt.addarg("filter",      group, pt.parray, None,
      "Waveband filter filenames.")
  # Retrieval options:
  group = parser.add_argument_group("Retrieval options")
  pt.addarg("bulk",        group, pt.parray, None,
      "Bulk-abundance atmospheric species.")
  pt.addarg("molscale",    group, pt.parray, None,
      "Variable-abundance atmospheric species.")
  pt.addarg("tmodel",      group, str,       None,
      "Temperature-profile model name.  Select from: isothermal or TCEA.")
  pt.addarg("params",      group, pt.parray, None,
      "Initial-guess for retrieval model-fitting parameter.")
  pt.addarg("stepsize",    group, pt.parray, None,
      "Stepsize for retrieval model-fitting parameter.")
  pt.addarg("tlow",        group, np.double, None,
      "Minimum valid temperature.")
  pt.addarg("thigh",        group, np.double, None,
      "Maximum valid temperature.")
  # System physical parameters:
  group = parser.add_argument_group("System physical variables")
  pt.addarg("starspec",    group, str,       None,
      "Stellar-spectrum model filename.")
  pt.addarg("kurucz",      group, str,       None,
      "Kurucz stellar-spectrum filename.")
  pt.addarg("marcs",       group, str,       None,
      "MARCS stellar-spectrum filename.")
  pt.addarg("phoenix",     group, str,       None,
      "PHOENIX stellar-spectrum filename.")
  pt.addarg("rstar",       group, str,       None,
      "Stellar radius (radunits).")
  pt.addarg("gstar",       group, np.double, None,
      "Stellar surface gravity (cm s-2).")
  pt.addarg("tstar",       group, np.double, None,
      "Stellar effective temperature (kelvin).")
  pt.addarg("mstar",       group, str,       None,
      "Stellar mass (default units: gram)")
  pt.addarg("rplanet",     group, str,       None,
      "Planetary radius (in radunits)")
  pt.addarg("refpressure", group, str,       None,
      "Pressure reference level corresponding to rplanet (in punits).")
  pt.addarg("mplanet",     group, str,       None,
      "Planetary mass (default units: gram)")
  pt.addarg("gplanet",     group, np.double, None,
      "Planetaty surface gravity (cm s-2).")
  pt.addarg("smaxis",     group, str,       None,
      "Orbital semi-major axis (default in radunits).")
  pt.addarg("tint",       group, np.double, None,
      "Planetary internal temperature (kelvin) [default: 100].")
  # Output file options:
  group = parser.add_argument_group("Output File's Options")
  pt.addarg("outspec",     group, str,       None,
      "Output spectrum file [default: 'outspec.dat']")
  pt.addarg("outsample",   group, str,       None,
      "Output samplings file [default: %(default)s]")
  pt.addarg("outmaxdepth", group, str,       None,
      "Filename to store the radius at maxdepth (per wavelength) "
      "[default: %(default)s]")
  pt.addarg("logfile",     group, str,       None,
      "Screen-output log filename.")

  # Set the defaults from the configuration file:
  parser.set_defaults(**defaults)
  # Set values from command line:
  user, unknown = parser.parse_known_args(remaining_argv)

  # Put user arguments into pyrat input:
  pyrat.inputs.configfile = args.configfile
  pyrat.inputs.runmode    = user.runmode
  pyrat.inputs.verb       = user.verb
  pyrat.inputs.nproc      = user.nproc
  # Input file:
  pyrat.inputs.atmfile    = user.atmfile
  pyrat.inputs.linedb     = user.linedb
  pyrat.inputs.csfile     = user.csfile
  pyrat.inputs.molfile    = user.molfile
  pyrat.inputs.extfile    = user.extfile
  # Wavelength:
  pyrat.inputs.wlunits    = user.wlunits
  pyrat.inputs.wllow      = user.wllow
  pyrat.inputs.wlhigh     = user.wlhigh
  # Wavenumber:
  pyrat.inputs.wnunits    = user.wnunits
  pyrat.inputs.wnlow      = user.wnlow
  pyrat.inputs.wnhigh     = user.wnhigh
  pyrat.inputs.wnstep     = user.wnstep
  pyrat.inputs.wnosamp    = user.wnosamp
  # Atmospheric radius:
  pyrat.inputs.radlow     = user.radlow
  pyrat.inputs.radhigh    = user.radhigh
  pyrat.inputs.radstep    = user.radstep
  pyrat.inputs.radunits   = user.radunits
  # Atmospheric pressure:
  pyrat.inputs.plow       = user.plow
  pyrat.inputs.phigh      = user.phigh
  pyrat.inputs.punits     = user.punits
  pyrat.inputs.nlayers    = user.nlayers
  # Hydrostatic-equilibrium base-level variables:
  pyrat.inputs.refpressure = user.refpressure
  # Extinction:
  pyrat.inputs.ethresh = user.ethresh
  pyrat.inputs.tmin    = user.tmin
  pyrat.inputs.tmax    = user.tmax
  pyrat.inputs.tstep   = user.tstep
  # Voigt-profile:
  pyrat.inputs.vextent    = user.vextent
  pyrat.inputs.Dmin       = user.Dmin
  pyrat.inputs.Dmax       = user.Dmax
  pyrat.inputs.nDop       = user.nDop
  pyrat.inputs.Lmin       = user.Lmin
  pyrat.inputs.Lmax       = user.Lmax
  pyrat.inputs.nLor       = user.nLor
  pyrat.inputs.DLratio    = user.DLratio
  # Hazes and clouds:
  pyrat.inputs.hazes      = user.hazes
  pyrat.inputs.hpars      = user.hpars
  # Alkali compounds:
  pyrat.inputs.alkali     = user.alkali
  # Optical depth:
  pyrat.inputs.path       = user.path
  pyrat.inputs.maxdepth   = user.maxdepth
  pyrat.inputs.raygrid    = user.raygrid
  pyrat.inputs.quadrature = user.quadrature
  # System physical variables:
  pyrat.inputs.rstar      = user.rstar
  pyrat.inputs.gstar      = user.gstar
  pyrat.inputs.tstar      = user.tstar
  pyrat.inputs.mstar      = user.mstar
  pyrat.inputs.rplanet    = user.rplanet
  pyrat.inputs.mplanet    = user.mplanet
  pyrat.inputs.gplanet    = user.gplanet
  pyrat.inputs.smaxis     = user.smaxis
  pyrat.inputs.tint       = user.tint
  pyrat.inputs.starspec   = user.starspec
  pyrat.inputs.kurucz     = user.kurucz
  pyrat.inputs.marcs      = user.marcs
  pyrat.inputs.phoenix    = user.phoenix
  # Observing variables:
  pyrat.inputs.data     = user.data
  pyrat.inputs.uncert   = user.uncert
  pyrat.inputs.filter   = user.filter
  # Retrieval variables:
  pyrat.inputs.bulk     = user.bulk
  pyrat.inputs.molscale = user.molscale
  pyrat.inputs.tmodel   = user.tmodel
  pyrat.inputs.params   = user.params
  pyrat.inputs.stepsize = user.stepsize
  pyrat.inputs.tlow     = user.tlow
  pyrat.inputs.thigh    = user.thigh
  # Output files:
  pyrat.inputs.outspec     = user.outspec
  pyrat.inputs.outsample   = user.outsample
  pyrat.inputs.outmaxdepth = user.outmaxdepth
  pyrat.inputs.logfile     = user.logfile

  # Open the Pyrat log file if requested:
  if   log is not None:  # Take pre-existing log
    pyrat.log = log
    pyrat.logfile = os.path.realpath(log.name)
  elif pyrat.inputs.logfile is not None:  # Start new log
    pyrat.logfile = os.path.realpath(pyrat.inputs.logfile)
    pyrat.log = open(pyrat.logfile, "w")

  # Welcome message:
  pt.msg(pyrat.inputs.verb-2,
         "{:s}\n  Python Radiative Transfer (PyRaT).\n"
         "  Version {:d}.{:d}.{:d}.\n"
         "  Copyright (c) 2016-2017 Patricio Cubillos and collaborators.\n"
         "  Pyrat is (temporarily) proprietaty software (see LICENSE).\n"
         "{:s}\n\n".format(pt.sep, ver.PYRAT_VER, ver.PYRAT_MIN,
                                   ver.PYRAT_REV, pt.sep), pyrat.log)

  pt.msg(pyrat.inputs.verb-3, "Read command-line arguments from "
         "configuration file: '{:s}'".format(cfile), pyrat.log)


def checkinputs(pyrat):
  """
  Check that user input arguments make sense.
  """
  # Shortcuts:
  inputs = pyrat.inputs
  phy    = pyrat.phy

  # Path to source parent's folder:
  pyratdir = os.path.dirname(os.path.realpath(__file__))

  # Verbose level:
  pyrat.verb = np.amax([0, inputs.verb])

  # Pyrat runmode:
  pyrat.runmode = inputs.runmode
  if pyrat.runmode is None:
    pt.warning(pyrat.verb-2, "Defaulted Pyrat's runmode to: spectrum.",
               pyrat.log, pyrat.wlog)
  if pyrat.runmode not in ['tli', 'pt', 'atmosphere', 'opacity', 'spectrum',
                           'mcmc']:
    pt.error("Invalid runmode ({:s}).  Select from: tli, pt, atmosphere, "
             "spectrum, opacity, mcmc.".format(pyrat.runmode), pyrat.log)

  # Check that input files exist:
  if inputs.atmfile is None:
    pt.error("Undefined atmospheric file (atmfile).", pyrat.log)
  elif not os.path.isfile(inputs.atmfile):
    pt.error("atmfile: '{:s}' does not exist.".format(inputs.atmfile),
             pyrat.log)
  pyrat.atmfile = inputs.atmfile

  if inputs.linedb is not None:
    for linedb in inputs.linedb:
      if not os.path.isfile(linedb):
        pt.error("linedb file: '{:s}' does not exist.".format(linedb),
                 pyrat.log)
  pyrat.linedb = pyrat.inputs.linedb

  if inputs.csfile is not None:
    for cs in pyrat.inputs.csfile:
      if not os.path.isfile(cs):
        pt.error("Cross-section file: '{:s}' does not exist.".format(cs),
                 pyrat.log)
  pyrat.cs.files = pyrat.inputs.csfile

  if inputs.molfile is None: # Set default
    inputs.molfile = os.path.realpath(pyratdir + "/../../inputs/molecules.dat")
  if not os.path.isfile(inputs.molfile):
    pt.error("Molecular-data file: '{:s}' does not exist.".
             format(inputs.molfile), pyrat.log)
  pyrat.molfile = os.path.realpath(inputs.molfile)

  if inputs.extfile is not None:
    if not os.path.exists(os.path.realpath(os.path.dirname(inputs.extfile))):
      pt.error("Directory for extinction-coefficient file '{:s}' does "
               "not exist.".format(inputs.extfile), pyrat.log)
    pyrat.ex.extfile = os.path.realpath(inputs.extfile)

  # Check spectrum arguments:
  pyrat.spec.wnunits = pt.defaultp(inputs.wnunits, "cm",
         "wnunits input variable defaulted to '{:s}'.", pyrat.wlog, pyrat.log)
  pyrat.spec.wlunits = pt.defaultp(inputs.wlunits, "um",
         "wlunits input variable defaulted to '{:s}'.", pyrat.wlog, pyrat.log)

  pyrat.spec.wllow = pt.getparam(inputs.wllow, pyrat.spec.wlunits)
  isgreater(pyrat.spec.wllow, "um", 0, False,
            "Low wavelength boundary ({:.2e} um) must be >= 0.", pyrat.log)

  pyrat.spec.wlhigh = pt.getparam(inputs.wlhigh, pyrat.spec.wlunits)
  isgreater(pyrat.spec.wlhigh, "um", 0, True,
            "High wavelength boundary ({:.2e} um) must be >= 0.", pyrat.log)

  # Wavenumber must be taken care differently (take inverse of units):
  if inputs.wnlow is not None:
    if len(inputs.wnlow.split()) == 2:
      wnunits = inputs.wnlow.split()[1]
    else:
      wnunits = pyrat.spec.wnunits
    wnlow = float(inputs.wnlow.split()[0])
    if wnlow < 0.0:
      pt.error("Low wavenumber boundary ({:.2e} {:s}-1) must be >= 0.".
               format(wnlow, wnunits), pyrat.log)
    pyrat.spec.wnlow = wnlow / pt.u(wnunits)

  if   inputs.wnhigh is not None:
    if len(inputs.wnhigh.split()) == 2:
      wnunits = inputs.wnhigh.split()[1]
    else:
      wnunits = pyrat.spec.wnunits
    wnhigh = float(inputs.wnhigh.split()[0])
    if wnhigh <= 0.0:
      pt.error("High wavenumber boundary ({:.2e} {:s}-1) must be > 0.".
               format(wnhigh, wnunits), pyrat.log)
    pyrat.spec.wnhigh = wnhigh / pt.u(wnunits)

  wnstep = pt.defaultp(inputs.wnstep, "1.0 cm",
     "Input wavenumber sampling step (wnstep) defaulted to '{:s}'.",
     pyrat.wlog, pyrat.log)
  if len(wnstep.split()) == 2:
    wnunits = wnstep.split()[1]
  else:
    wnunits = pyrat.spec.wnunits
  wnstep = float(wnstep.split()[0])
  if wnstep <= 0:
    pt.error("Wavenumber sampling step ({:.2e} {:s}-1) must be be > 0.".
             format(wnstep, wnunits), pyrat.log)
  pyrat.spec.wnstep = wnstep / pt.u(wnunits)

  pyrat.spec.wnosamp = pt.defaultp(inputs.wnosamp, 2160,
     "Input wavenumber oversampling factor (wnosamp) defaulted to {:d}.",
     pyrat.wlog, pyrat.log)
  isgreater(pyrat.spec.wnosamp, "none", 1, False,
            "Wavenumber oversampling factor ({:d}) must be >= 1.", pyrat.log)

  # Check atmospheric layers arguments:
  pyrat.punits = pt.defaultp(inputs.punits, "bar",
     "Input pressure units (punits) defaulted to '{:s}'.", pyrat.wlog,pyrat.log)
  pyrat.radunits = pt.defaultp(inputs.radunits, "km",
     "Input radius units (punits) defaulted to '{:s}'.", pyrat.wlog, pyrat.log)

  # Pressure boundaries:
  pyrat.phigh = pt.getparam(inputs.phigh, pyrat.punits)
  isgreater(pyrat.phigh, "bar", 0, True,
            "High atm pressure boundary ({:.2e} bar) must be > 0.0", pyrat.log)
  pyrat.plow  = pt.getparam(inputs.plow,    pyrat.punits)
  isgreater(pyrat.plow, "bar",  0, True,
            "Low atm pressure boundary ({:.2e} bar) must be > 0.0", pyrat.log)
  # Radius boundaries:
  pyrat.radlow  = pt.getparam(inputs.radlow,  pyrat.radunits)
  isgreater(pyrat.radlow, "cm", 0, False,
            "Low atm radius boundary ({:.2e} cm) must be >= 0.0", pyrat.log)
  pyrat.radhigh = pt.getparam(inputs.radhigh, pyrat.radunits)
  isgreater(pyrat.radhigh, "cm", 0, True,
            "High atm radius boundary ({:.2e} cm) must be > 0.0", pyrat.log)
  pyrat.radstep = pt.getparam(inputs.radstep, pyrat.radunits)
  isgreater(pyrat.radstep, "cm", 0, True,
            "Radius step size ({:.2f} cm) must be > 0.", pyrat.log)
  # System physical parameters:
  phy.rplanet = pt.getparam(inputs.rplanet, pyrat.radunits)
  isgreater(phy.rplanet, "cm",   0, True,
            "Planetary radius ({:.3e} cm) must be > 0.", pyrat.log)

  pyrat.refpressure = pt.getparam(inputs.refpressure, pyrat.punits)
  isgreater(pyrat.refpressure, "bar", 0, True,
      "Planetary reference pressure level ({:8g} bar) must be > 0.", pyrat.log)

  phy.gplanet  = pt.getparam(inputs.gplanet,  "none")
  isgreater(phy.gplanet, "none", 0, True,
            "Planetary surface gravity ({:.2f} cm s-2) must be > 0.", pyrat.log)

  phy.mplanet  = pt.getparam(inputs.mplanet,  "gram")
  isgreater(phy.mplanet, "mearth", 0, True,
            "Planetary mass ({:.2f} Mearth) must be > 0.", pyrat.log)

  # Check planetary surface gravity:
  if phy.mplanet is not None:
    pyrat.hydrom = True  # Use mass value for hydrostatic equilibrium
    if phy.rplanet is None and phy.gplanet is not None:
      phy.rplanet = np.sqrt(pc.G * phy.mplanet / phy.gplanet)
    if phy.rplanet is not None:
      gplanet = pc.G * phy.mplanet / phy.rplanet**2
      if phy.gplanet is None:
        phy.gplanet = gplanet
      elif np.abs(gplanet-phy.gplanet)/phy.gplanet > 0.05:
        pt.error("Both mplanet and gplanet were provided, but values are "
          "inconsistent (>5%): g(mplanet) = {:7.1f} cm s-2 and gplanet = "
          "{:7.1f} cm s-2.".format(gplanet, phy.gplanet))
  elif phy.gplanet is not None and phy.rplanet is not None:
    phy.mplanet = phy.gplanet * phy.rplanet**2 / pc.G

  pyrat.phy.rstar = pt.getparam(inputs.rstar, pyrat.radunits)
  isgreater(pyrat.phy.rstar, "cm",   0, True,
            "Stellar radius ({:.3e} cm) must be > 0.", pyrat.log)

  pyrat.phy.gstar  = pt.getparam(inputs.gstar,  "none")
  isgreater(pyrat.phy.gstar, "none", 0, True,
            "Stellar surface gravity ({:.2f} cm s-2) must be > 0.", pyrat.log)

  pyrat.phy.tstar  = pt.getparam(inputs.tstar,  "none")
  isgreater(pyrat.phy.tstar, "none", 0, True,
            "Stellar effective temperature ({:.1f} K) must be > 0.", pyrat.log)

  pyrat.phy.smaxis = pt.getparam(inputs.smaxis, pyrat.radunits)
  isgreater(pyrat.phy.smaxis, "cm",   0, True,
            "Planetary radius ({:.3e} cm) must be > 0.", pyrat.log)

  phy.mstar  = pt.getparam(inputs.mstar, "gram")
  isgreater(phy.mstar, "msun", 0, True,
            "Stellar mass ({:.2f} Msun) must be > 0.", pyrat.log)

  pyrat.phy.tint = pt.defaultp(inputs.tint, 100.0,
            "Planetary internal temperature (tint) defaulted to {:.1f} K.",
            pyrat.wlog, pyrat.log)
  isgreater(phy.tint, "none", 0, True,
            "Planetary internal temperature ({:.1f} K) must be > 0.", pyrat.log)

  # Compute the Hill radius for the planet:
  if (phy.mstar is not None and phy.mplanet is not None and
      phy.smaxis is not None):
    phy.rhill = phy.smaxis * (phy.mplanet/(3*phy.mstar))**(1.0/3.0)

  pyrat.atm.nlayers = pt.getparam(inputs.nlayers, "none", integer=True)
  isgreater(pyrat.atm.nlayers, "none", 0, True,
            "The number of atmospheric layers ({:d}) must be > 0.", pyrat.log)

  # Check Voigt-profile arguments:
  pyrat.voigt.extent = pt.defaultp(inputs.vextent, 20.0,
     "Input Voigt extent (vextent) defaulted to {:g}.", pyrat.wlog, pyrat.log)
  isgreater(pyrat.voigt.extent, "none", 1, False,
            "Voigt extent ({:g}) must be >= 1.0", pyrat.log)

  # Doppler width:
  pyrat.voigt.nDop = pt.defaultp(inputs.nDop, 40,
       "Number of Doppler-width samples (nDop) defaulted to {:d}.",
       pyrat.wlog, pyrat.log)
  isgreater(pyrat.voigt.nDop, "none", 1, False,
       "The number of Doppler-width samples ({:d}) must be >= 1", pyrat.log)

  pyrat.voigt.Dmin = pt.getparam(inputs.Dmin, "none")
  isgreater(pyrat.voigt.Dmin, "none", 0, True,
            "Dmin ({:g} cm-1) must be > 0.", pyrat.log)

  pyrat.voigt.Dmax = pt.getparam(inputs.Dmax, "none")
  isgreater(pyrat.voigt.Dmax, "none", 0, True,
            "Dmax ({:g} cm-1) must be > 0.", pyrat.log)

  if (pyrat.voigt.Dmin is not None and pyrat.voigt.Dmax is not None and
      pyrat.voigt.Dmax <= pyrat.voigt.Dmin):
    pt.error("Dmax ({:g} cm-1) must be > Dmin ({:g} cm-1).".
             format(pyrat.voigt.Dmax, pyrat.voigt.Dmin), pyrat.log)

  # Lorentz width:
  pyrat.voigt.nLor = pt.defaultp(inputs.nLor, 40,
       "Number of Lorentz-width samples (nLor) defaulted to {:d}.",
       pyrat.wlog, pyrat.log)
  isgreater(pyrat.voigt.nLor, "none", 1, False,
       "The number of Lorentz-width samples ({:d}) must be >= 1", pyrat.log)

  pyrat.voigt.Lmin = pt.getparam(inputs.Lmin, "none")
  isgreater(pyrat.voigt.Lmin, "none", 0, True,
            "Lmin ({:g} cm-1) must be > 0.", pyrat.log)

  pyrat.voigt.Lmax = pt.getparam(inputs.Lmax, "none")
  isgreater(pyrat.voigt.Lmax, "none", 0, True,
            "Lmax ({:g} cm-1) must be > 0.", pyrat.log)

  if (pyrat.voigt.Lmin is not None and pyrat.voigt.Lmax is not None and
      pyrat.voigt.Lmax <= pyrat.voigt.Lmin):
    pt.error("Lmax ({:g} cm-1) must be > Lmin ({:g} cm-1).".
             format(pyrat.voigt.Lmax, pyrat.voigt.Lmin), pyrat.log)

  pyrat.voigt.DLratio = pt.defaultp(inputs.DLratio, 0.1,
     "Doppler/Lorentz-width ratio threshold (DLratio) defaulted to {:g}.",
     pyrat.wlog, pyrat.log)
  isgreater(pyrat.voigt.DLratio, "none", 0, True,
     "Doppler/Lorentz-width ratio threshold ({:g}) must be > 0.", pyrat.log)

  # Check extinction-coefficient arguments:
  pyrat.ex.ethresh = pt.getparam(inputs.ethresh, "none")
  isgreater(pyrat.ex.ethresh, "none", 0, True,
        "Extinction-coefficient threshold ({:g}) must be positive.", pyrat.log)
  if pyrat.ex.extfile is not None:
    if inputs.tmin is None:
      pt.error("Undefined lower boundary (tmin) of temperature grid for "
               "extinction-coefficient grid.", pyrat.log)
    else:
      pyrat.ex.tmin = pt.getparam(inputs.tmin, "kelvin")
      isgreater(pyrat.ex.tmin,  "kelvin", 0, True,
            "Minimum temperature sample ({:g} K) must be positive.", pyrat.log)
    if inputs.tmax is None:
      pt.error("Undefined upper boundary (tmax) of temperature grid for "
               "extinction-coefficient grid.", pyrat.log)
    else:
      pyrat.ex.tmax  = pt.getparam(inputs.tmax, "kelvin")
      isgreater(pyrat.ex.tmax,  "kelvin", 0, True,
            "Maximum temperature sample ({:g} K) must be positive.", pyrat.log)

    pyrat.ex.tstep = pt.defaultp(inputs.tstep, 100,
      "Extinction-coefficient grid's temperature sampling interval (tstep) "
      "defaulted to {:g} K.", pyrat.wlog, pyrat.log)

    isgreater(pyrat.ex.tstep, "kelvin", 0, True,
      "Temperature sample step interval ({:g} K) must be positive.", pyrat.log)

    if pyrat.ex.tmax <= pyrat.ex.tmin:
      pt.error("Extinction-coefficient grid's maximum temperature ({:g} K) "
               "must be > minimum temperature ({:g} K).".
                format(pyrat.ex.tmax, pyrat.ex.tmin), pyrat.log)

  # Check haze models:
  if inputs.hazes is not None:
    nhpars = 0
    for hmodel in inputs.hazes:
      if hmodel not in hz.hnames:
        pt.error("Haze model '{:s}' is not in the list of available models:"
                 "\n{:s}".format(hmodel, hz.hnames), pyrat.log)
      else:
        ihaze = np.where(hz.hnames == hmodel)[0][0]
        pyrat.haze.model.append(hz.hmodels[ihaze])
        pyrat.haze.nmodels += 1
        nhpars += pyrat.haze.model[-1].npars
    # Process the haze parameters
    if inputs.hpars is not None:
      if nhpars != len(inputs.hpars):
        pt.error("Number of input haze params ({:d}) does not match the"
                 "number of required haze params({:d}).".
                 format(inputs.hpars, nhpars), pyrat.log)
      j = 0
      for i in np.arange(pyrat.haze.nmodels):
        pyrat.haze.model[i].pars = inputs.hpars[j:j+pyrat.haze.model[i].npars]
        j += pyrat.haze.model[i].npars

  # Check alkali arguments:
  if inputs.alkali is not None:
    nalkali = 0
    for amodel in inputs.alkali:
      if amodel not in al.mnames:
        pt.error("Alkali model '{:s}' is not in the list of available models:"
                 "\n{:s}".format(amodel, al.mnames), pyrat.log)
      ialkali = np.where(al.mnames == amodel)[0][0]
      pyrat.alkali.model.append(al.models[ialkali])
      pyrat.alkali.nmodels += 1

  # Check optical-depth arguments:
  pyrat.od.maxdepth = pt.defaultp(inputs.maxdepth, 10.0,
   "Maximum optical-depth (maxdepth) defaulted to {:g}.", pyrat.wlog, pyrat.log)
  isgreater(pyrat.od.maxdepth, "none", 0, False,
            "Maximum optical-depth limit ({:g}) must be >= 0.0", pyrat.log)

  # Accept ray-path argument:
  pyrat.od.path  = inputs.path
  if pyrat.od.path is None:
    pt.error("Undefined observing geometry (path).  Select between 'transit' "
             "or 'eclipse'.", pyrat.log)
  elif pyrat.od.path not in ['transit', 'eclipse']:
    pt.error("Unknown observing geometry (path = {:s}).  Select between "
             "'transit' or 'eclipse'.".format(pyrat.od.path), pyrat.log)

  # Accept output files:
  pyrat.outspec = pt.defaultp(inputs.outspec, 'outpsec.dat',
     "Output spectrum filename (outspec) defaulted to '{:s}'.",
     pyrat.wlog, pyrat.log)

  pyrat.outsample   = inputs.outsample
  pyrat.outmaxdepth = inputs.outmaxdepth

  # Check system arguments:
  if pyrat.od.path == "transit" and pyrat.phy.rstar is None:
    pt.error("Undefined stellar radius (rstar), required for transmission "
             "spectrum calculation.",  pyrat.log)
  # Stellar-spectrum models:
  pyrat.phy.starspec = inputs.starspec
  pyrat.phy.kurucz   = inputs.kurucz
  pyrat.phy.marcs    = inputs.marcs
  pyrat.phy.phoenix  = inputs.phoenix
  if inputs.starspec is not None and not os.path.isfile(inputs.starspec):
    pt.error("Stellar-spectrum model file: '{:s}' does not exist.".
             format(inputs.starspec), pyrat.log)
  if inputs.kurucz  is not None and not os.path.isfile(inputs.kurucz):
    pt.error("Stellar Kurucz model file: '{:s}' does not exist.".
             format(inputs.kurucz), pyrat.log)
  if inputs.marcs   is not None and not os.path.isfile(inputs.marcs):
    pt.error("Stellar MARCS model file: '{:s}' does not exist.".
             format(inputs.marcs), pyrat.log)
  if inputs.phoenix is not None and not os.path.isfile(inputs.phoenix):
    pt.error("Stellar PHOENIX model file: '{:s}' does not exist.".
             format(inputs.phoenix), pyrat.log)

  # Check raygrid:
  if pyrat.od.path == "eclipse" and inputs.raygrid is None:
    raygrid = pt.defaultp(inputs.raygrid, np.array([0, 20, 40, 60, 80.]),
       "Defaulted raygrid to {:s}.", pyrat.wlog, pyrat.log)
    if raygrid[0] != 0:
      pt.error("First angle in raygrid must be 0.0 (normal to surface).",
               pyrat.log)
    if np.any(raygrid < 0):
      pt.error("raygrid angles must lie between 0 and 90 deg.", pyrat.log)
    if np.any(np.ediff1d(raygrid) <= 0):
      pt.error("raygrid angles must be monotonically increasing.", pyrat.log)
    # Store raygrid values in radians:
    pyrat.raygrid = raygrid * sc.degree

  # Gauss quadrature integration variables:
  pyrat.quadrature = inputs.quadrature
  if inputs.quadrature is not None:
    qnodes, qweights = ss.p_roots(inputs.quadrature)
    pyrat.qnodes   = 0.5*(qnodes + 1.0)
    pyrat.qweights = 0.5 * qweights

  # Observational parameters:
  pyrat.obs.data   = inputs.data
  pyrat.obs.uncert = inputs.uncert
  pyrat.obs.filter = inputs.filter
  # Number of datapoints and filters:
  if inputs.data is not None:
    pyrat.obs.ndata = len(inputs.data)
  if inputs.filter is not None:
    pyrat.obs.nfilters = len(inputs.filter)
  # Number checks:
  if pyrat.obs.uncert is not None and pyrat.obs.ndata != len(pyrat.obs.uncert):
    pt.error("The number of data uncertainty values ({:d}) does not match "
       "the number of data points ({:d}).".
        format(len(pyrat.obs.uncert)), pyrat.obs.ndata)
  if pyrat.obs.filter is not None:
    for f in pyrat.obs.filter:
      if not os.path.isfile(f):
        pt.error("Filter file: '{:s}' does not exist.".format(f), pyrat.log)
    if pyrat.obs.ndata > 0  and  pyrat.obs.ndata != pyrat.obs.nfilters:
      pt.error("The number of filter bands ({:d}) does not match the number "
          "of data points ({:d}).".format(pyrat.obs.nfilters), pyrat.obs.ndata)

  # Retrieval variables:
  # Accept species lists, check after we load the atmospheric model:
  pyrat.ret.bulk     = inputs.bulk
  pyrat.ret.molscale = inputs.molscale
  pyrat.ret.params   = inputs.params
  if pyrat.ret.params is not None:
    pyrat.ret.nparams = len(pyrat.ret.params)
  pyrat.ret.stepsize = inputs.stepsize # FINDME checks
  pyrat.ret.tlow     = pt.getparam(inputs.tlow,  "kelvin")
  pyrat.ret.thigh    = pt.getparam(inputs.thigh, "kelvin")
  if inputs.tmodel is not None and inputs.tmodel not in ["TCEA", "isothermal"]:
    pt.error("Invalid temperature model '{:s}'.  Select from: TCEA or "
             "isothermal".format(inputs.tmodel), pyrat.log)
  pyrat.ret.tmodelname = inputs.tmodel
  if pyrat.ret.tmodelname == "TCEA":
    if pyrat.phy.rstar is None:
      pt.error("Undefined stellar radius (rstar), required for temperature "
               "model.", pyrat.log)
    if pyrat.phy.tstar is None:
      pt.error("Undefined stellar temperature (tstar), required for "
               "temperature model.", pyrat.log)
    if pyrat.phy.smaxis is None:
      pt.error("Undefined orbital semi-major axis (smaxis), required for "
               "temperature model.", pyrat.log)
    if pyrat.phy.gplanet is None:
      pt.error("Undefined planetary surface gravity (gplanet), required for "
               "temperature model.", pyrat.log)

  # Number of processors:
  pyrat.nproc = pt.getparam(inputs.nproc, "none", integer=True)
  isgreater(pyrat.nproc, "none", 1, False,
            "The number of processors ({:d}) must be >= 1.", pyrat.log)
  if pyrat.nproc >= mpr.cpu_count():
    pt.warning(pyrat.verb-2, "The number of requested CPUs ({:d}) is >= "
       "than the number of available CPUs ({:d}).  Enforced nproc to {:d}.".
       format(pyrat.nproc, mpr.cpu_count(), mpr.cpu_count()-1),
       pyrat.log, pyrat.wlog)
    pyrat.nproc = mpr.cpu_count() - 1
  pt.msg(pyrat.verb-3, "Done.", pyrat.log)


def isgreater(value, units, thresh, equal=False, text="", log=None):
  """
  Check that value (if not None) is greater than thresh.
  Throw error if not.

  Parameters
  ----------
  value: Scalar
    The value being tested.
  units: String
    The units of the value.
  thresh: Scalar
    Threshold against which value is being compared.
  equal: Boolean
    If True, strictly require greater than.
  text: String
    Text to show if condition is not satisfied.
  log: File
    Pyrat screen-output log file.

  Returns
  -------
  The value in the pyrat units.
  """
  # Set comparison command:
  if equal:
    compare = np.less_equal
  else:
    compare = np.less

  # Check value:
  if value is None:
    return

  if compare(value, thresh):
    pt.error(text.format(value/pt.u(units)), log, -3)


def setup(pyrat):
  """
  Process retrieval variables: bulk, molscale.
  Process stellar spectrum.
  Process the oberving filter bands.
  """
  # Shortcuts:
  obs = pyrat.obs
  phy = pyrat.phy
  ret = pyrat.ret

  # Setup bulk and variable-abundance species:
  species = pyrat.mol.name
  # Non-overlapping species:
  if ret.bulk is not None  and  len(np.setdiff1d(ret.bulk, species)) > 0:
    pt.error("These bulk species are not present in the atmosphere: {:s}.".
      format(str(np.setdiff1d(ret.bulk, species))), pyrat.log)
  if ret.molscale is not None and len(np.setdiff1d(ret.molscale, species)) > 0:
    pt.error("These variable-abundance species are not present in the "
             "atmosphere: {:s}.".
              format(str(np.setdiff1d(ret.molscale, species))), pyrat.log)
  # Overlapping species:
  if (ret.bulk is not None  and  ret.molscale is not None  and
      len(np.intersect1d(ret.bulk, ret.molscale)) > 0):
    pt.error("These species were marked as both bulk and variable-abundance: "
             "{:s}.".format(np.intersect1d(ret.bulk, ret.molscale)), pyrat.log)

  if pyrat.runmode == "mcmc":
    if ret.bulk is None:
      pt.error("Undefined bulk species list (bulk).", pyrat.log)
    if ret.molscale is None:
      pt.warning(pyrat.verb-2, "There are no variable-abundance species "
                               "(molscale).", pyrat.log, pyrat.wlog)

  # Obtain abundance ratios between the bulk species:
  if ret.bulk is not None:
    ret.ibulk = []
    for mol in ret.bulk:
      ret.ibulk  += list(np.where(species==mol)[0])
    ret.bulkratio, ret.invsrat = atm.ratio(pyrat.atm.q, ret.ibulk)
  if ret.molscale is not None:
    ret.iscale = []
    for mol in ret.molscale:
      ret.iscale += list(np.where(species==mol)[0])
    nabund = len(ret.iscale)
  else:
    nabund = 0


  # Read stellar spectrum model:
  if phy.starspec is not None:
    starwn, starflux = sspec.readpyrat(phy.starspec)
  # Kurucz stellar model:
  elif phy.kurucz is not None:
    if phy.tstar is None:
      pt.error("Undefined stellar temperature (tstar), required for Kurucz "
               "model.", pyrat.log)
    if phy.gstar is None:
      pt.error("Undefined stellar gravity (tstar), required for Kurucz "
               "model.", pyrat.log)
    starflux, starwn, kuruczt, kuruczg = sspec.readkurucz(phy.kurucz,
                                           phy.tstar, np.log10(phy.gstar))
    pt.msg(pyrat.verb-4, "Input stellar params: T={:7.1f} K, log(g)={:4.2f}\n"
                         "Best Kurucz match:    T={:7.1f} K, log(g)={:4.2f}".
          format(phy.tstar, np.log10(phy.gstar), kuruczt, kuruczg, pyrat.log))
  # MARCS stellar model:
  elif phy.marcs:
    pass
  # PHOENIX stellar model:
  elif phy.phoenix:
    pass
  else:
    starflux, starwn = None, None

  # Store input stellar spectrum into pyrat:
  phy.starflux  = starflux
  phy.starwn    = starwn
  # Store interpolated stellar spectrum:
  if phy.starflux is not None:
    sinterp = si.interp1d(phy.starwn, phy.starflux)
    pyrat.spec.starflux = sinterp(pyrat.spec.wn)


  # Skip if there are no filter bands:
  if obs.filter is not None:
    # Load filters:
    bandidx   = []  # Filter wavenumber indices
    starflux  = []  # Interpolated stellar flux
    bandtrans = []  # Normalized interpolated filter transmission
    bandwn    = []  # Band's mean wavenumber
    for i in np.arange(obs.nfilters):
      # Read filter wavenumber and transmission curves:
      filterwn, filtertr = w.readfilter(obs.filter[i])
      # Resample the filters into the stellar wavenumber array:
      btr, wni, isf = w.resample(pyrat.spec.wn, filterwn,   filtertr,
                                                phy.starwn, phy.starflux)
      bandidx.append(wni)
      bandtrans.append(btr)
      starflux.append(isf)
      bandwn.append(np.sum(filterwn*filtertr)/sum(filtertr))

    # Per-band variables:
    obs.bandidx   = bandidx
    obs.bandtrans = bandtrans
    obs.starflux  = starflux
    obs.bandwn    = np.asarray(bandwn)
    obs.bandflux  = np.zeros(obs.nfilters, np.double)

  # Planet-to-star radius ratio:
  if phy.rplanet is not None and phy.rstar is not None:
    phy.rprs = phy.rplanet/phy.rstar

  # Temperature model and arguments:
  if ret.tmodelname == "TCEA":
    ret.tmodel = PT.TCEA
    ret.ntpars = 5
    ret.targs  = [pyrat.atm.press, phy.rstar, phy.tstar, phy.tint,
                  phy.smaxis, phy.gplanet]
    ret.parname = [r"$\log_{10}(\kappa)$", r"$\log_{10}(\gamma_1)$",
                   r"$\log_{10}(\gamma2)$", r"$\alpha$", r"$\beta$"]
  elif ret.tmodelname == "isothermal":
    ret.tmodel = PT.isothermal
    ret.ntpars = 1
    ret.targs  = [pyrat.atm.nlayers]
    ret.parname = [r"$T\ ({\rm K})$"]
  else:
    ret.ntpars = 0

  # Number of free parameters:
  ntemp  = ret.ntpars
  nrad   = int(pyrat.od.path == "transit")
  if nrad == 1:
    ret.parname += [r"${\rm Radius\ (km)}$"]
  for i in np.arange(nabund):
    ret.parname += [r"$\log_{{10}}(f_{{\rm {:s}}})$".format(ret.molscale[i])]

  nhaze  = 0
  if pyrat.ret.nparams > ntemp+nrad+nabund:
    for i in np.arange(pyrat.haze.nmodels):
      nhaze += pyrat.haze.model[i].npars
      ret.parname += pyrat.haze.model[i].parname

  # Indices to parse the array of fitting parameters:
  pyrat.ret.itemp  = np.arange(0,          ntemp)
  pyrat.ret.irad   = np.arange(ntemp,      ntemp+nrad)
  pyrat.ret.iabund = np.arange(ntemp+nrad, ntemp+nrad+nabund)
  pyrat.ret.ihaze  = np.arange(ntemp+nrad+nabund, ntemp+nrad+nabund+nhaze)

  if pyrat.runmode == "mcmc":
    if pyrat.ret.nparams != ntemp+nrad+nabund+nhaze:
      pt.error("The input number of fitting parameters ({:d}) does not "
               "match the number of model parameters ({:d}).".
                format(len(pyrat.ret.params), ntemp+nrad+nabund+nhaze))
