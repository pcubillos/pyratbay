# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
if sys.version_info.major == 3:
    import configparser
else:
    import ConfigParser as configparser
import argparse
import multiprocessing as mp
from datetime import date

import numpy as np
import scipy.constants   as sc
import scipy.interpolate as si
import scipy.special     as ss

from .. import tools      as pt
from .. import constants  as pc
from .. import wine       as pw
from .. import starspec   as ps
from .. import atmosphere as pa
from .. import io         as io
from .. import VERSION    as ver

from .  import haze      as hz
from .  import rayleigh  as ray
from .  import alkali    as al

sys.path.append(pc.ROOT + "/pyratbay/atmosphere/")
import MadhuTP

sys.path.append(pc.ROOT + "/modules/MCcubed")
import MCcubed.utils as mu


def parse(pyrat, cfile, log=None):
  """
  Parse the command-line arguments into the pyrat object.

  Parameters
  ----------
  pyrat: Pyrat object
      A Pyrat instance where to store the command-line arguments.
  cfile: String
      A Pyrat Bay configuration file.
  log: Log object
      An MCcubed.utils.Log instance to log screen outputs to file.
      If None, start a new log from logfile argument in the cfile.
  """
  if not os.path.isfile(cfile):
      print("Configuration file: '{:s}' not found.".format(cfile))
      sys.exit(0)

  # Use argparse only to define the data types and defaults:
  config = configparser.ConfigParser()
  config.optionxform = str  # Enable case-sensitive variable names
  config.read([cfile])
  defaults = dict(config.items("pyrat"))

  parser = argparse.ArgumentParser(
               formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("--verb", type=int, default=None)
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
  pt.addarg("resolution",  group, np.double, None,
      "Output resolving power")
  # Atmospheric sampling options:
  group = parser.add_argument_group("Atmosphere Sampling Options")
  pt.addarg("tmodel",      group, str,       None,
      "Temperature-profile model name.  Select from: isothermal or TCEA.")
  pt.addarg("tpars",       group, pt.parray, None,
      "Temperature model fitting parameters.")
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
  pt.addarg("ethresh",     group, np.double, 1e-15,
      "Extinction-coefficient threshold [default: %(default)s]")
  pt.addarg("nproc",       group, int,       1,
      "Number of processors [default: %(default)s]")
  # Voigt-profile options:
  group = parser.add_argument_group("Voigt-profile  Options")
  pt.addarg("vextent",     group, np.double, None,
      "Extent of Voigt profile in number of Voigt widths [default: 40]")
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
  pt.addarg("rayleigh",    group, pt.parray, None,
      "Rayleigh models [default: %(default)s].")
  pt.addarg("rpars",       group, pt.parray, None,
      "Rayleigh model fitting parameters.")
  pt.addarg("fpatchy",     group, np.double, None,
      "Patchy-clouds factor [default: None].")
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
  pt.addarg("retflag",     group, pt.parray, None,
      "The list of retrieval models, select from: pt mol rad ray haze "
      "cloud patchy.")
  pt.addarg("bulk",        group, pt.parray, None,
      "Bulk-abundance atmospheric species.")
  pt.addarg("molmodel",    group, pt.parray, None,
      "Model to vary species abundance profile. Select from: vert or scale.")
  pt.addarg("molfree",    group, pt.parray, None,
      "Species to vary (one species per each molmodel).")
  pt.addarg("molpars",    group, pt.parray, None,
      "molmodel parameters.")
  pt.addarg("qcap",        group, np.double, 0.99,
      "Maximum acceptable cumulative abundance fraction of traces.")
  pt.addarg("params",      group, pt.parray, None,
      "Initial-guess for retrieval model-fitting parameter.")
  pt.addarg("stepsize",    group, pt.parray, None,
      "Stepsize for retrieval model-fitting parameter.")
  pt.addarg("tlow",        group, np.double, -np.inf,
      "Minimum valid temperature.")
  pt.addarg("thigh",        group, np.double, np.inf,
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
  user, unknown = parser.parse_known_args()

  # Put user arguments into pyrat input:
  for key, value in vars(user).items():
      setattr(pyrat.inputs, key, value)
  pyrat.inputs.configfile = cfile

  # Verbosity level:
  pyrat.verb = pyrat.inputs.verb

  # Open the Pyrat log file if requested:
  if log is not None:  # Take pre-existing log
      pyrat.log = log
      pyrat.logfile = log.logname
  elif pyrat.inputs.logfile is not None:  # Start new log
      pyrat.logfile = pt.path(pyrat.inputs.logfile)
      log = pyrat.log = mu.Log(pyrat.logfile, verb=pyrat.verb, width=80)

  # Welcome message:
  log.msg("{:s}\n"
          "  Python Radiative Transfer (Pyrat).\n"
          "  Version {:d}.{:d}.{:d}.\n"
          "  Copyright (c) 2016-{:d} Patricio Cubillos and collaborators.\n"
          "  Pyrat is (temporarily) proprietaty software (see LICENSE).\n"
          "{:s}\n\n".format(log.sep, ver.PYRAT_VER, ver.PYRAT_MIN,
                            ver.PYRAT_REV, date.today().year, log.sep), verb=0)

  log.msg("Read command-line arguments from configuration file: '{:s}'".
          format(cfile))


def checkinputs(pyrat):
  """
  Check that user input arguments make sense.
  """
  # Shortcuts:
  inputs = pyrat.inputs
  phy    = pyrat.phy
  log    = pyrat.log

  # Pyrat runmode:
  pyrat.runmode = inputs.runmode
  if pyrat.runmode is None:
    log.warning("Defaulted Pyrat's runmode to: spectrum.")
  if pyrat.runmode not in pc.rmodes:
    log.error("Invalid runmode ({:s}).  Select from: {}.".
              format(pyrat.runmode, pc.rmodes))

  # Check that input files exist:
  if inputs.atmfile is None:
    log.error("Undefined atmospheric file (atmfile).")
  elif not os.path.isfile(inputs.atmfile):
    log.error("atmfile: '{:s}' does not exist.".format(inputs.atmfile))
  pyrat.atmfile = inputs.atmfile

  if inputs.linedb is not None:
    for linedb in inputs.linedb:
      if not os.path.isfile(linedb):
        log.error("linedb file: '{:s}' does not exist.".format(linedb))
  pyrat.linedb = pyrat.inputs.linedb

  if inputs.csfile is not None:
    for cs in pyrat.inputs.csfile:
      if not os.path.isfile(cs):
        log.error("Cross-section file: '{:s}' does not exist.".format(cs))
  pyrat.cs.files = pyrat.inputs.csfile

  if inputs.molfile is None: # Set default
    inputs.molfile = os.path.realpath(pc.ROOT + "/inputs/molecules.dat")
  if not os.path.isfile(inputs.molfile):
    log.error("Molecular-data file: '{:s}' does not exist.".
              format(inputs.molfile))
  pyrat.molfile = os.path.realpath(inputs.molfile)

  if inputs.extfile is not None:
    if not os.path.exists(os.path.realpath(os.path.dirname(inputs.extfile))):
      log.error("Directory for extinction-coefficient file '{:s}' does "
                "not exist.".format(inputs.extfile))
    pyrat.ex.extfile = os.path.realpath(inputs.extfile)

  # Check spectrum arguments:
  pyrat.spec.wnunits = pt.defaultp(inputs.wnunits, "cm",
         "wnunits input variable defaulted to '{:s}'.", log)
  pyrat.spec.wlunits = pt.defaultp(inputs.wlunits, "um",
         "wlunits input variable defaulted to '{:s}'.", log)

  pyrat.spec.wllow = pt.getparam(inputs.wllow, pyrat.spec.wlunits, log)
  isgreater(pyrat.spec.wllow, "um", 0, False,
            "Low wavelength boundary ({:.2e} um) must be >= 0.", log)

  pyrat.spec.wlhigh = pt.getparam(inputs.wlhigh, pyrat.spec.wlunits, log)
  isgreater(pyrat.spec.wlhigh, "um", 0, True,
            "High wavelength boundary ({:.2e} um) must be >= 0.", log)

  # Wavenumber must be taken care differently (take inverse of units):
  if inputs.wnlow is not None:
    if len(inputs.wnlow.split()) == 2:
      wnunits = inputs.wnlow.split()[1]
    else:
      wnunits = pyrat.spec.wnunits
    wnlow = float(inputs.wnlow.split()[0])
    if wnlow < 0.0:
      log.error("Low wavenumber boundary ({:.2e} {:s}-1) must be >= 0.".
               format(wnlow, wnunits))
    pyrat.spec.wnlow = wnlow / pt.u(wnunits)

  if   inputs.wnhigh is not None:
    if len(inputs.wnhigh.split()) == 2:
      wnunits = inputs.wnhigh.split()[1]
    else:
      wnunits = pyrat.spec.wnunits
    wnhigh = float(inputs.wnhigh.split()[0])
    if wnhigh <= 0.0:
      log.error("High wavenumber boundary ({:.2e} {:s}-1) must be > 0.".
               format(wnhigh, wnunits))
    pyrat.spec.wnhigh = wnhigh / pt.u(wnunits)

  wnstep = pt.defaultp(inputs.wnstep, "1.0 cm",
     "Input wavenumber sampling step (wnstep) defaulted to '{:s}'.", log)
  if len(wnstep.split()) == 2:
    wnunits = wnstep.split()[1]
  else:
    wnunits = pyrat.spec.wnunits
  wnstep = float(wnstep.split()[0])
  if wnstep <= 0:
    log.error("Wavenumber sampling step ({:.2e} {:s}-1) must be be > 0.".
             format(wnstep, wnunits))
  pyrat.spec.wnstep = wnstep / pt.u(wnunits)

  pyrat.spec.wnosamp = pt.defaultp(inputs.wnosamp, 2160,
     "Input wavenumber oversampling factor (wnosamp) defaulted to {:d}.", log)
  isgreater(pyrat.spec.wnosamp, "none", 1, False,
     "Wavenumber oversampling factor ({:d}) must be >= 1.", log)

  pyrat.spec.resolution = inputs.resolution

  # Check atmospheric layers arguments:
  pyrat.punits = pt.defaultp(inputs.punits, "bar",
     "Input pressure units (punits) defaulted to '{:s}'.", log)
  pyrat.radunits = pt.defaultp(inputs.radunits, "km",
     "Input radius units (punits) defaulted to '{:s}'.", log)

  # Pressure boundaries:
  pyrat.phigh = pt.getparam(inputs.phigh, pyrat.punits, log)
  isgreater(pyrat.phigh, "bar", 0, True,
            "High atm pressure boundary ({:.2e} bar) must be > 0.0", log)
  pyrat.plow  = pt.getparam(inputs.plow,  pyrat.punits, log)
  isgreater(pyrat.plow, "bar",  0, True,
            "Low atm pressure boundary ({:.2e} bar) must be > 0.0", log)
  # Radius boundaries:
  pyrat.radlow  = pt.getparam(inputs.radlow,  pyrat.radunits, log)
  isgreater(pyrat.radlow, "cm", 0, False,
            "Low atm radius boundary ({:.2e} cm) must be >= 0.0", log)
  pyrat.radhigh = pt.getparam(inputs.radhigh, pyrat.radunits, log)
  isgreater(pyrat.radhigh, "cm", 0, True,
            "High atm radius boundary ({:.2e} cm) must be > 0.0", log)
  pyrat.radstep = pt.getparam(inputs.radstep, pyrat.radunits, log)
  isgreater(pyrat.radstep, "cm", 0, True,
            "Radius step size ({:.2f} cm) must be > 0.", log)
  # System physical parameters:
  phy.rplanet = pt.getparam(inputs.rplanet, pyrat.radunits, log)
  isgreater(phy.rplanet, "cm",   0, True,
            "Planetary radius ({:.3e} cm) must be > 0.", log)

  pyrat.refpressure = pt.getparam(inputs.refpressure, pyrat.punits, log)
  isgreater(pyrat.refpressure, "bar", 0, True,
      "Planetary reference pressure level ({:8g} bar) must be > 0.", log)

  phy.gplanet  = pt.getparam(inputs.gplanet, "none", log)
  isgreater(phy.gplanet, "none", 0, True,
            "Planetary surface gravity ({:.2f} cm s-2) must be > 0.", log)

  phy.mplanet  = pt.getparam(inputs.mplanet, "gram", log)
  isgreater(phy.mplanet, "mearth", 0, True,
            "Planetary mass ({:.2f} Mearth) must be > 0.", log)

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
        log.error("Both mplanet and gplanet were provided, but values are "
          "inconsistent (>5%): g(mplanet) = {:7.1f} cm s-2 and gplanet = "
          "{:7.1f} cm s-2.".format(gplanet, phy.gplanet))
  elif phy.gplanet is not None and phy.rplanet is not None:
    phy.mplanet = phy.gplanet * phy.rplanet**2 / pc.G

  pyrat.phy.rstar = pt.getparam(inputs.rstar, pyrat.radunits, log)
  isgreater(pyrat.phy.rstar, "cm",   0, True,
            "Stellar radius ({:.3e} cm) must be > 0.", log)

  pyrat.phy.gstar  = pt.getparam(inputs.gstar, "none", log)
  isgreater(pyrat.phy.gstar, "none", 0, True,
            "Stellar surface gravity ({:.2f} cm s-2) must be > 0.", log)

  pyrat.phy.tstar  = pt.getparam(inputs.tstar, "none", log)
  isgreater(pyrat.phy.tstar, "none", 0, True,
            "Stellar effective temperature ({:.1f} K) must be > 0.", log)

  pyrat.phy.smaxis = pt.getparam(inputs.smaxis, pyrat.radunits, log)
  isgreater(pyrat.phy.smaxis, "cm",   0, True,
            "Planetary radius ({:.3e} cm) must be > 0.", log)

  phy.mstar  = pt.getparam(inputs.mstar, "gram", log)
  isgreater(phy.mstar, "msun", 0, True,
            "Stellar mass ({:.2f} Msun) must be > 0.", log)

  pyrat.phy.tint = pt.defaultp(inputs.tint, 100.0,
            "Planetary internal temperature (tint) defaulted to {:.1f} K.", log)
  isgreater(phy.tint, "none", 0, True,
            "Planetary internal temperature ({:.1f} K) must be > 0.", log)

  # Compute the Hill radius for the planet:
  if (phy.mstar is not None and phy.mplanet is not None and
      phy.smaxis is not None):
    phy.rhill = phy.smaxis * (phy.mplanet/(3*phy.mstar))**(1.0/3.0)

  pyrat.atm.nlayers = pt.getparam(inputs.nlayers, "none", log, integer=True)
  isgreater(pyrat.atm.nlayers, "none", 0, True,
            "The number of atmospheric layers ({:d}) must be > 0.", log)

  # Check Voigt-profile arguments:
  pyrat.voigt.extent = pt.defaultp(inputs.vextent, 20.0,
     "Input Voigt extent (vextent) defaulted to {:g}.", log)
  isgreater(pyrat.voigt.extent, "none", 1, False,
            "Voigt extent ({:g}) must be >= 1.0", log)

  # Doppler width:
  pyrat.voigt.nDop = pt.defaultp(inputs.nDop, 40,
       "Number of Doppler-width samples (nDop) defaulted to {:d}.", log)
  isgreater(pyrat.voigt.nDop, "none", 1, False,
       "The number of Doppler-width samples ({:d}) must be >= 1", log)

  pyrat.voigt.Dmin = pt.getparam(inputs.Dmin, "none", log)
  isgreater(pyrat.voigt.Dmin, "none", 0, True,
            "Dmin ({:g} cm-1) must be > 0.", log)

  pyrat.voigt.Dmax = pt.getparam(inputs.Dmax, "none", log)
  isgreater(pyrat.voigt.Dmax, "none", 0, True,
            "Dmax ({:g} cm-1) must be > 0.", log)

  if (pyrat.voigt.Dmin is not None and pyrat.voigt.Dmax is not None and
      pyrat.voigt.Dmax <= pyrat.voigt.Dmin):
    log.error("Dmax ({:g} cm-1) must be > Dmin ({:g} cm-1).".
             format(pyrat.voigt.Dmax, pyrat.voigt.Dmin))

  # Lorentz width:
  pyrat.voigt.nLor = pt.defaultp(inputs.nLor, 40,
       "Number of Lorentz-width samples (nLor) defaulted to {:d}.", log)
  isgreater(pyrat.voigt.nLor, "none", 1, False,
       "The number of Lorentz-width samples ({:d}) must be >= 1", log)

  pyrat.voigt.Lmin = pt.getparam(inputs.Lmin, "none", log)
  isgreater(pyrat.voigt.Lmin, "none", 0, True,
            "Lmin ({:g} cm-1) must be > 0.", log)

  pyrat.voigt.Lmax = pt.getparam(inputs.Lmax, "none", log)
  isgreater(pyrat.voigt.Lmax, "none", 0, True,
            "Lmax ({:g} cm-1) must be > 0.", log)

  if (pyrat.voigt.Lmin is not None and pyrat.voigt.Lmax is not None and
      pyrat.voigt.Lmax <= pyrat.voigt.Lmin):
    log.error("Lmax ({:g} cm-1) must be > Lmin ({:g} cm-1).".
             format(pyrat.voigt.Lmax, pyrat.voigt.Lmin))

  pyrat.voigt.DLratio = pt.defaultp(inputs.DLratio, 0.1,
     "Doppler/Lorentz-width ratio threshold (DLratio) defaulted to {:g}.", log)
  isgreater(pyrat.voigt.DLratio, "none", 0, True,
     "Doppler/Lorentz-width ratio threshold ({:g}) must be > 0.", log)

  # Check extinction-coefficient arguments:
  pyrat.ex.ethresh = pt.getparam(inputs.ethresh, "none", log)
  isgreater(pyrat.ex.ethresh, "none", 0, True,
        "Extinction-coefficient threshold ({:g}) must be positive.", log)
  # Require tmin, tmax:
  if (pyrat.runmode == "opacity" or
      (pyrat.runmode in ["spectrum", "mcmc"] and
       pyrat.ex.extfile is not None and
       not os.path.isfile(pyrat.ex.extfile))):
    if inputs.tmin is None:
      log.error("Undefined lower boundary (tmin) of temperature grid for "
               "extinction-coefficient grid.")
    if inputs.tmax is None:
      log.error("Undefined upper boundary (tmax) of temperature grid for "
               "extinction-coefficient grid.")

  if inputs.tmin is not None:
    pyrat.ex.tmin = pt.getparam(inputs.tmin, "kelvin", log)
    isgreater(pyrat.ex.tmin,  "kelvin", 0, True,
          "Minimum temperature sample ({:g} K) must be positive.", log)
  if inputs.tmax is not None:
    pyrat.ex.tmax  = pt.getparam(inputs.tmax, "kelvin", log)
    isgreater(pyrat.ex.tmax,  "kelvin", 0, True,
          "Maximum temperature sample ({:g} K) must be positive.", log)

    pyrat.ex.tstep = pt.defaultp(inputs.tstep, 100,
      "Extinction-coefficient grid's temperature sampling interval (tstep) "
      "defaulted to {:g} K.", log)

    isgreater(pyrat.ex.tstep, "kelvin", 0, True,
      "Temperature sample step interval ({:g} K) must be positive.", log)

    if pyrat.ex.tmax <= pyrat.ex.tmin:
      log.error("Extinction-coefficient grid's maximum temperature ({:g} K) "
               "must be > minimum temperature ({:g} K).".
                format(pyrat.ex.tmax, pyrat.ex.tmin))

  # Check haze models:
  if inputs.hazes is not None:
    nhpars = 0
    for hmodel in inputs.hazes:
      if hmodel not in hz.hnames:
        log.error("Haze model '{:s}' is not in the list of available models:"
                 "\n{:s}".format(hmodel, hz.hnames))
      else:
        ihaze = np.where(hz.hnames == hmodel)[0][0]
        pyrat.haze.model.append(hz.hmodels[ihaze])
        pyrat.haze.nmodels += 1
        nhpars += pyrat.haze.model[-1].npars
    # Process the haze parameters:
    pyrat.haze.pars = inputs.hpars
    if pyrat.haze.pars is not None:
      if nhpars != len(pyrat.haze.pars):
        log.error("The number of input haze parameters ({:d}) does not match "
                 "the number of required haze parameters ({:d}).".
                 format(len(pyrat.haze.pars), nhpars))
      j = 0
      for i in np.arange(pyrat.haze.nmodels):
        npars = pyrat.haze.model[i].npars
        pyrat.haze.model[i].pars = pyrat.haze.pars[j:j+npars]
        j += npars

  if inputs.fpatchy is not None:
    if inputs.fpatchy < 0 or inputs.fpatchy > 1:
      log.error("Invalid patchy-cloud fraction ({:g}).  fpatchy must be "
               "in the range 0--1.".format(inputs.fpatchy))
    pyrat.haze.fpatchy = inputs.fpatchy

  # Check Rayleigh models:
  if inputs.rayleigh is not None:
    nrpars = 0
    for rmodel in inputs.rayleigh:
      if rmodel not in ray.rnames:
        log.error("Rayleigh model '{:s}' is not in the list of available "
                  "models:\n{:s}".format(rmodel, ray.rnames))
      j = np.where(ray.rnames == rmodel)[0][0]
      pyrat.rayleigh.model.append(ray.rmodels[j])
      pyrat.rayleigh.nmodels += 1
      nrpars += pyrat.rayleigh.model[-1].npars
    # Process the Rayleigh parameters:
    pyrat.rayleigh.pars = inputs.rpars
    if pyrat.rayleigh.pars is not None:
      if nrpars != len(pyrat.rayleigh.pars):
        log.error("The number of input Rayleigh parameters ({:d}) does not "
                 "match the number of required parameters ({:d}).".
                 format(len(pyrat.rayleigh.pars), nrpars))
      j = 0
      for i in np.arange(pyrat.rayleigh.nmodels):
        npars = pyrat.rayleigh.model[i].npars
        pyrat.rayleigh.model[i].pars = pyrat.rayleigh.pars[j:j+npars]
        j += npars

  # Check alkali arguments:
  if inputs.alkali is not None:
    for amodel in inputs.alkali:
      if amodel not in al.mnames:
        log.error("Alkali model '{:s}' is not in the list of available models:"
                 "\n{:s}.".format(amodel, al.mnames))
      ialkali = np.where(al.mnames == amodel)[0][0]
      pyrat.alkali.model.append(al.models[ialkali])
      pyrat.alkali.nmodels += 1

  # Check optical-depth arguments:
  pyrat.od.maxdepth = pt.defaultp(inputs.maxdepth, 10.0,
   "Maximum optical-depth (maxdepth) defaulted to {:g}.", log)
  isgreater(pyrat.od.maxdepth, "none", 0, False,
            "Maximum optical-depth limit ({:g}) must be >= 0.0", log)

  # Accept ray-path argument:
  pyrat.od.path  = inputs.path
  if pyrat.runmode in ["spectrum", "mcmc"]: # Check only if computing spectrum
    if pyrat.od.path is None:
      log.error("Undefined observing geometry (path).  Select between "
               "'transit' or 'eclipse'.")
    elif pyrat.od.path not in ['transit', 'eclipse']:
      log.error("Unknown observing geometry (path = {:s}).  Select between "
               "'transit' or 'eclipse'.".format(pyrat.od.path))

  # Accept output files:
  pyrat.outspec = pt.defaultp(inputs.outspec, 'outspec.dat',
     "Output spectrum filename (outspec) defaulted to '{:s}'.", log)

  pyrat.outsample   = inputs.outsample
  pyrat.outmaxdepth = inputs.outmaxdepth

  # Check system arguments:
  if pyrat.od.path == "transit" and pyrat.phy.rstar is None:
    log.error("Undefined stellar radius (rstar), required for transmission "
             "spectrum calculation.")
  # Stellar-spectrum models:
  pyrat.phy.starspec = inputs.starspec
  pyrat.phy.kurucz   = inputs.kurucz
  pyrat.phy.marcs    = inputs.marcs
  pyrat.phy.phoenix  = inputs.phoenix
  if inputs.starspec is not None and not os.path.isfile(inputs.starspec):
    log.error("Stellar-spectrum model file: '{:s}' does not exist.".
             format(inputs.starspec))
  if inputs.kurucz  is not None and not os.path.isfile(inputs.kurucz):
    log.error("Stellar Kurucz model file: '{:s}' does not exist.".
             format(inputs.kurucz))
  if inputs.marcs   is not None and not os.path.isfile(inputs.marcs):
    log.error("Stellar MARCS model file: '{:s}' does not exist.".
             format(inputs.marcs))
  if inputs.phoenix is not None and not os.path.isfile(inputs.phoenix):
    log.error("Stellar PHOENIX model file: '{:s}' does not exist.".
             format(inputs.phoenix))

  # Check raygrid:
  if inputs.raygrid is None:
    raygrid = pt.defaultp(inputs.raygrid, np.array([0, 20, 40, 60, 80.]),
        "Defaulted emission raygrid to {}.", log)
  else:
    raygrid = inputs.raygrid
    if raygrid[0] != 0:
      log.error("First angle in raygrid must be 0.0 (normal to surface).")
    if np.any(raygrid < 0) or np.any(raygrid > 90):
      log.error("raygrid angles must lie between 0 and 90 deg.")
    if np.any(np.ediff1d(raygrid) <= 0):
      log.error("raygrid angles must be monotonically increasing.")
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
    log.error("The number of data uncertainty values ({:d}) does not match "
       "the number of data points ({:d}).".
        format(len(pyrat.obs.uncert), pyrat.obs.ndata))
  if pyrat.obs.filter is not None:
    for f in pyrat.obs.filter:
      if not os.path.isfile(f):
        log.error("Filter file: '{:s}' does not exist.".format(f))
    if pyrat.obs.ndata > 0  and  pyrat.obs.ndata != pyrat.obs.nfilters:
      log.error("The number of filter bands ({:d}) does not match the number "
          "of data points ({:d}).".format(pyrat.obs.nfilters, pyrat.obs.ndata))

  # Retrieval variables:
  # Accept species lists, check after we load the atmospheric model:
  pyrat.ret.retflag  = inputs.retflag
  pyrat.ret.qcap     = inputs.qcap
  pyrat.ret.params   = inputs.params
  if pyrat.ret.params is not None:
      pyrat.ret.nparams = len(pyrat.ret.params)
  pyrat.ret.stepsize = inputs.stepsize # FINDME checks
  pyrat.ret.tlow     = pt.getparam(inputs.tlow,  "kelvin", log)
  pyrat.ret.thigh    = pt.getparam(inputs.thigh, "kelvin", log)

  # Atmospheric model:
  pyrat.atm.molmodel = inputs.molmodel
  pyrat.atm.molfree  = inputs.molfree
  pyrat.atm.molpars  = inputs.molpars
  pyrat.atm.bulk     = inputs.bulk
  if inputs.tmodel is not None and inputs.tmodel not in \
          ["TCEA", "isothermal", "MadhuInv", "MadhuNoInv"]:
    log.error("Invalid temperature model '{:s}'.  Select from: "
              "TCEA, MadhuInv, MadhuNoInv or isothermal.".format(inputs.tmodel))
  pyrat.atm.tmodelname = inputs.tmodel

  pyrat.atm.tpars = inputs.tpars

  if np.abs(pyrat.ret.qcap-0.5) > 0.5:
    log.error("Trace abundances cap (qcap={:.3f}) must lie in the range "
             "between 0.0 and 1.0.".format(pyrat.ret.qcap))
  if pyrat.atm.tmodelname == "TCEA":
    if pyrat.phy.rstar is None:
      log.error("Undefined stellar radius (rstar), required for temperature "
               "model.")
    if pyrat.phy.tstar is None:
      log.error("Undefined stellar temperature (tstar), required for "
               "temperature model.")
    if pyrat.phy.smaxis is None:
      log.error("Undefined orbital semi-major axis (smaxis), required for "
               "temperature model.")
    if pyrat.phy.gplanet is None:
      log.error("Undefined planetary surface gravity (gplanet), required for "
               "temperature model.")

  # Number of processors:
  pyrat.nproc = pt.getparam(inputs.nproc, "none", log, integer=True)
  isgreater(pyrat.nproc, "none", 1, False,
            "The number of processors ({:d}) must be >= 1.", log)
  if pyrat.nproc >= mp.cpu_count():
    log.warning("The number of requested CPUs ({:d}) is >= than the number "
       "of available CPUs ({:d}).  Enforced nproc to {:d}.".
       format(pyrat.nproc, mp.cpu_count(), mp.cpu_count()-1))
    pyrat.nproc = mp.cpu_count() - 1
  log.msg("Check inputs done.")


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
    log.error(text.format(value/pt.u(units)), tracklev=-3)


def setup(pyrat):
  """
  Process retrieval variables: bulk, molmodel.
  Process stellar spectrum.
  Process the oberving filter bands.
  """
  # Shortcuts:
  phy = pyrat.phy
  ret = pyrat.ret
  atm = pyrat.atm
  log = pyrat.log

  # Setup bulk and variable-abundance species:
  species = pyrat.mol.name
  # Non-overlapping species:
  if atm.bulk is not None  and  len(np.setdiff1d(atm.bulk, species)) > 0:
    log.error("These bulk species are not present in the atmosphere: {:s}.".
      format(str(np.setdiff1d(atm.bulk, species))))

  if atm.molmodel is not None:
      if atm.molfree is None:
          log.error("molmodel is set, but there are no molfree.")
      if len(atm.molmodel) != len(atm.molfree):
          log.error("There should be one molfree for each molmodel:\n"
              "molmodel: {}\nmolfree: {}".format(atm.molmodel, atm.molfree))
      if len(np.setdiff1d(atm.molfree, species)) > 0:
          log.error("These species are not present in the atmosphere: {:s}.".
                    format(str(np.setdiff1d(atm.molfree, species))))

  # Overlapping species:
  if (atm.bulk is not None  and  atm.molfree is not None  and
      len(np.intersect1d(atm.bulk, atm.molfree)) > 0):
    log.error("These species were marked as both bulk and variable-abundance: "
             "{:s}.".format(np.intersect1d(atm.bulk, atm.molfree)))

  if pyrat.runmode == "mcmc":
    if ret.retflag is None:
      log.error("Unspecified retrieval model flags.  Set the retflag list "
               "of models selecting from: {:s}.".format(ret.rmodels))
    elif not np.all(np.in1d(ret.retflag, ret.rmodels)):
      log.error("Invalid retrieval model flags in retflag={}.  Available "
               "options are: {}.".format(ret.retflag, ret.rmodels))
    if atm.bulk is None and "mol" in ret.retflag:
      log.error("Undefined bulk species list (bulk).")
    if atm.molmodel is None and "mol" in ret.retflag:
      log.error("Species abundances included for retrieval (retflag contains "
               "'mol') but there are no abundance model (molmodel).")

  # Obtain abundance ratios between the bulk species:
  spec = list(species)
  if atm.bulk is not None:
      atm.ibulk = [spec.index(mol) for mol in atm.bulk]
      atm.bulkratio, atm.invsrat = pa.ratio(pyrat.atm.q, atm.ibulk)
  if atm.molmodel is not None:
      atm.ifree = [spec.index(mol) for mol in atm.molfree]
      nabund = len(atm.ifree)
      # Abundance free-parameter names:
      mpnames   = ["log({:s})".format(mol) for mol in atm.molfree]
      mtexnames = [r"$\log_{{10}}(f_{{\rm {:s}}})$".format(mol)
                   for mol in atm.molfree]
  else:
      nabund = 0
      mpnames, mtexnames = [], []

  # Read stellar spectrum model:
  if phy.starspec is not None:
    starwn, starflux = io.read_spectrum(phy.starspec)
  # Kurucz stellar model:
  elif phy.kurucz is not None:
    if phy.tstar is None:
      log.error("Undefined stellar temperature (tstar), required for Kurucz "
               "model.")
    if phy.gstar is None:
      log.error("Undefined stellar gravity (gstar), required for Kurucz model.")
    starflux, starwn, kuruczt, kuruczg = ps.readkurucz(phy.kurucz,
                                           phy.tstar, np.log10(phy.gstar))
    log.msg("Input stellar params: T={:7.1f} K, log(g)={:4.2f}\n"
            "Best Kurucz match:    T={:7.1f} K, log(g)={:4.2f}".
            format(phy.tstar, np.log10(phy.gstar), kuruczt, kuruczg), verb=2)
  # MARCS stellar model:
  elif phy.marcs is not None:
    pass
  # PHOENIX stellar model:
  elif phy.phoenix is not None:
    pass
  # Blackbody stellar model:
  elif phy.tstar is not None:
    starwn   = pyrat.spec.wn
    starflux = ps.bbflux(starwn, phy.tstar)
  else:
    starflux, starwn = None, None

  # Store input stellar spectrum into pyrat:
  phy.starflux  = starflux
  phy.starwn    = starwn
  # Store interpolated stellar spectrum:
  if phy.starflux is not None:
    sinterp = si.interp1d(phy.starwn, phy.starflux)
    pyrat.spec.starflux = sinterp(pyrat.spec.wn)

  # Set observational variables (for given filters and other parameters):
  setfilters(pyrat.obs, pyrat.spec, pyrat.phy)

  # Planet-to-star radius ratio:
  if phy.rplanet is not None and phy.rstar is not None:
    phy.rprs = phy.rplanet/phy.rstar

  # Temperature models and arguments:
  if atm.tmodelname == "TCEA":
    ntemp = 5
    atm.tmodel = pa.temp_TCEA
    atm.targs  = [pyrat.atm.press, phy.rstar, phy.tstar, phy.tint,
                  phy.gplanet, phy.smaxis]
    tpnames   = ["log(kappa)", "log(gamma1)", "log(gamma2)", "alpha", "beta"]
    ttexnames = [r"$\log_{10}(\kappa)$", r"$\log_{10}(\gamma_1)$",
                 r"$\log_{10}(\gamma2)$", r"$\alpha$", r"$\beta$"]
  elif atm.tmodelname == "isothermal":
    ntemp = 1
    atm.tmodel = pa.temp_isothermal
    atm.targs = [pyrat.atm.nlayers]
    tpnames   = ["T (K)"]
    ttexnames = [r"$T\ ({\rm K})$"]
  elif atm.tmodelname == "MadhuNoInv":
    ntemp = 5
    atm.tmodel = MadhuTP.no_inversion
    atm.targs  = [pyrat.atm.press*1e-6]
    tpnames    = ["a1", "a2", "p1", "p3", "T3"]
    ttexnames  = [r"$a_1$", r"$a_2$", r"$p_1$", r"$p_3$", r"$T_3$"]
  elif atm.tmodelname == "MadhuInv":
    ntemp = 6
    atm.tmodel = MadhuTP.inversion
    atm.targs  = [pyrat.atm.press*1e-6]
    tpnames    = ["a1", "a2", "p1", "p2", "p3", "T3"]
    ttexnames  = [r"$a_1$", r"$a_2$", r"$p_1$", r"$p_2$", r"$p_3$", r"$T_3$"]
  else:
    ntemp = 0
    tpnames, ttexnames = [], []

  # Rayleigh models:
  nray     = 0
  rpnames, rtexnames = [], []
  for i in np.arange(pyrat.rayleigh.nmodels):
    rpnames   += pyrat.rayleigh.model[i].pnames
    rtexnames += pyrat.rayleigh.model[i].texnames
    nray += pyrat.rayleigh.model[i].npars

  # Haze models:
  nhaze    = 0
  hpnames, htexnames = [], []
  for i in np.arange(pyrat.haze.nmodels):
    hpnames   += pyrat.haze.model[i].pnames
    htexnames += pyrat.haze.model[i].texnames
    nhaze += pyrat.haze.model[i].npars

  # Indices to parse the array of fitting parameters:
  if ret.retflag is None:
    ret.retflag = []
  nparams = 0
  ret.pnames, ret.texnames = [], []
  if "pt" in ret.retflag:
    ret.itemp  = np.arange(nparams, nparams + ntemp)
    ret.pnames    += tpnames
    ret.texnames += ttexnames
    nparams += ntemp
  if "rad" in ret.retflag:
    ret.irad   = np.arange(nparams, nparams + 1)  # nrad is always 1
    ret.pnames   += ["Radius (km)"]
    ret.texnames += [r"${\rm Radius\ (km)}$"]
    nparams += 1
  if "mol" in ret.retflag:
    ret.iabund = np.arange(nparams, nparams + nabund)
    ret.pnames   += mpnames
    ret.texnames += mtexnames
    nparams += nabund
  if "ray" in ret.retflag:
    ret.iray   = np.arange(nparams, nparams + nray)
    ret.pnames    += rpnames
    ret.texnames += rtexnames
    nparams += nray
  if "haze" in ret.retflag:
    ret.ihaze  = np.arange(nparams, nparams + nhaze)
    ret.pnames   += hpnames
    ret.texnames += htexnames
    nparams += nhaze
  #if "cloud" in ret.retflag:
  #  ret.icloud  = np.arange(nparams, nparams + ncloud)
  #  ret.pnames    += cpnames
  #  ret.texnames += ctexnames
  #  nparams += ncloud
  if "patchy" in ret.retflag:
    ret.ipatchy = np.arange(nparams, nparams + 1)  # npatchy is always 1
    ret.pnames   += ["f_patchy"]
    ret.texnames += [r"$f_{\rm patchy}$"]
    nparams += 1

  if pyrat.runmode == "mcmc":
    if ret.nparams != nparams:
      log.error("The input number of fitting parameters ({:d}) does not "
                "match the number of model parameters ({:d}).".
                 format(ret.nparams, nparams))

  # Check for non-retrieval model/parameters:
  if (pyrat.rayleigh.nmodels > 0 and
      (pyrat.runmode != "mcmc" or "ray" not in ret.retflag)):
    if pyrat.rayleigh.pars is None:
      log.error("Rayleigh parameters (rpars) have not been specified.")
  if (pyrat.haze.nmodels > 0 and
      (pyrat.runmode != "mcmc" or "haze" not in ret.retflag)):
    if pyrat.haze.pars is None:
      log.error("Haze parameters (hpars) have not been specified.")


def setfilters(obs, spec, phy):
  """
  Set observational variables (pyrat.obs) based on given parameters.
  """
  # Skip if there are no filter bands:
  if obs.filter is None:
    return
  # Load filters:
  bandidx   = []  # Filter wavenumber indices
  starflux  = []  # Interpolated stellar flux
  bandtrans = []  # Normalized interpolated filter transmission
  bandwn    = []  # Band's mean wavenumber
  for i in np.arange(obs.nfilters):
    # Read filter wavenumber and transmission curves:
    filterwn, filtertr = io.read_filter(obs.filter[i])
    # Resample the filters into the stellar wavenumber array:
    btr, wni, isf = pw.resample(spec.wn, filterwn,   filtertr,
                                         phy.starwn, phy.starflux)
    bandidx.append(wni)
    bandtrans.append(btr)
    starflux.append(isf)
    bandwn.append(np.sum(filterwn*filtertr)/np.sum(filtertr))

  # Per-band variables:
  obs.bandidx   = bandidx
  obs.bandtrans = bandtrans
  obs.starflux  = starflux
  obs.bandwn    = np.asarray(bandwn)
  obs.bandflux  = np.zeros(obs.nfilters, np.double)
