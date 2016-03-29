import sys, os
import argparse, ConfigParser
import numpy as np
import scipy.constants as sc

from .. import tools     as pt
from .. import constants as pc
from .. import VERSION   as ver
from .  import haze      as hz
from .  import alkali    as al


def parse(pyrat):
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
  group.add_argument("--atmfile",    dest="atmfile",
                     help="Atmospheric file [default: %(default)s]",
                     action="store", type=str, default=None) 
  group.add_argument("--linedb",     dest="linedb",
                     help="Line database files [default: %(default)s]",
                     action="store", type=pt.parray, default=None) 
  group.add_argument("--csfile",     dest="csfile",
                     help="Cross-section files [default: %(default)s]",
                     action="store", type=pt.parray, default=None)
  group.add_argument("--molfile",    dest="molfile",
                     help="Molecular info file [default: "
                          "'pyrat/inputs/molecules.dat']",
                     action="store", type=str, default=None) 
  group.add_argument("--extfile",        dest="extfile",
                     help="Extinction-coefficient table file [default: "
                          "%(default)s]",
                     action="store", type=str, default=None)
  # Spectrum sampling options:
  group = parser.add_argument_group("Spectrum Sampling Options")
  group.add_argument("--wlunits",    dest="wlunits",
                     help="Wavelength (input) units [default: %(default)s]",
                     action="store", type=str, default='um',
                     choices=('A','nm','um','mm','cm','m'))
  group.add_argument("--wllow",      dest="wllow",
                     help="Wavelength low boundary [default: %(default)s]",
                     action="store", type=str, default=None)
  group.add_argument("--wlhigh",     dest="wlhigh",
                     help="Wavelength high boundary [default: %(default)s]",
                     action="store", type=str, default=None)

  group.add_argument("--wnunits",    dest="wnunits",
                     help="Wavenumber (input) inverse units [default: "
                          "%(default)s]",
                     action="store", type=str, default='cm',
                     choices=('A','nm','um','mm','cm','m'))
  group.add_argument("--wnlow",      dest="wnlow",
                     help="Wavenumber low boundary [default: %(default)s]",
                     action="store", type=str, default=None)
  group.add_argument("--wnhigh",     dest="wnhigh",
                     help="Wavenumber high boundary [default: %(default)s]",
                     action="store", type=str, default=None)
  group.add_argument("--wnstep",     dest="wnstep",
                     help="Wavenumber sampling step [default: %(default)s]",
                     action="store", type=str, default="1.0 cm")
  group.add_argument("--wnosamp",       dest="wnosamp",
                     help="Wavenumber oversampling factor "
                          "[default: %(default)s]",
                     action="store", type=int, default=2160)
  # Atmospheric sampling options:
  group = parser.add_argument_group("Atmosphere Sampling Options")
  group.add_argument("--radlow",     dest="radlow",
                     help="Atmospheric radius low boundary [default: "
                          "Use atmospheric file value]",
                     action="store", type=str, default=None)
  group.add_argument("--radhigh",    dest="radhigh",
                     help="Atmospheric radius high boundary [default: "
                          "Use atmospheric file value]",
                     action="store", type=str, default=None)
  group.add_argument("--radstep",        dest="radstep",
                     help="Atmospheric radius sampling step [default: "
                          "Use atmospheric file value]",
                     action="store", type=str, default=None)
  group.add_argument("--radunits",      dest="radunits",
                     help="Radius (user) units [default: %(default)s]",
                     action="store", type=str, default='km')
  group.add_argument("--plow",          dest="plow",
                     help="Atmospheric pressure low boundary (overrides "
                          "radius high boundary) [default: %(default)s]",
                     action="store", type=str, default=None)
  group.add_argument("--phigh",         dest="phigh",
                     help="Atmospheric pressure high boundary (overrides "
                          "radius  low boundary) [default: %(default)s]",
                     action="store", type=str, default=None)
  group.add_argument("--nlayers",       dest="nlayers",
                     help="Number of atmospheric pressure samples [default: "
                          "Use atmospheric file value]",
                     action="store", type=np.int, default=None)
  group.add_argument("--punits",        dest="punits",
                     help="Pressure (user) units [default: %(default)s]",
                     action="store", type=str, default='bar',
                     choices=('bar',))
  group.add_argument("--radiusbase",    dest="radiusbase",
                     help="Planetary radius base level (in radunits)",
                     action="store", type=str, default=None)
  group.add_argument("--pressurebase",  dest="pressurebase",
                     help="Planetary pressure base level (in punits)",
                     action="store", type=str, default=None)
  group.add_argument("--surfgravity",   dest="surfgravity",
                     help="Planet's surface gravity in cm s-2",
                     action="store", type=np.double, default=None)
  # Extinction options:
  group = parser.add_argument_group("Extinction Calculations Options")
  group.add_argument("--tmin",          dest="tmin",
                     help="Minimum temperature to sample/consider "
                     " in Kelvin [default: %(default)s]",
                     action="store", type=np.double, default=None)
  group.add_argument("--tmax",          dest="tmax",
                     help="Maximum temperature to sample/consider "
                     "in Kelvin [default: %(default)s]",
                     action="store", type=np.double, default=None)
  group.add_argument("--tstep",          dest="tstep",
                     help="Temperature sample step interval "
                     "in Kelvin [default: %(default)s]",
                     action="store", type=np.double, default=100.0)
  group.add_argument("--ethresh",       dest="ethresh",
                     help="Extinction-coefficient threshold "
                          "[default: %(default)s]",  # FINDME: Explain better
                     action="store", type=np.double, default=1e-6)
  # Voigt-profile options:
  group = parser.add_argument_group("Voigt-profile  Options")
  group.add_argument(      "--vextent",    dest="vextent",
                     help="Extent of Voigt profile in number of Voigt  "
                          "widths [default: %(default)s]",
                     action="store", type=np.double, default=20)
  group.add_argument(      "--Dmin",          dest="Dmin",
                     help="Minimum Doppler width to sample in cm-1 "
                          "[default: %(default)s]",
                     action="store", type=np.double, default=1e-3)
  group.add_argument(      "--Dmax",          dest="Dmax",
                     help="Maximum Doppler width to sample in cm-1 "
                          "[default: %(default)s]",
                     action="store", type=np.double, default=0.25)
  group.add_argument(      "--nDop",          dest="nDop",
                     help="Number of Doppler-width samples"
                          "[default: %(default)s]",
                     action="store", type=np.int, default=40)
  group.add_argument(      "--Lmin",          dest="Lmin",
                     help="Minimum Lorentz width to sample in cm-1 "
                          "[default: %(default)s]",
                     action="store", type=np.double, default=1e-4)
  group.add_argument(      "--Lmax",          dest="Lmax",
                     help="Maximum Lorentz width to sample in cm-1 "
                          "[default: %(default)s]",
                     action="store", type=np.double, default=10)
  group.add_argument(      "--nLor",          dest="nLor",
                     help="Number of Lorentz-width samples"
                          "[default: %(default)s]",
                     action="store", type=np.int, default=40)
  group.add_argument(      "--DLratio",      dest="DLratio",
                     help="Minimum Doppler/Lorentz-width ratio to re-calculate"
                          "a Voigt profile [default: %(default)s]",
                     action="store", type=np.double, default=0.1)
  # Hazes and clouds options:
  group = parser.add_argument_group("Hazes and Clouds Options")
  group.add_argument("--hazes",   dest="hazes",
                     help="Haze models [default: %(default)s].",
                     action="store", type=pt.parray, default=None)
  group.add_argument("--hpars",   dest="hpars",
                     help="Haze model fitting parameters.",
                     action="store", type=pt.parray, default=None)
  # Alkali opacity options:
  group = parser.add_argument_group("Alkali Options")
  group.add_argument("--alkali",   dest="alkali",
                     help="Alkali absorption models [default: %(default)s].",
                     action="store", type=pt.parray, default=None)
  # Optical depth options:
  group = parser.add_argument_group("Optical Depth Options")
  group.add_argument("--path",          dest="path",
                     help="Lightray-path geometry [default: %(default)s]",
                     action="store", type=str, default='transit',
                     choices=('transit', 'eclipse'))
  group.add_argument("--maxdepth",       dest="maxdepth",
                     help="Maximum optical depth to calculate [default: "
                          "%(default)s]",
                     action="store", type=np.double, default=10)
  group.add_argument("--raygrid",   dest="raygrid",
                     help="Grid of incident angles (degrees).",
                     action="store", type=pt.parray, default="0 20 40 60 80")
  # System options:
  group = parser.add_argument_group("System Options")
  group.add_argument(      "--rstar",       dest="rstar",
                     help="Stellar radius [default: %(default)s]",
                     action="store", type=str, default=None)
  # Output file options:
  group = parser.add_argument_group("Output File's Options")
  group.add_argument("-o", "--outspec",       dest="outspec",
                     help="Output spectrum file [default: %(default)s]",
                     action="store", type=str, default='output.dat') 
  group.add_argument(      "--outsample",     dest="outsample",
                     help="Output samplings file [default: %(default)s]",
                     action="store", type=str, default=None) 
  group.add_argument(      "--outmaxdepth",    dest="outmaxdepth",
                     help="Filename to store the radius at maxdepth "
                          "(per wavelength) [default: %(default)s]",
                     action="store", type=str, default=None) 
  group.add_argument(      "--logfile", dest="logfile",
                     help="Log file.",
                     action="store", default=None)

  # Set the defaults from the configuration file:
  parser.set_defaults(**defaults)
  # Set values from command line:
  user, unknown = parser.parse_known_args(remaining_argv)

  # Put user arguments into pyrat input:
  pyrat.inputs.configfile = args.configfile
  pyrat.inputs.verb       = user.verb
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
  pyrat.inputs.radiusbase   = user.radiusbase
  pyrat.inputs.pressurebase = user.pressurebase
  pyrat.inputs.surfgravity  = user.surfgravity
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
  # System:
  pyrat.inputs.rstar = user.rstar
  # Output files:
  pyrat.inputs.outspec     = user.outspec
  pyrat.inputs.outsample   = user.outsample
  pyrat.inputs.outmaxdepth = user.outmaxdepth
  pyrat.inputs.logfile     = user.logfile

  # Open the Pyrat log file if requested:
  if pyrat.inputs.logfile is not None:
    pyrat.logfile = os.path.realpath(pyrat.inputs.logfile)
    pyrat.log = open(pyrat.logfile, "w")

  # Welcome message:
  pt.msg(1, "{:s}\n  Python Radiative Transfer (PyRaT).\n"
            "  Version {:d}.{:d}.{:d}.\n"
            "  Copyright (c) 2016 Patricio Cubillos and collaborators.\n"
            "  Pyrat is open-source software under the FINDME license.\n"
            "{:s}\n\n".format(pt.sep, ver.PYRAT_VER, ver.PYRAT_MIN,
                                      ver.PYRAT_REV, pt.sep), pyrat.log)

  pt.msg(1, "Read command-line arguments from configuration file: '{:s}'".
            format(cfile), pyrat.log)

def checkinputs(pyrat):
  """
  Check that user input arguments make sense.
  """
  # User-inputs object:
  inputs = pyrat.inputs

  # Path to source parent's folder:
  pyratdir = os.path.dirname(os.path.realpath(__file__))

  # Check that input files exist:
  if not os.path.isfile(inputs.atmfile):
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
    inputs.molfile = pyratdir + "/../../inputs/molecules.dat"
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
  pyrat.spec.wnunits = inputs.wnunits  # Accept units
  pyrat.spec.wlunits = inputs.wlunits

  pyrat.spec.wllow = pt.getparam(inputs.wllow, pyrat.spec.wlunits)
  isgreater(pyrat.spec.wllow, "um", 0, False,
            "Low wavelength boundary ({:.2e} um) must be >= 0.")

  pyrat.spec.wlhigh = pt.getparam(inputs.wlhigh, pyrat.spec.wlunits)
  isgreater(pyrat.spec.wlhigh, "um", 0, True,
            "High wavelength boundary ({:.2e} um) must be >= 0.")

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

  if inputs.wnstep is None:
    pt.error("Wavenumber sampling step must (wnstep) be defined.", pyrat.log)
  if len(inputs.wnstep.split()) == 2:
    wnunits = inputs.wnstep.split()[1]
  else:
    wnunits = pyrat.spec.wnunits
  wnstep = float(inputs.wnstep.split()[0])
  if inputs.wnstep is None or wnstep <= 0:
    pt.error("Wavenumber sampling step ({:.2e} {:s}-1) must be be > 0.".
             format(wnstep, wnunits), pyrat.log)
  pyrat.spec.wnstep = wnstep / pt.u(wnunits)

  pyrat.spec.wnosamp = pt.getparam(inputs.wnosamp, "none", integer=True)
  isgreater(pyrat.spec.wnosamp, "none", 1, False,
            "Wavenumber oversampling factor ({:d}) must be >= 1.")

  # Check atmospheric layers arguments:
  pyrat.radunits = inputs.radunits
  pyrat.punits   = inputs.punits

  # Pressure boundaries:
  pyrat.phigh = pt.getparam(inputs.phigh, pyrat.punits)
  isgreater(pyrat.phigh, "bar", 0, True,
            "High atm pressure boundary ({:.2e} bar) must be > 0.0")
  pyrat.plow  = pt.getparam(inputs.plow,    pyrat.punits)
  isgreater(pyrat.plow, "bar",  0, True,
            "Low atm pressure boundary ({:.2e} bar) must be > 0.0")
  # Radius boundaries:
  pyrat.radlow  = pt.getparam(inputs.radlow,  pyrat.radunits)
  isgreater(pyrat.radlow, "cm", 0, False,
            "Low atm radius boundary ({:.2e} cm) must be >= 0.0")
  pyrat.radhigh = pt.getparam(inputs.radhigh, pyrat.radunits)
  isgreater(pyrat.radhigh, "cm", 0, True,
            "High atm radius boundary ({:.2e} cm) must be > 0.0")
  pyrat.radstep = pt.getparam(inputs.radstep, pyrat.radunits)
  isgreater(pyrat.radstep, "cm", 0, True,
            "Radius step size ({:.2f} cm) must be > 0.")
  # Pressure-radius reference level:
  pyrat.radiusbase   = pt.getparam(inputs.radiusbase, pyrat.radunits)
  isgreater(pyrat.radiusbase, "cm",   0, True,
            "Planetary radius base ({:.3e} cm) must be > 0.")
  pyrat.pressurebase = pt.getparam(inputs.pressurebase, pyrat.punits)
  isgreater(pyrat.pressurebase, "bar", 0, True,
            "Planetary pressure base ({:8g} bar) must be > 0.")
  pyrat.surfgravity  = pt.getparam(inputs.surfgravity,  "none")
  isgreater(pyrat.surfgravity, "none", 0, True,
            "Planetary surface gravity ({:.2f} cm s-2) must be > 0.")

  pyrat.atm.nlayers = pt.getparam(inputs.nlayers, "none", integer=True)
  isgreater(pyrat.atm.nlayers, "none", 0, True,
            "The number of atmospheric layers ({:d}) must be > 0.")

  # Check Voigt-profile arguments:
  pyrat.voigt.extent = pt.getparam(inputs.vextent, "none")
  isgreater(pyrat.voigt.extent, "none", 1, False,
            "Voigt extent ({:g}) must be >= 1.0")

  # Doppler width:
  pyrat.voigt.nDop = pt.getparam(inputs.nDop, "none", integer=True)
  isgreater(pyrat.voigt.nDop, "none", 1, False,
            "The number of Doppler samples ({:d}) must be >= 1")

  pyrat.voigt.Dmin = pt.getparam(inputs.Dmin, "none")
  isgreater(pyrat.voigt.Dmin, "none", 0, True, "Dmin ({:g} cm-1) must be > 0.")

  pyrat.voigt.Dmax = pt.getparam(inputs.Dmax, "none")
  isgreater(pyrat.voigt.Dmax, "none", 0, True, "Dmax ({:g} cm-1) must be > 0.")

  if pyrat.voigt.Dmax <= pyrat.voigt.Dmin:
    pt.error("Dmax ({:g} cm-1) must be > Dmin ({:g} cm-1).".format(
             pyrat.voigt.Dmax, pyrat.voigt.Dmin), pyrat.log)

  # Lorentz width:
  pyrat.voigt.nLor = pt.getparam(inputs.nLor, "none", integer=True)
  isgreater(pyrat.voigt.nLor, "none", 1, False,
            "The number of Lorentz samples ({:d}) must be >= 1")

  pyrat.voigt.Lmin = pt.getparam(inputs.Lmin, "none")
  isgreater(pyrat.voigt.Lmin, "none", 0, True, "Lmin ({:g} cm-1) must be > 0.")

  pyrat.voigt.Lmax = pt.getparam(inputs.Lmax, "none")
  isgreater(pyrat.voigt.Lmax, "none", 0, True, "Lmax ({:g} cm-1) must be > 0.")

  if pyrat.voigt.Lmax <= pyrat.voigt.Lmin:
    pt.error("Lmax ({:g} cm-1) must be > Lmin ({:g} cm-1).".format(
             pyrat.voigt.Lmax, pyrat.voigt.Lmin), pyrat.log)

  pyrat.voigt.DLratio = pt.getparam(inputs.DLratio, "none")
  isgreater(pyrat.voigt.DLratio, "none", 0, True,
            "Doppler/Lorentz width ratio threshold ({:g}) must be > 0.")

  # Check extinction-coefficient arguments:
  pyrat.ex.ethresh = pt.getparam(inputs.ethresh, "none")
  isgreater(pyrat.ex.ethresh, "none", 0, True,
               "Extinction-coefficient threshold ({:g}) must be positive.")
  if inputs.tmin is not None:
    pyrat.ex.tmin  = pt.getparam(inputs.tmin, "kelvin")
    isgreater(pyrat.ex.tmin,  "kelvin", 0, True,
               "Minimum temperature sample ({:g} K) must be positive.")
  if inputs.tmax is not None:
    pyrat.ex.tmax  = pt.getparam(inputs.tmax, "kelvin")
    isgreater(pyrat.ex.tmax,  "kelvin", 0, True,
               "Maximum temperature sample ({:g} K) must be positive.")
  if inputs.tstep is not None:
    pyrat.ex.tstep = pt.getparam(inputs.tstep, 'kelvin')
    isgreater(pyrat.ex.tstep, "kelvin", 0, True,
               "Temperature sample step interval ({:g} K) must be positive.")

  if pyrat.ex.tmax is not None and pyrat.ex.tmin is not None:
    if pyrat.ex.tmax <= pyrat.ex.tmin:
      pt.error("Maximum temperature limit ({:g} K) must be > minimum "
       "temperature ({:g} K).".format(pyrat.ex.tmax, pyrat.ex.tmin), pyrat.log)

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
  pyrat.od.maxdepth = pt.getparam(inputs.maxdepth, "none")
  isgreater(pyrat.od.maxdepth, "none", 0, False,
            "Maximum optical-depth limit ({:g}) must be >= 0.0")

  # Check system arguments:
  pyrat.rstar = pt.getparam(inputs.rstar, pyrat.radunits)
  isgreater(pyrat.rstar, "cm", 0, True,
            "Stellar radius ({:.4e} cm) must be > 0.")

  # Accept ray-path argument:
  pyrat.od.path  = inputs.path
  # Accept output files:
  pyrat.outspec     = inputs.outspec    
  pyrat.outsample   = inputs.outsample  
  pyrat.outmaxdepth = inputs.outmaxdepth

  # Check raygrid:
  if inputs.raygrid[0] != 0:
    pt.error("First angle in raygrid must be 0.0 (normal to surface).",
             pyrat.log)
  if np.any(inputs.raygrid < 0):
    pt.error("raygrid angles must lie between 0 and 90 deg.", pyrat.log)
  if np.any(np.ediff1d(inputs.raygrid) <= 0):
    pt.error("raygrid angles must be monotonically increasing.", pyrat.log)
  # Store raygrid values in radians:
  pyrat.raygrid = inputs.raygrid * sc.degree

  # Verbose level:
  pyrat.verb = np.amax([0, inputs.verb])
  pt.msg(pyrat.verb, "Done.", pyrat.log)


def isgreater(value, units, thresh, equal=False, text=""):
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
    pt.error(text.format(value/pt.u(units)), pyrat.log, -3)

