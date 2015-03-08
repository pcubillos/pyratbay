import sys, os
import argparse, ConfigParser
import numpy as np

import ptools     as pt
import pconstants as pc

def parse(pyrat):
  """
  Parse the command-line arguments into the pyrat object

  Parameters:
  -----------
  pyrat: Object
     A pyrat object where to store the CLA.

  Modification History:
  ---------------------
  2014-04-26  patricio  Initial implementation.
  2014-06-29  patricio  Added radius/pressure base levels and surface gravity.
  2014-08-15  patricio  Added Voigt-profile section and arguments.
  """

  pt.msg(1, "Processing command-line arguments:")
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
  pt.msg(1, "Configuration file: '{:s}'".format(cfile), 2) 
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
  group.add_argument("--cia",        dest="cia",
                     help="Collision Induced Absorption files [default: "
                          "%(default)s]",
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
                     action="store", type=np.double, default=None)
  group.add_argument("--wlhigh",     dest="wlhigh",
                     help="Wavelength high boundary [default: %(default)s]",
                     action="store", type=np.double, default=None)

  group.add_argument("--wnunits",    dest="wnunits",
                     help="Wavenumber (input) inverse units [default: "
                          "%(default)s]",
                     action="store", type=str, default='cm',
                     choices=('A','nm','um','mm','cm','m'))
  group.add_argument("--wnlow",      dest="wnlow",
                     help="Wavenumber low boundary [default: %(default)s]",
                     action="store", type=np.double, default=None)
  group.add_argument("--wnhigh",     dest="wnhigh",
                     help="Wavenumber high boundary [default: %(default)s]",
                     action="store", type=np.double, default=None)
  group.add_argument("--wnstep",     dest="wnstep",
                     help="Wavenumber sampling step [default: %(default)s]",
                     action="store", type=np.double, default=1.0)
  group.add_argument("--wnosamp",       dest="wnosamp",
                     help="Wavenumber oversampling factor "
                          "[default: %(default)s]",
                     action="store", type=int, default=2160)
  # Atmospheric sampling options:
  group = parser.add_argument_group("Atmosphere Sampling Options")
  group.add_argument("--radlow",     dest="radlow",
                     help="Atmospheric radius low boundary [default: "
                          "Use atmospheric file value]",
                     action="store", type=np.double, default=None)
  group.add_argument("--radhigh",    dest="radhigh",
                     help="Atmospheric radius high boundary [default: "
                          "Use atmospheric file value]",
                     action="store", type=np.double, default=None)
  group.add_argument("--radstep",        dest="radstep",
                     help="Atmospheric radius sampling step [default: "
                          "Use atmospheric file value]",
                     action="store", type=np.double, default=None)
  group.add_argument("--radunits",      dest="radunits",
                     help="Radius (user) units [default: %(default)s]",
                     action="store", type=str, default='km',
                     choices=('cm','m', 'km', 'atmfile'))
  group.add_argument("--plow",          dest="plow",
                     help="Atmospheric pressure low boundary (overrides "
                          "radius high boundary) [default: %(default)s]",
                     action="store", type=np.double, default=None)
  group.add_argument("--phigh",         dest="phigh",
                     help="Atmospheric pressure high boundary (overrides "
                          "radius  low boundary) [default: %(default)s]",
                     action="store", type=np.double, default=None)
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
                     action="store", type=np.double, default=None)
  group.add_argument("--pressurebase",  dest="pressurebase",
                     help="Planetary pressure base level (in punits)",
                     action="store", type=np.double, default=None)
  group.add_argument("--surfgravity",   dest="surfgravity",
                     help="Planet's surface gravity in cm s-2",
                     action="store", type=np.double, default=None)
  # Extinction options:
  group = parser.add_argument_group("Extinction Calculations Options")
  group.add_argument("--tmin",          dest="tmin",
                     help="Minimum temperature to sample/consider "
                     " in Kelvin [default: %(default)s]",
                     action="store", type=np.double, default=500.0)
  group.add_argument("--tmax",          dest="tmax",
                     help="Maximum temperature to sample/consider "
                     "in Kelvin [default: %(default)s]",
                     action="store", type=np.double, default=3000.0)
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
  group.add_argument(      "--voigtwidth",    dest="voigtwidth",
                     help="Width of Voigt profile in number of max(Doppler "
                          "width, Lorentz width) [default: %(default)s]",
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
  # Optical depth options:
  group = parser.add_argument_group("Optical Depth Options")
  group.add_argument("-p", "--path",          dest="path",
                     help="Lightray-path geometry [default: %(default)s]",
                     action="store", type=str, default='transit',
                     choices=('transit', 'eclipse'))
  group.add_argument("--maxdepth",       dest="maxdepth",
                     help="Maximum optical depth to calculate [default: "
                          "%(default)s]",
                     action="store", type=np.double, default=10)
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
  pyrat.inputs.cia        = user.cia
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
  pyrat.inputs.voigtwidth = user.voigtwidth
  pyrat.inputs.Dmin       = user.Dmin
  pyrat.inputs.Dmax       = user.Dmax
  pyrat.inputs.nDop       = user.nDop
  pyrat.inputs.Lmin       = user.Lmin
  pyrat.inputs.Lmax       = user.Lmax
  pyrat.inputs.nLor       = user.nLor
  pyrat.inputs.DLratio    = user.DLratio
  # Optical depth:
  pyrat.inputs.path       = user.path
  pyrat.inputs.maxdepth   = user.maxdepth
  # Output files:
  pyrat.inputs.outspec    = user.outspec
  pyrat.inputs.outsample  = user.outsample
  pyrat.inputs.outmaxdepth = user.outmaxdepth


def checkinputs(pyrat):
  """
  Check that user input arguments make sense.

  Modification History:
  ---------------------
  2014-04-26  patricio  Initial python implementation.
  2014-06-29  patricio  Added radius/pressure base levels and surface gravity.
  2014-08-15  patricio  Added Voigt variables check. Put extinction variables 
                        in pyrat.ex object.
  """
  # User-inputs object:
  inputs = pyrat.inputs

  # Path to source parent's folder:
  pyratdir = os.path.dirname(os.path.realpath(__file__))

  # Check that input files exist:
  if not os.path.isfile(inputs.atmfile):
    pt.error("atmfile: '{:s}' does not exist.".format(inputs.atmfile))
  pyrat.atmfile = inputs.atmfile

  if inputs.linedb is not None:
    for linedb in inputs.linedb:
      if not os.path.isfile(linedb):
        pt.error("linedb file: '{:s}' does not exist.".format(linedb))
  pyrat.linedb = pyrat.inputs.linedb

  if inputs.cia is not None:
    for cia in pyrat.inputs.cia:
      if not os.path.isfile(cia):
        pt.error("CIA file: '{:s}' does not exist.".format(cia))
  pyrat.cia.files = pyrat.inputs.cia

  if inputs.molfile is None: # Set default
    inputs.molfile = pyratdir + "/../inputs/molecules.dat"
  if not os.path.isfile(inputs.molfile):
    pt.error("Molecular-data file: '{:s}' does not exist.".
             format(inputs.molfile))
  pyrat.molfile = inputs.molfile

  if inputs.extfile is not None:
    if not os.path.exists(os.path.realpath(os.path.dirname(inputs.extfile))):
      pt.error("Directory for extinction-coefficient file '{:s}' does "
               "not exist.".format(inputs.extfile))
  pyrat.ex.extfile = os.path.realpath(inputs.extfile)

  # Check spectrum arguments:
  pyrat.wnunits = inputs.wnunits  # Accept units
  pyrat.wlunits = inputs.wlunits

  pyrat.wllow  = isgreater(inputs.wllow,  pyrat.wlunits,   0, False,
                  "Low wavelength boundary ({:.2e} {:s}) must be >= 0.")

  pyrat.wlhigh = isgreater(inputs.wlhigh, pyrat.wlunits,   0, True,
                  "Low wavelength boundary ({:.2e} {:s}) must be >= 0.")

  if   inputs.wnlow is not None:
    if inputs.wnlow < 0.0:
      pt.error("Low wavenumber boundary ({:.2e} {:s}-1) must be "
               ">= 0.".format(inputs.wnlow, inputs.wnunits))
    pyrat.wnlow = inputs.wnlow / pc.units[pyrat.wnunits]

  if   inputs.wnhigh is not None:
    if inputs.wnhigh <= 0.0:
      pt.error("High wavenumber boundary ({:.2e} {:s}-1) must be "
               "> 0.".format(inputs.wnhigh, inputs.wnunits))
    pyrat.wnhigh = inputs.wnhigh / pc.units[pyrat.wnunits]

  if inputs.wnstep is None or inputs.wnstep <= 0:
    pt.error("Wavenumber sampling step ({:.2e} {:s}-1) must be defined and "
             "be > 0.".format(inputs.wnstep, inputs.wnunits))
  pyrat.wnstep = inputs.wnstep / pc.units[pyrat.wnunits]

  pyrat.wnosamp = isgreater(inputs.wnosamp, None, 1, False,
                     "Wavenumber oversampling factor ({:d}) must be >= 1.")

  # Check atmospheric layers arguments:
  pyrat.radunits = inputs.radunits
  pyrat.punits   = inputs.punits

  pyrat.phigh   = isgreater(inputs.phigh,   pyrat.punits,   0, True,
                     "High atm pressure boundary ({:.2e} {:s}) must be > 0.0")
  pyrat.plow    = isgreater(inputs.plow,    pyrat.punits,   0, True,
                     "Low atm pressure boundary ({:.2e} {:s}) must be > 0.0")

  pyrat.radlow  = isgreater(inputs.radlow,  pyrat.radunits, 0, False,
                     "Low atm radius boundary ({:.2e} {:s}) must be >= 0.0")
  pyrat.radhigh = isgreater(inputs.radhigh, pyrat.radunits, 0, True,
                     "High atm radius boundary ({:.2e} {:s}) must be > 0.0")
  pyrat.radstep = isgreater(inputs.radstep, pyrat.radunits, 0, True,
                     "Radius step size ({:.2f} {:s}) must be > 0.")

  pyrat.radiusbase   = isgreater(inputs.radiusbase,   pyrat.radunits, 0, True,
                     "Planetary radius base ({:.3e} {:s}) must be > 0.")
  pyrat.pressurebase = isgreater(inputs.pressurebase, pyrat.punits,   0, True,
                     "Planetary pressure base ({:8g} {:s}) must be > 0.")
  pyrat.surfgravity  = isgreater(inputs.surfgravity,  None,           0, True,
                     "Planetary surface gravity ({:.2f} cm s-2) must be > 0.")

  pyrat.atm.nlayers = isgreater(inputs.nlayers, None, 0, True,
                     "The number of atmospheric layers ({:d}) must be > 0.")

  # Check Voigt-profile arguments:
  pyrat.voigt.width = isgreater(inputs.voigtwidth, None, 1, False,
                                "Voigt width ({:g}) must be >= 1.0")

  # Doppler width:
  pyrat.voigt.nDop = isgreater(inputs.nDop, None, 1, False,
                         "The number of Doppler samples ({:d}) must be >= 1")

  pyrat.voigt.Dmin = isgreater(inputs.Dmin, None, 0, True,
                                   "Dmin ({:g} cm-1) must be > 0.")

  pyrat.voigt.Dmax = isgreater(inputs.Dmax, None, 0, True,
                                   "Dmax ({:g} cm-1) must be > 0.")

  if pyrat.voigt.Dmax <= pyrat.voigt.Dmin:
    pt.error("Dmax ({:g} cm-1) must be > Dmin ({:g} cm-1).".format(
             pyrat.voigt.Dmax, pyrat.voigt.Dmin))

  # Lorentz width:
  pyrat.voigt.nLor = isgreater(inputs.nLor, None, 1, False,
             "The number of Lorentz samples ({:d}) must be >= 1")

  pyrat.voigt.Lmin = isgreater(inputs.Lmin, None, 0, True,
             "Lmin ({:g} cm-1) must be > 0.")

  pyrat.voigt.Lmax = isgreater(inputs.Lmax, None, 0, True,
             "Lmax ({:g} cm-1) must be > 0.")

  if pyrat.voigt.Lmax <= pyrat.voigt.Lmin:
    pt.error("Lmax ({:g} cm-1) must be > Lmin ({:g} cm-1).".format(
             pyrat.voigt.Lmax, pyrat.voigt.Lmin))

  pyrat.voigt.DLratio = isgreater(inputs.DLratio, None, 0, True,
             "Doppler/Lorentz width ratio threshold ({:g}) must be > 0.")

  # Check extinction-coefficient arguments:
  pyrat.ex.ethresh = isgreater(inputs.ethresh,  None, 0, True,
             "Extinction-coefficient threshold ({:g} K) must be positive.")
  pyrat.ex.tmin  = isgreater(inputs.tmin,  None, 0, True,
             "Minimum temperature sample ({:g} K) must be positive.")
  pyrat.ex.tmax  = isgreater(inputs.tmax,  None, 0, True,
             "Maximum temperature sample ({:g} K) must be positive.")
  pyrat.ex.tstep = isgreater(inputs.tstep, None, 0, True,
             "Temperature sample step interval ({:g} K) must be positive.")

  if pyrat.ex.tmax is not None and pyrat.ex.tmin is not None:
    if pyrat.ex.tmax <= pyrat.ex.tmin:
      pt.error("Maximum temperature limit ({:g} K) must be > minimum "
               "temperature ({:g} K).".format(pyrat.ex.tmax, pyrat.ex.tmin))

  # Check optical-depth arguments:
  pyrat.maxdepth = isgreater(inputs.maxdepth, None, 0, False,
                        "Maximum-optical-depth limit ({:g}) must be >= 0.0")

  # Accept ray-path argument:
  pyrat.path  = inputs.path
  # Accept output files:
  pyrat.outspec     = inputs.outspec    
  pyrat.outsample   = inputs.outsample  
  pyrat.outmaxdepth = inputs.outmaxdepth

  # Verbose level:
  pyrat.verb = np.amax([0, inputs.verb])
  pt.msg(pyrat.verb, "Done.", 0)

  # FINDME: set system geometry variables


def isgreater(value, units, thresh, equal=False, text=""):
  """
  Check that value (if not None) is greater than thresh.
  Throw error if not.

  Parameters:
  -----------
  value: Scalar
    The value being tested.
  units: String
    The units of value as given in pyrat.units
  thresh: Scalar
    Threshold against which value is being compared.
  equal: Boolean
    If True, strictly require greater than.
  text: String
    Text to show if condition is not satisfied.

  Returns:
  --------
  The value in the pyrat units.

  Modification History:
  ---------------------
  2015-01-18  patricio  initial implementation.
  """
  # Set comparison command:
  if equal:
    compare = np.less_equal
  else:
    compare = np.less

  # Get the units factor value:
  if units is None:
    unitsval = 1
  else:
    unitsval = pc.units[units]

  # Check value:
  if value is not None:
    if compare(value, thresh):
      pt.error(text.format(value, units), -3)
    return value * unitsval
  return None
