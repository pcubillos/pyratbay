import sys, os
import argparse, ConfigParser
import numpy as np
import ptools as pt

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
  group.add_argument("-a", "--atmfile",    dest="atmfile",
                     help="Atmospheric file [default: %(default)s]",
                     action="store", type=str, default=None) 
  group.add_argument("-l", "--linedb",     dest="linedb",
                     help="Line database files [default: %(default)s]",
                     action="store", type=pt.parray, default=None) 
  group.add_argument(      "--cia",        dest="cia",
                     help="Collision Induced Absorption files [default: "
                          "%(default)s]",
                     action="store", type=pt.parray, default=None)
  # Spectrum sampling options:
  group = parser.add_argument_group("Spectrum Sampling Options")
  group.add_argument("-i", "--wllow",      dest="wllow",
                     help="Wavelength low boundary [default: %(default)s]",
                     action="store", type=np.double, default=1.0)
  group.add_argument("-f", "--wlhigh",     dest="wlhigh",
                     help="Wavelength high boundary [default: %(default)s]",
                     action="store", type=np.double, default=10.0)
  group.add_argument("-s", "--wlstep",     dest="wlstep",
                     help="Wavelength sampling step (the wavenumber array will in fact be equispaced in wavenumber, with a number of samples equal to the number of samples of the produced by this step in the wavelength range) [default: %(default)s]",
                     action="store", type=np.double, default=0.001)
  group.add_argument("-u", "--wlunits",    dest="wlunits",
                     help="Wavelength (user) units [default: %(default)s]",
                     action="store", type=str, default='cm',
                     choices=('A','nm','um','mm','cm','m'))
  group.add_argument(      "--wnstep",     dest="wnstep",
                     help="Wavenumber sampling step (overrides wlstep) "
                          "[default: %(default)s]",
                     action="store", type=np.double, default=None)
  group.add_argument(      "--wnunits",    dest="wnunits",
                     help="Wavenumber (user) inverse units [default: "
                          "%(default)s]",
                     action="store", type=str, default='cm',
                     choices=('A','nm','um','mm','cm','m'))
  group.add_argument(     "--nspec",       dest="nspec",
                     help="Number of spectral samples [default: %(default)s]",
                     action="store", type=int, default=None)
  # Atmospheric sampling options:
  group = parser.add_argument_group("Atmosphere Sampling Options")
  group.add_argument(      "--radlow",     dest="radlow",
                     help="Atmospheric radius low boundary [default: "
                          "Use atmospheric file value]",
                     action="store", type=np.double, default=None)
  group.add_argument(      "--radhigh",    dest="radhigh",
                     help="Atmospheric radius high boundary [default: "
                          "Use atmospheric file value]",
                     action="store", type=np.double, default=None)
  group.add_argument(      "--radstep",        dest="radstep",
                     help="Atmospheric radius sampling step [default: "
                          "Use atmospheric file value]",
                     action="store", type=np.double, default=None)
  group.add_argument(      "--radunits",      dest="radunits",
                     help="Radius (user) units [default: %(default)s]",
                     action="store", type=str, default='km',
                     choices=('cm','m', 'km', 'atmfile'))
  group.add_argument(      "--plow",          dest="plow",
                     help="Atmospheric pressure low boundary (overrides "
                          "radius high boundary) [default: %(default)s]",
                     action="store", type=np.double, default=None)
  group.add_argument(      "--phigh",         dest="phigh",
                     help="Atmospheric pressure high boundary (overrides "
                          "radius  low boundary) [default: %(default)s]",
                     action="store", type=np.double, default=None)
  group.add_argument(      "--nlayers",       dest="nlayers",
                     help="Number of atmospheric pressure samples [default: "
                          "Use atmospheric file value]",
                     action="store", type=np.int, default=-1)
  group.add_argument(      "--punits",        dest="punits",
                     help="Pressure (user) units [default: %(default)s]",
                     action="store", type=str, default='bar',
                     choices=('bar',))
  group.add_argument(      "--radiusbase",    dest="radiusbase",
                     help="Planetary radius base level (in radunits)",
                     action="store", type=np.double, default=None)
  group.add_argument(      "--pressurebase",  dest="pressurebase",
                     help="Planetary pressure base level (in punits)",
                     action="store", type=np.double, default=None)
  group.add_argument(      "--surfgravity",   dest="surfgravity",
                     help="Planet's surface gravity in cm/s^2",
                     action="store", type=np.double, default=None)
  # Extinction options:
  group = parser.add_argument_group("Extinction Calculations Options")
  group.add_argument(      "--voigtbin",      dest="voigtbin",
                     help="Oversampling bin factor for Voigt profile "
                          "radius low boundary) [default: %(default)s]",
                     action="store", type=int, default=5)
  group.add_argument(      "--voigtwidth",    dest="voigtwidth",
                     help="Width of Voigt profile in number of max(Doppler "
                          "width, Lorentz width) [default: %(default)s]",
                     action="store", type=np.double, default=50)
  group.add_argument(      "--Tmin",          dest="tmin",
                     help="Minimum temperature to sample/consider "
                     " in Kelvin [default: %(default)s]",
                     action="store", type=np.double, default=500.0)
  group.add_argument(      "--Tmax",          dest="tmax",
                     help="Maximum temperature to sample/consider "
                     "in Kelvin [default: %(default)s]",
                     action="store", type=np.double, default=3000.0)
  group.add_argument(      "--minelow",       dest="minelow",
                     help="Minimum Elow to consider line transition "
                          "[default: %(default)s]",
                     action="store", type=np.double, default=0)
  # Voigt-profile options:
  group = parser.add_argument_group("Voigt-profile  Options")
  group.add_argument(      "--Dmin",          dest="Dmin",
                     help="Minimum Doppler width to sample in cm-1 "
                          "[default: %(default)s]",
                     action="store", type=np.double, default=None)
  group.add_argument(      "--Dmax",          dest="Dmax",
                     help="Maximum Doppler width to sample in cm-1 "
                          "[default: %(default)s]",
                     action="store", type=np.double, default=None)
  group.add_argument(      "--nDop",          dest="nDop",
                     help="Number of Doppler-width samples"
                          "[default: %(default)s]",
                     action="store", type=np.int, default=100)
  group.add_argument(      "--Lmin",          dest="Lmin",
                     help="Minimum Lorentz width to sample in cm-1 "
                          "[default: %(default)s]",
                     action="store", type=np.double, default=None)
  group.add_argument(      "--Lmax",          dest="Lmax",
                     help="Maximum Lorentz width to sample in cm-1 "
                          "[default: %(default)s]",
                     action="store", type=np.double, default=None)
  group.add_argument(      "--nLor",          dest="nLor",
                     help="Number of Lorentz-width samples"
                          "[default: %(default)s]",
                     action="store", type=np.int, default=100)
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
  group.add_argument("-t", "--toomuch",       dest="toomuch",
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
  group.add_argument(      "--outtoomuch",    dest="outtoomuch",
                     help="Output toomuch-level file [default: %(default)s]",
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
  # Wavelength:
  pyrat.inputs.wllow      = user.wllow
  pyrat.inputs.wlhigh     = user.wlhigh
  pyrat.inputs.wlstep     = user.wlstep
  pyrat.inputs.wlunits    = user.wlunits
  # Wavenumber:
  #pyrat.inputs.wnlow      = user.wnlow
  #pyrat.inputs.wnhigh     = user.wnhigh
  pyrat.inputs.wnstep     = user.wnstep
  pyrat.inputs.wnunits    = user.wnunits
  pyrat.inputs.nspec      = user.nspec
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
  pyrat.inputs.voigtbin   = user.voigtbin
  pyrat.inputs.voigtwidth = user.voigtwidth
  pyrat.inputs.minelow    = user.minelow
  pyrat.inputs.tmin       = user.tmin
  pyrat.inputs.tmax       = user.tmax
  # Voigt-profile:
  pyrat.inputs.Dmin       = user.Dmin
  pyrat.inputs.Dmax       = user.Dmax
  pyrat.inputs.nDop       = user.nDop
  pyrat.inputs.Lmin       = user.Lmin
  pyrat.inputs.Lmax       = user.Lmax
  pyrat.inputs.nLor       = user.nLor
  pyrat.inputs.DLratio    = user.DLratio
  # Optical depth:
  pyrat.inputs.path       = user.path
  pyrat.inputs.toomuch    = user.toomuch
  # Output files:
  pyrat.inputs.outspec    = user.outspec
  pyrat.inputs.outsample  = user.outsample
  pyrat.inputs.outtoomuch = user.outtoomuch


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

  # Check that input files exist:
  if not os.path.isfile(pyrat.inputs.atmfile):
    pt.error("atmfile: '{:s}' does not exist.".format(pyrat.inputs.atmfile))
  pyrat.atmfile = pyrat.inputs.atmfile

  if pyrat.inputs.linedb is not None:
    for linedb in pyrat.inputs.linedb:
      if not os.path.isfile(linedb):
        pt.error("linedb file: '{:s}' does not exist.".format(linedb))
  pyrat.linedb = pyrat.inputs.linedb

  if pyrat.inputs.cia is not None:
    for cia in pyrat.inputs.cia:
      if not os.path.isfile(cia):
        pt.error("CIA file: '{:s}' does not exist.".format(cia))
  pyrat.cia.files = pyrat.inputs.cia

  # Check Voigt-profile arguments:
  if pyrat.inputs.voigtbin < 1:
    pt.error("Voigt bin oversampling ({:d}) factor must be >= 1".format(
                                                       pyrat.inputs.voigtbin))
  pyrat.voigt.osamp = pyrat.inputs.voigtbin

  if pyrat.inputs.voigtwidth < 1:
    pt.error("Voigt width ({:g}) must be >= 1.0".format(
                                                     pyrat.inputs.voigtwidth))
  pyrat.voigt.width = pyrat.inputs.voigtwidth

  # Doppler width:
  if pyrat.inputs.Dmin is not None and pyrat.inputs.Dmin <= 0:
    pt.error("The minimum Doppler width ({:g} cm-1) to sample must be "
             "positive.".format(pyrat.inputs.Dmin))
  pyrat.voigt.Dmin = pyrat.inputs.Dmin

  if pyrat.inputs.Dmax is not None and pyrat.inputs.Dmax <= 0:
    pt.error("The maximum Doppler width ({:g} cm-1) to sample must be "
             "positive.".format(pyrat.inputs.Dmax))
  pyrat.voigt.Dmax = pyrat.inputs.Dmax

  if pyrat.voigt.Dmax is not None and pyrat.voigt.Dmin is not None:
    if pyrat.voigt.Dmax <= pyrat.voigt.Dmin:
      pt.error("Maximum Doppler width ({:g} cm-1) must be > minimum Doppler "
               "width ({:g} cm-1).".format(pyrat.voigt.Dmax, pyrat.voigt.Dmin))

  if pyrat.inputs.nDop < 1:
    pt.error("The number of Doppler samples ({:d}) must be "
             ">= 1".format(pyrat.inputs.nDop))
  pyrat.voigt.nDop = pyrat.inputs.nDop

  # Lorentz width:
  if pyrat.inputs.Lmin is not None and pyrat.inputs.Lmin <= 0:
    pt.error("the minimum Lorentz width ({:g} cm-1) to sample must be "
             "positive.".format(pyrat.inputs.Lmin))
  pyrat.voigt.Lmin = pyrat.inputs.Lmin

  if pyrat.inputs.Lmax is not None and pyrat.inputs.Lmax <= 0:
    pt.error("The maximum Lorentz width ({:g} cm-1) to sample must be "
             "positive.".format(pyrat.inputs.Lmax))
  pyrat.voigt.Lmax = pyrat.inputs.Lmax

  if pyrat.voigt.Lmax is not None and pyrat.voigt.Lmin is not None:
    if pyrat.voigt.Lmax <= pyrat.voigt.Lmin:
      pt.error("Maximum Lorentz width ({:g} cm-1) must be > minimum Lorentz "
               "width ({:g} cm-1).".format(pyrat.voigt.Lmax, pyrat.voigt.Lmin))

  if pyrat.inputs.nLor < 1:
    pt.error("The number of Lorentz samples ({:d}) must be "
             ">= 1".format(pyrat.inputs.nLor))
  pyrat.voigt.nLor = pyrat.inputs.nLor

  if pyrat.inputs.DLratio <= 0:
    pt.error("Lorentz/Doppler width ratio threshold ({:g}) must be "
             "positive.".format(pyrat.inputs.DLratio))
  pyrat.voigt.DLratio = pyrat.inputs.DLratio

  # Check extinction arguments:
  if pyrat.inputs.minelow < 0.0:
    pt.error("Minimum Elow ({:g}) must be >= 0.0".format(pyrat.inputs.minelow))
  pyrat.ex.minelow = pyrat.inputs.minelow

  if pyrat.inputs.tmin is not None and pyrat.inputs.tmin <= 0.0:
    pt.error("Minimum temperature sample ({:g} K) must be positive.".format(
                                                            pyrat.inputs.tmin))
  pyrat.ex.tmin = pyrat.inputs.tmin

  if pyrat.inputs.tmax is not None and pyrat.inputs.tmax <= 0.0:
    pt.error("Maximum temperature sample ({:g} K) must be positive.".format(
                                                            pyrat.inputs.tmax))
  pyrat.ex.tmax = pyrat.inputs.tmax

  if pyrat.ex.tmax is not None and pyrat.ex.tmin is not None:
    if pyrat.ex.tmax <= pyrat.ex.tmin:
      pt.error("Maximum temperature limit ({:g} K) must be > minimum "
               "temperature ({:g} K).".format(pyrat.ex.tmax, pyrat.ex.tmin))

  # Check opacity arguments:
  if pyrat.inputs.toomuch < 0.0:
    pt.error("Tau max limit ({:g}) must be >= 0.0".format(pyrat.inputs.toomuch))
  pyrat.toomuch = pyrat.inputs.toomuch

  # Accept ray-path argument:
  pyrat.path  = pyrat.inputs.path
  # Accept output files:
  pyrat.outspec    = pyrat.inputs.outspec    
  pyrat.outsample  = pyrat.inputs.outsample  
  pyrat.outtoomuch = pyrat.inputs.outtoomuch 

  # Verbose level:
  pyrat.verb = np.amax([0, pyrat.inputs.verb])
  pt.msg(pyrat.verb, "Done.", 0)

  # FINDME: set system geometry variables
