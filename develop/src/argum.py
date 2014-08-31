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
  pyrat.user.configfile = args.configfile
  pyrat.user.verb       = user.verb
  # Input file:
  pyrat.user.atmfile    = user.atmfile
  pyrat.user.linedb     = user.linedb
  pyrat.user.cia        = user.cia
  # Wavelength:
  pyrat.user.wllow      = user.wllow
  pyrat.user.wlhigh     = user.wlhigh
  pyrat.user.wlstep     = user.wlstep
  pyrat.user.wlunits    = user.wlunits
  # Wavenumber:
  #pyrat.user.wnlow      = user.wnlow
  #pyrat.user.wnhigh     = user.wnhigh
  pyrat.user.wnstep     = user.wnstep
  pyrat.user.wnunits    = user.wnunits
  pyrat.user.nspec      = user.nspec
  # Atmospheric radius:
  pyrat.user.radlow     = user.radlow
  pyrat.user.radhigh    = user.radhigh
  pyrat.user.radstep    = user.radstep
  pyrat.user.radunits   = user.radunits
  # Atmospheric pressure:
  pyrat.user.plow       = user.plow
  pyrat.user.phigh      = user.phigh
  pyrat.user.punits     = user.punits
  pyrat.user.nlayers    = user.nlayers
  # Hydrostatic-equilibrium base-level variables:
  pyrat.user.radiusbase   = user.radiusbase
  pyrat.user.pressurebase = user.pressurebase
  pyrat.user.surfgravity  = user.surfgravity
  # Extinction:
  pyrat.user.voigtbin   = user.voigtbin
  pyrat.user.voigtwidth = user.voigtwidth
  pyrat.user.minelow    = user.minelow
  pyrat.user.tmin       = user.tmin
  pyrat.user.tmax       = user.tmax
  # Voigt-profile:
  pyrat.user.Dmin       = user.Dmin
  pyrat.user.Dmax       = user.Dmax
  pyrat.user.nDop       = user.nDop
  pyrat.user.Lmin       = user.Lmin
  pyrat.user.Lmax       = user.Lmax
  pyrat.user.nLor       = user.nLor
  pyrat.user.DLratio    = user.DLratio
  # Optical depth:
  pyrat.user.path       = user.path
  pyrat.user.toomuch    = user.toomuch
  # Output files:
  pyrat.user.outspec    = user.outspec
  pyrat.user.outsample  = user.outsample
  pyrat.user.outtoomuch = user.outtoomuch


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
  if not os.path.isfile(pyrat.user.atmfile):
    pt.error("atmfile: '{:s}' does not exist.".format(pyrat.user.atmfile))
  pyrat.atmfile = pyrat.user.atmfile

  if pyrat.user.linedb is not None:
    for linedb in pyrat.user.linedb:
      if not os.path.isfile(linedb):
        pt.error("linedb file: '{:s}' does not exist.".format(linedb))
  pyrat.linedb = pyrat.user.linedb

  if pyrat.user.cia is not None:
    for cia in pyrat.user.cia:
      if not os.path.isfile(cia):
        pt.error("CIA file: '{:s}' does not exist.".format(cia))
  pyrat.cia = pyrat.user.cia

  # Check Voigt-profile arguments:
  if pyrat.user.voigtbin < 1:
    pt.error("Voigt bin oversampling ({:d}) factor must be >= 1".format(
                                                           pyrat.user.voigtbin))
  pyrat.voigt.osamp = pyrat.user.voigtbin

  if pyrat.user.voigtwidth < 1:
    pt.error("Voigt width ({:g}) must be >= 1.0".format(pyrat.user.voigtwidth))
  pyrat.voigt.width = pyrat.user.voigtwidth

  # Doppler width:
  if pyrat.user.Dmin is not None and pyrat.user.Dmin <= 0:
    pt.error("The minimum Doppler width ({:g} cm-1) to sample must be "
             "positive.".format(pyrat.user.Dmin))
  pyrat.voigt.Dmin = pyrat.user.Dmin

  if pyrat.user.Dmax is not None and pyrat.user.Dmax <= 0:
    pt.error("The maximum Doppler width ({:g} cm-1) to sample must be "
             "positive.".format(pyrat.user.Dmax))
  pyrat.voigt.Dmax = pyrat.user.Dmax

  if pyrat.voigt.Dmax is not None and pyrat.voigt.Dmin is not None:
    if pyrat.voigt.Dmax <= pyrat.voigt.Dmin:
      pt.error("Maximum Doppler width ({:g} cm-1) must be > minimum Doppler "
               "width ({:g} cm-1).".format(pyrat.voigt.Dmax, pyrat.voigt.Dmin))

  if pyrat.user.nDop < 1:
    pt.error("The number of Doppler samples ({:d}) must be "
             ">= 1".format(pyrat.user.nDop))
  pyrat.voigt.nDop = pyrat.user.nDop

  # Lorentz width:
  if pyrat.user.Lmin is not None and pyrat.user.Lmin <= 0:
    pt.error("the minimum Lorentz width ({:g} cm-1) to sample must be "
             "positive.".format(pyrat.user.Lmin))
  pyrat.voigt.Lmin = pyrat.user.Lmin

  if pyrat.user.Lmax is not None and pyrat.user.Lmax <= 0:
    pt.error("The maximum Lorentz width ({:g} cm-1) to sample must be "
             "positive.".format(pyrat.user.Lmax))
  pyrat.voigt.Lmax = pyrat.user.Lmax

  if pyrat.voigt.Lmax is not None and pyrat.voigt.Lmin is not None:
    if pyrat.voigt.Lmax <= pyrat.voigt.Lmin:
      pt.error("Maximum Lorentz width ({:g} cm-1) must be > minimum Lorentz "
               "width ({:g} cm-1).".format(pyrat.voigt.Lmax, pyrat.voigt.Lmin))

  if pyrat.user.nLor < 1:
    pt.error("The number of Lorentz samples ({:d}) must be "
             ">= 1".format(pyrat.user.nLor))
  pyrat.voigt.nLor = pyrat.user.nLor

  if pyrat.user.DLratio <= 0:
    pt.error("Lorentz/Doppler width ratio threshold ({:g}) must be "
             "positive.".format(pyrat.user.DLratio))
  pyrat.voigt.DLratio = pyrat.user.DLratio

  # Check extinction arguments:
  if pyrat.user.minelow < 0.0:
    pt.error("Minimum Elow ({:g}) must be >= 0.0".format(pyrat.user.minelow))
  pyrat.ex.minelow = pyrat.user.minelow

  if pyrat.user.tmin is not None and pyrat.user.tmin <= 0.0:
    pt.error("Minimum temperature sample ({:g} K) must be positive.".format(
                                                            pyrat.user.tmin))
  pyrat.ex.tmin = pyrat.user.tmin

  if pyrat.user.tmax is not None and pyrat.user.tmax <= 0.0:
    pt.error("Maximum temperature sample ({:g} K) must be positive.".format(
                                                            pyrat.user.tmax))
  pyrat.ex.tmax = pyrat.user.tmax

  if pyrat.ex.tmax is not None and pyrat.ex.tmin is not None:
    if pyrat.ex.tmax <= pyrat.ex.tmin:
      pt.error("Maximum temperature limit ({:g} K) must be > minimum "
               "temperature ({:g} K).".format(pyrat.ex.tmax, pyrat.ex.tmin))

  # Check opacity arguments:
  if pyrat.user.toomuch < 0.0:
    pt.error("Too-much limit ({:g}) must be >= 0.0".format(pyrat.user.toomuch))
  pyrat.toomuch = pyrat.user.toomuch

  # Accept ray-path argument:
  pyrat.path  = pyrat.user.path
  # Accept output files:
  pyrat.outspec    = pyrat.user.outspec    
  pyrat.outsample  = pyrat.user.outsample  
  pyrat.outtoomuch = pyrat.user.outtoomuch 

  # Verbose level:
  pyrat.verb = np.amax([0, pyrat.user.verb])
  pt.msg(pyrat.verb, "Done.", 0)

  # FINDME: set system geometry variables
