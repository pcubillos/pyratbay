import numpy as np

class pyrat(object):
  """
  Main PyRaT object.
  """
  def __init__(self):
    # User inputs:
    self.inputs = inputs()         # User inputs
    self.atmf   = atm()            # Input-file atmosphere
    self.atm    = atm()            # Modeling atmosphere
    self.lt     = linetransition() # Line-transition data
    self.mol    = molecules()      # Molecules data
    self.iso    = isotopes()       # Isotopes data
    self.voigt  = voigt()          # Voigt profile
    self.ex     = extinction()     # Extinction
    self.cia    = cia()            # CIA
    self.od     = optdepth()       # Optical depth

    self.wnlow   = None
    self.wnhigh  = None

class inputs(object):
  """
  This is a holder class to store user-input arguments.

  Modification History:
  ---------------------
  2014-04-26  patricio  Initial implementation.
  2014-06-29  patricio  Added radius/pressure base levels and surface gravity.
  2014-08-31  patricio  Added CIA variables.
  """
  def __init__(self):
    # General arguments:
    self.configfile = None
    self.verb       = None
    # Input/output files arguments:
    self.atmfile  = None
    self.linedb   = None
    self.ciafiles = None
    self.molfile  = None
    self.extfile  = None
    # Wavelength arguments:
    self.wlunits = None
    self.wllow   = None
    self.wlhigh  = None
    self.wlstep  = None
    # Wavenumber arguments:
    self.wnunits = None
    self.wnlow   = None
    self.wnhigh  = None
    self.wnstep  = None
    self.wnosamp = None
    # Atmospheric radius arguments:
    self.radlow   = None
    self.radhigh  = None 
    self.radstep  = None
    self.radunits = None
    # Atmospheric pressure arguments:
    self.plow   = None 
    self.phigh  = None
    self.pstep  = None
    self.punits = None
    # Base radius-pressure level:
    self.zeroradius  = None
    self.zerpress    = None
    self.surfgravity = None
    # Voigt profile arguments:
    self.voigtwidth = None
    self.DLratio    = None
    self.Dmin       = None
    self.Dmax       = None
    self.nDop       = None
    self.Lmin       = None
    self.Lmax       = None
    self.nLor       = None
    # Extinction calculation arguments:
    self.exthresh = None
    self.tmin     = None
    self.tmax     = None
    # Optical depth arguments:
    self.path     = None
    self.maxdepth = None
    # System arguments:
    self.rstar = None
    # Output files arguments:
    self.outspec     = None
    self.outsample   = None
    self.outmaxdepth = None


class atm(object):
  def __init__(self):
    self.abundance = None      # Abundance by mass (True) or number (False)
    self.info      = None      # General info from atmfile
    self.runits    = 'km'      # Input radius units
    self.punits    = 'bar'     # Input pressure units
    self.tunits    = 'kelvin'  # Input temperature units
    self.nlayers   = None      # Number of layers
    self.radius    = None      # Radius array (cm)            [layers]
    self.press     = None      # Pressure array (barye)       [layers]
    self.temp      = None      # Temperature array (K)        [layers]
    self.mm        = None      # Mean molecular mass (gr/mol) [layers]
    self.q         = None      # Molecular abundances         [layers, nmol]
    self.d         = None      # Molecular densities          [layers, nmol]


class molecules(object):
  def __init__(self):
    self.nmol   = 0     # Number of molecules
    self.name   = None  # Molecule's name               [nmol]
    self.symbol = None  # Molecule's symbol             [nmol]
    self.mass   = None  # Molecule's mass  (gr/mol)     [nmol]
    self.radius = None  # Molecule's radius (Angstroms) [nmol]
    self.ID     = None  # Molecule's universal ID       [nmol]


class linetransition(object):
  def __init__(self):
    self.nTLI    = 0      # Number of TLI files
    self.ndb     = 0      # Number of data bases
    self.db      = []     # Data base objects
    self.ntransitions = 0 # Number of line transitions
    self.wn      = np.array([]) # Line wavenumber
    self.elow    = np.array([]) # Line lower energy level
    self.gf      = np.array([]) # Line gf value
    self.isoid   = np.array([], np.int) # Line isotope index


class database(object):
  def __init__(self):
    self.name    = None  # Data base name
    self.molname = None  # Molecule name
    self.niso    = None  # Number of isotopes in database
    self.iiso    = None  # Isotope correlative index
    self.ntemp   = None  # Number of temperature samples
    self.temp    = None  # Temperature array
    self.z       = None  # Isotopes' partition function array [niso, ntemp]


class isotopes(object):
  def __init__(self):
    self.niso    = 0            # Number of isotopes
    self.name    = np.array([]) # Isotope's name
    self.mass    = np.array([]) # Isotope's mass
    self.ratio   = np.array([]) # Isotopic abundance ratio
    self.dbindex = np.array([], np.int) # Isotope's data base index
    self.imol    = np.array([]) # Isotope's molecule index
    self.iext    = None         # Molecule index in ext. coef. table
    self.ntemp   = None         # Number of temperature samples
    self.temp    = None         # Temperature array
    self.z       = None         # Isotopes' partition function [niso, ntemp]


class voigt(object):
  def __init__(self):
    self.Dmin     = None  # Minimum Doppler width sampled
    self.Dmax     = None  # Maximum Doppler width sampled
    self.nDop     = None  # Number of Doppler-width samples
    self.Lmin     = None  # Minimum Lorentz width sampled
    self.Lmax     = None  # Maximum Lorentz width sampled
    self.nLor     = None  # Number of Lorentz-width samples
    self.doppler  = None  # Doppler-width sample array [nDop]
    self.lorentz  = None  # Lorentz-width sample array [nLor]
    self.DLratio  = None  # Doppler-Lorentz ratio threshold
    self.profile  = None  # Voigt profile [sum(2*size+1)]
    self.size     = None  # Profile wavenumber half-size [nDop, nLor]
    self.index    = None  # Index where each profile starts [nDop, nLor]


class extinction(object):
  def __init__(self):
    self.ec      = None # Molecular line-transition extinction-coefficient
                        #  in cm-1 [nlayers, nspec]
    self.extfile = None # Extinction-coefficient table filename
    self.etable  = None # Table of ext. coefficient [nmol, nlayer, nTemp, nspec]
    self.ethresh = None # Extinction-coefficient threshold
    self.tmin    = None # Minimum temperature to sample
    self.tmax    = None # Maximum temperature to sample
    self.tstep   = None # Temperature-sample step interval
    self.ntemp   = None # Number of temperature samples
    self.molID   = None # Tabulated species ID
    self.temp    = None # Tabulated temperatures
    self.press   = None # Tabulated pressures
    self.wn      = None # Tabulated wavenumber
    self.z       = None # Partition function at tabulated temperatures
    self.ciaext  = None # CIA extinction [nlayer, nwave]


class cia(object):
  def __init__(self):
    self.files      = None # CIA file names
    self.nfiles     = None # Number of files read
    self.molecules  = None # Molecules involved for each file
    self.ntemp      = None # Number of temperature samples per file
    self.nwave      = None # Number of wavenumber samples per file
    self.temp       = []   # Temperature sampling (in Kelvin)
    self.wavenumber = []   # Wavenumber sampling (in cm-1)
    self.absorption = []   # CIA extinction (in cm-1 amagat-2)
    self.ec         = None # Interpolated CIA extinction coefficient
                           #  in cm-1 [nlayer, nspec]

class optdepth(object):
  def __init__(self):
    self.ec      = None  # Total extinction coefficient [nlayers, nspec]
    self.raypath = []    # Distance along ray path  [nlayers]
    self.depth   = None  # Optical depth at raypath [nlayers, nspec]
    self.ideep   = None  # Layer index where depth reached maxdepth [nspec]

