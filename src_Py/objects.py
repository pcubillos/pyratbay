import numpy  as np

import ptools     as pt
import pconstants as pc

class Pyrat(object):
  """
  Main Pyrat object.
  """
  def __init__(self):
    # Sub-classes:
    self.inputs = Inputs()          # User inputs
    self.spec   = Spectrum()        # Spectrum data
    self.atmf   = Atm()             # Input-file atmospheric model
    self.atm    = Atm()             # Modeling atmospheric model
    self.lt     = Linetransition()  # Line-transition data
    self.mol    = Molecules()       # Molecules data
    self.iso    = Isotopes()        # Isotopes data
    self.voigt  = Voigt()           # Voigt profile
    self.ex     = Extinction()      # Extinction-coefficient
    self.cia    = Cia()             # Collision-induced absorption
    self.od     = Optdepth()        # Optical depth
    # Files:
    self.atmfile     = None  # Atmopheric-model file
    self.linedb      = None  # Line-transition data file
    self.molfile     = None  # Molecular-properties file
    self.outsample   = None  # 
    self.outmaxdepth = None  # 
    self.outspec     = None  # Modulation/Flux spectrum file
    # Photometric surface:
    self.pressurebase = None  # Pressure reference level
    self.radiusbase   = None  # Radius reference level
    # Atmosphere:
    self.radunits = None  # Radius physical units
    self.radstep  = None  # Radius sampling interval
    self.radlow   = None  # Lowest radius boundary
    self.radhigh  = None  # Highest radius boundary
    self.punits   = None  # Pressure physical units
    self.plow     = None  # Lowest pressure boundary
    self.phigh    = None  # Highest pressure boundary
    # Geometry:
    self.path     = None  # Observing geometry
    self.raygrid  = None  # Array of incident ray-angles
    self.maxdepth = None  # Maximum optical depth to calculate
    # Physical parameters:
    self.rstar       = None  # Stellar radius
    self.surfgravity = None  # Planetary surface gravity
    # Other:
    self.verb       = None  # Verbosity level
    self.timestamps = None  # Time stamps


class Inputs(object):
  """
  This is a holder class to store user-input arguments.
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
    self.raygrid  = None
    # System arguments:
    self.rstar = None
    # Output files arguments:
    self.outspec     = None
    self.outsample   = None
    self.outmaxdepth = None


class Spectrum(object):
  def __init__(self):
    self.nspec     = None  # Number of spectral samples
    # Wavenumber:
    self.wnunits   = None  # Wavenumber physical units
    self.wn        = None  # Wavenumber array
    self.wnlow     = None  # Lowest wavenumber boundary
    self.wnhigh    = None  # Highest wavenumber boundary
    self.wnstep    = None  # Wavenumber sampling interval
    # Oversampled-wavenumber:
    self.wnosamp   = None  # Wavenumber oversampling factor
    self.own       = None  # Oversampled wavenumber array
    self.ownstep   = None  # Oversampled wavenumber sampling interval
    self.onspec    = None  # Number of oversampled-wavenumber samples
    self.odivisors = None  # Oversampling-factor integer divisors
    # Wavelength:
    self.wlunits   = None  # Wavelength physical units
    self.wllow     = None  # Lowest wavelength boundary
    self.wlhigh    = None  # Highest wavelength boundary
    # Spectrum:
    self.intensity = None  # Intensity spectrum array
    self.spectrum  = None  # Modulation/Flux spectrum array

  def info(self):
    """
    Print the Spectral info.
    """
    pt.msg(1, "Spectral info:", 0)
    pt.msg(1, "Wavenumber:", 2)
    pt.msg(1, "Number of samples:      {:d}".format(self.nspec), 4)
    pt.msg(1, "User-input units:       {:s}-1".format(self.wnunits), 4)
    pt.msg(1, "Pyrat (internal) units: cm-1", 4)
    pt.msg(1, "Low  boundary:     {:9.3f} cm-1".format(self.wnlow),  4)
    pt.msg(1, "High boundary:     {:9.3f} cm-1".format(self.wnhigh), 4)
    pt.msg(1, "Sampling interval: {:9.3f} cm-1".format(self.wnstep), 4)
    pt.msg(1, "Wavenumber array (cm-1):\n  [{:.3f}, {:.3f}, {:.3f}, ..., "
              "{:.3f}, {:.3f}]".format(self.wn[ 0], self.wn[ 1], self.wn[2],
                                       self.wn[-2], self.wn[-1]), 4)
    pt.msg(1, "Oversampled wavenumber:", 2)
    pt.msg(1, "Oversampling factor:    {:d}".format(self.wnosamp), 4)
    pt.msg(1, "Number of samples:      {:d}".format(self.onspec), 4)
    pt.msg(1, "Sampling interval: {:.3e} cm-1".format(self.ownstep), 4)
    pt.msg(1, "Integer divisors for oversampling factor:\n{:s}".
                          format(str(self.odivisors).replace("\n", "")), 4)
    pt.msg(1, "Wavenumber:", 2)
    pt.msg(1, "User-input units: {:s}".format(self.wlunits), 4)
    pt.msg(1, "Low  boundary: {:7.3f} {:s}".
                  format(self.wllow/pc.units[self.wlunits], self.wlunits),  4)
    pt.msg(1, "High boundary: {:7.3f} {:s}".
                  format(self.wlhigh/pc.units[self.wlunits], self.wlunits), 4)
    pt.msg(1, "Spectrum:", 2)
    if self.intensity is not None:
      pt.msg(1, "Intensity spectrum array (erg/s/cm/sr): [{:.3f}, {:.3f}, "
                "{:.3f}, ..., {:.3f}, {:.3f}]".format(self.intensity[ 0],
                                  self.intensity[ 1], self.intensity[ 2],
                                  self.intensity[-2], self.intensity[-1]), 4)
    if self.spectrum is None:
      pt.msg(1, "Modulation/Flux spectrum array: None", 4)
    else:
      # FINDME: how to get the transit/eclipse geometry?
      pt.msg(1, "Modulation/Flux spectrum array: [{:.3f}, {:.3f}, {:.3f}, ..., "
                "{:.3f}, {:.3f}]".format(self.spectrum[ 0], self.spectrum[ 1],
                     self.spectrum[2], self.spectrum[-2], self.spectrum[-1]), 4)


class Atm(object):
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

  def info(self):
    pass


class Molecules(object):
  def __init__(self):
    self.nmol   = 0     # Number of species
    self.name   = None  # Species' name               [nmol]
    self.symbol = None  # Species' symbol             [nmol]
    self.mass   = None  # Species' mass  (gr/mol)     [nmol]
    self.radius = None  # Species' radius (Angstroms) [nmol]
    self.ID     = None  # Species' universal ID       [nmol]

  def info(self):
    pt.msg(1, "Atmospheric species info:", 0)
    pt.msg(1, "Number of species: {:d}".format(self.nmol), 2)
    pt.msg(1, "Species:   ID   Mass      Radius\n"
              "                (gr/mol)  (Angstrom)", 2)
    for i in np.arange(self.nmol):
      pt.msg(1, "{:>7s}:  {:3d}  {:8.4f}  {:.3f}".
             format(self.symbol[i], self.ID[i],
                    self.mass[i], self.radius[i]/pc.units["A"]), 2)


class Linetransition(object):
  def __init__(self):
    self.nTLI    = 0      # Number of TLI files
    self.ndb     = 0      # Number of data bases
    self.db      = []     # Data base objects
    self.ntransitions = 0 # Number of line transitions
    self.wn      = np.array([], np.double)  # Line wavenumber
    self.elow    = np.array([], np.double)  # Line lower energy level
    self.gf      = np.array([], np.double)  # Line gf value
    self.isoid   = np.array([], np.int)     # Line isotope index

  def info(self):
    pt.msg(1, "Line-transition info:", 0)
    pt.msg(1, "Number of TLI files:           {:d}".format(self.nTLI), 2)
    pt.msg(1, "Number of databases (species): {:d}".format(self.ndb),  2)
    for i in np.arange(self.ndb):
      self.db[i].info(2)
    pt.msg(1, "Number of line transitions:    {:d}".format(self.ntransitions),2)


class Database(object):
  def __init__(self):
    self.name    = None  # Data base name
    self.molname = None  # Molecule name
    self.niso    = None  # Number of isotopes in database
    self.iiso    = None  # Isotope correlative index
    self.ntemp   = None  # Number of temperature samples
    self.temp    = None  # Temperature array
    self.z       = None  # Isotopes' partition function array [niso, ntemp]

  def info(self, indent=0):
    pt.msg(1, "Database info:", 0+indent)
    pt.msg(1, "Database name: {:s}".format(self.name), 2+indent)
    pt.msg(1, "Species' name: {:s}".format(self.molname), 2+indent)
    pt.msg(1, "Number of isotopes: {:d}".format(self.niso), 2+indent)
    pt.msg(1, "Isotope correlative index: {:d}".format(self.iiso),  2+indent)
    pt.msg(1, "Number of temperature samples: {:d} (for partition function)".
                format(self.ntemp), 2+indent)
    pt.msg(1, "Temperature boundaries (K): [{:.1f}, {:.1f}]".
                format(self.temp[0], self.temp[-1]), 2+indent)


class Isotopes(object):
  def __init__(self):
    self.niso    = 0                     # Number of isotopes
    self.name    = np.array([])          # Isotope's name [niso]
    self.mass    = np.array([])          # Isotope's mass [niso]
    self.ratio   = np.array([])          # Isotopic abundance ratio  [niso]
    self.dbindex = np.array([], np.int)  # Isotope's data base index [niso]
    self.imol    = np.array([], np.int)  # Isotope's molecule index  [niso]
    self.iext    = None                  # Molecule index in ext-coef. table
    self.z       = None                  # Isotopes' partition function at
                                         #   atmospheric layer [niso, nlayers]

  def info(self, pyrat):
    # Patch self.iext if Pyrat doesn't use an EC table:
    if pyrat.ex.extfile is None:
      iext = [None]*self.niso
    else:
      iext = self.iext

    pt.msg(1, "Isotopes info:", 0)
    pt.msg(1, "Number of isotopes: {:d}".format(self.niso), 2)
    pt.msg(1,
          "Isotope:  Species  Mass      Isotopic   Database  Ext-coefficient\n"
          "                   (gr/mol)  ratio      index     table index", 2)
    for i in np.arange(self.niso):
      pt.msg(1, "{:>7s}:  {:>7s}  {:8.4f}  {:.3e}       {:3d}  {}".
             format(self.name[i], pyrat.mol.name[self.imol[i]],
                    self.mass[i], self.ratio[i], self.dbindex[i], iext[i]), 2)
    # FINDME: Partition function?
    #pt.msg(1, "Partition Function:", 2)
    #for i in np.arange(self.niso):
    #  pt.msg(1, "{:>7s}: [{:.2e}, {:.2e}, ..., {:.2e}]".
    #             format(db.z[j,0], db.z[j,1], db.z[j,-1]), 4)


class Voigt(object):
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

  def info(self):
    pass


class Extinction(object):
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

  def info(self):
    pass


class Cia(object):
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

  def info(self):
    pass

class Optdepth(object):
  def __init__(self):
    self.ec      = None  # Total extinction coefficient [nlayers, nspec]
    self.raypath = []    # Distance along ray path  [nlayers]
    self.depth   = None  # Optical depth at raypath [nlayers, nspec]
    self.ideep   = None  # Layer index where depth reached maxdepth [nspec]

  def info(self):
    pass

