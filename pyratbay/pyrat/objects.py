# Copyright (c) 2016-2017 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import numpy  as np

from .. import tools      as pt
from .. import constants  as pc
from .. import atmosphere as atm

class Pyrat(object):
  """
  Main Pyrat object.
  """
  def __init__(self):
    # Sub-classes:
    self.inputs = Inputs()          # User inputs
    self.spec   = Spectrum()        # Spectrum data
    self.atm    = Atm()             # Modeling atmospheric model
    self.lt     = Linetransition()  # Line-transition data
    self.mol    = Molecules()       # Molecules data
    self.iso    = Isotopes()        # Isotopes data
    self.voigt  = Voigt()           # Voigt profile
    self.ex     = Extinction()      # Extinction-coefficient
    self.cs     = Cross()           # Cross-section extinction
    self.od     = Optdepth()        # Optical depth
    self.haze   = Haze()            # Hazes
    self.rayleigh = Rayleigh()      # Rayleigh models
    self.alkali = Alkali()          # Alkali opacity models
    self.obs    = Observation()     # Observational data
    self.phy    = Physics()         # System physical parameters
    self.ret    = Retrieval()       # Retrieval variables
    # Files:
    self.atmfile     = None  # Atmopheric-model file
    self.linedb      = None  # Line-transition data file
    self.molfile     = None  # Molecular-properties file
    self.outsample   = None  # 
    self.outmaxdepth = None  # 
    self.outspec     = None  # Modulation/Flux spectrum file
    # Photometric surface:
    self.refpressure = None  # Pressure reference level
    # Atmosphere:
    self.radunits = None  # Radius physical units
    self.radstep  = None  # Radius sampling interval
    self.radlow   = None  # Lowest radius boundary
    self.radhigh  = None  # Highest radius boundary
    self.punits   = None  # Pressure physical units
    self.plow     = None  # Lowest pressure boundary
    self.phigh    = None  # Highest pressure boundary
    self.hydrom   = False # Variable/constant-g flag for hydrostatic equilib.
    # Geometry:
    self.raygrid  = None  # Array of incident ray-angles
    # Other:
    self.verb       = None  # Verbosity level
    self.logfile    = None  # Pyrat log filename
    self.log        = None  # Pyrat log file
    self.wlog       = []    # List of raised warnings
    self.timestamps = None  # Time stamps

  def hydro(self, pressure, temperature, mu, g, mass, p0, r0):
    """
    Hydrostatic-equilibrium driver.
    Depending on the self.hydrom flag, select between the g=GM/r**2
    (hydrom=True) or constant-g (hydrom=False) formula to compute
    the hydrostatic-equilibrium radii of the planet layers.

    Parameters
    ----------
    pressure: 1D float ndarray
       Atmospheric pressure for each layer (in barye).
    temperature: 1D float ndarray
       Atmospheric temperature for each layer (in K).
    mu: 1D float ndarray
       Mean molecular mass for each layer (in g mol-1).
    g: Float
       Atmospheric gravity (in cm s-2).
    mass: Float
       Planetary mass (in g).
    p0: Float
       Reference pressure level (in barye) where radius(p0) = r0.
    r0: Float
       Reference radius level (in cm) corresponding to p0.
    """
    # H.E. with  g=GM/r**2:
    if self.hydrom:
      return atm.hydro_m(pressure, temperature, mu, mass, p0, r0)
    # H.E. with constant g:
    return atm.hydro_g(pressure, temperature, mu, g, p0, r0)


class Inputs(object):
  """
  This is a holder class to store user-input arguments.
  """
  def __init__(self):
    self.atm = Atm()  # Input-file atmospheric model


class Spectrum(object):
  def __init__(self):
    self.nwave     = None  # Number of wavenumber spectral samples
    # Wavenumber:
    self.wnunits   = None  # User-input wavenumber physical units
    self.wn        = None  # Wavenumber array
    self.wnlow     = None  # Lowest wavenumber boundary
    self.wnhigh    = None  # Highest wavenumber boundary
    self.wnstep    = None  # Wavenumber sampling interval
    # Oversampled-wavenumber:
    self.wnosamp   = None  # Wavenumber oversampling factor
    self.own       = None  # Oversampled wavenumber array
    self.ownstep   = None  # Oversampled wavenumber sampling interval
    self.onwave    = None  # Number of oversampled-wavenumber samples
    self.odivisors = None  # Oversampling-factor integer divisors
    # Wavelength:
    self.wlunits   = None  # User-input wavelength physical units
    self.wllow     = None  # Lowest wavelength boundary
    self.wlhigh    = None  # Highest wavelength boundary
    # Spectrum:
    self.intensity = None  # Intensity spectrum array
    self.spectrum  = None  # Modulation/Flux spectrum
    self.clear     = None  # Clear modulation spectrum for patchy model
    self.cloudy    = None  # Cloudy modulation spectrum for patchy model
    self.starflux  = None  # Stellar flux spectrum

  def info(self):
    """
    Print the Spectral info.
    """
    pt.msg(1, "Spectral info:")
    pt.msg(1, "Wavenumber:", indent=2)
    pt.msg(1, "Number of samples:      {:d}".format(self.nwave), indent=4)
    pt.msg(1, "User-input units:       {:s}-1".format(self.wnunits), indent=4)
    pt.msg(1, "Pyrat (internal) units: cm-1", indent=4)
    pt.msg(1, "Low  boundary:     {:9.3f} cm-1".format(self.wnlow),  indent=4)
    pt.msg(1, "High boundary:     {:9.3f} cm-1".format(self.wnhigh), indent=4)
    pt.msg(1, "Sampling interval: {:9.3f} cm-1".format(self.wnstep), indent=4)
    pt.msg(1, "Wavenumber array (cm-1):\n  [{:.3f}, {:.3f}, {:.3f}, ..., "
              "{:.3f}, {:.3f}]".format(self.wn[ 0], self.wn[ 1], self.wn[2],
                                       self.wn[-2], self.wn[-1]), indent=4)
    pt.msg(1, "Oversampled wavenumber:", indent=2)
    pt.msg(1, "Oversampling factor:    {:d}".format(self.wnosamp),   indent=4)
    pt.msg(1, "Number of samples:      {:d}".format(self.onwave),    indent=4)
    pt.msg(1, "Sampling interval: {:.3e} cm-1".format(self.ownstep), indent=4)
    pt.msg(1, "Integer divisors for oversampling factor:\n{:s}".
                      format(str(self.odivisors).replace("\n", "")), indent=4)
    pt.msg(1, "Wavenumber:", indent=2)
    pt.msg(1, "User-input units: {:s}".format(self.wlunits), indent=4)
    pt.msg(1, "Low  boundary: {:7.3f} {:s}".
                format(self.wllow/pt.u(self.wlunits), self.wlunits),  indent=4)
    pt.msg(1, "High boundary: {:7.3f} {:s}".
                format(self.wlhigh/pt.u(self.wlunits), self.wlunits), indent=4)
    pt.msg(1, "Spectrum:", indent=2)
    if self.intensity is not None:
      pt.msg(1, "Intensity spectrum array (erg/s/cm/sr): [{:.3f}, {:.3f}, "
                "{:.3f}, ..., {:.3f}, {:.3f}]".format(self.intensity[ 0],
                             self.intensity[ 1], self.intensity[ 2],
                             self.intensity[-2], self.intensity[-1]), indent=4)
    if self.spectrum is None:
      pt.msg(1, "Modulation/Flux spectrum array: None", indent=4)
    else:
      # FINDME: how to get the transit/eclipse geometry?
      pt.msg(1, "Modulation/Flux spectrum array: [{:.3f}, {:.3f}, {:.3f}, ..., "
              "{:.3f}, {:.3f}]".format(self.spectrum[ 0], self.spectrum[ 1],
              self.spectrum[2], self.spectrum[-2], self.spectrum[-1]), indent=4)


class Atm(object):
  def __init__(self):
    self.qunits    = None      # Input abundance units ('mass' or 'number')
    self.runits    = None      # Input radius units
    self.punits    = None     # Input pressure units
    self.tunits    = 'kelvin'  # Input temperature units
    self.nlayers   = None      # Number of layers
    self.radius    = None      # Radius array (cm)            [layers]
    self.press     = None      # Pressure array (barye)       [layers]
    self.temp      = None      # Temperature array (K)        [layers]
    self.mm        = None      # Mean molecular mass (gr/mol) [layers]
    self.q         = None      # Molecular abundances         [layers, nmol]
    self.d         = None      # Molecular densities          [layers, nmol]
    self.rtop      = 0         # Index of topmost layer (within Hill radius)

  def info(self):
    pt.msg(1, "Atmospheric model info:")
    pt.msg(1, "Abundance input units:   {:s}.".format(self.qunits),  indent=2)
    pt.msg(1, "Radius input units:      {:s}.".format(self.runits),  indent=2)
    pt.msg(1, "Pressure input units:    {:s}.".format(self.punits),  indent=2)
    pt.msg(1, "Temperature input units: {:s}.".format(self.tunits),  indent=2)
    pt.msg(1, "Number of layers: {:d}".        format(self.nlayers), indent=2)
    pt.msg(1, "Radius (km):        [{:8.1f}, {:8.1f}, ..., {:8.1f}].".
              format(self.radius[0]/pc.km,
                     self.radius[1]/pc.km, self.radius[-1]/pc.km),   indent=4)
    pt.msg(1, "Pressure (bar):     [{:.2e}, {:.2e}, ..., {:.2e}].".
              format(self.press[0]/pc.bar,
                     self.press[1]/pc.bar, self.press[-1]/pc.bar),   indent=4)
    pt.msg(1, "Temperature (K):    [{:8.2f}, {:8.2f}, ..., {:8.2f}].".
              format(self.temp[0],   self.temp[1],   self.temp[-1]), indent=4)
    pt.msg(1, "Mean M. Mass (amu): [{:8.4f}, {:8.4f}, ..., {:8.4f}].".
              format(self.mm[0],     self.mm[1],     self.mm[-1]),   indent=4)
    pt.msg(1, "Number of species: {:d}".format(len(self.q[0])),      indent=2)
    pt.msg(1, "Abundances (mole mixing ratio):", indent=2)
    for i in np.arange(len(self.q[0])):
      pt.msg(1, "Species [{: 2d}]:       [{:.2e}, {:.2e}, ..., {:.2e}].".
                format(i, self.q[0,i], self.q[1,i], self.q[-1,i]), indent=4)
    pt.msg(1, "Density (gr/cm3):", indent=2)
    for i in np.arange(len(self.q[0])):
      pt.msg(1, "Species [{: 2d}]:       [{:.2e}, {:.2e}, ..., {:.2e}].".
                format(i, self.d[0,i], self.d[1,i], self.d[-1,i]), indent=4)


class Molecules(object):
  def __init__(self):
    self.nmol   = 0     # Number of species
    self.name   = None  # Species' name               [nmol]
    self.symbol = None  # Species' symbol             [nmol]
    self.mass   = None  # Species' mass  (gr/mol)     [nmol]
    self.radius = None  # Species' radius (Angstroms) [nmol]
    self.ID     = None  # Species' universal ID       [nmol]

  def info(self):
    pt.msg(1, "Atmospheric species info:")
    pt.msg(1, "Number of species: {:d}".format(self.nmol), indent=2)
    pt.msg(1, "Species:   ID   Mass      Radius\n"
              "                (gr/mol)  (Angstrom)",      indent=2)
    for i in np.arange(self.nmol):
      pt.msg(1, "{:>7s}:  {:3d}  {:8.4f}  {:.3f}".
             format(self.symbol[i], self.ID[i],
                    self.mass[i], self.radius[i]/pc.A),    indent=2)


class Linetransition(object):
  def __init__(self):
    self.nTLI    = 0        # Number of TLI files
    self.ndb     = 0        # Number of data bases
    self.db      = []       # Data base objects
    self.ntransitions = 0   # Number of line transitions
    self.tmin    = -np.inf  # Minimum temperature sampled by all TLI files
    self.tmax    =  np.inf  # Maximum temperature sampled by all TLI files
    self.wn      = np.array([], np.double)  # Line wavenumber
    self.elow    = np.array([], np.double)  # Line lower energy level
    self.gf      = np.array([], np.double)  # Line gf value
    self.isoid   = np.array([], np.int)     # Line isotope index

  def info(self):
    pt.msg(1, "Line-transition info:")
    pt.msg(1, "Number of TLI files:           {:d}".format(self.nTLI), indent=2)
    pt.msg(1, "Number of databases (species): {:d}".format(self.ndb),  indent=2)
    for i in np.arange(self.ndb):
      self.db[i].info(2)
    pt.msg(1, "Number of line transitions:    {:d}".format(self.ntransitions),
                                                                       indent=2)
    pt.msg(1, "Minimum and maximum covered temperatures: [{:.1f}, {:.1f}] K".
              format(self.tmin, self.tmax), indent=2)


class Database(object):
  def __init__(self):
    self.name    = None  # Data base name
    self.molname = None  # Molecule name
    self.niso    = None  # Number of isotopes in database
    self.iiso    = None  # Isotope correlative index
    self.ntemp   = None  # Number of temperature samples
    self.temp    = None  # Temperature array
    self.z       = None  # Isotopes' partition function array [niso, ntemp]

  def info(self, idt=0):
    pt.msg(1, "Database info:", indent=idt)
    pt.msg(1, "Database name: {:s}".format(self.name),      indent=2+idt)
    pt.msg(1, "Species' name: {:s}".format(self.molname),   indent=2+idt)
    pt.msg(1, "Number of isotopes: {:d}".format(self.niso), indent=2+idt)
    pt.msg(1, "Isotope correlative index: {:d}".format(self.iiso), indent=2+idt)
    pt.msg(1, "Number of temperature samples: {:d} (for partition function)".
                format(self.ntemp), indent=2+idt)
    pt.msg(1, "Temperature boundaries (K): [{:.1f}, {:.1f}]".
                format(self.temp[0], self.temp[-1]), indent=2+idt)


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
    # Print info to screen:
    pt.msg(1, "Isotopes info:")
    pt.msg(1, "Number of isotopes: {:d}".format(self.niso), indent=2)
    pt.msg(1,
      "Isotope:  Species  Mass      Isotopic   Database  Ext-coefficient\n"
      "                   (gr/mol)  ratio      index     table index", indent=2)
    for i in np.arange(self.niso):
      pt.msg(1, "{:>7s}:  {:>7s}  {:8.4f}  {:.3e}       {:3d}  {}".
             format(self.name[i], pyrat.mol.name[self.imol[i]],
               self.mass[i], self.ratio[i], self.dbindex[i], iext[i]), indent=2)
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
    self.extent   = None  # Extent covered by the profile (in number of HWHM)
    self.profile  = None  # Voigt profile [sum(2*size+1)]
    self.size     = None  # Profile wavenumber half-size [nDop, nLor]
    self.index    = None  # Index where each profile starts [nDop, nLor]

  def info(self):
    pt.msg(1, "Voigt profile info:", 0)
    pt.msg(1, "Number of Doppler-width samples:  {:d}".format(self.nDop),
           indent=2)
    pt.msg(1, "Number of Lorentz-width samples:  {:d}".format(self.nLor),
           indent=2)
    pt.msg(1, "Doppler-width array (cm-1):  [{:.3e}, {:.3e}, ..., {:.3e}]".
               format(self.doppler[0], self.doppler[1], self.doppler[-1]),
           indent=2)
    pt.msg(1, "Lorentz-width array (cm-1):  [{:.3e}, {:.3e}, ..., {:.3e}]".
               format(self.lorentz[0], self.lorentz[1], self.lorentz[-1]),
           indent=2)
    pt.msg(1, "Doppler--Lorentz ratio threshold:  {:.3e}".
               format(self.DLratio), indent=2)
    pt.msg(1, "Extent covered by a profile in units of Voigt half widths:  "
              "{:.2f}".format(self.extent), indent=2)
    pt.msg(1, "Total size of all Voigt profiles (~sum of: 2*size+1):  {:d}".
               format(np.size(self.profile)), indent=2)
    pt.msg(1, "Voigt-profile half-sizes [Ndop, Nlor]:", indent=2)
    pt.msg(1, "{}".format(self.size), indent=4)
    pt.msg(1, "Voigt-profile indices [Ndop, Nlor]:", indent=2)
    pt.msg(1, "{}".format(self.index), indent=4)


class Extinction(object):
  def __init__(self):
    self.ec      = None # Molecular line-transition extinction coefficient
                        #  in cm-1 [nlayers, nwave]
    self.ethresh = None # Extinction-coefficient threshold

    self.extfile = None # Extinction-coefficient table filename
    self.etable  = None # Table of ext. coefficient [nmol, nlayer, ntemp, nwave]

    self.tmin    = None # Minimum temperature to sample
    self.tmax    = None # Maximum temperature to sample
    self.tstep   = None # Temperature-sample step interval
    self.z       = None # Partition function at tabulated temperatures
                        #   [niso, ntemp]
    self.nmol    = None # Number of species
    self.ntemp   = None # Number of temperature samples
    self.nlayers = None # Number of pressure layers
    self.nwave   = None # Number of wavenumber spectral samples

    self.molID   = None # Tabulated species ID
    self.temp    = None # Tabulated temperatures
    self.press   = None # Tabulated pressures
    self.wn      = None # Tabulated wavenumber


  def info(self):
    pt.msg(1, "Extinction coefficient info:")
    pt.msg(1, "Line-transition strength threshold: {:.3e}".
               format(self.ethresh), indent=2)
    if self.extfile is None:
      pt.msg("No extinction-coefficient table defined.",   indent=2)
    else:
      pt.msg(1, "Extinction-coefficient table filename:",  indent=2)
      pt.msg(1, "'{:s}'".format(self.extfile), 4)
      pt.msg(1, "Minimum temperature:           {:6.1f} K".format(self.tmin),
                                                           indent=4)
      pt.msg(1, "Maximum temperature:           {:6.1f} K".format(self.tmax),
                                                           indent=4)
      pt.msg(1, "Temperature sampling interval: {:6.1f} K".format(self.tstep),
                                                           indent=4)
      pt.msg(1, "Number of tabulated species:          {:5d}".
                 format(self.nmol),    indent=4)
      pt.msg(1, "Number of tabulated temperatures:     {:5d}".
                 format(self.ntemp),   indent=4)
      pt.msg(1, "Number of tabulated layers:           {:5d}".
                 format(self.nlayers), indent=4)
      pt.msg(1, "Number of tabulated spectral samples: {:5d}".
                 format(self.nwave),   indent=4)
      pt.msg(1, "Temperature array (K):   [{:8.1f}, {:8.1f}, ..., {:8.1f}]".
                 format(self.temp[0], self.temp[1], self.temp[-1]), indent=4)
      pt.msg(1, "Partition function at tabulated temperatures:", 4)
      pt.msg(1, "{}".format(self.z), indent=6)
      pt.msg(1, "Species ID array: {:s}".
                 format(str(self.molID).replace("\n", "")), indent=4)
      pt.msg(1, "Pressure array: (bar)    [{:.2e}, {:.2e}, ..., {:.2e}]".
                 format(self.press[0]/pc.bar,
                        self.press[1]/pc.bar, self.press[-1]/pc.bar), indent=4)
      pt.msg(1, "Wavenumber array (cm-1): [{:8.3f}, {:8.3f}, ..., {:8.3f}]".
                 format(self.wn[0], self.wn[1], self.wn[-1]), indent=4)
      np.set_printoptions(formatter={'float': '{: .1e}'.format})
      pt.msg(1, "Tabulated extinction coefficient (cm2 gr-1)\n"
                "                       [spec, temp, layer, wave]:", indent=4)
      pt.msg(1, "{}".format((self.etable)), indent=4)
    if self.ec is not None:
      np.set_printoptions(formatter={'float': '{: .2e}'.format})
      pt.msg(1, "\nLine-transition extinction coefficient for the "
                   "atmospheric model (cm-1) [layer, wave]:", indent=2)
      pt.msg(1, "{}".format((self.ec)), indent=2)
    np.set_printoptions(formatter=None)



class Cross(object):
  def __init__(self):
    self.files      = None    # CS file names
    self.nfiles     = None    # Number of files read
    self.nmol       = None    # Number of species per CS file
    self.molecules  = None    # Species involved for each file
    self.ntemp      = None    # Number of temperature samples per file
    self.nwave      = None    # Number of wavenumber samples per file
    self.tmin       =     0.0 # Minimum temperature sampled by all CS files
    self.tmax       = 70000.0 # Maximum temperature sampled by all CS files
    self.temp       = []      # Temperature sampling (in Kelvin)
    self.wavenumber = []      # Wavenumber sampling (in cm-1)
    self.absorption = []      # CS extinction (in cm-1 amagat-2)
    self.iabsorp    = []      # wn-interpolated CS extinction (in cm-1 amagat-2)
    self.iz         = []      # Second derivatives of iabsorp
    self.iwnlo      = []      # Lower-wavenumber index for interpolation
    self.iwnhi      = []      # Upper-wavenumber index for interpolation
    self.ec         = None    # Interpolated CS extinction coefficient
                              #  in cm-1 [nlayer, nwave]

  def info(self):
    pt.msg(1, "Cross-section extinction info:")
    pt.msg(1, "Number of CS files: {:d}".format(self.nfiles), indent=2)
    for i in np.arange(self.nfiles):
      pt.msg(1, "CS file: '{:s}':".format(self.files[i]), indent=2)
      pt.msg(1, "Species: {:s}".
                 format("-".join(self.molecules[i,0:self.nmol[i]])), indent=4)
      pt.msg(1, "Number of temperatures:       {:4d}".format(self.ntemp[i]),
                                                                     indent=4)
      pt.msg(1, "Number of wavenumber samples: {:4d}".format(self.nwave[i]),
                                                                     indent=4)
      pt.msg(1, "Temperature array (K):  {}".
                 format(str(self.temp[i]).replace("\n", "")), indent=4, si=6)
      pt.msg(1, "Wavenumber array (cm-1): [{:7.1f}, {:7.1f}, ..., {:7.1f}]".
                 format(self.wavenumber[i][0], self.wavenumber[i][1],
                        self.wavenumber[i][-1]), indent=4)
      np.set_printoptions(formatter={'float': '{: .1e}'.format})
      pt.msg(1, "Tabulated CS extinction coefficient (cm-1 amagat-{:d}) "
                "[layer, wave]:".format(self.nmol[i]), indent=4)
      pt.msg(1, "{}".format((self.ec)),                indent=6)
      np.set_printoptions(formatter=None)
    pt.msg(1, "\nMinimum and maximum covered temperatures (K): "
              "[{:.1f}, {:.1f}]".format(self.tmin, self.tmax), indent=2)
    # FINDME iabsorp ?
    if self.ec is not None:
      np.set_printoptions(formatter={'float': '{: .2e}'.format})
      pt.msg(1, "CS extinction coefficient for the "
                   "atmospheric model (cm-1) [layer, wave]:", indent=2)
      pt.msg(1, "{}".format((self.ec)), indent=2)
    np.set_printoptions(formatter=None)


class Haze(object):
  def __init__(self):
    self.nmodels = 0     # Number of haze models
    self.model   = []    # List of haze models
    self.ec      = None  # Haze extinction coefficient
    self.fpatchy = None  # Pyatchy-cloud fraction
  def info(self, pyrat):
    # FINDME
    pass


class Rayleigh(object):
  def __init__(self):
    self.nmodels = 0     # Number of Rayleigh models
    self.model   = []    # List of Rayleigh models
    self.ec      = None  # Rayleigh extinction coefficient
  def info(self, pyrat):
    # FINDME
    pass


class Alkali(object):
  def __init__(self):
    self.nmodels = 0     # Number of alkali models
    self.model   = []    # List of alkali models
    self.ec      = None  # Alkali extinction coefficient
  def info(self, pyrat):
    # FINDME
    pass


class Optdepth(object):
  def __init__(self):
    self.maxdepth = None  # Maximum optical depth to calculate
    self.path     = None  # Observing geometry
    self.ec       = None  # Total extinction coefficient [nlayers, nwave]
    self.epatchy  = None  # Cloudy extinction coefficient for patchy model
    self.raypath  = []    # Distance along ray path  [nlayers]
    self.depth    = None  # Optical depth at raypath [nlayers, nwave]
    self.pdepth   = None  # Cloudy optical depth for patchy model
    self.B        = None  # Blackbody Planck emission [nlayers, nwave]
    self.ideep    = None  # Layer index where depth reached maxdepth [nwave]

  def info(self, pyrat):
    pt.msg(1, "Optical depth info:")
    pt.msg(1, "Ray-path geometry:  {:s}".format(self.path), indent=2)
    pt.msg(1, "Maximum optical depth to calculate:  {:.2f}".
               format(self.maxdepth), indent=2)
    if self.ec is not None:
      np.set_printoptions(formatter={'float': '{: .2e}'.format})
      pt.msg(1, "Total atmospheric-model extinction coefficient (cm-1) "
                "[layer, wave]:",       indent=2)
      pt.msg(1, "{}".format((self.ec)), indent=2)
      np.set_printoptions(formatter=None)
    if self.depth is not None:
      pt.msg(1, "Layer index where the optical depth reached maxdepth:",
                                         indent=2)
      pt.msg(1, "{}".format(self.ideep), indent=4)
      np.set_printoptions(formatter={'float': '{: .1f}'.format})
      # Raypath for transit geometry:
      if self.path == "transit":
        pt.msg(1, "\nDistance (km) along the raypath over each layer "
                  "(outside-in) for each impact parameter:", indent=2)
        pt.msg(1, "IP[  1] ({:.1f} km): {}".format(
          pyrat.atm.radius[1]/pc.km, self.raypath[1]/pc.km), indent=4)
        pt.msg(1, "IP[  2] ({:.1f} km): {}".format(
          pyrat.atm.radius[2]/pc.km, self.raypath[2]/pc.km), indent=4)
        pt.msg(1, "IP[  3] ({:.1f} km): {}".format(
          pyrat.atm.radius[3]/pc.km, self.raypath[3]/pc.km), indent=4)
        pt.msg(1, "...", indent=4)
        pt.msg(1, "IP[{:3d}] ({:.1f} km): {}".format(len(pyrat.atm.radius),
               pyrat.atm.radius[-1]/pc.km,
               str(self.raypath[-1]/pc.km)).replace("\n", ""), indent=4, si=6)
        pt.msg(1, "\nOptical depth for each impact parameter (outside-in) for "
                "each wavenumber:", indent=2)
      # Raypath for eclipse geometry:
      elif self.path == "eclipse":
        pt.msg(1, "\nDistance over each layer along a normal-incident "
                  "raypath (km):  {}".
          format(str(self.raypath/pc.km).replace("\n", "")), indent=2, si=4)
        pt.msg(1, "\nOptical depth over each layer (outside-in) along a "
                  "normal-incident raypath for each wavenumber:", indent=2)
      # Print optical depth:
      np.set_printoptions(formatter={'float': '{: .1e}'.format})
      pt.msg(1, "At {:7.1f} cm-1:  {}".format(pyrat.spec.wn[0],
             str(self.depth[0:self.ideep[0]+1,0]).replace("\n","")),
                       indent=4, si=6)
      pt.msg(1, "...", indent=4)
      index = pyrat.spec.nwave/2
      pt.msg(1, "At {:7.1f} cm-1:  {}".format(pyrat.spec.wn[index],
        str(self.depth[0:self.ideep[index]+1,index]).replace("\n","")),
                       indent=4, si=6)
      pt.msg(1, "...", indent=4)
      index = pyrat.spec.nwave-1
      pt.msg(1, "At {:7.1f} cm-1:  {}".format(pyrat.spec.wn[index],
        str(self.depth[0:self.ideep[index]+1,index]).replace("\n","")),
                       indent=4, si=6)


class Observation(object):
  def __init__(self):
    self.ndata    = 0     # Number of data points
    self.data     = None  # Transit or eclipse data point
    self.uncert   = None  # Data's 1-sigma uncertainty
    self.nfilters  = 0     # Number of filter bands
    self.filter    = None  # Observing filter filename
    self.bandidx   = None  # Band wavenumber indices
    self.bandtrans = None  # Band-interpolated transmission function
    self.starflux  = None  # Band-interpolated stellar flux
    self.bandflux  = None  # Band-integrated flux
    self.bandwn    = None  # Filter mean wavenumber
  def info(self, pyrat):
    # FINDME
    pass

# Retrieval variables:
class Retrieval(object):
  def __init__(self):
    # Available model types for retrieval:
    self.rmodels    = ["pt", "rad", "mol", "ray", "haze", "cloud", "patchy"]
    self.retflag    = None  # Flags for models to be included for retrieval
    self.nparams    = 0     # Number of free parameters
    self.tmodelname = None  # Temperature-model name
    self.tmodel     = None  # Temperature model
    self.targs      = None  # Temperature-model arguments
    self.tlow       = None  # Lower-temperature retrieval boundary
    self.thigh      = None  # Higher-temperature retrieval boundary
    self.bulk       = None  # Bulk species name list
    self.molscale   = None  # Variable-abundance species name list
    self.ibulk      = None  # Indices of bulk species in pyrat.mol.name
    self.iscale     = None  # Indices of variable-abundance species
    self.bulkratio  = None  # Abundance ratio among bulk species
    self.invsrat    = None  # Inverse of the sum of the bulk ratios/layer
    self.itemp  = None  # Temperature-model parameter indices
    self.irad   = None  # Reference-radius model parameter index
    self.iabund = None  # Abundance-model parameter indices
    self.iray   = None  # Haze-model parameter indices
    self.ihaze  = None  # Haze-model parameter indices
    self.icloud = None  # Haze-model parameter indices
    self.ipatchy = None  # Patchy-model parameter index
    self.parname = []   # Model parameter names
  def info(self, pyrat):
    # FINDME
    pass


# System physical variables:
class Physics(object):
  def __init__(self):
    self.tstar    = None  # Stellar effective temperature
    self.rstar    = None  # Stellar radius
    self.mstar    = None  # Stellar mass
    self.gstar    = None  # Stellar surface gravity
    self.rplanet  = None  # Planetary radius
    self.mplanet  = None  # Planetary mass
    self.gplanet  = None  # Planetary surface gravity
    self.rprs     = None  # Planet-to-star radius ratio
    self.smaxis   = None  # Orbital semi-major axis
    self.rhill    = np.inf  # Planetary Hill radius
    self.starspec = None  # Stellar spectrum filename
    self.kurucz   = None  # Kurucz stellar spectrum
    self.marcs    = None  # MARCS stellar spectrum
    self.phoenix  = None  # PHOENIX stellar spectrum
    self.starwn   = None  # Input stellar wavenumber array
    self.starflux = None  # Input stellar flux spectrum in  FINDME units
  def info(self, pyrat):
    # FINDME
    pass

