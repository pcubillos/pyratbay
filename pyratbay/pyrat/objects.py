# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import textwrap

import numpy  as np

from .. import tools      as pt
from .. import constants  as pc
from .. import atmosphere as atm
from .  import extinction as ex
from .  import crosssec   as cs
from .  import rayleigh   as ray
from .  import haze       as hz
from .  import alkali     as al


def wrap(outlist, text, indent=0):
    """
    Wrap input text, store it into outlist list.

    Parameters
    ----------
    outlist: List
    text: String
    indent: Integer
    """
    # TDB: Document and move to tools.
    indspace = " "*indent
    lines = text.splitlines()
    for line in lines:
        outlist.append(textwrap.fill(line, break_long_words=False,
                                     initial_indent=indspace,
                                     subsequent_indent=indspace, width=80))


class Pyrat(object):
  """
  Main Pyrat object.
  """
  def __init__(self):
    # Sub-classes:
    self.inputs   = Inputs()          # User inputs
    self.spec     = Spectrum()        # Spectrum data
    self.atm      = Atm()             # Modeling atmospheric model
    self.lt       = Linetransition()  # Line-transition data
    self.mol      = Molecules()       # Molecules data
    self.iso      = Isotopes()        # Isotopes data
    self.voigt    = Voigt()           # Voigt profile
    self.ex       = Extinction()      # Extinction-coefficient
    self.cs       = Cross()           # Cross-section extinction
    self.od       = Optdepth()        # Optical depth
    self.haze     = Haze()            # Hazes
    self.rayleigh = Rayleigh()        # Rayleigh models
    self.alkali   = Alkali()          # Alkali opacity models
    self.obs      = Observation()     # Observational data
    self.phy      = Physics()         # System physical parameters
    self.ret      = Retrieval()       # Retrieval variables
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
    self.timestamps = None  # Time stamps


  def run(self):
      """Evaluate a model"""
      # TBD: Bring code from driver.run() in here.
      pass


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


  def get_ec(self, layer):
    """
    Extract extinction-coefficient contribution (in cm-1) from each
    component of the atmosphere at the requested layer.

    Parameters
    ----------
    layer: Integer
       The index of the atmospheric layer where to extract the EC.

    Returns
    -------
    ec: 2D float ndarray
       An array of shape [ncomponents, nwave] with the EC spectra
       (in cm-1) from each component of the atmosphere.
    label: List of strings
       The names of each atmospheric component that contributed to EC.
    """
    # Allocate outputs:
    ec = np.empty((0, self.spec.nwave))
    label = []
    # Line-by-line extinction coefficient:
    if self.ex.nmol != 0:
      e, lab = ex.get_ec(self, layer)
      ec = np.vstack((ec, e))
      label += lab
    # Cross-section extinction coefficient:
    if self.cs.nfiles != 0:
      e, lab = cs.interpolate(self, layer)
      ec = np.vstack((ec, e))
      label += lab
    # Rayleigh scattering extinction coefficient:
    if self.rayleigh.nmodels != 0:
      e, lab = ray.get_ec(self, layer)
      ec = np.vstack((ec, e))
      label += lab
    # Haze/clouds extinction coefficient:
    if self.haze.nmodels != 0:
      e, lab = hz.get_ec(self, layer)
      ec = np.vstack((ec, e))
      label += lab
    # Alkali resonant lines extinction coefficient:
    if self.alkali.nmodels != 0:
      e, lab = al.get_ec(self, layer)
      ec = np.vstack((ec, e))
      label += lab
    return ec, label


  def __repr__(self):
      if self.spec.resolution is not None:
         wave = "R={.0f}".format(self.spec.resolution)
      else:
         wave = "dwn={:.3f} cm-1".format(self.spec.wnstep)

      opacities = []
      if self.ex.nmol != 0:
          for molID in self.ex.molID:
              imol = np.where(self.mol.ID == molID)[0][0]
              opacities.append(self.mol.name[imol])
      if self.cs.nfiles != 0:
          for molecs in self.cs.molecules:
              opacities.append("-".join(molecs))
      if self.rayleigh.nmodels != 0:
          for ray in self.rayleigh.model:
              opacities.append(ray.name)
      for haze in self.haze.model:
          opacities.append(haze.name)
      for alkali in self.alkali.model:
          opacities.append(self.alkali.mol)

      return ("Pyrat atmospheric model\n"
          "configuration file:  '{:s}'\n"
          "Pressure profile (bar):  {:.2e} -- {:.2e} ({:d} layers)\n"
          "Wavelength range (um):  {:.2f} -- {:.2f} ({:d} samples, {:s})\n"
          "Composition:  {}\n"
          "Opacity sources:  {}".format(
          self.inputs.configfile,
          self.atm.press[ 0]/pc.bar,
          self.atm.press[-1]/pc.bar,
          self.atm.nlayers,
          1.0/(self.spec.wn[ 0]*pc.um),
          1.0/(self.spec.wn[-1]*pc.um),
          self.spec.nwave,
          wave,
          self.mol.name,
          opacities))


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

  def __repr__(self):
    """
    Print the Spectral info.
    """
    info = []
    wrap(info, "Spectral info:")
    wrap(info, "Wavenumber:", 2)
    wrap(info, "Number of samples:      {:d}".format(self.nwave), 4)
    wrap(info, "User-input units:       {:s}-1".format(self.wnunits), 4)
    wrap(info, "Pyrat (internal) units: cm-1", 4)
    wrap(info, "Low  boundary:     {:9.3f} cm-1".format(self.wnlow),  4)
    wrap(info, "High boundary:     {:9.3f} cm-1".format(self.wnhigh), 4)
    wrap(info, "Sampling interval: {:9.3f} cm-1".format(self.wnstep), 4)
    wrap(info, "Wavenumber array (cm-1):\n  [{:.3f}, {:.3f}, {:.3f}, ..., "
           "{:.3f}, {:.3f}]".format(self.wn[ 0], self.wn[ 1], self.wn[2],
                                    self.wn[-2], self.wn[-1]), 4)
    wrap(info, "Oversampled wavenumber:", 2)
    wrap(info, "Oversampling factor:    {:d}".format(self.wnosamp),   4)
    wrap(info, "Number of samples:      {:d}".format(self.onwave),    4)
    wrap(info, "Sampling interval: {:.3e} cm-1".format(self.ownstep), 4)
    wrap(info, "Integer divisors for oversampling factor:\n{:s}".
                   format(str(self.odivisors).replace("\n", "")), 4)
    wrap(info, "Wavenumber:", 2)
    wrap(info, "User-input units: {:s}".format(self.wlunits), 4)
    wrap(info, "Low  boundary: {:7.3f} {:s}".
             format(self.wllow/pt.u(self.wlunits), self.wlunits),  4)
    wrap(info, "High boundary: {:7.3f} {:s}".
             format(self.wlhigh/pt.u(self.wlunits), self.wlunits), 4)
    wrap(info, "Spectrum:", 2)
    if self.intensity is not None:
      wrap(info, "Intensity spectrum array (erg/s/cm/sr): [{:.3f}, {:.3f}, "
             "{:.3f}, ..., {:.3f}, {:.3f}]".format(self.intensity[ 0],
                             self.intensity[ 1], self.intensity[ 2],
                             self.intensity[-2], self.intensity[-1]), 4)
    if self.spectrum is None:
      wrap(info, "Modulation/Flux spectrum array: None", 4)
    else:
      wrap(info, "Modulation/Flux spectrum array: [{:.3f}, {:.3f}, {:.3f}, ..., "
              "{:.3f}, {:.3f}]".format(self.spectrum[ 0], self.spectrum[ 1],
              self.spectrum[2], self.spectrum[-2], self.spectrum[-1]), 4)


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
    self.tmodel    = None
    self.tpars     = None
    self.rtop      = 0         # Index of topmost layer (within Hill radius)

  def __repr__(self):
    info = []
    wrap(info, "Atmospheric model info:")
    wrap(info, "Abundance input units:   {:s}.".format(self.qunits),  2)
    wrap(info, "Radius input units:      {:s}.".format(self.runits),  2)
    wrap(info, "Pressure input units:    {:s}.".format(self.punits),  2)
    wrap(info, "Temperature input units: {:s}.".format(self.tunits),  2)
    wrap(info, "Number of layers: {:d}".        format(self.nlayers), 2)
    wrap(info, "Radius (km):        [{:8.1f}, {:8.1f}, ..., {:8.1f}].".
              format(self.radius[0]/pc.km,
                     self.radius[1]/pc.km, self.radius[-1]/pc.km),   4)
    wrap(info, "Pressure (bar):     [{:.2e}, {:.2e}, ..., {:.2e}].".
              format(self.press[0]/pc.bar,
                     self.press[1]/pc.bar, self.press[-1]/pc.bar),   4)
    wrap(info, "Temperature (K):    [{:8.2f}, {:8.2f}, ..., {:8.2f}].".
              format(self.temp[0],   self.temp[1],   self.temp[-1]), 4)
    wrap(info, "Mean M. Mass (amu): [{:8.4f}, {:8.4f}, ..., {:8.4f}].".
              format(self.mm[0],     self.mm[1],     self.mm[-1]),   4)
    wrap(info, "Number of species: {:d}".format(len(self.q[0])),      2)
    wrap(info, "Abundances (mole mixing ratio):", 2)
    for i in np.arange(len(self.q[0])):
      wrap(info, "Species [{: 2d}]:       [{:.2e}, {:.2e}, ..., {:.2e}].".
                format(i, self.q[0,i], self.q[1,i], self.q[-1,i]), 4)
    wrap(info, "Density (gr/cm3):", 2)
    for i in np.arange(len(self.q[0])):
      wrap(info, "Species [{: 2d}]:       [{:.2e}, {:.2e}, ..., {:.2e}].".
                format(i, self.d[0,i], self.d[1,i], self.d[-1,i]), 4)
    return "\n".join(info)


class Molecules(object):
  def __init__(self):
    self.nmol   = 0     # Number of species
    self.name   = None  # Species' name               [nmol]
    self.symbol = None  # Species' symbol             [nmol]
    self.mass   = None  # Species' mass  (gr/mol)     [nmol]
    self.radius = None  # Species' radius (Angstroms) [nmol]
    self.ID     = None  # Species' universal ID       [nmol]

  def __repr__(self):
    info = []
    wrap(info, "Atmospheric species info:")
    wrap(info, "Number of species: {:d}\n"
               "Species:   ID   Mass      Radius\n"
               "                (gr/mol)  (Angstrom)".format(self.nmol), 2)
    for i in np.arange(self.nmol):
      wrap(info, "{:>7s}:  {:3d}  {:8.4f}  {:.3f}".
             format(self.symbol[i], self.ID[i],
                    self.mass[i], self.radius[i]/pc.A), 2)
    return "\n".join(info)


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

  def __repr__(self):
    info = []
    wrap(info, "Line-transition info:")
    wrap(info, "Number of TLI files:           {:d}".format(self.nTLI), 2)
    wrap(info, "Number of databases (species): {:d}".format(self.ndb),  2)
    for i in np.arange(self.ndb):
      self.db[i].info(2)
    wrap(info, "Number of line transitions:    {:d}\n"
               "Minimum and maximum covered temperatures: [{:.1f}, {:.1f}] K".
               format(self.ntransitions, self.tmin, self.tmax), 2)
    return "\n".join(info)


class Database(object):
  def __init__(self):
    self.name    = None  # Data base name
    self.molname = None  # Molecule name
    self.niso    = None  # Number of isotopes in database
    self.iiso    = None  # Isotope correlative index
    self.ntemp   = None  # Number of temperature samples
    self.temp    = None  # Temperature array
    self.z       = None  # Isotopes' partition function array [niso, ntemp]

  def __repr__(self):
    info = []
    wrap(info, "Database info:")
    wrap(info, "Database name: {:s}".format(self.name),      2)
    wrap(info, "Species' name: {:s}".format(self.molname),   2)
    wrap(info, "Number of isotopes: {:d}".format(self.niso), 2)
    wrap(info, "Isotope correlative index: {:d}".format(self.iiso), 2)
    wrap(info, "Number of temperature samples: {:d} (for partition function)".
                format(self.ntemp), 2)
    wrap(info, "Temperature boundaries (K): [{:.1f}, {:.1f}]".
                format(self.temp[0], self.temp[-1]), 2)
    return "\n".join(info)


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

  def __repr__(self):
    info = []
    # Patch self.iext if Pyrat doesn't use an EC table:
    if self.ex.extfile is None:
      iext = [None]*self.niso
    else:
      iext = self.iext
    # Print info to screen:
    wrap(info, "Isotopes info:")
    wrap(info, "Number of isotopes: {:d}".format(self.niso), 2)
    wrap(info,
      "Isotope:  Species  Mass      Isotopic   Database  Ext-coefficient\n"
      "                   (gr/mol)  ratio      index     table index", 2)
    for i in np.arange(self.niso):
      wrap(info, "{:>7s}:  {:>7s}  {:8.4f}  {:.3e}       {:3d}  {}".
             format(self.name[i], self.mol.name[self.imol[i]],
               self.mass[i], self.ratio[i], self.dbindex[i], iext[i]), 2)
    # FINDME: Partition function?
    #wrap(info, "Partition Function:", 2)
    #for i in np.arange(self.niso):
    #  wrap(info, "{:>7s}: [{:.2e}, {:.2e}, ..., {:.2e}]".
    #             format(db.z[j,0], db.z[j,1], db.z[j,-1]), 4)
    return "\n".join(info)


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

  def __repr__(self):
    info = []
    wrap(info, "Voigt profile info:", 0)
    wrap(info, "Number of Doppler-width samples:  {:d}".format(self.nDop),
           2)
    wrap(info, "Number of Lorentz-width samples:  {:d}".format(self.nLor),
           2)
    wrap(info, "Doppler-width array (cm-1):  [{:.3e}, {:.3e}, ..., {:.3e}]".
               format(self.doppler[0], self.doppler[1], self.doppler[-1]),
           2)
    wrap(info, "Lorentz-width array (cm-1):  [{:.3e}, {:.3e}, ..., {:.3e}]".
               format(self.lorentz[0], self.lorentz[1], self.lorentz[-1]),
           2)
    wrap(info, "Doppler--Lorentz ratio threshold:  {:.3e}".
               format(self.DLratio), 2)
    wrap(info, "Extent covered by a profile in units of Voigt half widths:  "
              "{:.2f}".format(self.extent), 2)
    wrap(info, "Total size of all Voigt profiles (~sum of: 2*size+1):  {:d}".
               format(np.size(self.profile)), 2)
    wrap(info, "Voigt-profile half-sizes [Ndop, Nlor]:", 2)
    wrap(info, "{}".format(self.size), 4)
    wrap(info, "Voigt-profile indices [Ndop, Nlor]:", 2)
    wrap(info, "{}".format(self.index), 4)
    return "\n".join(info)


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

  def __repr__(self):
    info = []
    wrap(info, "Extinction coefficient info:")
    wrap(info, "Line-transition strength threshold: {:.3e}".
               format(self.ethresh), 2)
    if self.extfile is None:
      wrap(info, "No extinction-coefficient table defined.", 2)
    else:
      wrap(info, "Extinction-coefficient table filename:",  2)
      wrap(info, "'{:s}'".format(self.extfile), 4)
      wrap(info, "Minimum temperature:           {:6.1f} K".format(self.tmin), 4)
      wrap(info, "Maximum temperature:           {:6.1f} K".format(self.tmax), 4)
      wrap(info, "Temperature sampling interval: {:6.1f} K".format(self.tstep), 4)
      wrap(info, "Number of tabulated species:          {:5d}".
                 format(self.nmol), 4)
      wrap(info, "Number of tabulated temperatures:     {:5d}".
                 format(self.ntemp), 4)
      wrap(info, "Number of tabulated layers:           {:5d}".
                 format(self.nlayers), 4)
      wrap(info, "Number of tabulated spectral samples: {:5d}".
                 format(self.nwave), 4)
      wrap(info, "Temperature array (K):   [{:8.1f}, {:8.1f}, ..., {:8.1f}]".
                 format(self.temp[0], self.temp[1], self.temp[-1]), 4)
      wrap(info, "Partition function at tabulated temperatures:", 4)
      wrap(info, "{}".format(self.z), 6)
      wrap(info, "Species ID array: {:s}".
                 format(str(self.molID).replace("\n", "")), 4)
      wrap(info, "Pressure array: (bar)    [{:.2e}, {:.2e}, ..., {:.2e}]".
                 format(self.press[0]/pc.bar,
                        self.press[1]/pc.bar, self.press[-1]/pc.bar), 4)
      wrap(info, "Wavenumber array (cm-1): [{:8.3f}, {:8.3f}, ..., {:8.3f}]".
                 format(self.wn[0], self.wn[1], self.wn[-1]), 4)
      np.set_printoptions(formatter={'float': '{: .1e}'.format})
      wrap(info, "Tabulated extinction coefficient (cm2 gr-1)\n"
                "                       [spec, temp, layer, wave]:", 4)
      wrap(info, "{}".format((self.etable)), 4)
    if self.ec is not None:
      np.set_printoptions(formatter={'float': '{: .2e}'.format})
      wrap(info, "\nLine-transition extinction coefficient for the "
                   "atmospheric model (cm-1) [layer, wave]:", 2)
      wrap(info, "{}".format((self.ec)), 2)
    np.set_printoptions(formatter=None)
    return "\n".join(info)


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
    info = []
    wrap(info, "Cross-section extinction info:")
    wrap(info, "Number of CS files: {:d}".format(self.nfiles), 2)
    for i in np.arange(self.nfiles):
      wrap(info, "CS file: '{:s}':".format(self.files[i]), 2)
      wrap(info, "Species: {:s}".
                 format("-".join(self.molecules[i,0:self.nmol[i]])), 4)
      wrap(info, "Number of temperatures:       {:4d}".format(self.ntemp[i]),
                                                                     4)
      wrap(info, "Number of wavenumber samples: {:4d}".format(self.nwave[i]),
                                                                     4)
      wrap(info, "Temperature array (K):  {}".
                 format(str(self.temp[i]).replace("\n", "")), 4, si=6)
      wrap(info, "Wavenumber array (cm-1): [{:7.1f}, {:7.1f}, ..., {:7.1f}]".
                 format(self.wavenumber[i][0], self.wavenumber[i][1],
                        self.wavenumber[i][-1]), 4)
      np.set_printoptions(formatter={'float': '{: .1e}'.format})
      wrap(info, "Tabulated CS extinction coefficient (cm-1 amagat-{:d}) "
                "[layer, wave]:".format(self.nmol[i]), 4)
      wrap(info, "{}".format((self.ec)),                6)
      np.set_printoptions(formatter=None)
    wrap(info, "\nMinimum and maximum covered temperatures (K): "
              "[{:.1f}, {:.1f}]".format(self.tmin, self.tmax), 2)
    # FINDME iabsorp ?
    if self.ec is not None:
      np.set_printoptions(formatter={'float': '{: .2e}'.format})
      wrap(info, "CS extinction coefficient for the "
                   "atmospheric model (cm-1) [layer, wave]:", 2)
      wrap(info, "{}".format((self.ec)), 2)
    np.set_printoptions(formatter=None)
    return "\n".join(info)


class Haze(object):
  def __init__(self):
    self.nmodels = 0     # Number of haze models
    self.model   = []    # List of haze models
    self.ec      = None  # Haze extinction coefficient
    self.fpatchy = None  # Patchy-cloud fraction
    self.pars    = None  # Input haze parameters

  def __repr__(self):
    info = []
    return "\n".join(info)


class Rayleigh(object):
  def __init__(self):
    self.nmodels = 0     # Number of Rayleigh models
    self.model   = []    # List of Rayleigh models
    self.ec      = None  # Rayleigh extinction coefficient
    self.pars    = None  # Input rayleigh parameters

  def __repr__(self):
    info = []
    return "\n".join(info)


class Alkali(object):
  def __init__(self):
    self.nmodels = 0     # Number of alkali models
    self.model   = []    # List of alkali models
    self.ec      = None  # Alkali extinction coefficient
    self.imol    = None  # Species indices in atmospheric file
    self.doppler = None  # Tabulated Doppler widths
    self.lorentz = None  # Tabulated Lorentz widths
    self.voigt   = None  # Tabulated alkali Voigt profiles
    self.vsize   = None  # Size of the Voigt profiles
    self.vindex  = None  # Starting indices of the Voigt profiles

  def __repr__(self):
    info = []
    return "\n".join(info)


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

  def __repr__(self):
    info = []
    wrap(info, "Optical depth info:")
    wrap(info, "Ray-path geometry:  {:s}".format(self.path), 2)
    wrap(info, "Maximum optical depth to calculate:  {:.2f}".
               format(self.maxdepth), 2)
    if self.ec is not None:
      np.set_printoptions(formatter={'float': '{: .2e}'.format})
      wrap(info, "Total atmospheric-model extinction coefficient (cm-1) "
                 "[layer, wave]:", 2)
      wrap(info, "{}".format((self.ec)), 2)
      np.set_printoptions(formatter=None)
    if self.depth is not None:
      wrap(info, "Layer index where the optical depth reached maxdepth:", 2)
      wrap(info, "{}".format(self.ideep), 4)
      np.set_printoptions(formatter={'float': '{: .1f}'.format})
      # Raypath for transit geometry:
      if self.path == "transit":
        wrap(info, "\nDistance (km) along the raypath over each layer "
                  "(outside-in) for each impact parameter:", 2)
        wrap(info, "IP[  1] ({:.1f} km): {}".format(
          pyrat.atm.radius[1]/pc.km, self.raypath[1]/pc.km), 4)
        wrap(info, "IP[  2] ({:.1f} km): {}".format(
          pyrat.atm.radius[2]/pc.km, self.raypath[2]/pc.km), 4)
        wrap(info, "IP[  3] ({:.1f} km): {}".format(
          pyrat.atm.radius[3]/pc.km, self.raypath[3]/pc.km), 4)
        wrap(info, "...", 4)
        wrap(info, "IP[{:3d}] ({:.1f} km): {}".format(len(pyrat.atm.radius),
               pyrat.atm.radius[-1]/pc.km,
               str(self.raypath[-1]/pc.km)).replace("\n", ""), 4, si=6)
        wrap(info, "\nOptical depth for each impact parameter (outside-in) for "
                "each wavenumber:", 2)
      # Raypath for eclipse geometry:
      elif self.path == "eclipse":
        wrap(info, "\nDistance over each layer along a normal-incident "
                  "raypath (km):  {}".
          format(str(self.raypath/pc.km).replace("\n", "")), 2, si=4)
        wrap(info, "\nOptical depth over each layer (outside-in) along a "
                  "normal-incident raypath for each wavenumber:", 2)
      # Print optical depth:
      np.set_printoptions(formatter={'float': '{: .1e}'.format})
      wrap(info, "At {:7.1f} cm-1:  {}".format(pyrat.spec.wn[0],
             str(self.depth[0:self.ideep[0]+1,0]).replace("\n","")), 4, si=6)
      wrap(info, "...", 4)
      index = pyrat.spec.nwave/2
      wrap(info, "At {:7.1f} cm-1:  {}".format(pyrat.spec.wn[index],
        str(self.depth[0:self.ideep[index]+1,index]).replace("\n","")), 4, si=6)
      wrap(info, "...", 4)
      index = pyrat.spec.nwave-1
      wrap(info, "At {:7.1f} cm-1:  {}".format(pyrat.spec.wn[index],
        str(self.depth[0:self.ideep[index]+1,index]).replace("\n","")), 4, si=6)
    return "\n".join(info)


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

  def __repr__(self):
    info = []
    return "\n".join(info)


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
    self.icloud = None  # Cloud-model parameter indices
    self.ipatchy = None  # Patchy-model parameter index
    self.pnames   = []   # Model parameter names (screen)
    self.texnames = []   # Model parameter names (figures)

  def __repr__(self):
    info = []
    return "\n".join(info)


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

  def __repr__(self):
    info = []
    return "\n".join(info)
