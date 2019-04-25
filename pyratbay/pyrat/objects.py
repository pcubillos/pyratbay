# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ['Spectrum', 'Atm', 'Linetransition', 'Molecules',
           'Isotopes', 'Voigt', 'Extinction', 'Cross', 'Optdepth',
           'Haze', 'Rayleigh', 'Alkali', 'Observation', 'Physics',
           'Retrieval']

import numpy  as np

from .. import tools      as pt
from .. import constants  as pc


class Spectrum(object):
  def __init__(self):
    self.outspec   = None  # Modulation/Flux spectrum file
    # Wavenumber:
    self.nwave     = None  # Number of wavenumber spectral samples
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
    self.raygrid   = None  # Array of incident ray-angles (emission)
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
    pt.wrap(info, "Spectral info:")
    pt.wrap(info, "Wavenumber:", 2)
    pt.wrap(info, "Number of samples:      {:d}".format(self.nwave), 4)
    pt.wrap(info, "User-input units:       {:s}-1".format(self.wnunits), 4)
    pt.wrap(info, "Pyrat (internal) units: cm-1", 4)
    pt.wrap(info, "Low  boundary:     {:9.3f} cm-1".format(self.wnlow),  4)
    pt.wrap(info, "High boundary:     {:9.3f} cm-1".format(self.wnhigh), 4)
    pt.wrap(info, "Sampling interval: {:9.3f} cm-1".format(self.wnstep), 4)
    pt.wrap(info, "Wavenumber array (cm-1):\n  [{:.3f}, {:.3f}, {:.3f}, ..., "
        "{:.3f}, {:.3f}]".format(self.wn[ 0], self.wn[ 1], self.wn[2],
                                 self.wn[-2], self.wn[-1]), 4)
    pt.wrap(info, "Oversampled wavenumber:", 2)
    pt.wrap(info, "Oversampling factor:    {:d}".format(self.wnosamp),   4)
    pt.wrap(info, "Number of samples:      {:d}".format(self.onwave),    4)
    pt.wrap(info, "Sampling interval: {:.3e} cm-1".format(self.ownstep), 4)
    pt.wrap(info, "Integer divisors for oversampling factor:\n{:s}".
                  format(str(self.odivisors).replace("\n", "")), 4)
    pt.wrap(info, "Wavenumber:", 2)
    pt.wrap(info, "User-input units: {:s}".format(self.wlunits), 4)
    pt.wrap(info, "Low  boundary: {:7.3f} {:s}".
                format(self.wllow/pt.u(self.wlunits), self.wlunits),  4)
    pt.wrap(info, "High boundary: {:7.3f} {:s}".
                format(self.wlhigh/pt.u(self.wlunits), self.wlunits), 4)
    pt.wrap(info, "Spectrum:", 2)
    if self.intensity is not None:
      pt.wrap(info, "Intensity spectrum array (erg/s/cm/sr): [{:.3f}, {:.3f}, "
              "{:.3f}, ..., {:.3f}, {:.3f}]".format(self.intensity[ 0],
                             self.intensity[ 1], self.intensity[ 2],
                             self.intensity[-2], self.intensity[-1]), 4)
    if self.spectrum is None:
        pt.wrap(info, "Modulation/Flux spectrum array: None", 4)
    else:
        pt.wrap(info, "Modulation/Flux spectrum array: [{:.3f}, {:.3f}, "
            "{:.3f}, ..., {:.3f}, {:.3f}]".
            format(self.spectrum[ 0], self.spectrum[ 1], self.spectrum[2],
                   self.spectrum[-2], self.spectrum[-1]), 4)


class Atm(object):
  def __init__(self):
    # From pyrat:
    self.refpressure = None  # Pressure reference level
    self.radstep  = None  # Radius sampling interval
    self.radlow   = None  # Lowest radius boundary
    self.radhigh  = None  # Highest radius boundary
    self.ptop     = None  # Lowest pressure boundary
    self.pbottom  = None  # Highest pressure boundary
    self.hydrom   = False # Variable/constant-g flag for hydrostatic equilib.

    self.atmfile   = None      # Atmopheric-model file
    self.qunits    = None      # Input abundance units ('mass' or 'number')
    self.runits    = None      # Input radius units
    self.punits    = None      # Input pressure units
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
    pt.wrap(info, "Atmospheric model info:")
    pt.wrap(info, "Abundance input units:   {:s}.".format(self.qunits),  2)
    pt.wrap(info, "Radius input units:      {:s}.".format(self.runits),  2)
    pt.wrap(info, "Pressure input units:    {:s}.".format(self.punits),  2)
    pt.wrap(info, "Temperature input units: {:s}.".format(self.tunits),  2)
    pt.wrap(info, "Number of layers: {:d}".        format(self.nlayers), 2)
    pt.wrap(info, "Radius (km):        [{:8.1f}, {:8.1f}, ..., {:8.1f}].".
              format(self.radius[0]/pc.km,
                     self.radius[1]/pc.km, self.radius[-1]/pc.km),   4)
    pt.wrap(info, "Pressure (bar):     [{:.2e}, {:.2e}, ..., {:.2e}].".
              format(self.press[0]/pc.bar,
                     self.press[1]/pc.bar, self.press[-1]/pc.bar),   4)
    pt.wrap(info, "Temperature (K):    [{:8.2f}, {:8.2f}, ..., {:8.2f}].".
              format(self.temp[0],   self.temp[1],   self.temp[-1]), 4)
    pt.wrap(info, "Mean M. Mass (amu): [{:8.4f}, {:8.4f}, ..., {:8.4f}].".
              format(self.mm[0],     self.mm[1],     self.mm[-1]),   4)
    pt.wrap(info, "Number of species: {:d}".format(len(self.q[0])),      2)
    pt.wrap(info, "Abundances (mole mixing ratio):", 2)
    for i in np.arange(len(self.q[0])):
        pt.wrap(info, "Species [{: 2d}]:       [{:.2e}, {:.2e}, ..., {:.2e}].".
                format(i, self.q[0,i], self.q[1,i], self.q[-1,i]), 4)
    pt.wrap(info, "Density (gr/cm3):", 2)
    for i in np.arange(len(self.q[0])):
        pt.wrap(info, "Species [{: 2d}]:       [{:.2e}, {:.2e}, ..., {:.2e}].".
                format(i, self.d[0,i], self.d[1,i], self.d[-1,i]), 4)
    return "\n".join(info)


class Molecules(object):
  def __init__(self):
    self.molfile = None  # Molecular-properties file
    self.nmol   = 0     # Number of species
    self.name   = None  # Species' name               [nmol]
    self.symbol = None  # Species' symbol             [nmol]
    self.mass   = None  # Species' mass  (gr/mol)     [nmol]
    self.radius = None  # Species' radius (Angstroms) [nmol]
    self.ID     = None  # Species' universal ID       [nmol]

  def __repr__(self):
    info = []
    pt.wrap(info, "Atmospheric species info:")
    pt.wrap(info, "Number of species: {:d}\n"
                  "Species:   ID   Mass      Radius\n"
                  "                (gr/mol)  (Angstrom)".format(self.nmol), 2)
    for i in np.arange(self.nmol):
        pt.wrap(info, "{:>7s}:  {:3d}  {:8.4f}  {:.3f}".
                format(self.symbol[i], self.ID[i],
                       self.mass[i], self.radius[i]/pc.A), 2)
    return "\n".join(info)


class Linetransition(object):
  def __init__(self):
    self.tlifile = None     # Line-transition data file
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
    pt.wrap(info, "Line-transition info:")
    pt.wrap(info, "Number of TLI files:           {:d}".format(self.nTLI), 2)
    pt.wrap(info, "Number of databases (species): {:d}".format(self.ndb),  2)
    for i in np.arange(self.ndb):
        self.db[i].info(2)
    pt.wrap(info, "Number of line transitions:    {:d}\n"
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
    pt.wrap(info, "Database info:")
    pt.wrap(info, "Database name: {:s}".format(self.name),      2)
    pt.wrap(info, "Species' name: {:s}".format(self.molname),   2)
    pt.wrap(info, "Number of isotopes: {:d}".format(self.niso), 2)
    pt.wrap(info, "Isotope correlative index: {:d}".format(self.iiso), 2)
    pt.wrap(info, "Number of temperature samples: {:d} (for partition "
                  "function)".format(self.ntemp), 2)
    pt.wrap(info, "Temperature boundaries (K): [{:.1f}, {:.1f}]".
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
    pt.wrap(info, "Isotopes info:")
    pt.wrap(info, "Number of isotopes: {:d}".format(self.niso), 2)
    pt.wrap(info,
      "Isotope:  Species  Mass      Isotopic   Database  Ext-coefficient\n"
      "                   (gr/mol)  ratio      index     table index", 2)
    for i in np.arange(self.niso):
      pt.wrap(info, "{:>7s}:  {:>7s}  {:8.4f}  {:.3e}       {:3d}  {}".
             format(self.name[i], self.mol.name[self.imol[i]],
               self.mass[i], self.ratio[i], self.dbindex[i], iext[i]), 2)
    # FINDME: Partition function?
    #pt.wrap(info, "Partition Function:", 2)
    #for i in np.arange(self.niso):
    #  pt.wrap(info, "{:>7s}: [{:.2e}, {:.2e}, ..., {:.2e}]".
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
    pt.wrap(info, "Voigt profile info:", 0)
    pt.wrap(info, "Number of Doppler-width samples:  {:d}".format(self.nDop),
           2)
    pt.wrap(info, "Number of Lorentz-width samples:  {:d}".format(self.nLor),
           2)
    pt.wrap(info, "Doppler-width array (cm-1):  [{:.3e}, {:.3e}, ..., {:.3e}]".
               format(self.doppler[0], self.doppler[1], self.doppler[-1]),
           2)
    pt.wrap(info, "Lorentz-width array (cm-1):  [{:.3e}, {:.3e}, ..., {:.3e}]".
               format(self.lorentz[0], self.lorentz[1], self.lorentz[-1]),
           2)
    pt.wrap(info, "Doppler--Lorentz ratio threshold:  {:.3e}".
               format(self.DLratio), 2)
    pt.wrap(info, "Extent covered by a profile in units of Voigt half widths:  "
              "{:.2f}".format(self.extent), 2)
    pt.wrap(info, "Total size of all Voigt profiles (~sum of: 2*size+1):  {:d}".
               format(np.size(self.profile)), 2)
    pt.wrap(info, "Voigt-profile half-sizes [Ndop, Nlor]:", 2)
    pt.wrap(info, "{}".format(self.size), 4)
    pt.wrap(info, "Voigt-profile indices [Ndop, Nlor]:", 2)
    pt.wrap(info, "{}".format(self.index), 4)
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
    pt.wrap(info, "Extinction coefficient info:")
    pt.wrap(info, "Line-transition strength threshold: {:.3e}".
               format(self.ethresh), 2)
    if self.extfile is None:
        pt.wrap(info, "No extinction-coefficient table defined.", 2)
    else:
        pt.wrap(info, "Extinction-coefficient table filename:",  2)
        pt.wrap(info, "'{:s}'".format(self.extfile), 4)
        pt.wrap(info, "Minimum temperature:           {:6.1f} K".
                      format(self.tmin), 4)
        pt.wrap(info, "Maximum temperature:           {:6.1f} K".
                      format(self.tmax), 4)
        pt.wrap(info, "Temperature sampling interval: {:6.1f} K".
                      format(self.tstep), 4)
        pt.wrap(info, "Number of tabulated species:          {:5d}".
                      format(self.nmol), 4)
        pt.wrap(info, "Number of tabulated temperatures:     {:5d}".
                      format(self.ntemp), 4)
        pt.wrap(info, "Number of tabulated layers:           {:5d}".
                      format(self.nlayers), 4)
        pt.wrap(info, "Number of tabulated spectral samples: {:5d}".
                      format(self.nwave), 4)
        pt.wrap(info, "Temperature array (K):   [{:8.1f}, {:8.1f}, ..., "
                "{:8.1f}]".format(self.temp[0], self.temp[1], self.temp[-1]), 4)
        pt.wrap(info, "Partition function at tabulated temperatures:", 4)
        pt.wrap(info, "{}".format(self.z), 6)
        pt.wrap(info, "Species ID array: {:s}".
                      format(str(self.molID).replace("\n", "")), 4)
        pt.wrap(info, "Pressure array: (bar)    [{:.2e}, {:.2e}, ..., {:.2e}]".
                      format(self.press[0]/pc.bar,
                             self.press[1]/pc.bar, self.press[-1]/pc.bar), 4)
        pt.wrap(info, "Wavenumber array (cm-1): [{:8.3f}, {:8.3f}, ..., "
                "{:8.3f}]".format(self.wn[0], self.wn[1], self.wn[-1]), 4)
        np.set_printoptions(formatter={'float': '{: .1e}'.format})
        pt.wrap(info, "Tabulated extinction coefficient (cm2 gr-1)\n"
                      "                       [spec, temp, layer, wave]:", 4)
        pt.wrap(info, "{}".format((self.etable)), 4)
    if self.ec is not None:
        np.set_printoptions(formatter={'float': '{: .2e}'.format})
        pt.wrap(info, "\nLine-transition extinction coefficient for the "
                      "atmospheric model (cm-1) [layer, wave]:", 2)
        pt.wrap(info, "{}".format((self.ec)), 2)
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
    pt.wrap(info, "Cross-section extinction info:")
    pt.wrap(info, "Number of CS files: {:d}".format(self.nfiles), 2)
    for i in np.arange(self.nfiles):
        pt.wrap(info, "CS file: '{:s}':".format(self.files[i]), 2)
        pt.wrap(info, "Species: {:s}".
                format("-".join(self.molecules[i,0:self.nmol[i]])), 4)
        pt.wrap(info, "Number of temperatures:       {:4d}".
                format(self.ntemp[i]), 4)
        pt.wrap(info, "Number of wavenumber samples: {:4d}".
                format(self.nwave[i]), 4)
        pt.wrap(info, "Temperature array (K):  {}".
                format(str(self.temp[i]).replace("\n", "")), 4, 6)
        pt.wrap(info, "Wavenumber array (cm-1): [{:7.1f}, {:7.1f}, ..., "
                "{:7.1f}]".format(self.wavenumber[i][0], self.wavenumber[i][1],
                                  self.wavenumber[i][-1]), 4)
        np.set_printoptions(formatter={'float': '{: .1e}'.format})
        pt.wrap(info, "Tabulated CS extinction coefficient (cm-1 amagat-{:d}) "
                      "[layer, wave]:".format(self.nmol[i]), 4)
        pt.wrap(info, "{}".format((self.ec)),                6)
        np.set_printoptions(formatter=None)
    pt.wrap(info, "\nMinimum and maximum covered temperatures (K): "
                  "[{:.1f}, {:.1f}]".format(self.tmin, self.tmax), 2)
    # FINDME iabsorp ?
    if self.ec is not None:
        np.set_printoptions(formatter={'float': '{: .2e}'.format})
        pt.wrap(info, "CS extinction coefficient for the "
                      "atmospheric model (cm-1) [layer, wave]:", 2)
        pt.wrap(info, "{}".format((self.ec)), 2)
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
    pt.wrap(info, "Optical depth info:")
    pt.wrap(info, "Ray-path geometry:  {:s}".format(self.path), 2)
    pt.wrap(info, "Maximum optical depth to calculate:  {:.2f}".
               format(self.maxdepth), 2)
    if self.ec is not None:
        np.set_printoptions(formatter={'float': '{: .2e}'.format})
        pt.wrap(info, "Total atmospheric-model extinction coefficient (cm-1) "
                      "[layer, wave]:", 2)
        pt.wrap(info, "{}".format((self.ec)), 2)
        np.set_printoptions(formatter=None)
    if self.depth is not None:
        pt.wrap(info, "Layer index where the optical depth reached maxdepth:",2)
        pt.wrap(info, "{}".format(self.ideep), 4)
        np.set_printoptions(formatter={'float': '{: .1f}'.format})
        # Raypath for transit geometry:
        if self.path == "transit":
            pt.wrap(info, "\nDistance (km) along the raypath over each layer "
                          "(outside-in) for each impact parameter:", 2)
            pt.wrap(info, "IP[  1] ({:.1f} km): {}".format(
                    self.atm.radius[1]/pc.km, self.raypath[1]/pc.km), 4)
            pt.wrap(info, "IP[  2] ({:.1f} km): {}".format(
                    self.atm.radius[2]/pc.km, self.raypath[2]/pc.km), 4)
            pt.wrap(info, "IP[  3] ({:.1f} km): {}".format(
                    self.atm.radius[3]/pc.km, self.raypath[3]/pc.km), 4)
            pt.wrap(info, "...", 4)
            pt.wrap(info, "IP[{:3d}] ({:.1f} km): {}".format(
                    len(self.atm.radius),
                    self.atm.radius[-1]/pc.km,
                    str(self.raypath[-1]/pc.km)).replace("\n", ""), 4, 6)
            pt.wrap(info, "\nOptical depth for each impact parameter "
                          "(outside-in) for each wavenumber:", 2)
        # Raypath for eclipse geometry:
        elif self.path == "eclipse":
            pt.wrap(info, "\nDistance over each layer along a normal-incident "
                          "raypath (km):  {}".format(
                          str(self.raypath/pc.km).replace("\n", "")), 2, 4)
            pt.wrap(info, "\nOptical depth over each layer (outside-in) along "
                          "a normal-incident raypath for each wavenumber:", 2)
        # Print optical depth:
        np.set_printoptions(formatter={'float': '{: .1e}'.format})
        pt.wrap(info, "At {:7.1f} cm-1:  {}".format(self.spec.wn[0],
               str(self.depth[0:self.ideep[0]+1,0]).replace("\n","")), 4, 6)
        pt.wrap(info, "...", 4)
        index = self.spec.nwave/2
        pt.wrap(info, "At {:7.1f} cm-1:  {}".format(self.spec.wn[index],
                str(self.depth[0:self.ideep[index]+1,index]).replace("\n","")),
                4, 6)
        pt.wrap(info, "...", 4)
        index = self.spec.nwave-1
        pt.wrap(info, "At {:7.1f} cm-1:  {}".format(self.spec.wn[index],
                str(self.depth[0:self.ideep[index]+1,index]).replace("\n","")),
                4, 6)
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
    self.rmodels = ["pt", "rad", "mol", "ray", "haze", "cloud", "patchy"]
    self.retflag = None  # Flags for models to be included for retrieval
    self.nparams = 0     # Number of free parameters
    self.tlow    = None  # Lower-temperature retrieval boundary
    self.thigh   = None  # Higher-temperature retrieval boundary
    self.itemp   = None  # Temperature-model parameter indices
    self.irad    = None  # Reference-radius model parameter index
    self.iabund  = None  # Abundance-model parameter indices
    self.iray    = None  # Haze-model parameter indices
    self.ihaze   = None  # Haze-model parameter indices
    self.icloud  = None  # Cloud-model parameter indices
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
