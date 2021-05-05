# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'Spectrum',
    'Atm',
    'Molecules',
    'Linetransition',
    'Isotopes',
    'Voigt',
    'Extinction',
    'Cross',
    'Optdepth',
    'Cloud',
    'Rayleigh',
    'Alkali',
    'Observation',
    'Physics',
    'Retrieval',
    ]

import os
import numpy as np
import scipy.constants as sc

from .. import tools     as pt
from .. import constants as pc


class Spectrum(object):
  def __init__(self):
      self.specfile  = None  # Transmission/Emission spectrum file
      # Wavenumber:
      self.nwave     = None  # Number of wavenumber spectral samples
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

  def __str__(self):
      fmt = {'float': '{: .3e}'.format}
      fw = pt.Formatted_Write()
      fw.write('Spectral information:')
      fw.write('Wavenumber internal units: cm-1')
      fw.write('Wavelength internal units: cm')
      fw.write('Wavelength display units (wlunits): {:s}', self.wlunits)
      fw.write('Low wavenumber boundary (wnlow):   {:10.3f} cm-1  '
               '(wlhigh = {:6.2f} {})', self.wnlow,
               self.wlhigh/pt.u(self.wlunits), self.wlunits)
      fw.write('High wavenumber boundary (wnhigh): {:10.3f} cm-1  '
               '(wllow  = {:6.2f} {})', self.wnhigh,
               self.wllow/pt.u(self.wlunits), self.wlunits)
      fw.write('Number of samples (nwave): {:d}', self.nwave)
      if self.resolution is None:
          fw.write('Sampling interval (wnstep): {:.3f} cm-1', self.wnstep)
      else:
          fw.write('Spectral resolving power (resolution): {:.1f}',
              self.resolution)
      fw.write('Wavenumber array (wn, cm-1):\n    {}', self.wn,
          fmt={'float': '{: .3f}'.format})
      fw.write('Oversampling factor (wnosamp): {:d}', self.wnosamp)

      if self._rt_path in pc.emission_rt:
          if self.quadrature is not None:
              fw.write('Number of Gaussian-quadrature points for intensity '
                  'integration into flux (quadrature): {}', self.quadrature)
          fw.write('\nIntensity zenithal angles (raygrid, degree): {}',
              self.raygrid/sc.degree, prec=3)
          fw.write('raygrid internal units: radian')
          if self.intensity is not None:
              fw.write('Intensity spectra (intensity, erg s-1 cm-2 sr-1 cm):')
              for intensity in self.intensity:
                  fw.write('    {}', intensity, fmt=fmt, edge=3)
          fw.write('Emission spectrum (spectrum, erg s-1 cm-2 cm):\n'
                   '    {}', self.spectrum, fmt=fmt, edge=3)
      elif self._rt_path in pc.transmission_rt:
          fw.write('\nModulation spectrum, (Rp/Rs)**2 (spectrum):\n'
                   '    {}', self.spectrum, fmt=fmt, edge=3)
      return fw.text


class Atm(object):
  def __init__(self):
      self.refpressure = None  # Pressure reference level
      self.radstep = None     # Radius sampling interval
      self.radlow  = None     # Lowest radius boundary
      self.radhigh = None     # Highest radius boundary
      self.ptop    = None     # Lowest pressure boundary
      self.pbottom = None     # Highest pressure boundary
      self.atmfile = None     # Atmopheric-model file
      self.qunits  = None     # Input abundance units ('volume' or 'mass')
      self.runits  = None     # Input radius units
      self.punits  = None     # Input pressure units
      self.tunits  = 'kelvin' # Input temperature units
      self.nlayers = None     # Number of layers
      self.radius  = None     # Radius array (cm)            [layers]
      self.press   = None     # Pressure array (barye)       [layers]
      self.temp    = None     # Temperature array (K)        [layers]
      self.mm      = None     # Mean molecular mass (gr/mol) [layers]
      self.q       = None     # Molecular abundances         [layers, nmol]
      self.d       = None     # Molecular densities          [layers, nmol]
      self.tmodel  = None
      self.tpars   = None
      self.rtop    = 0        # Index of topmost layer (within Hill radius)

  def __str__(self):
      fmt = {'float': '{: .3e}'.format}
      fw = pt.Formatted_Write()
      press  = self.press/pt.u(self.punits)
      radius = self.radius/pt.u(self.runits)
      fw.write('Atmospheric model information:')
      fw.write("Atmospheric file name (atmfile): '{}'", self.atmfile)
      fw.write('Number of layers (nlayers): {:d}', self.nlayers)

      fw.write('\nPressure display units (punits): {}', self.punits)
      fw.write('Pressure internal units: barye')
      fw.write('Pressure at top of atmosphere (ptop):        {:.2e} {}',
          self.ptop/pt.u(self.punits), self.punits)
      fw.write('Pressure at bottom of atmosphere (pbottom):  {:.2e} {}',
          self.pbottom/pt.u(self.punits), self.punits)
      fw.write('Reference pressure at rplanet (refpressure): {:.2e} {}',
          self.refpressure/pt.u(self.punits), self.punits)
      fw.write('Pressure profile (press, {}):\n    {}', self.punits, press,
          fmt=fmt, edge=3)

      fw.write('\nRadius display units (runits): {}', self.runits)
      fw.write('Radius internal units: cm', self.runits)
      fw.write('Radius model name (rmodelname): {}', self.rmodelname)
      if self.radstep is not None:
          fw.write('Radius step size (radstep): {} {}',
              self.radstep/pt.u(self.runits), self.runits, prec=3, edge=3)
      if self.radhigh is not None:
          fw.write('Radius at top of atmosphere (radhigh): {} {}',
              self.radhigh/pt.u(self.runits), self.runits, prec=3, edge=3)
          fw.write('Radius at bottom of atmosphere (radlow): {} {}',
              self.radlow/pt.u(self.runits), self.runits, prec=3, edge=3)
      fw.write('Radius profile (radius, {}):\n    {}', self.runits, radius,
          prec=4, edge=3, lw=800)

      fw.write('\nTemperature units (tunits): {}', self.tunits)
      fw.write('Temperature model name (tmodelname): {}', self.tmodelname)
      if self.tmodel is not None:
          fw.write('  tmodel parameters (tpars): {}', self.tpars)
      fw.write('Temperature profile (temp, K):\n    {}', self.temp,
          fmt={'float': '{:9.3f}'.format}, edge=3)

      fw.write('\nMean molecular mass (mm, amu):\n    {}', self.mm,
          fmt={'float': '{:8.4f}'.format}, edge=3)
      fw.write('\nAbundance units (qunits): {}', self.qunits)
      fw.write('Abundance internal units: mole mixing fraction')
      fw.write('Number of atmospheric species: {:d}', len(self.q[0]))
      if self.molmodel is not None:
          molpars = [None for _ in self.molmodel] if self.molpars is None \
                    else self.molpars
          fw.write('Abundance models:\n'
                   '  ifree  molfree     molmodel    molpars')
          for molvals in zip(self.ifree, self.molfree, self.molmodel, molpars):
              fw.write('     {:2d}  {:10s}  {:10s}  {}', *molvals)
          fw.write('Bulk species:\n'
                   '  ibulk  bulk')
          for ibulk, bulk in zip(self.ibulk, self.bulk):
              fw.write('     {:2d}  {:10s}', ibulk, bulk)

      fw.write('Abundance profiles (q, mole mixing fraction):')
      for i, q in enumerate(self.q.T):
          fw.write('    species [{:2d}]:   {}', i, q,    fmt=fmt, edge=2)
      fw.write('Density profiles (d, molecules cm-3):')
      for i, dens in enumerate(self.d.T):
          fw.write('    species [{:2d}]:   {}', i, dens, fmt=fmt, edge=2)

      return fw.text


class Molecules(object):
  def __init__(self):
      self.molfile = None  # Molecular-properties file
      self.nmol   = 0     # Number of species
      self.name   = None  # Species' name               [nmol]
      self.symbol = None  # Species' symbol             [nmol]
      self.mass   = None  # Species' mass  (gr/mol)     [nmol]
      self.radius = None  # Species' radius (Angstroms) [nmol]

  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write('Atmospheric species information:')
      fw.write('Number of species (nmol): {:d}\n', self.nmol)

      fw.write('\nMolecule    Mass       Radius\n'
                 '            g/mol      Angstrom\n'
                 '(name)      (mass)     (radius)  ')
      for i in range(self.nmol):
          fw.write('  {:8s}  {:8.4f}  {:10.3f}',
              self.name[i], self.mass[i], self.radius[i]/pc.A)
      fw.write("Molecular data taken from (molfile): '{}'", self.molfile)
      return fw.text


class Linetransition(object):
  def __init__(self):
      self.tlifile = None     # Line-transition data file
      self.dblist  = None
      self.ndb     = 0        # Number of data bases
      self.db      = []       # Data base objects
      self.ntransitions = 0   # Number of line transitions
      self.tmin    = -np.inf  # Minimum temperature sampled by all TLI files
      self.tmax    =  np.inf  # Maximum temperature sampled by all TLI files
      self.wn      = np.array([], np.double)  # Line wavenumber
      self.elow    = np.array([], np.double)  # Line lower energy level
      self.gf      = np.array([], np.double)  # Line gf value
      self.isoid   = np.array([], np.int)     # Line isotope index

  def clone_new(self, pyrat):
      """Return a new LT instance (as returned by Pyrat.__init__)."""
      lt = Linetransition()
      lt.tlifile = pyrat.lt.tlifile
      lt.dblist  = pyrat.lt.dblist
      return lt

  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write('Line-transition information:')
      if self.tlifile is None:
          fw.write('No input TLI files.')
          return fw.text
      fw.write('Input TLI files (tlifile): {}', self.tlifile)
      fw.write('Number of databases (ndb): {:d}', self.ndb)
      for db in self.db:
          fw.write('\n'+str(db))
      fw.write('\nTotal number of line transitions (ntransitions): {:,d}\n'
          'Minimum and maximum temperatures (tmin, tmax): [{:.1f}, {:.1f}] K',
          self.ntransitions, self.tmin, self.tmax)
      fw.write('Line-transition isotope IDs (isoid):\n    {}', self.isoid,
          edge=7)
      fw.write('Line-transition wavenumbers (wn, cm-1):\n    {}',
          self.wn,   fmt={'float':'{:.3f}'.format}, edge=3)
      fw.write('Line-transition lower-state energy (elow, cm-1):\n    {}',
          self.elow, fmt={'float':'{: .3e}'.format}, edge=3)
      fw.write('Line-transition gf (gf, cm-1):\n    {}',
          self.gf,   fmt={'float':'{: .3e}'.format}, edge=3)
      return fw.text


class Database(object):
  def __init__(self):
      self.name    = None  # Data base name
      self.molname = None  # Molecule name
      self.niso    = None  # Number of isotopes in database
      self.iiso    = None  # Isotope correlative index
      self.ntemp   = None  # Number of temperature samples
      self.temp    = None  # Temperature array
      self.z       = None  # Isotopes' partition function array [niso, ntemp]

  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write('Database name (name): {:s}', self.name)
      fw.write('Species name (molname):  {:s}', self.molname)
      fw.write('Number of isotopes (niso): {:d}', self.niso)
      fw.write('Isotope correlative index (iiso): {:d}', self.iiso)
      fw.write('Number of temperature samples (ntemp): {:d}', self.ntemp)
      fw.write('Temperature (temp, K):\n    {}', self.temp, prec=3, edge=3)
      fw.write('Partition function for each isotope (z):')
      for z in self.z:
          fw.write('    {}', z, fmt={'float':'{: .3e}'.format}, edge=3)
      return fw.text


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

  def __str__(self):
      fw = pt.Formatted_Write()
      if self.iext is None:
          iext = [None for _ in self.name]
      else:
          iext = self.iext
      fw.write('Isotopes information:')
      fw.write('Number of isotopes (niso): {:d}', self.niso)
      fw.write(
          '\nIsotope  Molecule      Mass    Isotopic   Database   Extinc-coeff'
          '\n            index     g/mol       ratio      index    index'
          '\n (name)    (imol)    (mass)     (ratio)   (dbindex)  (iext)')
      for i in np.arange(self.niso):
          fw.write('{:>7s}  {:8d}  {:8.4f}   {:.3e}   {:8d}   {}',
              self.name[i], self.imol[i], self.mass[i], self.ratio[i],
              self.dbindex[i], iext[i])
      fw.write('Partition function evaluated at atmosperic layers (z):')
      for z in self.z:
          fw.write('    {}', z, fmt={'float':'{: .3e}'.format}, edge=3)
      return fw.text


class Voigt(object):
  def __init__(self):
      self.dmin     = None  # Minimum Doppler width sampled
      self.dmax     = None  # Maximum Doppler width sampled
      self.ndop     = None  # Number of Doppler-width samples
      self.lmin     = None  # Minimum Lorentz width sampled
      self.lmax     = None  # Maximum Lorentz width sampled
      self.nlor     = None  # Number of Lorentz-width samples
      self.doppler  = None  # Doppler-width sample array [ndop]
      self.lorentz  = None  # Lorentz-width sample array [nlor]
      self.dlratio  = None  # Doppler-Lorentz ratio threshold
      self.extent   = None  # Extent covered by the profile (in number of HWHM)
      self.cutoff   = None  # Max cutoff extent (in cm-1)
      self.profile  = None  # Voigt profile [sum(2*size+1)]
      self.size     = None  # Profile wavenumber half-size [ndop, nlor]
      self.index    = None  # Index where each profile starts [ndop, nlor]

  def __str__(self):
      fw = pt.Formatted_Write(fmt={'float':'{: .3e}'.format}, edge=3)
      fw.write('Voigt-profile information:')
      fw.write('\nNumber of Doppler-width samples (ndop): {:d}', self.ndop)
      fw.write('Number of Lorentz-width samples (nlor): {:d}', self.nlor)
      fw.write('Doppler HWHM (doppler, cm-1):\n    {}', self.doppler)
      fw.write('Lorentz HWMH (lorentz, cm-1):\n    {}', self.lorentz)
      fw.write('Doppler--Lorentz ratio threshold (dlratio): {:.3e}',
          self.dlratio)
      fw.write("\nVoigt-profiles' extent (extent, in HWHMs): {:.1f}",
          self.extent)
      fw.write("Voigt-profiles' cutoff extent (cutoff in cm-1): {:.1f}",
          self.cutoff)
      fw.write('Voigt-profile half-sizes (size) of shape [ndop, nlor]:\n{}',
          self.size, edge=2)
      fw.write('Voigt-profile indices (index) of shape [ndop, nlor]:\n{}',
          self.index, edge=2)

      index, size = self.index[0,0], 2*self.size[0,0]+1
      fw.write('\nVoigt profiles:\n  profile[ 0, 0]: {}',
          self.profile[index:index+size],
          fmt={'float':'{: .5e}'.format}, edge=2)
      index =  self.index[self.ndop-1,self.nlor-1]
      size  = 2*self.size[self.ndop-1,self.nlor-1] + 1
      fw.write('  ...\n  profile[{:2d},{:2d}]: {}', self.ndop-1, self.nlor-1,
          self.profile[index:index+size],
          fmt={'float':'{: .5e}'.format}, edge=2)
      return fw.text


class Extinction(object):
  def __init__(self):
      self.ec      = None # Molecular line-transition extinction coefficient
                          #  in cm-1 [nlayers, nwave]
      self.ethresh = None # Extinction-coefficient threshold
      self.extfile = None # Extinction-coefficient table filename
      self.etable  = None # Tabulated extinction coefficient (cm-2 molecule-1)
                          # with shape [nmol, nlayer, ntemp, nwave]
      self.tmin    = None # Minimum temperature to sample
      self.tmax    = None # Maximum temperature to sample
      self.tstep   = None # Temperature-sample step interval
      self.z       = None # Partition function at tabulated temperatures
                          #   [niso, ntemp]
      self.nspec   = 0    # Number of species
      self.ntemp   = 0    # Number of temperature samples
      self.nlayers = 0    # Number of pressure layers
      self.nwave   = 0    # Number of wavenumber spectral samples

      self.temp    = None # Tabulated temperatures
      self.press   = None # Tabulated pressures
      self.wn      = None # Tabulated wavenumber

  def __str__(self):
      fmt = {'float': '{:.2e}'.format}
      fw = pt.Formatted_Write()
      fw.write('Extinction-coefficient information:')
      fw.write('Line-transition strength threshold (ethresh): {:.2e}',
          self.ethresh)
      if self.ec is not None:
          fw.write('\nLBL extinction coefficient for the atmospheric model '
              '(ec, cm-1) [layer, wave]:\n{}', self.ec, fmt=fmt)
      extfile = ['None'] if self.extfile is None else self.extfile
      fw.write("Extinction-coefficient table filename(s) (extfile): {}",
          '\n    '.join(extfile))
      if self.extfile is None:
          return fw.text
      fw.write('Minimum temperature (tmin, K): {:6.1f}', self.tmin)
      fw.write('Maximum temperature (tmax, K): {:6.1f}', self.tmax)
      fw.write('Temperature sampling interval (tstep, K): {:6.1f}', self.tstep)
      fw.write('\nNumber of species (nspec):          {:5d}', self.nspec)
      fw.write('Number of temperatures (ntemp):     {:5d}', self.ntemp)
      fw.write('Number of layers (nlayers):         {:5d}', self.nlayers)
      fw.write('Number of spectral samples (nwave): {:5d}', self.nwave)
      fw.write('\nSpecies array (species): {}', self.species)
      fw.write('Temperature array (temp, K):\n   {}', self.temp)
      fw.write('Partition function (z): {}', self.z)
      fw.write('Pressure array (press, bar):\n   {}', self.press/pc.bar,
               fmt=fmt, lw=8000)
      fw.write('Wavenumber array (wn, cm-1):\n    {}', self.wn, prec=4)
      fw.write('Tabulated extinction coefficient (etable, cm2 molecule-1) '
               'of shape\n    [nmol, ntemp, nlayers, nwave]:\n{}', self.etable,
               fmt=fmt, edge=2)
      return fw.text


class Cross(object):
  def __init__(self):
      self.files      = None  # CS file names
      self.nfiles     = 0     # Number of files read
      self.tmin       = 0.0   # Minimum temperature sampled by all CS files
      self.tmax       = 1e6   # Maximum temperature sampled by all CS files
      self.molecules  = []    # Species involved for each file
      self.temp       = []    # Temperature sampling (in Kelvin)
      self.wavenumber = []    # Wavenumber sampling (in cm-1)
      self.absorption = []    # CS extinction (in cm-1 amagat-2)
      self.iabsorp    = []    # wn-interpolated CS extinction (in cm-1 amagat-2)
      self.iz         = []    # Second derivatives of iabsorp
      self.iwnlo      = []    # Lower-wavenumber index for interpolation
      self.iwnhi      = []    # Upper-wavenumber index for interpolation
      self.ec         = None  # Interpolated CS extinction coefficient
                              #  in cm-1 [nlayer, nwave]

  def clone_new(self, pyrat):
      """Return a new Cross instance (as returned by Pyrat.__init__)."""
      cs = Cross()
      cs.files = pyrat.cs.files
      return cs

  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write('Cross-section extinction information:')
      fw.write('Number of cross-section files (nfiles): {:d}', self.nfiles)
      for i in range(self.nfiles):
          fw.write("\nCross-section file name (files[{}]): '{:s}'", i,
              self.files[i])
          fw.write('Cross-section species (molecules): {:s}',
              '-'.join(self.molecules[i]))
          fw.write('Number of temperature samples: {}', len(self.temp[i]))
          fw.write('Number of wavenumber samples: {}', len(self.wavenumber[i]))
          fw.write('Temperature array (temp, K):\n  {}', self.temp[i],
              prec=1, lw=800, edge=50)
          fw.write('Wavenumber array (wavenumber, cm-1):\n  {}',
              self.wavenumber[i], fmt={'float':'{:.1f}'.format}, lw=80, edge=3)
          fw.write('Input extinction coefficient (absorption, cm-1 '
              'amagat-{:d}):\n{}', len(self.molecules[i]), self.absorption[i],
              fmt={'float': '{: .2e}'.format}, edge=2)
      fw.write('\nMinimum and maximum temperatures (tmin, tmax) in K: '
                '[{:.1f}, {:.1f}]', self.tmin, self.tmax)
      if self.ec is not None:
          fw.write('Atmospheric-model extinction coefficient (ec, cm-1):\n{}',
                   self.ec, fmt={'float':'{: .2e}'.format})
      return fw.text


class Cloud(object):
  def __init__(self):
      self.models  = []    # List of cloud models
      self.ec      = None  # Cloud extinction coefficient
      self.fpatchy = None  # Patchy-cloud fraction
      self.pars    = None  # Input cloud parameters

  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write('Cloud-opacity models (models):')
      for model in self.models:
          fw.write('\n' + str(model))
      fw.write('\nPatchiness fraction (fpatchy): {:.3f}', self.fpatchy)
      fw.write('Total atmospheric cloud extinction-coefficient '
               '(ec, cm-1):\n{}', self.ec, fmt={'float':' {:.3e}'.format})
      return fw.text


class Rayleigh(object):
  def __init__(self):
      self.models  = []    # List of Rayleigh models
      self.ec      = None  # Rayleigh extinction coefficient
      self.pars    = None  # Input rayleigh parameters

  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write('Rayleigh-opacity models (models):')
      for model in self.models:
          fw.write('\n' + str(model))
      fw.write('\nTotal atmospheric Rayleigh extinction-coefficient '
               '(ec, cm-1):\n{}', self.ec, fmt={'float': '{: .3e}'.format})
      return fw.text


class Alkali(object):
  def __init__(self):
      self.models  = []    # List of alkali models
      self.ec      = None  # Alkali extinction coefficient
      self.cutoff  = 4500  # Profiles cutoff from line center (cm)

  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write('Alkali-opacity models (models):')
      for model in self.models:
          fw.write('\n' + str(model))
      fw.write('\nTotal atmospheric alkali extinction-coefficient '
               '(ec, cm-1):\n{}', self.ec, fmt={'float': '{: .3e}'.format})
      return fw.text


class Optdepth(object):
  def __init__(self):
      self.maxdepth = None  # Maximum optical depth to calculate
      self.rt_path  = None  # Radiative=transfer observing geometry
      self.ec       = None  # Total extinction coefficient [nlayers, nwave]
      self.epatchy  = None  # Cloudy extinction coefficient for patchy model
      self.raypath  = []    # Distance along ray path  [nlayers]
      self.depth    = None  # Optical depth at raypath [nlayers, nwave]
      self.pdepth   = None  # Cloudy optical depth for patchy model
      self.B        = None  # Blackbody Planck emission [nlayers, nwave]
      self.ideep    = None  # Layer index where depth reached maxdepth [nwave]

  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write('Optical depth information:')
      fw.write('Observing geometry (rt_path): {}', self.rt_path)
      if self.ec is not None:
          fw.write('Total atmospheric extinction coefficient (ec, cm-1) [layer'
                   ', wave]:\n{}', self.ec, fmt={'float':'{: .3e}'.format})
      if self.depth is None:
          fw.write('\nMaximum optical depth to calculate (maxdepth): {:.2f}',
              self.maxdepth)
          return fw.text

      ideepest = np.amax(self.ideep)
      if self.rt_path in pc.transmission_rt:
          fw.write('\nDistance along the ray path across each layer '
                   '(outside-in) at each impact parameter (raypath, km):')
          with np.printoptions(formatter={'float':'{:.1f}'.format},threshold=6):
              fw.write('    IP[{:3d}]: {}', 1, self.raypath[1]/pc.km)
              fw.write('    IP[{:3d}]: {}', 2, self.raypath[2]/pc.km)
              fw.write('    IP[{:3d}]: {}', 3, self.raypath[3]/pc.km)
              fw.write('    ...')
              fw.write('    IP[{:3d}]: {}', len(self.raypath)-1,
                       self.raypath[len(self.raypath)-1]/pc.km)
          od_text = ('\nOptical depth at each impact parameter, down to '
                     'max(ideep) (depth):')
      elif self.rt_path in pc.emission_rt:
          fw.write('\nDistance across each layer along a normal ray path '
              '(raypath, km):\n    {}', self.raypath/pc.km,
              fmt={'float':'{:.1f}'.format}, edge=4)
          od_text = ('\nOptical depth at each layer along a normal ray '
                     'path into the planet, down to max(ideep) (depth):')

      fw.write('\nMaximum optical depth to calculate (maxdepth): {:.2f}',
               self.maxdepth)
      fw.write('Layer index where the optical depth reaches maxdepth (ideep):'
               '\n    {}', self.ideep, fmt={'int': '{:3d}'.format}, edge=7)
      fw.write('Maximum ideep (deepest layer reaching maxdepth): {}', ideepest)

      if self.rt_path in pc.emission_rt:
          fw.write('\nPlanck emission down to max(ideep) (B, erg s-1 cm-2 '
                   'sr-1 cm):\n{}', self.B[0:ideepest+1],
                   fmt={'float':'{: .3e}'.format})

      fw.write('{}\n{}', od_text, self.depth[0:ideepest+1],
               fmt={'float':'{: .3e}'.format})
      return fw.text


class Observation(object):
  def __init__(self):
      self.ndata     = 0     # Number of data points
      self.data      = None  # Transit or eclipse data point
      self.uncert    = None  # Data's 1-sigma uncertainty
      self.nfilters  = 0     # Number of filter bands
      self.filters   = None  # Observing filter filename
      self.bandidx   = None  # Band wavenumber indices
      self.bandtrans = None  # Band-interpolated transmission function
      self.starflux  = None  # Band-interpolated stellar flux
      self.bandflux  = None  # Band-integrated flux
      self.bandwn    = None  # Filter mean wavenumber
      self.units     = 'none' # Data units

  def __str__(self):
      units = pt.u(self.units)
      fw = pt.Formatted_Write()
      fw.write('Observing information:')
      if self.data is not None or self.filters is not None:
          fw.write('Data/bandflux display units (units): {}', self.units)
          fw.write('Data/bandflux internal units: none')
      fw.write('Number of data points (ndata): {}', self.ndata)
      if self.data is not None:
          fw.write('        Data  Uncertainty   Wavenumber  Wavelength\n'
                   '     {:>7s}      {:>7s}         cm-1          um\n'
                   '      (data)     (uncert)     (bandwn)',
                   self.units, self.units)
          for data, uncert, bandwn in zip(self.data, self.uncert, self.bandwn):
              fw.write('  {:10.5f}   {:10.5f}    {:9.2f}  {:10.3f}',
              data/units, uncert/units, bandwn, 1.0/(bandwn*pc.um))

      fw.write('\nNumber of filter pass bands (nfilters): {}', self.nfilters)
      if self.filters is None:
          return fw.text
      fw.write('Wavenumber  Wavelength    Bandflux  Filter name\n'
               '      cm-1          um     {:>7s}\n'
               '  (bandwn)              (bandflux)  (filters)', self.units)
      for filter,bandwn,bflux in zip(self.filters, self.bandwn, self.bandflux):
          fw.write(' {:9.2f}  {:10.3f}  {:10.5f}  {:s}', bandwn,
              1.0/(bandwn*pc.um), bflux/units, os.path.basename(filter))
      # TBD: Do I want to show bandidx, bandtrans, and starflux?
      return fw.text


class Retrieval(object):
  def __init__(self):
      self.retflag = None  # Flags for models to be included for retrieval
      self.nparams = 0     # Number of free parameters
      self.tlow    = None  # Lower-temperature retrieval boundary
      self.thigh   = None  # Higher-temperature retrieval boundary
      self.itemp   = None  # Temperature-model parameter indices
      self.irad    = None  # Reference-radius model parameter index
      self.imol    = None  # Abundance-model parameter indices
      self.iray    = None  # Rayleigh-model parameter indices
      self.icloud  = None  # Cloud-model parameter indices
      self.ipatchy = None  # Patchy-model parameter index
      self.imass   = None
      self.posterior = None
      self.bestp     = None
      self.spec_best = None
      self.spec_low1 = None
      self.spec_low2 = None
      self.spec_high1 = None
      self.spec_high2 = None
      self.pnames   = []   # Model parameter names (screen)
      self.texnames = []   # Model parameter names (figures)

  def __str__(self):
      {'float':'{: .3e}'.format}
      flags = []
      for flag in self.retflag:
          flags += [flag for _ in getattr(self, 'i'+flag)]

      fw = pt.Formatted_Write()
      fw.write('Retrieval information:')
      if self.params is None:
          fw.write('No retrieval parameters set.')
          return fw.text

      pmin = [None for _ in self.params] if self.pmin is None else self.pmin
      pmax = [None for _ in self.params] if self.pmax is None else self.pmax
      psteps = [None for _ in self.params] if self.pstep is None \
                else self.pstep

      fw.write('  Parameter name        value        pmin        pmax'
               '       pstep  Model type')
      fw.write('  {:15}  {:>10}  {:>10}  {:>10}  {:>10}  {}',
          '(pnames)', '(params)', '(pmin)', '(pmax)', '(pstep)', '(retflag)')
      for pname, param, min, max, pstep, flag in zip(self.pnames,
              self.params, pmin, pmax, psteps, flags):
          fw.write('  {:15s}  {:10.3e}  {:10.3e}  {:10.3e}  {:10.3e}  {}',
              pname, param, min, max, pstep, flag)

      if self.prior is not None:
          fw.write('\nParameter name     Prior')
          for i, pname in enumerate(self.pnames):
              if psteps[i] == 0.0:
                  fw.write('  {:15s}  Fixed at  {:10.3e}', pname,self.params[i])
              elif psteps[i] < 0.0:
                  j = -int(psteps[i]) - 1
                  fw.write('  {:15s}  Shared with  {}', pname, self.pnames[j])
              elif self.priorlow[i]==0 and self.priorup[i]==0:
                  fw.write('  {:15s}  Uniform between     [{:10.3e}, {:10.3e}]',
                      pname, pmin[i], pmax[i])
              else:
                  fw.write('  {:15s}  Gaussian  {:10.3e} -{:.3e}  {:+.3e}',
                      pname, self.prior[i], self.priorlow[i], self.priorup[i])

      fw.write('\nRetrieval algorithm (sampler): {}', self.sampler)
      if self.sampler is None:
          return fw.text
      fw.write('Number of retrieval samples (nsamples): {:,}', self.nsamples)
      # if self.sampler == 'snooker':
      fw.write('Number of parallel chains (nchains):   {}', self.nchains)
      fw.write('Number of burned-in samples (burnin):  {:,}', self.burnin)
      fw.write('Thinning factor (thinning): {}', self.thinning)
      if self.grbreak > 0.0:
          fw.write('Gelman-Rubin convergence criterion to stop (grbreak): {}',
              self.grbreak)
          if self.grnmin > 1:
              fw.write('Minimum number of samples before GR stop (grnmin): {}',
                  int(self.grnmin))
          else:
              fw.write('Minimum fraction of samples before GR stop (grnmin): '
                  ' {}', self.grnmin)

      fw.write('\nUpper boundary for sum of metal abundances (qcap): {}',
          self.qcap)
      fw.write('Temperature upper boundary (tlow, K):  {:6.1f}', self.tlow)
      fw.write('Temperature lower boundary (thigh, K): {:6.1f}', self.thigh)

      fw.write('\nRetrieval posterior file (mcmcfile): {}', self.mcmcfile)
      if self.posterior is not None:
          fw.write('\nParameter name     Best-fit   Posterior distribution '
                   'of shape [{},{}]\n'
                   '                   (bestp)    (posterior)',
                   *self.posterior.shape)
          post = iter(self.posterior.T)
          for pname, bestp, pstep in zip(self.pnames, self.bestp, psteps):
              if pstep > 0:
                  fw.write('  {:15} {:10.3e}  {}', pname, bestp, next(post),
                      fmt={'float':'{: .3e}'.format}, edge=2)
              else:
                  fw.write('  {:15} {:10.3e}', pname, bestp)
          fw.write('\nBest-fit spectrum (spec_best):\n    {}', self.spec_best,
              fmt={'float':'{: .3e}'.format})
      return fw.text


def none_div(a, b):
    if a is None:
        return None
    return a/b


class Physics(object):
  def __init__(self):
      self.tstar    = None  # Stellar effective temperature
      self.rstar    = None  # Stellar radius
      self.mstar    = None  # Stellar mass
      self.gstar    = None  # Stellar surface gravity
      self.rplanet  = None  # Planetary radius
      self.mplanet  = None  # Planetary mass
      self.gplanet  = None  # Planetary surface gravity
      self.smaxis   = None  # Orbital semi-major axis
      self.rhill    = np.inf  # Planetary Hill radius
      self.starspec = None  # Stellar spectrum filename
      self.kurucz   = None  # Kurucz stellar spectrum
      self.marcs    = None  # MARCS stellar spectrum
      self.phoenix  = None  # PHOENIX stellar spectrum
      self.starwn   = None  # Input stellar wavenumber array
      self.starflux = None  # Input stellar flux spectrum in  FINDME units

  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write('Physical properties information:')
      fw.write('\nStellar effective temperature (tstar, K): {:.1f}', self.tstar)
      fw.write('Stellar radius (rstar, Rsun): {:.3f}',
          none_div(self.rstar, pc.rsun))
      fw.write('Stellar mass (mstar, Msun):   {:.3f}',
          none_div(self.mstar, pc.msun))
      fw.write('Stellar surface gravity (gstar, cm s-2): {:.1f}', self.gstar)

      fw.write('\nPlanetary radius (rplanet, Rjup): {:.3f}',
          none_div(self.rplanet, pc.rjup))
      fw.write('Planetary mass (mplanet, Mjup):   {:.3f}',
          none_div(self.mplanet, pc.mjup))
      fw.write('Planetary surface gravity (gplanet, cm s-2): {:.1f}',
          self.gplanet)
      fw.write('Planetary internal temperature (tint, K):  {:.1f}', self.tint)
      fw.write('Orbital semi-major axis (smaxis, AU): {:.4f}',
          none_div(self.smaxis, pc.au))
      fw.write('Planet-to-star radius ratio (rprs):   {:.5f}',
          self.rplanet/self.rstar)
      fw.write('Planetary Hill radius (rhill, Rjup):  {:.3f}',
          none_div(self.rhill, pc.rjup))

      if self.starspec is not None:
          fw.write("\nInput stellar spectrum (starspec): '{}'", self.starspec)
      elif self.kurucz is not None:
          fw.write("Input Kurucz stellar spectrum (kurucz): '{}'", self.kurucz)
      elif self.marcs is not None:
          fw.write("Input MARCS stellar spectrum (marcs): '{}'", self.marcs)
      elif self.phoenix is not None:
          fw.write("Input PHOENIX stellar spectrum (phoenix): '{}'",
              self.phoenix)
      elif self.starflux is not None:
          fw.write("Input stellar spectrum is a blackbody at Teff = {:.1f} K.",
              self.tstar)
      fw.write('Stellar spectrum wavenumber (starwn, cm-1):\n    {}',
          self.starwn, fmt={'float': '{:10.3f}'.format})
      fw.write('Stellar flux spectrum (starflux, erg s-1 cm-2 cm):\n    {}',
          self.starflux, fmt={'float': '{: .3e}'.format})
      return fw.text


