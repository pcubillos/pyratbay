# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ['SodiumVdW', 'PotassiumVdW']

import numpy as np

from ... import constants  as pc
from ... import tools      as pt
from ... import broadening as broad


class VanderWaals(object):
  """
  Base class for Van der Waals plus statistical theory model
  Burrows et al. (2000), ApJ, 531, 438.
  """
  def __init__(self):
      self.ec    = None   # Opacity cross section (cm2 molecule-1)
      self.imol  = -1     # Index of mol in the atmosphere
      self.voigt = broad.Voigt()

  def setup(self, molecules, mass, dwave):
      if self.mol in molecules:
          self.imol = np.where(molecules==self.mol)[0][0]
          self.mass = mass[self.imol]
      self.dwave = dwave

  def absorption(self, press, temp, wn):
      """Evaluate alkali model's opacity (cm2 molecule-1)."""
      nlayers = len(press)
      nwave   = len(wn)
      self.ec = np.zeros((nlayers, nwave), np.double)

      # Species is not in atmosphere:
      if self.imol < 0:
          return np.zeros((nlayers, nwave), np.double)

      # Detuning frequency (cm-1):
      dsigma = self.detuning * (temp/500.0)**0.6
      # Doppler half width (cm-1):
      doppler = (np.sqrt(2*pc.k*temp / (self.mass*pc.amu))
                 * np.expand_dims(self.wn, axis=1)/pc.c)
      # Lorentz half width (cm-1):
      lor = self.lpar * (temp/2000.0)**(-0.7) * press/pc.atm

      # Calculate cross section:
      for wn0, gf, dwave, dop in zip(self.wn, self.gf, self.dwave, doppler):
          ec = np.zeros((nlayers, nwave), np.double)
          # Update Voigt model:
          self.voigt.x0 = wn0
          fwidth = 2*(0.5346*lor + np.sqrt(0.2166*lor**2 + dop**2))
          for j in range(nlayers):
              # Profile ranges:
              det = np.abs(wn - wn0) < dsigma[j]
              wlo = pt.ifirst(det, default_ret=-1)
              whi = pt.ilast( det, default_ret=-2) + 1
              self.voigt.hwhmL = lor[j]
              self.voigt.hwhmG = dop[j]
              wndet = wn[wlo:whi]
              # EC at the detuning boundary:
              edet = self.voigt(wn0+dsigma[j])
              # Extinction outside the detuning region (power law):
              ec[j] += edet * (np.abs(wn-wn0)/dsigma[j])**-1.5
              # Extinction in the detuning region (Voigt profile):
              profile = self.voigt(wndet)
              # Correction for undersampled line:
              if whi > wlo and fwidth[j] < 2.0*dwave:
                  i0 = np.argmin(np.abs(wn0-wndet))
                  profile[i0] = 0.0
                  profile[i0] = 1.0 - np.trapz(profile, wndet)
              ec[j, wlo:whi] = profile
          # Add up contribution (include exponential cutoff):
          self.ec += (pc.C3 * ec * gf/self.Z
              * np.exp(-pc.C2*np.abs(wn-wn0) / np.expand_dims(temp, axis=1)))
          # Note this equation neglects the exp(-Elow/T)*(1-exp(-wn0/T))
          # terms because they are approximately 1.0 at T=[100K--4000K]

  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write("Model name (name): '{}'", self.name)
      fw.write('Model species (mol): {}', self.mol)
      fw.write('Species index in atmosphere (imol): {}', self.imol)
      fw.write('Detuning parameter (detuning): {}', self.detuning)
      fw.write('Lorentz-width parameter (lpar): {}', self.lpar)
      fw.write('Partition function (Z): {}', self.Z)
      fw.write('Wavenumber  Wavelength          gf   Lower-state energy\n'
               '      cm-1          um               cm-1\n'
               '      (wn)                    (gf)   (elow)')
      for wn, gf, elow in zip(self.wn, self.gf, self.elow):
          fw.write('  {:8.2f}  {:10.6f}   {:.3e}   {:.3e}',
              wn, 1.0/(wn*pc.um), gf, elow)
      fw.write('Extinction-coefficient (ec, cm2 molecule-1):\n{}', self.ec,
          fmt={'float': '{: .3e}'.format}, edge=2)
      return fw.text


class SodiumVdW(VanderWaals):
  """Sodium Van der Waals model (Burrows et al. 2000, ApJ, 531)."""
  def __init__(self):
      super(SodiumVdW, self).__init__()
      self.name  = 'sodium_vdw'          # Model name
      self.mol   = 'Na'                  # Species causing the extinction
      # Line-transition properties (from VALD, Piskunov 1995):
      self.wn    = [16960.87, 16978.07]  # Wavenumber (cm-1)
      self.elow  = [0.0,      0.0]       # Lower-state energy (cm-1)
      self.gf    = [0.65464,  1.30918]   # gf (unitless)
      self.lpar  = 0.071     # Lorentz width parameter (Iro et al. 2005)
      self.Z     = 2.0       # Partition function (valid for temp < 4000 K)
      self.detuning = 30.0   # Detuning parameter


class PotassiumVdW(VanderWaals):
  """Potassium Van der Waals model (Burrows et al. 2000, ApJ, 531)."""
  def __init__(self):
      super(PotassiumVdW, self).__init__()
      self.name  = 'potassium_vdw'        # Model name
      self.mol   = 'K'                    # Species causing the extinction
      # Line-transition properties (from VALD, Piskunov 1995):
      self.wn    = [12988.76, 13046.486]  # Wavenumber (cm-1)
      self.elow  = [0.0,      0.0]        # Lower-state energy (cm-1)
      self.gf    = [0.701455, 1.40929]    # gf (unitless)
      self.lpar  = 0.14      # Lorentz width parameter (Iro et al. 2005)
      self.Z     = 2.0       # Partition function (temp < 4000 K, Barklem 2016)
      self.detuning = 20.0   # Detuning parameter

