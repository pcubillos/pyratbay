# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import numpy as np

from .. import broadening as broad
from .. import constants  as pc
from .. import tools      as pt


def init(pyrat):
    """Set species index in atmosphere."""
    if pyrat.alkali.models != []:
        pyrat.log.msg("\nSetup Alkali opacity models.")
    for alkali in pyrat.alkali.models:
        alkali.set_imol(pyrat.mol.name)


def absorption(pyrat):
  """
  Evaluate the total alkali absorption in the atmosphere.
  """
  # Initialize extinction coefficient variable:
  pyrat.alkali.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  for alkali in pyrat.alkali.models:
      imol = alkali.imol

      if alkali.imol < 0:
          pyrat.log.warning("Alkali species '{:s}' is not present in the "
                            "atmospheric file.".format(alkali.mol))
          continue
      dens = np.expand_dims(pyrat.atm.d[:,imol], axis=1)
      temp = pyrat.atm.temp
      wn   = pyrat.spec.wn

      # Detuning frequency (cm-1):
      dsigma = alkali.detuning * (temp/500.0)**0.6
      # Doppler half width (cm-1):
      dop = (np.sqrt(2*pc.k*temp / (pyrat.mol.mass[imol]*pc.amu)) *
             np.expand_dims(alkali.wn, axis=1)/pc.c)
      # Lorentz half width (cm-1):
      lor = alkali.lpar * (temp/2000.0)**(-0.7) * pyrat.atm.press/pc.atm

      # Initialize Voigt model:
      voigt = broad.Voigt()
      # Calculate cross section:
      alkali.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave), np.double)
      for k in range(len(alkali.wn)):
          ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave), np.double)
          # Update Voigt model:
          voigt.x0 = alkali.wn[k]
          fwidth = 2*(0.5346*lor + np.sqrt(0.2166*lor**2 + dop[k]**2))
          # Model spectral sampling rate at alkali.wn:
          dwave = pyrat.spec.wnstep if pyrat.spec.resolution is None else \
                  alkali.wn[k]/pyrat.spec.resolution
          for j in range(pyrat.atm.nlayers):
              # Profile ranges:
              det = np.abs(wn - alkali.wn[k]) < dsigma[j]
              wlo = pt.ifirst(det, default_ret=-1)
              whi = pt.ilast( det, default_ret=-2) + 1
              voigt.hwhmL = lor[j]
              voigt.hwhmG = dop[k,j]
              wndet = wn[wlo:whi]
              # EC at the detuning boundary:
              edet = voigt(alkali.wn[k]+dsigma[j])
              # Extinction outside the detuning region (power law):
              ec[j] += edet * (np.abs(wn-alkali.wn[k])/dsigma[j])**-1.5
              # Extinction in the detuning region (Voigt profile):
              profile = voigt(wndet)
              # Correction for undersampled line:
              if whi > wlo and fwidth[j] < 2.0*dwave:
                  i0 = np.argmin(np.abs(alkali.wn[k]-wndet))
                  profile[i0] = 0.0
                  profile[i0] = 1.0-np.trapz(profile, wndet)
              ec[j, wlo:whi] = profile
          # Add up contribution (include exponential cutoff):
          alkali.ec += (pc.C3 * ec * alkali.gf[k]/alkali.Z * dens
                        * np.exp(-pc.C2*np.abs(wn-alkali.wn[k])
                                 / np.expand_dims(temp,axis=1)))
          # Note this equation neglects the exp(-Elow/T)*(1-exp(-nu0/T))
          # terms because they are approximately 1.0 at T=[100K--4000K]

      # Sum alkali extinction coefficient (cm-1):
      pyrat.alkali.ec += alkali.ec


def get_ec(pyrat, layer):
  """
  Extract per-species extinction coefficient at requested layer.
  """
  absorption(pyrat)
  ec, label = [], []
  for alkali in pyrat.alkali.models:
      if alkali.imol >= 0:
          ec.append(alkali.ec[layer])
          label.append(alkali.mol)
  return ec, label


class VanderWaals():
  """
  Van der Waals plus statistical theory model
  Burrows et al. (2000), ApJ, 531, 438.
  """
  def __init__(self):
      self.name  = None     # Model name
      self.ec    = None     # Opacity cross section (cm2 molec-1)
      self.mol   = None     # Species causing the extinction
      self.imol  = -1       # Index of mol in the atmosphere
      self.wn    = []       # Line-transition wavenumber (cm-1)
      self.elow  = []       # Line-transition lower-state energy (cm-1)
      self.gf    = []       # Line-transition gf
      self.lpar  = None     # Lorentz width parameter
      self.Z     = None     # Partition function
      self.detuning = None  # Detuning parameter

  def set_imol(self, molecules):
      if self.mol in molecules:
          self.imol = np.where(molecules==self.mol)[0][0]

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
      fw.write('Extinction-coefficient (ec, cm-1):\n{}', self.ec,
          fmt={'float': '{: .3e}'.format}, edge=2)
      return fw.text


class SodiumVdWst(VanderWaals):
  """Sodium Van der Waals model (Burrows et al. 2000, ApJ, 531)."""
  def __init__(self):
      self.name  = 'SodiumVdWst'
      self.ec    = None           # Opacity cross section (cm2 molec-1)
      self.mol   = 'Na'           # Species causing the extinction
      self.imol  = -1
      # Line properties (from VALD, Psikunov ):
      self.wn    = [16960.87, 16978.07] # Wavenumber in cm-1
      self.elow  = [0.0,      0.0]
      self.gf    = [0.65464,  1.30918]
      self.detuning = 30.0   # Detuning parameter
      self.lpar = 0.071  # Lorentz width parameter (From Iro et al. 2005)
      self.Z    = 2.0    # Partition function (valid for temp < 4000 K)


class PotassiumVdWst(VanderWaals):
  """Potassium Van der Waals model (Burrows et al. 2000, ApJ, 531)."""
  def __init__(self):
      self.name  = 'PotassiumVdWst' # Model name
      self.ec    = None             # Opacity cross section (cm2 molec-1)
      self.mol   = 'K'
      self.imol  = -1
      # Line properties (VALD):
      self.wn    = [12988.76, 13046.486]  # Wavenumber in cm-1
      self.elow  = [0.0,      0.0]
      self.gf    = [0.701455, 1.40929]
      self.detuning = 20.0  # Detuning parameter
      self.lpar = 0.14      # Lorentz width parameter (From Iro et al. 2005)
      self.Z    = 2.0       # Partition function (temp < 4000 K) # Barklem 2016


# List of available alkali models:
models = [SodiumVdWst(),
          PotassiumVdWst()]

# Compile list of alkali-model names:
mnames = []
for model in models:
    mnames.append(model.name)
mnames = np.asarray(mnames)

