# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import numpy as np

from .. import broadening as broad
from .. import constants  as pc
from .. import tools      as pt


def init(pyrat):
    if pyrat.alkali.models != []:
        pyrat.log.msg("\nSetup Alkali opacity models.")
    # Species index in atmosphere:
    for alkali in pyrat.alkali.models:
        if alkali.mol in pyrat.mol.name:
            alkali.imol = np.where(pyrat.mol.name==alkali.mol)[0][0]


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

      # Detuning frequency (cm-1):
      dsigma = alkali.detuning * (temp/500.0)**0.6
      # Doppler half width (cm-1):
      dop = (np.sqrt(2*pc.k*temp / (pyrat.mol.mass[imol]*pc.amu)) *
             np.expand_dims(alkali.wn, axis=1)/pc.c  )
      # Lorentz half width (cm-1):
      lor = alkali.lpar * (temp/2000.0)**(-0.7) * pyrat.atm.press/pc.atm


      # Initialize Voigt model:
      voigt = broad.Voigt()
      # Calculate cross section:
      alkali.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave), np.double)
      for k in np.arange(len(alkali.wn)):
          ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave), np.double)
          # Update Voigt model:
          voigt.x0 = alkali.wn[k]
          fwidth = 0.5346*(2*lor) + np.sqrt(0.2166*(2*lor)**2 + (2*dop[k])**2)
          # Model spectral sampling rate at alkali.wn:
          dwave = pyrat.spec.wnstep if pyrat.spec.resolution is None else \
                  alkali.wn[k]/pyrat.spec.resolution
          for j in np.arange(pyrat.atm.nlayers):
              # Profile ranges:
              det = np.abs(pyrat.spec.wn - alkali.wn[k]) < dsigma[j]
              wlo = pt.ifirst(det, default_ret=-1)
              whi = pt.ilast( det, default_ret=-2) + 1
              voigt.hwhmL = lor[j]
              voigt.hwhmG = dop[k,j]
              wndet = pyrat.spec.wn[wlo:whi]
              # EC at the detuning boundary:
              edet = voigt(alkali.wn[k]+dsigma[j])
              # Extinction outside the detuning region (power law):
              ec[j] += edet / (dsigma[j] / np.abs(pyrat.spec.wn-alkali.wn[k]))**-1.5
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
                        * np.exp(-pc.C2*np.abs(pyrat.spec.wn-alkali.wn[k])
                                 / np.expand_dims(temp,axis=1)))
          # Note this equation is neglecting the exp(-Elow/T)*(1-exp(-nu0/T))
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


class SodiumVdWst():
  """
  Sodium Van der Waals plus statistical theory opacity model from
  Burrows et al. (2000), ApJ, 531, 438.
  """
  def __init__(self):
    self.name  = "SodiumVdWst"  # Model name
    self.ec    = None           # Opacity cross section (cm2 molec-1)
    self.mol   = "Na"           # Species causing the extinction
    self.imol  = -1
    # Line properties (from VALD):
    self.wn    = [16960.87, 16978.07] # Wavenumber in cm-1
    self.elow  = [0.0,      0.0]
    self.gf    = [0.65464,  1.30918]
    self.detuning = 30.0   # Detuning parameter
    self.lpar = 0.071  # Lorentz width parameter (From Iro et al. 2005)
    self.Z    = 2.0    # Partition function (valid for temp < 4000 K)


class PotassiumVdWst():
  """
  Potassium Van der Waals plus statistical theory opacity model from
  Burrows et al. (2000), ApJ, 531, 438.
  """
  def __init__(self):
    self.name  = "PotassiumVdWst" # Model name
    self.ec    = None             # Opacity cross section (cm2 molec-1)
    self.mol   = "K"
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

