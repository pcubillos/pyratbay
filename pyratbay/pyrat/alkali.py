# Copyright (c) 2016-2017 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys, os
import numpy as np

from .. import tools as pt
from .. import constants as pc
from .. import broadening as broad


def init(pyrat):
  if pyrat.alkali.nmodels > 0:

    # Species index in atmosphere:
    pyrat.alkali.imol = -np.ones(pyrat.alkali.nmodels, int)
    for i in np.arange(pyrat.alkali.nmodels):
      imol = np.where(pyrat.mol.name == pyrat.alkali.model[i].mol)[0]
      if np.size(imol) != 0:
        pyrat.alkali.imol[i] = imol[0]
    pt.msg(pyrat.verb-3, "Done.", pyrat.log)


def absorption(pyrat):
  """
  Evaluate the total alkali absorption in the atmosphere.
  """
  # Initialize extinction coefficient variable:
  pyrat.alkali.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  for i in np.arange(pyrat.alkali.nmodels):
    alkali = pyrat.alkali.model[i]

    if pyrat.alkali.imol[i] < 0:
      pt.warning(pyrat.verb-2, "Alkali species '{:s}' is not present in "
        "the atmospheric file.".format(alkali.mol), pyrat.log, pyrat.wlog)
      continue
    imol = pyrat.alkali.imol[i]
    dens = np.expand_dims(pyrat.atm.d[:,imol], axis=1)
    temp = pyrat.atm.temp

    # Detuning frequency (cm-1):
    dsigma = alkali.detuning * (temp/500.0)**0.6
    # Doppler half width (cm-1):
    dop = (np.sqrt(2*pc.k*temp / (pyrat.mol.mass[imol]*pc.amu)) *
           np.expand_dims(alkali.wn, axis=1)/pc.c  )
    # Lorentz half width (cm-1):
    lor = alkali.lpar * (temp/2000.0)**(-0.7) * pyrat.atm.press/pc.atm

    idet = np.asarray(dsigma/pyrat.spec.wnstep, np.int)

    # Initialize Voigt model:
    voigt = broad.voigt()
    # Calculate cross section:
    alkali.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave), np.double)
    for k in np.arange(len(alkali.wn)):
      ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave), np.double)
      # Profile ranges:
      offset = int((alkali.wn[k] - pyrat.spec.wn[0])/pyrat.spec.wnstep) + 1
      wlo = np.clip(offset-idet,   0, pyrat.spec.nwave)
      whi = np.clip(offset+idet+1, 0, pyrat.spec.nwave)
      # Update Voigt model:
      voigt.x0 = alkali.wn[k]
      fwidth = 0.5346*(2*lor) + np.sqrt(0.2166*(2*lor)**2 + (2*dop[k])**2)
      for j in np.arange(pyrat.atm.nlayers):
        voigt.hwhmL = lor[j]
        voigt.hwhmG = dop[k,j]
        wndet = pyrat.spec.wn[wlo[j]:whi[j]]
        # EC at the detuning boundary:
        edet = voigt(alkali.wn[k]+dsigma[j])
        # Extinction outside the detuning region (power law):
        ec[j] += edet / (dsigma[j] / np.abs(pyrat.spec.wn-alkali.wn[k]))**-1.5
        # Extinction in the detuning region (Voigt profile):
        profile = voigt(wndet)
        # Correction for undersampled line:
        if fwidth[j] < 2.0*pyrat.spec.wnstep:
          i0 = np.argmin(np.abs(alkali.wn[k]-wndet))
          profile[i0] = 0.0
          profile[i0] = 1.0-np.trapz(profile, wndet)
        ec[j, wlo[j]:whi[j]] = profile
      # Add up contribution (include exponential cutoff):
      alkali.ec += (pc.C3 * ec * alkali.gf[k] /alkali.Z * dens *
                    np.exp(-pc.C2*np.abs(pyrat.spec.wn-alkali.wn[k])/temp[j]))

    # Sum alkali extinction coefficient (cm-1):
    pyrat.alkali.ec += alkali.ec


def get_ec(pyrat, layer):
  """
  Extract per-species extinction coefficient at requested layer.
  """
  absorption(pyrat)
  label = []
  ec = np.zeros((0, pyrat.spec.nwave))
  for i in np.arange(pyrat.alkali.nmodels):
    if pyrat.alkali.imol[i] >= 0:
      e = pyrat.alkali.model[i].ec[layer]
      ec = np.vstack((ec, e))
      label.append(pyrat.alkali.model[i].mol)
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
    # Line properties (VALD):
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

