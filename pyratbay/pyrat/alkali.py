# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys, os
import numpy as np

from .. import tools as pt
from .. import constants as pc

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../lib')
import vprofile as vp
import cutils   as cu


def init(pyrat):
  if pyrat.alkali.nmodels > 0:
    pressure = pyrat.atm.press

    # Minimum-maximum Doppler widths (cm-1) for Na & K @50--5000 K:
    doppler = np.logspace(np.log10(0.006), np.log10(0.11), 15)
    # Minimum-maximum Lorentz widths (cm-1) for Na & K @50--5000 K:
    lorentz = np.logspace(np.log10(0.035*np.amin(pressure/pc.atm)),
                          np.log10(1.900*np.amax(pressure/pc.atm)), 30)

    # Maximum number of samples in detuning region for Na & K @3000K:
    ndetuning = int(90.0/pyrat.spec.wnstep) + 1
    vindex  = np.zeros((len(lorentz),len(doppler)), np.int)
    vsize   = np.tile(ndetuning, (len(lorentz),len(doppler)))
    for i in np.arange(len(lorentz)):
      iskip = np.where(doppler/lorentz[i] < 0.1)[0][1:]
      vsize[i][iskip] = 0

    profile = np.zeros(np.sum(2*vsize+1), np.double)
    logtext = " "*800
    # Calculate Voigt profiles:
    pt.msg(pyrat.verb-3, "\nComputing grid of Voigt profiles for "
           "alkali species.", pyrat.log)
    vp.alkali(profile, vsize, vindex, lorentz, doppler, pyrat.spec.wnstep,
              logtext, pyrat.verb)
    pyrat.alkali.voigt   = profile
    pyrat.alkali.lorentz = lorentz
    pyrat.alkali.doppler = doppler
    pyrat.alkali.vindex  = vindex
    pyrat.alkali.vsize   = vsize[0,0]
    pt.msg(pyrat.verb-3, "Done.", pyrat.log)


def absorption(pyrat):
  """
  Evaluate the total alkali absorption in the atmosphere.
  """
  # Initialize extinction coefficient variable:
  pyrat.alkali.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  for i in np.arange(pyrat.alkali.nmodels):
    alkali = pyrat.alkali.model[i]

    imol = np.where(pyrat.mol.name == alkali.mol)[0]
    if np.size(imol) == 0:
      pt.warning(pyrat.verb-2, "Alkali species '{:s}' is not present in "
        "the atmospheric file.".format(alkali.mol), pyrat.log, pyrat.wlog)
      continue
    else:
      imol = imol[0]
    dens = np.expand_dims(pyrat.atm.d[:,imol], axis=1)
    temp     = pyrat.atm.temp
    pressure = pyrat.atm.press

    # Detuning frequency (cm-1):
    dsigma = alkali.detuning * (temp/500.0)**0.6
    # Doppler width (cm-1):
    dop = (np.sqrt(2*pc.k*temp / (pyrat.mol.mass[imol]*pc.amu)) *
           np.expand_dims(alkali.wn, axis=1)/pc.c  )
    # Lorentz width (cm-1):
    lor = alkali.lpar * (temp/2000.0)**(-0.7) * pressure/pc.atm

    # Index of closest doppler/lorentz width for each layer:
    idop = cu.arrbinsearch(dop[0], pyrat.alkali.doppler)
    ilor = cu.arrbinsearch(lor,    pyrat.alkali.lorentz)
    # Number of detuning spectral samples for each layer:
    idet = np.asarray(dsigma/pyrat.spec.wnstep, np.int)

    # Calculate cross section:
    alkali.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave), np.double)
    for k in np.arange(len(alkali.wn)):
      ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave), np.double)

      # Profile ranges:
      offset = int((alkali.wn[k] - pyrat.spec.wn[0])/pyrat.spec.wnstep)+1
      wlo = np.clip(offset-idet,   0, pyrat.spec.nwave)
      whi = np.clip(offset+idet+1, 0, pyrat.spec.nwave)
      pran = np.vstack((wlo - offset + idet, whi - offset + idet)).T
      for j in np.arange(pyrat.atm.nlayers):
        # Voigt profile inside the detuning region:
        profile = pyrat.alkali.voigt[
           pyrat.alkali.vindex[ilor[j],idop[j]]+pyrat.alkali.vsize-idet[j]:
           pyrat.alkali.vindex[ilor[j],idop[j]]+pyrat.alkali.vsize+idet[j]+1]
        # Extinction outside the detuning region (power law):
        ec[j] += (profile[0] / dsigma[j]**-1.5 *
                  np.abs(pyrat.spec.wn-alkali.wn[k])**-1.5)
        # Extinction in the detuning region (Voigt profile):
        ec[j, wlo[j]:whi[j]] = profile[pran[j,0]:pran[j,1]]
      # Add up contribution (include exponential cutoff):
      alkali.ec += (pc.C3 * ec * alkali.gf[k] /alkali.Z * dens *
                    np.exp(-pc.C2*np.abs(pyrat.spec.wn-alkali.wn[k])/temp[j]))

    # Sum alkali extinction coefficient (cm-1):
    pyrat.alkali.ec += alkali.ec


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

