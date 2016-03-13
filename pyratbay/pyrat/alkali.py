import sys, os
import numpy as np

from .. import tools as pt
from .. import constants as pc

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../lib')
import vprofile as vp


def absorption(pyrat):
  """
  Evaluate the total alkali absorption in the atmosphere.
  """
  # Initialize extinction coefficient variable:
  pyrat.alkali.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  for i in np.arange(pyrat.alkali.nmodels):
    alkali = pyrat.alkali.model[i]

    imol = np.where(pyrat.mol.name == alkali.mol)[0][0]
    iH2  = np.where(pyrat.mol.name == "H2"      )[0][0]
    iHe  = np.where(pyrat.mol.name == "He"      )[0][0]
    dens = np.expand_dims(pyrat.atm.d[:,imol], axis=1)
    temp     = pyrat.atm.temp
    pressure = pyrat.atm.press

    H2diam = (pyrat.mol.radius[imol] + pyrat.mol.radius[iH2])
    Hediam = (pyrat.mol.radius[imol] + pyrat.mol.radius[iHe])

    # Detuning frequency (cm-1):
    dsigma = alkali.detuning * (temp/500.0)**0.6
    # Doppler width (cm-1):
    dop = (np.sqrt(2*pc.k*temp / (pyrat.mol.mass[imol]*pc.amu)) *
           np.expand_dims(alkali.wn, axis=1)/pc.c  )
    # Lorentz width (cm-1):
    lor = alkali.lpar * (temp/2000.0)**(-0.7) * pressure/pc.atm
    lor2 = (np.sqrt(2/(np.pi*pc.k*temp*pc.amu)) * pressure/pc.c *
            ((H2diam**2.0 * np.sqrt(1.0/pyrat.mol.mass[imol] +
                                    1.0/pyrat.mol.mass[iH2])) +
             (Hediam**2.0 * np.sqrt(1.0/pyrat.mol.mass[imol] +
                                    1.0/pyrat.mol.mass[iHe]))))
    #print(lor/lor2)

    # Number of spectral samples in detuning region:
    vsize   = np.asarray(dsigma/pyrat.spec.wnstep+1, np.int)
    vindex  = np.zeros(pyrat.atm.nlayers, np.int)
    profile = np.zeros(np.sum(2*vsize+1), np.double)
    # Calculate Voigt profiles:
    vp.voigt(profile, vsize, vindex, lor, dop[0], pyrat.spec.wnstep, 50)

    # Calculate cross section:
    alkali.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave), np.double)
    for k in np.arange(len(alkali.wn)):
      ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave), np.double)
      wlo = np.where(pyrat.spec.wn > alkali.wn[k])[0][0] - vsize
      whi = np.where(pyrat.spec.wn > alkali.wn[k])[0][0] + vsize + 1
      for j in np.arange(pyrat.atm.nlayers):
        ec[j] += (profile[vindex[j]] *
                  (np.abs(pyrat.spec.wn-alkali.wn[k])/dsigma[j])**(-1.5))
        ec[j, wlo[j]:whi[j]] = profile[vindex[j]:vindex[j]+2*vsize[j]+1]
      alkali.ec += (pc.C3 * ec * alkali.gf[k] /alkali.Z * dens *
         np.exp(-pc.h*pc.c*np.abs(pyrat.spec.wn-alkali.wn[k])/(pc.k*temp[j])))
    #dsigma, dop, lor, lor2, profile, vsize, vindex

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

