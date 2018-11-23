# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys, os
import numpy as np

from .. import tools as pt
from .. import constants as pc


def absorption(pyrat):
  """
  Evaluate the total absorption in the atmosphere.
  """
  # Initialize extinction coefficient variable:
  pyrat.rayleigh.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  for i in np.arange(pyrat.rayleigh.nmodels):
    # Calculate the extinction coefficient (in cm2 molecule-1):
    pyrat.rayleigh.model[i].extinction(pyrat.spec.wn, pyrat.atm.press)

    # Get molecule index:
    imol = np.where(pyrat.mol.name == pyrat.rayleigh.model[i].mol)[0][0]
    # Densities in molecules cm-3:
    dens = pyrat.atm.d[:,imol]
    # Absorption (cm-1):
    pyrat.rayleigh.ec += pyrat.rayleigh.model[i].ec * \
                          np.expand_dims(dens, axis=1)


def get_ec(pyrat, layer):
  """
  Extract per-model extinction coefficient at requested layer.
  """
  ec = np.zeros((pyrat.rayleigh.nmodels, pyrat.spec.nwave))
  label = []
  for i in np.arange(pyrat.rayleigh.nmodels):
    imol = np.where(pyrat.mol.name == pyrat.rayleigh.model[i].mol)[0][0]
    pyrat.rayleigh.model[i].extinction(pyrat.spec.wn, pyrat.atm.press)
    ec[i] = pyrat.rayleigh.model[i].ec * pyrat.atm.d[layer,imol]
    label.append(pyrat.rayleigh.model[i].name)
  return ec, label


class Dalgarno():
  """
  Rayleigh-scattering model from Dalgarno (1962), Kurucz (1970), and
  Dalgarno & Williams (1962).
  """
  def __init__(self, mol):
    """
    Parameters
    ----------
    mol: String
       The species, which can be H, He, or H2.
    """
    self.name  = "dalgarno_{:s}".format(mol)  # Model name
    self.npars = 0              # Number of model fitting parameters
    self.pars  = None           # Model fitting parameters
    self.ec    = None           # Model extinction coefficient (cm2 molec-1)
    self.mol   = mol            # Species causing the extinction
    self.pnames    = []         # Fitting-parameter names
    self.figpnames = []         # Fitting-parameter names

    if   self.mol == "H":
      self.coef = np.array([5.799e-45, 1.422e-54, 2.784e-64])
      self.extinction = self._extH
    elif self.mol == "He":
      self.coef = np.array([5.484e-46, 2.440e-11, 5.940e-42, 2.900e-11])
      self.extinction = self._extHe
    elif self.mol == "H2":
      self.coef = np.array([8.140e-45, 1.280e-54, 1.610e-64])
      self.extinction = self._extH

  def _extH(self, wn, pressure):
    """
    Calculate the opacity cross-section in cm2 molec-1 units.

    Parameters
    ----------
    wn: 1D float ndarray
       Wavenumber in cm-1.
    """
    self.ec = self.coef[0]*wn**4.0 + self.coef[1]*wn**6.0 + self.coef[2]*wn**8.0

  def _extHe(self, wn, pressure):
    """
    Calculate the opacity cross-section in cm2 molec-1 units.

    Parameters
    ----------
    wn: 1D float ndarray
       Wavenumber in cm-1.
    """
    self.ec = self.coef[0]*wn**4 * (1 + self.coef[1]*wn**2 +
                           self.coef[2]*wn**4/(1 - self.coef[2]*wn**2))**2


class Lecavelier():
  """
  Rayleigh-scattering model from Lecavelier des Etangs et al. (2008).
  AA, 485, 865.
  """
  def __init__(self):
    self.name  = "lecavelier"     # Model name
    self.pars  = [ 0.0,           # Cross-section scale factor (unitless)
                  -4.0]           # Power-law exponent
    self.npars = len(self.pars)   # Number of model fitting parameters
    self.ec    = None             # Model extinction coefficient
    self.mol   = "H2"             # Species causing the extinction
    self.pnames    = ["log(f_Ray)", "alpha_Ray"]
    self.figpnames = [r"$\log_{10}(f_{\rm Ray})$", r"$\alpha_{\rm Ray}$"]
    self.s0    = 5.31e-27         # Cross section (cm-2 molec-1) at l0
    self.l0    = 3.5e-5           # Nominal wavelength (cm)

  def extinction(self, wn, pressure):
    """
    Calculate the H2 Rayleigh cross section in cm2 molec-1:
       cross section = 10**pars[0] * s0 * (lambda/l0)**(pars[1])
    With lambda the wavelength = 1/wavenumber.

    Parameters
    ----------
    wn:  1D float ndarray
       Wavenumber array in cm-1.
    """
    # Rayleigh opacity cross section in cm2 molec-1 (aka. extinction coef.):
    self.ec = 10.0**self.pars[0] * self.s0 * (wn * self.l0)**(-self.pars[1])



# List of available Rayleigh models:
rmodels = [Dalgarno("H"),
           Dalgarno("He"),
           Dalgarno("H2"),
           Lecavelier()]

# Compile list of Rayleigh-model names:
rnames = []
for rmodel in rmodels:
  rnames.append(rmodel.name)
rnames = np.asarray(rnames)

