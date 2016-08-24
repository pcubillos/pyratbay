# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys, os
import numpy as np

from .. import tools as pt
from .. import constants as pc


def absorption(pyrat):
  """
  Evaluate the total haze absorption in the atmosphere.
  """
  # Initialize haze extinction coefficient variable:
  pyrat.haze.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  for i in np.arange(pyrat.haze.nmodels):
    # Calculate the extinction coefficient (in cm2 molecule-1):
    pyrat.haze.model[i].extinction(pyrat.spec.wn, pyrat.atm.press)

    # Get molecule index:
    imol = np.where(pyrat.mol.name == pyrat.haze.model[i].mol)[0][0]
    # Densities in molecules cm-3:
    dens = pyrat.atm.d[:,imol]
    # Haze absorption (cm-1):
    pyrat.haze.ec += pyrat.haze.model[i].ec * np.expand_dims(dens, axis=1)


class rayleighH2():
  """
  Rayleigh-scattering model from Dalgarno & Williams (1962).
  """
  def __init__(self):
    self.name  = "rayleigh_H2"  # Model name
    self.npars = 0              # Number of model fitting parameters
    self.pars  = None           # Model fitting parameters
    self.ec    = None           # Model extinction coefficient (cm2 molec-1)
    self.mol   = "H2"           # Species causing the extinction
    self.parname = []           # Fitting-parameter names
    self.coef  = np.array([8.14e-45, 1.28e-54, 1.61e-64])

  def extinction(self, wn, pressure):
    """
    Calculate the opacity cross-section in cm2 molec-1 units.

    Parameters
    ----------
    wn: 1D float ndarray
       Wavenumber in cm-1.
    """
    self.ec = self.coef[0]*wn**4.0 + self.coef[1]*wn**6.0 + self.coef[2]*wn**8.0


class rayleighLdE():
  """
  Rayleigh-scattering model from Lecavelier des Etangs et al. (2008).
  AA, 485, 865.
  """
  def __init__(self):
    self.name  = "rayleigh_LdE"   # Model name
    self.pars  = [ 1.0,           # Cross-section scale factor (unitless)
                  -4.0]           # Power-law exponent
    self.npars = len(self.pars)   # Number of model fitting parameters
    self.ec    = None             # Model extinction coefficient
    self.mol   = "H2"             # Species causing the extinction
    self.parname = [r"$f_{\rm Rayleigh}$", # Fitting-parameter names
                    r"$\alpha$"]
    self.s0    = 5.31e-27         # Cross section (cm-2 molec-1) at l0
    self.l0    = 3.5e-5           # Nominal wavelength (cm)

  def extinction(self, wn, pressure):
    """
    Calculate the H2 Rayleigh cross section in cm2 molec-1:
       cross section = pars[0] * s0 * (lambda/l0)**(pars[1])
    With lambda the wavelength = 1/wavenumber.

    Parameters
    ----------
    wn:  1D float ndarray
       Wavenumber array in cm-1.
    """
    # Rayleigh opacity cross section in cm2 molec-1 (aka. extinction coef.):
    self.ec = (self.pars[0] * self.s0) * (wn * self.l0)**(-self.pars[1])


class grey():
  """
  Constant cross-section cloud model.
  """
  def __init__(self):
    self.name  = "grey"           # Model name
    self.pars  = [ 0.0,           # log10 of cross-section scale factor, top,
                   -4, 2]         #  and bottom pressure (bar) boundaries
    self.npars = len(self.pars)   # Number of model fitting parameters
    self.ec    = None             # Model extinction coefficient
    self.mol   = "H2"             # Species causing the extinction
    self.parname = [r"$\log_{10}(f_{\rm gray})$",  # Fitting-parameter names
                    r"$\log_{10}(p_{\rm top})\ ({\rm bar})$",
                    r"$\log_{10}(p_{\rm bot})\ ({\rm bar})$"]
    self.s0    = 5.31e-27         # Default coss-section (cm-2 molec-1)

  def extinction(self, wn, pressure):
    """
    Calculate a uniform grey-cloud cross section in cm2 molec-1:
       cross section = s0 * 10**pars[0],
    between layers with pressure 10**pars[1] -- 10**pars[2] bar
    (top and bottom layers, respectively).
    s0 is the H2 Rayleigh cross section at 0.35 um.

    Parameters
    ----------
    wn:  1D float ndarray
       Wavenumber array in cm-1.
    """
    nlayers = len(pressure)
    nwave   = len(wn)
    # Get indices for cloud layer boundaries:
    itop    = np.where(pressure >= 10**self.pars[1]*pc.bar)[0][ 0]
    ibottom = np.where(pressure <= 10**self.pars[2]*pc.bar)[0][-1]
    # Rayleigh opacity cross section in cm2 molec-1 (aka. extinction coef.):
    self.ec = np.zeros((nlayers, nwave))
    self.ec[itop:ibottom,:] = 10**self.pars[0] * self.s0


# List of available haze models:
hmodels = [rayleighH2(),
           rayleighLdE(),
           grey()]

# Compile list of haze-model names:
hnames = []
for hmodel in hmodels:
  hnames.append(hmodel.name)
hnames = np.asarray(hnames)

