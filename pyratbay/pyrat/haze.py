# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
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
    # Little hack (clean up later):
    if pyrat.haze.model[i].name == "deck":
      pyrat.haze.model[i].extinction(pyrat.spec.wn, pyrat.atm.press,
                                     pyrat.atm.radius)
      pyrat.haze.ec += pyrat.haze.model[i].ec
      continue

    # Calculate the extinction coefficient (in cm2 molecule-1):
    pyrat.haze.model[i].extinction(pyrat.spec.wn, pyrat.atm.press)

    # Get molecule index:
    imol = np.where(pyrat.mol.name == pyrat.haze.model[i].mol)[0][0]
    # Densities in molecules cm-3:
    dens = pyrat.atm.d[:,imol]
    # Haze absorption (cm-1):
    pyrat.haze.ec += pyrat.haze.model[i].ec * np.expand_dims(dens, axis=1)


def get_ec(pyrat, layer):
  """
  Extract per-model extinction coefficient at requested layer.
  """
  ec = np.zeros((pyrat.haze.nmodels, pyrat.spec.nwave))
  label = []
  for i in np.arange(pyrat.haze.nmodels):
    if pyrat.haze.model[i].name == "deck":
      pyrat.haze.model[i].extinction(pyrat.spec.wn, pyrat.atm.press,
                                     pyrat.atm.radius)
      ec[i] = pyrat.haze.model[i].ec[layer]
    else:
      imol = np.where(pyrat.mol.name == pyrat.haze.model[i].mol)[0][0]
      ec[i] = pyrat.haze.model[i].ec[layer] * pyrat.atm.d[layer,imol]
    label.append(pyrat.haze.model[i].name)
  return ec, label


class CCSgray():
  """
  Constant cross-section gray cloud model.
  """
  def __init__(self):
    self.name  = "ccsgray"        # Model name
    self.pars  = [ 0.0,          # log10 of cross-section scale factor, top,
                   -4, 2]        #  and bottom pressure (bar) boundaries
    self.npars = len(self.pars)  # Number of model fitting parameters
    self.ec    = None            # Model extinction coefficient
    self.mol   = "H2"            # Species causing the extinction
    self.parname = [r"$\log_{10}(f_{\rm gray})$",  # Fitting-parameter names
                    r"$\log_{10}(p_{\rm top})\ ({\rm bar})$",
                    r"$\log_{10}(p_{\rm bot})\ ({\rm bar})$"]
    self.s0    = 5.31e-27         # Default coss-section (cm-2 molec-1)

  def extinction(self, wn, pressure):
    """
    Calculate a uniform gray-cloud cross section in cm2 molec-1:
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
    # Gray opacity cross section in cm2 molec-1 (aka. extinction coef.):
    self.ec = np.zeros((nlayers, nwave))
    self.ec[itop:ibottom,:] = 10**self.pars[0] * self.s0


class Deck():
  """
  Instantly opaque gray cloud deck at given pressure.
  """
  def __init__(self):
    self.name  = "deck"         # Model name
    self.pars  = [-1.0]          # log10(Pressure[bar]) of cloud top
    self.npars = len(self.pars)  # Number of model fitting parameters
    self.ec    = None            # Model extinction coefficient
    self.mol   = None            # Species causing the extinction
    self.parname = [r"$\log_{10}(p_{\rm top})$"]  # Fitting-parameter names

  def extinction(self, wn, pressure, radius):
    """
    Calculate gray-cloud absorption (in cm-1) that's optically thin
    above ptop, and becomes nearly-instantly opaque at ptop, with
       ptop (bar) = 10**pars[0].

    Parameters
    ----------
    wn: 1D float ndarray
       Wavenumber array (in cm-1).
    pressure: 1D float ndarray
       Atmospheric pressure profile (in barye).
    radius: 1D float ndarray
       Atmospheric radius profile (in cm).
    """
    nlayers = len(pressure)
    nwave   = len(wn)
    ptop = 10**self.pars[0]*pc.bar
    # Index of layer directly below cloud top:
    if ptop >= pressure[-1]:  # Atmosphere boundary cases
      itop = nlayers-1
    elif ptop < pressure[0]:
      itop = 1
    else:
      itop = np.where(pressure>ptop)[0][0]

    # Set the extinction at itop such that vertical tau[itop] ~ 2/3:
    alpha = (2.0/3)/ptop
    ec = np.zeros(nlayers)
    ec[itop] = 2*alpha*(pressure[itop])/(radius[itop-1]-radius[itop])
    # Below itop, scale the extinction proportional to the pressure:
    ec[itop:] = ec[itop] * pressure[itop:]/pressure[itop]

    # Gray absorption in cm-1:
    self.ec = np.zeros((nlayers, nwave))
    self.ec[itop:,:] = np.expand_dims(ec[itop:], axis=1)


# List of available haze models:
hmodels = [CCSgray(), Deck()]

# Compile list of haze-model names:
hnames = []
for hmodel in hmodels:
  hnames.append(hmodel.name)
hnames = np.asarray(hnames)

