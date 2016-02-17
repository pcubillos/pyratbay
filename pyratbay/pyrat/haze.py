import sys, os
import numpy as np

from .. import tools as pt
from .. import constants as pc

def extinction(pyrat):
  """
  Calculate the hazes extinction coefficient.
  """
  # Get list of haze models
  for i in np.arange(pyrat.haze.nmodels):
    # Feed params into the model:
    pass
    # Calculate the extinction coefficient (in cm2 molecule-1)
    pyrat.haze.model[i].extinction(pyrat.spec.wn)


def evaluate(pyrat):
  """
  Evaluate the total haze extinction coefficient in the atmosphere.
  """
  # Initialize haze extinction coefficient variable:
  pyrat.haze.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  for i in np.arange(pyrat.haze.nmodels):
    # Get molecule index:
    imol = np.where(pyrat.mol.name == pyrat.haze.model[i].mol)[0][0]
    # Densities in molecules cm-3:
    dens = pyrat.atm.d[:,imol]
    # Haze absorption (cm-1):
    pyrat.haze.ec += pyrat.haze.model[i].ec * np.expand_dims(dens, axis=1)



# List of available haze models:
hmodels = []

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
    self.coef  = np.array([8.14e-45, 1.28e-54, 1.61e-64])

  def extinction(self, wn):
    """
    Calculate the extinction coefficient in FINDME units.

    Parameters
    ----------
    wn: 1D float ndarray
       Wavenumber in cm-1.
    """
    self.ec = self.coef[0]*wn**4.0 + self.coef[1]*wn**6.0 + self.coef[2]*wn**8.0

# Add model to list:
hmodels.append(rayleighH2())


class rayleighLdE():
  """
  Rayleigh-scattering model from Lecavelier des Etangs et al. (2008).
  """
  def __init__(self):
    self.name  = "rayleigh_LdE"  # Model name
    self.npars = 3               # Number of model fitting parameters
    self.pars  = None            # Model fitting parameters
    self.ec    = None            # Model extinction coefficient

# Add model to list:
hmodels.append(rayleighLdE())

# Compile list of haze-model names:
hnames = []
for hmodel in hmodels:
  hnames.append(hmodel.name)
hnames = np.asarray(hnames)

