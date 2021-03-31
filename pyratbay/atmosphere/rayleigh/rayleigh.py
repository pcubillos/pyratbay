# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'Dalgarno',
    'Lecavelier',
    ]

import numpy as np

from ... import tools as pt


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
      self.name  = 'dalgarno_{:s}'.format(mol)  # Model name
      self.mol   = mol     # Species causing the extinction
      self.npars = 0       # Number of model fitting parameters
      self.pars  = []      # Model fitting parameters
      self.ec    = None    # Model extinction coefficient (cm2 molec-1)
      self.pnames   = []   # Fitting-parameter names
      self.texnames = []   # Fitting-parameter names

      if self.mol == 'H':
          self.coef = np.array([5.799e-45, 1.422e-54, 2.784e-64])
          self.extinction = self._extH
      elif self.mol == 'He':
          self.coef = np.array([5.484e-46, 2.440e-11, 5.940e-42, 2.900e-11])
          self.extinction = self._extHe
      elif self.mol == 'H2':
          self.coef = np.array([8.140e-45, 1.280e-54, 1.610e-64])
          self.extinction = self._extH


  def _extH(self, wn):
      """
      Calculate the opacity cross-section in cm2 molec-1 units.

      Parameters
      ----------
      wn: 1D float ndarray
         Wavenumber in cm-1.
      """
      self.ec = (self.coef[0]*wn**4.0 + self.coef[1]*wn**6.0
                                      + self.coef[2]*wn**8.0)

  def _extHe(self, wn):
      """
      Calculate the opacity cross-section in cm2 molec-1 units.

      Parameters
      ----------
      wn: 1D float ndarray
         Wavenumber in cm-1.
      """
      self.ec = self.coef[0]*wn**4 * (1 + self.coef[1]*wn**2 +
                             self.coef[2]*wn**4/(1 - self.coef[3]*wn**2))**2

  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write("Model name (name): '{}'", self.name)
      fw.write('Model species (mol): {}', self.mol)
      fw.write('Number of model parameters (npars): {}', self.npars)
      fw.write('Extinction-coefficient (ec, cm2 molec-1):\n    {}', self.ec,
          fmt={'float':'{: .3e}'.format}, edge=3)
      return fw.text


class Lecavelier():
  """
  Rayleigh-scattering model from Lecavelier des Etangs et al. (2008).
  AA, 485, 865.
  """
  def __init__(self):
      self.name = 'lecavelier'  # Model name
      self.mol = 'H2'           # Species causing the extinction
      self.pars = [ 0.0,        # Cross-section scale factor (unitless)
                   -4.0]        # Power-law exponent
      self.npars = len(self.pars)  # Number of model fitting parameters
      self.ec = None               # Model extinction coefficient
      self.pnames = ['log(f_ray)', 'alpha_ray']
      self.texnames = [r'$\log_{10}(f_{\rm ray})$', r'$\alpha_{\rm ray}$']
      self.s0 = 5.31e-27  # Cross section (cm2 molec-1) at l0
      self.l0 = 3.5e-5    # Nominal wavelength (cm)

  def extinction(self, wn):
      """
      Calculate the Rayleigh cross section in cm2 molec-1:
          cross section = f_ray * s0 * (lambda/l0)**alpha_ray,
      parameterized as params = [log10(f_ray), alpha_ray).

      Parameters
      ----------
      wn:  1D float ndarray
          Wavenumber array in cm-1.
      """
      # Rayleigh cross section opacity in cm2 molec-1:
      self.ec = 10.0**self.pars[0] * self.s0 * (wn*self.l0)**(-self.pars[1])

  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write("Model name (name): '{}'", self.name)
      fw.write('Model species (mol): {}', self.mol)
      fw.write('Number of model parameters (npars): {}', self.npars)
      fw.write('Parameter name     Value\n'
               '  (pnames)         (pars)\n')
      for pname, param in zip(self.pnames, self.pars):
          fw.write('  {:15s}  {: .3e}', pname, param)
      fw.write('Opacity cross section (ec, cm2 molec-1):\n    {}', self.ec,
          fmt={'float':'{: .3e}'.format}, edge=3)
      return fw.text
