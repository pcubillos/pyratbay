# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["lorentz", "gauss", "voigt",
           "min_widths", "max_widths"]

import numpy as np
import scipy.special as ss

import pyratbay.constants as pc


class lorentz():
  """
  1D Lorentz profile model.

  Parameters
  ----------
  x0: Float
     Profile center location.
  hwhm: Float
     Profile's half-width at half maximum.
  scale: Float
     Scale of the profile (scale=1 returns a profile with integral=1.0).
  """
  def __init__(self, x0=0.0, hwhm=1.0, scale=1.0):
    self.x0    = x0
    self.hwhm  = hwhm
    self.scale = scale


  def __call__(self, x):
    return self.eval(x)


  def eval(self, x):
    """
    Compute Lorentz profile over the specified coordinates range.

    Parameters
    ----------
    x: 1D float ndarray
       Input coordinates where to evaluate the profile.

    Returns
    -------
    l: 1D float ndarray
       The line profile at the x locations.
    """
    return self.scale * self.hwhm/np.pi / (self.hwhm**2 + (x-self.x0)**2)


class gauss():
  """
  1D Gaussian profile model.

  Parameters
  ----------
  x0: Float
     Profile center location.
  hwhm: Float
     Profile's half-width at half maximum.
  scale: Float
     Scale of the profile (scale=1 returns a profile with integral=1.0).
  """
  def __init__(self, x0=0.0, hwhm=1.0, scale=1.0):
    self.x0    = x0
    self.hwhm  = hwhm
    self.scale = scale
    self._c1 = 1.0/np.sqrt(2*np.log(2))
    self._c2 = 1.0/np.sqrt(2*np.pi)


  def __call__(self, x):
    return self.eval(x)


  def eval(self, x):
    """
    Compute Gaussian profile over the specified coordinates range.

    Parameters
    ----------
    x: 1D float ndarray
       Input coordinates where to evaluate the profile.

    Returns
    -------
    g: 1D float ndarray
       The line profile at the x locations.
    """
    sigma = self.hwhm * self._c1
    return self.scale * np.exp(-0.5*((x-self.x0)/sigma)**2) * self._c2 / sigma


class voigt():
  """
  1D Voigt profile model.

  Parameters
  ----------
  x0: Float
     Line center location.
  hwhmL: Float
     Half-width at half maximum of the Lorentz distribution.
  hwhmG: Float
     Half-width at half maximum of the Gaussian distribution.
  scale: Float
     Scale of the profile (scale=1 returns a profile with integral=1.0).

  Example
  -------
  >>> import broadening as b

  >>> Nl = 5
  >>> Nw = 10.0

  >>> hG = 1.0
  >>> HL = np.logspace(-2, 2, Nl)
  >>> l = b.lorentz(x0=0.0)
  >>> d = b.gauss  (x0=0.0, hwhm=hG)
  >>> v = b.voigt  (x0=0.0, hwhmG=hG)
  >>>
  >>> plt.figure(11, (6,6))
  >>> plt.clf()
  >>> plt.subplots_adjust(0.15, 0.1, 0.95, 0.95, wspace=0, hspace=0)
  >>> for i in np.arange(Nl):
  >>>   hL = HL[i]
  >>>   ax = plt.subplot(Nl, 1, 1+i)
  >>>   v.hwhmL = hL
  >>>   l.hwhm  = hL
  >>>   width = 0.5346*hL + np.sqrt(0.2166*hL**2+hG**2)
  >>>   x = np.arange(-Nw*width, Nw*width, width/1000.0)
  >>>   plt.plot(x/width, l(x), lw=1.5, color="b",         label="Lorentz")
  >>>   plt.plot(x/width, d(x), lw=1.5, color="limegreen", label="Doppler")
  >>>   plt.plot(x/width, v(x), lw=1.5, color="orange",    label="Voigt")
  >>>   plt.ylim(np.amin([l(x), v(x)]), 3*np.amax([l(x), v(x), d(x)]))
  >>>   ax.set_yscale("log")
  >>>   plt.text(0.025, 0.75, r"$\rm HW_L/HW_G={:4g}$".format(hL/hG),
  >>>            transform=ax.transAxes)
  >>>   plt.xlim(-Nw, Nw)
  >>>   plt.xlabel(r"$\rm X/HW_V$", fontsize=12)
  >>>   plt.ylabel(r"$\rm Profile$")
  >>>   if i != Nl-1:
  >>>     ax.set_xticklabels([""])
  >>>   if i == 0:
  >>>     plt.legend(loc="upper right", fontsize=11)
  """
  def __init__(self, x0=0.0, hwhmL=1.0, hwhmG=1.0, scale=1.0):
    # Profile parameters:
    self.x0    = x0
    self.hwhmL = hwhmL
    self.hwhmG = hwhmG
    self.scale = scale
    # Constants:
    self._A = np.array([-1.2150, -1.3509, -1.2150, -1.3509])
    self._B = np.array([ 1.2359,  0.3786, -1.2359, -0.3786])
    self._C = np.array([-0.3085,  0.5906, -0.3085,  0.5906])
    self._D = np.array([ 0.0210, -1.1858, -0.0210,  1.1858])
    self._sqrtln2 = np.sqrt(np.log(2.0))
    self._sqrtpi  = np.sqrt(np.pi)


  def __call__(self, x):
    return self.eval(x)


  def eval(self, x):
    """
    Compute Voigt profile over the specified coordinates range.

    Parameters
    ----------
    x: 1D float ndarray
       Input coordinates where to evaluate the profile.

    Returns
    -------
    v: 1D float ndarray
       The line profile at the x locations.
    """
    if self.hwhmL/self.hwhmG < 0.1:
      # sigma * sqrt(2):
      sigmaroot2 = self.hwhmG / (self._sqrtln2 * np.sqrt(2))
      z = (x + 1j * self.hwhmL - self.x0) / sigmaroot2
      return self.scale * ss.wofz(z).real / (sigmaroot2 * self._sqrtpi)

    # This is faster than the previous script (but it fails for HWl/HWg > 1.0):
    X = (x-self.x0) * self._sqrtln2 / self.hwhmG
    Y = self.hwhmL * self._sqrtln2 / self.hwhmG

    V = 0.0
    for i in np.arange(4):
      V += (self._C[i]*(Y-self._A[i]) + self._D[i]*(X-self._B[i])) /\
           ((Y-self._A[i])**2 + (X-self._B[i])**2)
    V /= np.pi * self.hwhmL
    return self.scale * self.hwhmL/self.hwhmG * self._sqrtpi*self._sqrtln2 * V


def min_widths(min_temp, min_wn, max_mass, DLratio=0.1):
  """
  Estimate the minimum Doppler and Lorentz half-widths at half maximum
  (cm-1) for a given atmosphere.

  Parameters
  ----------
  min_temp: Float
      Minimum atmospheric tmperature (Kelvin degrees).
  min_wn: Float
      Minimum spectral wavenumber (cm-1).
  max_mass: Float
      Maximum mass of molecule/isotope (amu).
  DLratio: Float
      Doppler--Lorentz width ratio.

  Returns
  -------
  dmin: Float
      Minimum Doppler HWHM (cm-1).
  lmin: Float
      Minimum Lorentz HWHM (cm-1).

  Examples
  --------
  >>> wn = np.amin(pyrat.spec.wn)
  >>> temperature = np.amin(pyrat.atm.temp)
  >>> mols = np.unique(pyrat.iso.imol)
  >>> mols = mols[np.where(mols>=0)]
  >>> molmass = np.amax(pyrat.mol.mass[mols])
  >>> # TBD
  """
  # Minimum Doppler and Lorenz widths (cm-1):
  dmin = np.sqrt(np.log(2)*2.0*pc.k*min_temp/(max_mass*pc.amu)) * min_wn / pc.c
  lmin = dmin * DLratio
  return dmin, lmin


def max_widths(min_temp, max_temp, max_wn, min_mass, max_rad, max_press):
  """
  Estimate the maximum Doppler and Lorentz half-widths at half maximum
  (cm-1) for a given atmosphere.

  Parameters
  ----------
  min_temp: Float
      Minimum atmospheric tmperature (Kelvin degrees).
  max_temp: Float
      Maximum atmospheric tmperature (Kelvin degrees).
  max_wn: Float
      Maximum spectral wavenumber (cm-1).
  min_mass: Float
      Minimum mass of molecule/isotope (amu).
  max_press: Float
      Maximum collision radius of molecule/isotope (Angstrom).
  max_press: Float
      Maximum atmospheric pressure (barye).

  Returns
  -------
  dmax: Float
      Maximum Doppler HWHM (cm-1).
  lmax: Float
      Maximum Lorentz HWHM (cm-1).

  Examples
  --------
  >>> pmax = np.amax(pyrat.atm.press)
  >>> wn = np.amax(pyrat.spec.wn)
  >>> temperature = np.amax(pyrat.atm.temp)
  >>> mols = np.unique(pyrat.iso.imol) # Molecules with transitions
  >>> mols = mols[np.where(mols>=0)]   # Remove -1's
  >>> mmin = np.amin(pyrat.mol.mass[mols])
  >>> rmax = np.amax(pyrat.mol.radius[mols])
  >>> # TBD
  """
  # TBD: Extract values from files instead
  H2_diam = 2.89    # Angstrom
  H2_mass = 2.01588 # amu

  # Get max collision diameter:
  max_diam = H2_diam*pc.A/2.0 + max_rad

  # Doppler widths (cm-1):
  dmax = np.sqrt(np.log(2)*2.0*pc.k*max_temp/(min_mass*pc.amu)) * max_wn / pc.c

  # Approximate Sum_a (n_a * d_a**2 ...) ~ n_H2 *d_H2 ...
  # (assuming H2-dominated atmosphere)
  lmax = (np.sqrt(2.0/(np.pi * pc.k * min_temp * pc.amu)) * max_press / pc.c *
                        max_diam**2.0 * np.sqrt(1.0/min_mass + 1.0/H2_mass))
  return dmax, lmax
