# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    # Classes:
    'Lorentz',
    'Gauss',
    'Voigt',
    # Functions:
    'doppler_hwhm',
    'lorentz_hwhm',
    'min_widths',
    'max_widths',
    ]

import numpy as np
import scipy.special as ss

from ... import constants as pc


class Lorentz(object):
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

  Examples
  --------
  >>> import numpy as np
  >>> import matplotlib.pyplot as plt
  >>> import pyratbay.opacity.broadening as b
  >>> lor = b.Lorentz(x0=0.0, hwhm=2.5, scale=1.0)
  >>> # Half-width at half maximum is ~2.5:
  >>> x = np.linspace(-10.0, 10.0, 100001)
  >>> print(0.5 * np.ptp(x[lor(x)>0.5*np.amax(lor(x))]))
  2.4998
  >>> # Integral is ~ 1.0:
  >>> x = np.linspace(-5000.0, 5000.0, 100001)
  >>> print(np.trapz(lor(x), x))
  0.999681690140321
  >>> # Take a look at a Lorenzt profile:
  >>> x = linspace(-10, 10, 101)
  >>> plt.plot(x, lor(x))
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


class Gauss(object):
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

  Examples
  --------
  >>> import numpy as np
  >>> import matplotlib.pyplot as plt
  >>> import pyratbay.opacity.broadening as b
  >>> gauss = b.Gauss(x0=0.0, hwhm=2.5, scale=1.0)
  >>> # Half-width at half maximum is ~2.5:
  >>> x = np.linspace(-10.0, 10.0, 100001)
  >>> print(0.5 * np.ptp(x[gauss(x)>0.5*np.amax(gauss(x))]))
  2.4998
  >>> # Integral is ~ 1.0:
  >>> x = np.linspace(-5000.0, 5000.0, 100001)
  >>> print(np.trapz(gauss(x), x))
  1.0
  >>> # Take a look at a Lorenzt profile:
  >>> x = linspace(-10, 10, 101)
  >>> plt.plot(x, gauss(x))
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


class Voigt(object):
  r"""
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
  >>> import numpy as np
  >>> import matplotlib.pyplot as plt
  >>> import pyratbay.opacity.broadening as b
  >>> Nl = 5
  >>> Nw = 10.0
  >>> hG = 1.0
  >>> HL = np.logspace(-2, 2, Nl)
  >>> l = b.Lorentz(x0=0.0)
  >>> d = b.Gauss  (x0=0.0, hwhm=hG)
  >>> v = b.Voigt  (x0=0.0, hwhmG=hG)

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
  >>>   plt.plot(x/width, l(x), lw=2.0, color="b",         label="Lorentz")
  >>>   plt.plot(x/width, d(x), lw=2.0, color="limegreen", label="Doppler")
  >>>   plt.plot(x/width, v(x), lw=2.0, color="orange",    label="Voigt",
  >>>            dashes=(8,2))
  >>>   plt.ylim(np.amin([l(x), v(x)]), 3*np.amax([l(x), v(x), d(x)]))
  >>>   ax.set_yscale("log")
  >>>   plt.text(0.025, 0.75, r"$\rm HW_L/HW_G={:4g}$".format(hL/hG),
  >>>            transform=ax.transAxes)
  >>>   plt.xlim(-Nw, Nw)
  >>>   plt.xlabel(r"$\rm x/HW_V$", fontsize=12)
  >>>   plt.ylabel(r"$\rm Profile$")
  >>>   if i != Nl-1:
  >>>       ax.set_xticklabels([""])
  >>>   if i == 0:
  >>>       plt.legend(loc="upper right", fontsize=11)
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

      # This is faster than the previous script (but fails for HWl/HWg > 1.0):
      X = (x-self.x0) * self._sqrtln2 / self.hwhmG
      Y = self.hwhmL * self._sqrtln2 / self.hwhmG

      V = 0.0
      for i in np.arange(4):
          V += (self._C[i]*(Y-self._A[i]) + self._D[i]*(X-self._B[i])) \
               / ((Y-self._A[i])**2 + (X-self._B[i])**2)
      V /= np.pi * self.hwhmL
      return self.scale * self.hwhmL/self.hwhmG * self._sqrtpi*self._sqrtln2 * V


def doppler_hwhm(temperature, mass, wn):
    r"""
    Get Doppler half-width at half maximum broadening.

    Parameters
    ----------
    temperature: Float scalar or ndarray
        Atmospheric temperature (Kelvin degree).
    mass: Float scalar or ndarray
        Mass of the species (AMU).
    wn: Float scalar or ndarray
        Wavenumber (cm-1).

    Returns
    -------
    dop_hwhm: Float scalar or ndarray
        The Doppler half-width at half maximum broadening (cm-1).

    Note
    ----
    All inputs must have compatible data shapes to be broadcastable.

    Examples
    --------
    >>> import pyratbay.opacity.broadening as b
    >>> # Doppler HWHM at 1000K and 1 micron, for H2O and CO2:
    >>> temperature = 1000.0
    >>> wn = 10000.0
    >>> mass = np.array([18.0, 44.0])
    >>> dop_hw = b.doppler_hwhm(temperature, mass, wn)
    >>> print(f'Doppler broadening:\n H2O        CO2\n{dop_hw}')
    Doppler broadening:
     H2O        CO2
    [0.02669241 0.01707253]
    """
    dop_hwhm = wn / pc.c  \
               * np.sqrt(2.0*np.log(2) * pc.k*temperature/(mass*pc.amu))
    return dop_hwhm


def lorentz_hwhm(temperature, pressure, masses, radii, vol_mix_ratio, imol):
    r"""
    Get Lorentz half-width at half maximum broadening.

    Parameters
    ----------
    temperature: Float scalar or ndarray
        Atmospheric temperature (Kelvin degree).
    pressure: Float scalar or ndarray
        Atmospheric pressure (barye).
    masses: 1D float ndarray
        Masses of atmospheric species (AMU).
    radii: 1D float ndarray
        Collision radius of atmospheric species (cm).
    vol_mix_ratio: 1D float ndarray
        Volume mixing ratio of atmospheric species.
    imol: Integer
        Index of species to calculate the HWHM (in masses/radii arrays).

    Returns
    -------
    lor_hwhm: Float scalar or ndarray
        The Lorentz half-width at half maximum broadening (cm-1).

    Note
    ----
    The temperature, pressure, and imol inputs must have compatible
    shapes to be broadcastable.

    Examples
    --------
    >>> import pyratbay.opacity.broadening as b
    >>> import pyratbay.constants as pc
    >>> # Lorenz HWHM at 1000K and 1 bar, for H2O and CO2:
    >>> temperature = 1000.0
    >>> pressure = 1.0 * pc.bar
    >>> #                  H2O   CO2   H2    He
    >>> masses = np.array([18.0, 44.0, 2.0,  4.0])
    >>> radii  = np.array([1.6,  1.9,  1.45, 1.4]) * pc.A
    >>> vmr    = np.array([1e-4, 1e-4, 0.85, 0.15])
    >>> imol = np.array([0, 1])
    >>> lor_hw = b.lorentz_hwhm(temperature, pressure, masses, radii, vmr, imol)
    >>> print(f'Lorentz broadening:\n H2O        CO2\n{lor_hw}')
    Lorentz broadening:
     H2O        CO2
    [0.03691111 0.04308068]
    """
    lor_hwhm = pressure / pc.c \
        * np.sqrt(2.0/(np.pi * pc.k * temperature * pc.amu)) \
        * sum(vmr * (radius+radii[imol])**2 * np.sqrt(1/mass + 1/masses[imol])
              for radius,mass,vmr in zip(radii, masses, vol_mix_ratio))
    return lor_hwhm


def min_widths(min_temp, max_temp, min_wn, max_mass, min_rad, min_press):
  """
  Estimate the minimum Doppler and Lorentz half-widths at half maximum
  (cm-1) for an H2-dominated atmosphere.

  Parameters
  ----------
  min_temp: Float
      Minimum atmospheric tmperature (Kelvin degrees).
  max_temp: Float
      Maximum atmospheric tmperature (Kelvin degrees).
  min_wn: Float
      Minimum spectral wavenumber (cm-1).
  max_mass: Float
      Maximum mass of molecule/isotope (amu).
  min_rad: Float
      Minimum collisional radius (cm).
  min_press: Float
      Minimum atmospheric pressure (barye).

  Returns
  -------
  dmin: Float
      Minimum Doppler HWHM (cm-1).
  lmin: Float
      Minimum Lorentz HWHM (cm-1).

  Examples
  --------
  >>> import pyratbay.opacity.broadening as b
  >>> import pyratbay.constants as pc
  >>> min_temp =  100.0
  >>> max_temp = 3000.0
  >>> min_wn   = 1.0/(10.0*pc.um)
  >>> max_mass = 18.015    # H2O molecule
  >>> min_rad  = 1.6*pc.A  # H2O molecule
  >>> min_press = 1e-5 * pc.bar
  >>> dmin, lmin = b.min_widths(min_temp, max_temp, min_wn, max_mass,
  >>>     min_rad, min_press)
  >>> print('Minimum Doppler half width: {:.2e} cm-1\n'
  >>>       'Minimum Lorentz half width: {:.2e} cm-1'.format(dmin,lmin))
  """
  # Minimum Doppler and Lorenz widths (cm-1):
  dmin = np.sqrt(np.log(2)*2.0*pc.k*min_temp/(max_mass*pc.amu)) * min_wn / pc.c
  # TBD: Extract values from atmosphere instead
  H2_radius = 1.445e-8  # cm
  H2_mass   = 2.01588   # amu

  # Get max collision diameter:
  min_diam = H2_radius + min_rad

  # Sum_a (n_a*d_a**2 ...) ~ n_H2*d_H2 ... (assuming H2-dominated atmosphere)
  lmin = (np.sqrt(2.0/(np.pi * pc.k * max_temp * pc.amu)) * min_press / pc.c *
                        min_diam**2.0 * np.sqrt(1.0/max_mass + 1.0/H2_mass))

  return dmin, lmin


def max_widths(min_temp, max_temp, max_wn, min_mass, max_rad, max_press):
  """
  Estimate the maximum Doppler and Lorentz half-widths at half maximum
  (cm-1) for an H2-dominated atmosphere.

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
  max_rad: Float
      Maximum collisional radius (cm).
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
  >>> import pyratbay.opacity.broadening as b
  >>> import pyratbay.constants as pc
  >>> min_temp =  100.0
  >>> max_temp = 3000.0
  >>> max_wn   = 1.0/(1.0*pc.um)
  >>> min_mass = 18.015    # H2O molecule
  >>> max_rad  = 1.6*pc.A  # H2O molecule
  >>> max_press = 100.0*pc.bar
  >>> dmax, lmax = b.max_widths(min_temp, max_temp, max_wn, min_mass,
  >>>     max_rad, max_press)
  >>> print('Maximum Doppler half width: {:.2e} cm-1\n'
  >>>       'Maximum Lorentz half width: {:.2e} cm-1'.format(dmax,lmax))
  """
  # TBD: Extract values from files instead
  H2_radius = 1.445e-8  # cm
  H2_mass   = 2.01588   # amu

  # Get max collision diameter:
  max_diam = H2_radius + max_rad

  # Doppler widths (cm-1):
  dmax = np.sqrt(np.log(2)*2.0*pc.k*max_temp/(min_mass*pc.amu)) * max_wn / pc.c

  # Approximate Sum_a (n_a * d_a**2 ...) ~ n_H2 *d_H2 ...
  # (assuming H2-dominated atmosphere)
  lmax = (np.sqrt(2.0/(np.pi * pc.k * min_temp * pc.amu)) * max_press / pc.c *
                        max_diam**2.0 * np.sqrt(1.0/min_mass + 1.0/H2_mass))
  return dmax, lmax
