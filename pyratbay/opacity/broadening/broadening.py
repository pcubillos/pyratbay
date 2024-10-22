# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

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


class Lorentz():
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
    >>> print(np.trapezoid(lor(x), x))
    0.999681690140321
    >>> # Take a look at a Lorenzt profile:
    >>> x = linspace(-10, 10, 101)
    >>> plt.plot(x, lor(x))
    """
    def __init__(self, x0=0.0, hwhm=1.0, scale=1.0):
        self.x0 = x0
        self.hwhm = hwhm
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


class Gauss():
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
    >>> print(np.trapezoid(gauss(x), x))
    1.0
    >>> # Take a look at a Lorenzt profile:
    >>> x = linspace(-10, 10, 101)
    >>> plt.plot(x, gauss(x))
    """
    def __init__(self, x0=0.0, hwhm=1.0, scale=1.0):
        self.x0 = x0
        self.hwhm = hwhm
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
        return (
            self.scale * self._c2 / sigma
            * np.exp(-0.5*((x-self.x0)/sigma)**2)
        )


class Voigt():
    r"""
    1D Voigt profile model.

    Parameters
    ----------
    x0: Float
        Line center location.
    hwhm_L: Float
        Half-width at half maximum of the Lorentz distribution.
    hwhm_G: Float
        Half-width at half maximum of the Gaussian distribution.
    scale: Float
        Scale of the profile (scale=1 returns a profile with integral=1.0).

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pyratbay.opacity.broadening as b

    >>> hwhm_G = 1.0
    >>> hwhm_L = 1.0
    >>> voigt = b.Voigt(x0=0.0, hwhm_L=hwhm_L, hwhm_G=hwhm_G)

    >>> plt.figure('A Voigt profile', (6,4))
    >>> plt.clf()
    >>> ax = plt.subplot(1, 1, 1)
    >>> x = np.linspace(-15.0, 15.0, 3001)
    >>> plt.plot(x, voigt(x), lw=2.0, color="orange")
    >>> plt.xlim(np.amin(x), np.amax(x))
    >>> plt.xlabel(r"x", fontsize=12)
    >>> plt.ylabel(r"Voigt profile", fontsize=12)

    >>> # Compare a range of Voigt, Lorentz, and Doppler profiles:
    >>> lorentz = b.Lorentz(x0=0.0, hwhm=hwhm_L)
    >>> doppler = b.Gauss(x0=0.0, hwhm=hwhm_G)
    >>> nplots = 5
    >>> HWHM_L = np.logspace(-2, 2, nplots)
    >>> nwidths = 10.0
    >>> plt.figure(11, (6,6))
    >>> plt.clf()
    >>> plt.subplots_adjust(0.15, 0.1, 0.95, 0.95, wspace=0, hspace=0)
    >>> for i,hwhm_L in enumerate(HWHM_L):
    >>>     ax = plt.subplot(nplots, 1, 1+i)
    >>>     voigt.hwhm_L = lorentz.hwhm = hwhm_L
    >>>     width = 0.5346*hwhm_L + np.sqrt(0.2166*hwhm_L**2+hwhm_G**2)
    >>>     x = np.arange(-nwidths*width, nwidths*width, width/1000.0)
    >>>     plt.plot(x/width, lorentz(x), lw=2.0, color="blue", label="Lorentz")
    >>>     plt.plot(
    >>>         x/width, doppler(x), lw=2.0, color="limegreen", label="Doppler")
    >>>     plt.plot(
    >>>         x/width, voigt(x), lw=2.0, color="orange", label="Voigt",
    >>>         dashes=(4,1))
    >>>     ymin = np.amin([lorentz(x), voigt(x)])
    >>>     ymax = np.amax([lorentz(x), voigt(x), doppler(x)])
    >>>     plt.ylim(ymin, 3*ymax)
    >>>     ax.set_yscale("log")
    >>>     plt.text(
    >>>         0.025, 0.75, rf"$\rm HWHM_L/HWHM_G={hwhm_L/hwhm_G:4g}$",
    >>>         transform=ax.transAxes)
    >>>     plt.xlim(-nwidths, nwidths)
    >>>     plt.xlabel(r"$\rm x/HWHM_V$", fontsize=12)
    >>>     plt.ylabel("Profile")
    >>>     if i != nplots-1:
    >>>         ax.set_xticklabels([])
    >>>     if i == 0:
    >>>         plt.legend(loc="upper right", fontsize=11)
    """
    def __init__(self, x0=0.0, hwhm_L=1.0, hwhm_G=1.0, scale=1.0):
        # Profile parameters:
        self.x0    = x0
        self.hwhm_L = hwhm_L
        self.hwhm_G = hwhm_G
        self.scale = scale
        # Constants:
        self._A = np.array([[-1.2150, -1.3509, -1.2150, -1.3509]]).T
        self._B = np.array([[ 1.2359,  0.3786, -1.2359, -0.3786]]).T
        self._C = np.array([[-0.3085,  0.5906, -0.3085,  0.5906]]).T
        self._D = np.array([[ 0.0210, -1.1858, -0.0210,  1.1858]]).T
        self._sqrt_pi_log2 = np.sqrt(np.pi*np.log(2.0))


    def __call__(self, x):
        return self.eval(x)


    def eval(self, x):
        """
        Evaluate the Voigt profile at the specified coordinates range.

        Parameters
        ----------
        x: 1D float ndarray
            Input coordinates where to evaluate the profile.

        Returns
        -------
        v: 1D float ndarray
            The line profile at the x locations.
        """
        if self.hwhm_L/self.hwhm_G < 0.1:
            sigma = self.hwhm_G / np.sqrt(np.log(2))
            z = (x + 1j * self.hwhm_L - self.x0) / sigma
            return self.scale * ss.wofz(z).real / (sigma * np.sqrt(np.pi))

        # This is faster than the previous algorithm but it gets inaccurate
        # and breaks down for HWHM_L / HWHM_G ~< 0.1:
        X = (x-self.x0) * np.sqrt(np.log(2)) / self.hwhm_G
        Y = self.hwhm_L * np.sqrt(np.log(2)) / self.hwhm_G

        V = np.sum(
            (self._C*(Y-self._A) + self._D*(X-self._B))
            / ((Y-self._A)**2 + (X-self._B)**2),
            axis=0,
        )
        if np.isscalar(x):
            V = np.squeeze(V)

        return V * self.scale * self._sqrt_pi_log2 / (np.pi * self.hwhm_G)


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
    dop_hwhm = (
        wn / pc.c
        * np.sqrt(2.0*np.log(2) * pc.k*temperature/(mass*pc.amu))
    )
    return dop_hwhm


def lorentz_hwhm(temperature, pressure, masses, radii, vol_mix_ratio, imol):
    r"""
    Get Lorentz half-width at half maximum broadening.

    Parameters
    ----------
    temperature: Float scalar or ndarray
        Atmospheric temperature (Kelvin degree).
    pressure: Float scalar or ndarray
        Atmospheric pressure (bar).
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
    >>> pressure = 1.0  # bar
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
    coll_factor = sum(
        vmr * (radius+radii[imol])**2 * np.sqrt(1/mass + 1/masses[imol])
        for radius,mass,vmr in zip(radii, masses, vol_mix_ratio)
    )
    lor_hwhm = (
        pressure*pc.bar / pc.c
        * np.sqrt(2.0/(np.pi * pc.k * temperature * pc.amu))
        * coll_factor
    )
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
        Minimum atmospheric pressure (bar).

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
    >>> min_press = 1e-5 # bar
    >>> dmin, lmin = b.min_widths(min_temp, max_temp, min_wn, max_mass,
    >>>     min_rad, min_press)
    >>> print('Minimum Doppler half width: {:.2e} cm-1\n'
    >>>       'Minimum Lorentz half width: {:.2e} cm-1'.format(dmin,lmin))
    Minimum Doppler half width: 8.44e-04 cm-1
    Minimum Lorentz half width: 2.21e-07 cm-1
    """
    # Doppler width (cm-1):
    dmin = (
        np.sqrt(2.0*np.log(2.0) * pc.k*min_temp/(max_mass*pc.amu))
        * min_wn / pc.c
    )

    # Lorentz width (cm-1)
    H2_radius = 1.445e-8  # cm
    H2_mass = 2.01588   # amu

    # Get max collision diameter:
    min_diam = H2_radius + min_rad

    # Sum_a (n_a*d_a**2 ...) ~ n_H2*d_H2 ... (assuming H2-dominated atmosphere)
    lmin = (
        np.sqrt(2.0/(np.pi * pc.k * max_temp * pc.amu))
        * min_press*pc.bar * min_diam**2.0 / pc.c
        * np.sqrt(1.0/max_mass + 1.0/H2_mass)
    )

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
        Maximum atmospheric pressure (bar).

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
    >>> max_press = 100.0 # bar
    >>> dmax, lmax = b.max_widths(min_temp, max_temp, max_wn, min_mass,
    >>>     max_rad, max_press)
    >>> print('Maximum Doppler half width: {:.2e} cm-1\n'
    >>>       'Maximum Lorentz half width: {:.2e} cm-1'.format(dmax,lmax))
    Maximum Doppler half width: 4.62e-02 cm-1
    Maximum Lorentz half width: 1.21e+01 cm-1
    """
    # Doppler width (cm-1):
    dmax = (
        np.sqrt(2.0*np.log(2.0) * pc.k*max_temp/(min_mass*pc.amu))
        * max_wn / pc.c
    )

    # Lorentz width (cm-1)
    H2_radius = 1.445e-8  # cm
    H2_mass = 2.01588  # amu

    # Get max collision diameter:
    max_diam = H2_radius + max_rad

    # Approximate Sum_a (n_a * d_a**2 ...) ~ n_H2 *d_H2 ...
    # (assuming H2-dominated atmosphere)
    lmax = (
        np.sqrt(2.0/(np.pi * pc.k * min_temp * pc.amu))
        * max_press*pc.bar * max_diam**2.0 / pc.c
        * np.sqrt(1.0/min_mass + 1.0/H2_mass)
    )
    return dmax, lmax
