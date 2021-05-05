# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'Isothermal',
    'TCEA',
    'Madhu',
    'get_model',
    ]

import functools

import numpy as np
from numpy.core.numeric import isscalar
from scipy.ndimage import gaussian_filter1d
from collections import Iterable

from ... import constants as pc
from ...lib import _pt


def check_params(func):
    """Decorator to check that the number of model parameters is correct."""
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        self = args[0]
        params = kwargs['params'] if 'params' in kwargs else args[1]
        if np.size(params) != self.npars:
            raise ValueError(
                f'Number of temperature parameters ({np.size(params)}) does '
                f'not match the required number of parameters ({self.npars}) '
                f'of the {self.name} model')
        return func(*args, **kwargs)
    return new_func


class Isothermal(object):
    """Isothermal temperature profile model."""
    def __init__(self, nlayers):
        """
        Parameters
        ----------
        nlayers: Integer
            Number of layers in temperature profile.
        """
        self.name = 'isothermal'
        self.pnames = ['T (K)']
        self.texnames = [r'$T\ ({\rm K})$']
        self.npars = len(self.pnames)
        self.temp = np.zeros(nlayers, np.double)

    @check_params
    def __call__(self, params):
        """
        Parameters
        ----------
        params: scalar or iterable
            Temperature of the isothermal profile (Kelvin degree).

        Returns
        -------
        temp: 1D float ndarray
            Temperature profile in K.

        Examples
        --------
        >>> import pyratbay.atmosphere as pa
        >>> nlayers = 8
        >>> iso = pa.tmodels.Isothermal(nlayers)
        >>> print(iso(1500.0))
        [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.]
        >>> print(iso([1500.0]))
        [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.]
        """
        if isinstance(params, Iterable):
            self.temp[:] = params[0]
        else:
            self.temp[:] = params
        return np.copy(self.temp)


class TCEA(object):
    """
    Three-channel Eddington Approximation (TCEA) temperature profile
    model from Line et al. (2013), which is based on Guillot (2010).
    """
    def __init__(self, pressure, gravity=None):
        """
        Parameters
        ----------
        pressure: 1D float ndarray
            Atmospheric pressure profile (barye).
        gravity: 1D float ndarray or scalar
            Atmospheric gravity profile (cm s-2).
            If None, assume a constant gravity of 1 cm s-2, in which
            case, one should regard the kappa parameter as
            kappa' = kappa/gravity.

        Note that the input gravity can be a scalar value used at all
        atmospheric pressures (as it has been used so far in the
        literature.  However, from a parametric point of view, this is
        redundant, as it only acts as a scaling factor for kappa.
        Ideally, one would wish to input a pressure-dependent gravity,
        but such profile would need to be derived from a hydrostatic
        equilibrium calculation, for example.  Unfortunately, HE cannot
        be solved without knowing the temperature, thus making this a
        circular problem (shrug emoji).
        """
        self.name = 'tcea'
        self.pnames = [
            "log(kappa')",
            'log(gamma1)',
            'log(gamma2)',
            'alpha',
            'T_irr (K)',
            'T_int (K)']
        self.texnames = [
            r"$\log_{10}(\kappa')$",
            r'$\log_{10}(\gamma_1)$',
            r'$\log_{10}(\gamma_2)$',
            r'$\alpha$',
            r'$T_{\rm irr} (K)$',
            r'$T_{\rm int} (K)$']
        self.npars = len(self.pnames)

        if gravity is None:
            gravity = np.tile(1.0, len(pressure))
        elif isscalar(gravity):
            gravity = np.tile(gravity, len(pressure))
        self.pressure = np.asarray(pressure, float)
        self.gravity = np.asarray(gravity, float)
        self.temp = np.zeros_like(pressure, float)

    @check_params
    def __call__(self, params):
        """
        Parameters
        ----------
        params: 1D iterable
            TCEA model parameters:
            log10(kappa'): Planck thermal IR opacity divided by gravity.
                This relates to the familiar kappa from Guillot (2010) as
                kappa' = kappa/g, unless gravity is defined on initialization,
                in which case the two kappa's are the same.
            log10(gamma1): Visible-to-thermal stream Planck mean opacity ratio
            log10(gamma2): Visible-to-thermal stream Planck mean opacity ratio
            alpha: Visible-stream partition (0.0--1.0)
            t_irr: Stellar irradiation temperature (Kelvin degrees)
                A good approximation is the planet's equilibrium temperature
                t_irr = t_star * sqrt(0.5*R_star/a) * ((1-A)/f)**0.25
            t_int: Planetary internal heat flux (in Kelvin degrees)

        Returns
        -------
        temp: 1D float ndarray
            Temperature profile in Kelvin degrees.

        Examples
        --------
        >>> import pyratbay.atmosphere as pa
        >>> import pyratbay.constants as pc
        >>> pressure = pa.pressure(1e-8, 1e2, 10, 'bar')
        >>> tcea = pa.tmodels.TCEA(pressure)
        >>> kappa, gamma1, gamma2, alpha = -4.8, -0.8, -0.8, 0.5
        >>> t_irr = pa.equilibrium_temp(5040, 0.756*pc.rsun, 0.031*pc.au)[0]
        >>> t_int = 100.0
        >>> print(tcea([kappa, gamma1, gamma2, alpha, t_irr, t_int]))
        [1047.04157666 1047.04205435 1047.04857799 1047.13739008 1048.34031821
         1064.09404968 1215.18944824 1608.78252538 1659.93776642 1665.89970977]
        """
        # Ensure Numpy array:
        if isinstance(params, (list, tuple)):
            params = np.array(params, np.double)
        self.temp[:] = _pt.tcea(params, self.pressure, self.gravity)
        return np.copy(self.temp)


class Madhu(object):
    """Temperature profile model by Madhusudhan & Seager (2009)"""
    def __init__(self, pressure):
        """
        Parameters
        ----------
        pressure: 1D float ndarray
            Pressure array in barye.
        """
        self.name = 'madhu'
        self.pnames = ['logp1', 'logp2', 'logp3', 'a1', 'a2', 'T0']
        self.texnames = [
            r'$\log_{10}(p_1)$', r'$\log_{10}(p_2)$', r'$\log_{10}(p_3)$',
            r'$a_1$', r'$a_2$', r'$T_0$']
        self.npars = len(self.pnames)

        self.logp = np.log10(pressure/pc.bar)
        self.temp = np.zeros_like(pressure)
        self.logp0 = np.amin(self.logp)
        # Standard deviation of smoothing kernel (~0.3 dex in pressure):
        self.fsmooth = 0.33/np.ediff1d(self.logp)[0]
        self.loge = np.log10(np.e)

    @check_params
    def __call__(self, params):
        """
        Parameters
        ----------
        params: 1D float ndarray
            Temperature model parameters: log10(p1), log10(p2), log10(p3),
            a1, a2, and T0; as defined in Madhusudhan & Seager (2009),
            with the pressure values given in bars.

        Returns
        -------
        temp: 1D float ndarray
            Temperature profile in K.

        Examples
        --------
        >>> import pyratbay.atmosphere as pa
        >>> import pyratbay.constants as pc
        >>> import matplotlib.pyplot as plt
        >>> # array of pressures, equally spaced in log space
        >>> pressure = pa.pressure(1e-6, 1e3, 100, 'bar')
        >>> madhu = pa.tmodels.Madhu(pressure)
        >>> # Thermally-inverted profile (p2 > p1):
        >>> #       [logp1, logp2, logp3, a1,  a2,  T0]
        >>> params = -4.0, -1.0, 1.5, 0.85, 0.67, 870.0
        >>> inv = madhu(params)
        >>> # Non thermally-inverted profile (p1 > p2):
        >>> params = -1.0, -4.0, 1.5, 0.85, 0.67, 870.0
        >>> non_inv = madhu(params)
        >>> # Plot the results:
        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.semilogy(inv,     pressure/pc.bar, color='orange')
        >>> plt.semilogy(non_inv, pressure/pc.bar, color='b')
        >>> plt.ylim(1e3, 1e-6)
        """
        logp1, logp2, logp3, a1, a2, T0 = params
        if logp1 > logp3:
            self.temp[:] = 0.0
            return np.copy(self.temp)

        # Calculate temperatures at layer boundaries:
        T1 = T0 + ((logp1 - self.logp0) / (a1*self.loge))**2
        T2 = T1 - ((logp1 - logp2) / (a2*self.loge))**2
        T3 = T2 + ((logp3 - logp2) / (a2*self.loge))**2

        layer1 =  self.logp <  logp1
        layer2 = (self.logp >= logp1) & (self.logp < logp3)
        layer3 =  self.logp >= logp3

        self.temp[layer1] = T0 + ((self.logp[layer1]-self.logp0)
                                  / (a1*self.loge))**2
        self.temp[layer2] = T2 + ((self.logp[layer2]-logp2) / (a2*self.loge))**2
        self.temp[layer3] = np.tile(T3, np.sum(layer3))

        # Smooth PT profile:
        self.temp = gaussian_filter1d(self.temp, sigma=self.fsmooth,
            mode='nearest')
        return np.copy(self.temp)


def get_model(name, *args, **kwargs):
    """Get a temperature-profile model by its name."""
    if name == 'isothermal':
        return Isothermal(*args, kwargs['nlayers'])
    if name == 'tcea':
        return TCEA(*args, kwargs['pressure'])
    if name == 'madhu':
        return Madhu(*args, kwargs['pressure'])

