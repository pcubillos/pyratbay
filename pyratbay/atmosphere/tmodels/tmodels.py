# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Isothermal',
    'Guillot',
    'TCEA',
    'Madhu',
    'get_model',
]

import functools
from collections.abc import Iterable

import numpy as np
from scipy.ndimage import gaussian_filter1d

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


class Isothermal():
    """Isothermal temperature profile model."""
    def __init__(self, pressure):
        """
        Parameters
        ----------
        pressure: 1D float iterable
            Pressure array (bar) where to evaluate the temperature profile.
        """
        self.name = 'isothermal'
        self.pnames = ['T_iso']
        self.texnames = [r'$T\ ({\rm K})$']
        self.npars = len(self.pnames)
        self.pressure = pressure
        self.temperature = np.zeros(len(pressure), np.double)

    @check_params
    def __call__(self, params):
        """
        Parameters
        ----------
        params: scalar or iterable
            Temperature of the isothermal profile (Kelvin degree).

        Returns
        -------
        temperature: 1D float ndarray
            Temperature profile in K.

        Examples
        --------
        >>> import pyratbay.atmosphere as pa

        >>> nlayers = 51
        >>> pressure = pa.pressure('1e-7 bar', '100 bar', nlayers)
        >>> iso = pa.tmodels.Isothermal(pressure)
        >>> print(iso(1500.0))
        [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.
         1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.
         1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.
         1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.
         1500. 1500. 1500.]
        >>> print(iso([1500.0]))
        [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.
         1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.
         1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.
         1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.
         1500. 1500. 1500.]
        """
        if isinstance(params, Iterable):
            self.temperature[:] = params[0]
        else:
            self.temperature[:] = params
        return np.copy(self.temperature)

    def __str__(self):
        with np.printoptions(formatter={'float':'{:.3e}'.format}):
            str_pressure = str(self.pressure)
        return (
            f'Model name: {self.name}\n'
            f'Number of parameters (npars): {self.npars}\n'
            f'Parameter names (pnames): {self.pnames}\n'
            f'Parameter Latex names (texnames): {self.texnames}\n'
            f'Pressure array (pressure, bar):\n {str_pressure}\n'
            f'Last evaluated profile (temperature, K):\n {self.temperature}\n'
        )


class Guillot():
    """
    Guillot (2010) temperature profile based on the three-channel
    Eddington approximation, as described Line et al. (2013)
    """
    def __init__(self, pressure, gravity=None):
        """
        Parameters
        ----------
        pressure: 1D float ndarray
            Atmospheric pressure profile (bar).
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
        be solved without knowing the temperature a priori, thus making
        this a circular problem (shrug emoji).
        """
        self.name = 'guillot'
        self.pnames = [
            "log_kappa'",
            'log_gamma1',
            'log_gamma2',
            'alpha',
            'T_irr',
            'T_int',
        ]
        self.texnames = [
            r"$\log\ \kappa'$",
            r'$\log\ \gamma_1$',
            r'$\log\ \gamma_2$',
            r'$\alpha$',
            r'$T_{\rm irr} (K)$',
            r'$T_{\rm int} (K)$',
        ]
        self.npars = len(self.pnames)

        if gravity is None:
            gravity = np.tile(1.0, len(pressure))
        elif np.isscalar(gravity):
            gravity = np.tile(gravity, len(pressure))
        self.pressure = np.asarray(pressure, float)
        self.gravity = np.asarray(gravity, float)
        self.temperature = np.zeros_like(pressure, float)

    @check_params
    def __call__(self, params):
        """
        Parameters
        ----------
        params: 1D iterable
            Guillot model parameters:
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

        >>> nlayers = 101
        >>> pressure = pa.pressure('1e-7 bar', '100 bar', nlayers)
        >>> guillot = pa.tmodels.Guillot(pressure)

        >>> log_k, log_g1 = -4.8, -0.6
        >>> log_g2, alpha = 0.0, 0.0
        >>> t_irr, t_int = 1200.0, 100.0
        >>> tp_inv = guillot([log_k, log_g1, log_g2, alpha, t_irr, t_int])
        >>> log_g1 = 0.2
        >>> tp_non_inv = guillot([log_k, log_g1, log_g2, alpha, t_irr, t_int])

        >>> # Plot the profiles:
        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.semilogy(tp_inv, pressure, color='darkorange')
        >>> plt.semilogy(tp_non_inv, pressure, color='red')
        >>> plt.ylim(1e2, 1e-7)
        """
        # Ensure Numpy array:
        if isinstance(params, (list, tuple)):
            params = np.array(params, np.double)
        self.temperature[:] = _pt.guillot(
            params, self.pressure*pc.bar, self.gravity,
        )
        return np.copy(self.temperature)

    def __str__(self):
        with np.printoptions(formatter={'float':'{:.3e}'.format}):
            str_pressure = str(self.pressure)
        return (
            f'Model name: {self.name}\n'
            f'Number of parameters (npars): {self.npars}\n'
            f'Parameter names (pnames): {self.pnames}\n'
            f'Parameter Latex names (texnames): {self.texnames}\n'
            f'Pressure array (pressure, bar):\n {str_pressure}\n'
            f'Last evaluated profile (temperature, K):\n {self.temperature}\n'
        )



# For backwards compatibility:
TCEA = Guillot


class Madhu():
    """Temperature profile model by Madhusudhan & Seager (2009)"""
    def __init__(self, pressure):
        """
        Parameters
        ----------
        pressure: 1D float ndarray
            Pressure array in bar.
        """
        self.name = 'madhu'
        self.pnames = ['log_p1', 'log_p2', 'log_p3', 'a1', 'a2', 'T0']
        self.texnames = [
            r'$\log\ p_1$',
            r'$\log\ p_2$',
            r'$\log\ p_3$',
            r'$a_1$',
            r'$a_2$',
            r'$T_0$',
        ]
        self.npars = len(self.pnames)

        self.pressure = pressure
        self.logp = np.log10(pressure)
        self.temperature = np.zeros_like(pressure)
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

        >>> nlayers = 101
        >>> pressure = pa.pressure('1e-7 bar', '100 bar', nlayers)

        >>> madhu = pa.tmodels.Madhu(pressure)
        >>> # Thermally-inverted profile (p2 > p1):
        >>> #       [logp1, logp2, logp3, a1,  a2,  T0]
        >>> params = -3.5, 0.0, 0.5, 3.0, 0.5, 1500.0
        >>> inv = madhu(params)
        >>> # Non thermally-inverted profile (p1 > p2):
        >>> params = -3.5, -4.0, 0.5, 3.0, 0.5, 1100.0
        >>> non_inv = madhu(params)

        >>> # Plot the profiles:
        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.semilogy(inv, pressure, color='darkgreen')
        >>> plt.semilogy(non_inv, pressure, color='limegreen')
        >>> plt.ylim(1e2, 1e-7)
        """
        logp1, logp2, logp3, a1, a2, T0 = params
        if logp1 > logp3:
            self.temperature[:] = 0.0
            return np.copy(self.temperature)

        # Calculate temperatures at layer boundaries:
        T1 = T0 + ((logp1 - self.logp0) / (a1*self.loge))**2
        T2 = T1 - ((logp1 - logp2) / (a2*self.loge))**2
        T3 = T2 + ((logp3 - logp2) / (a2*self.loge))**2

        layer1 =  self.logp <  logp1
        layer2 = (self.logp >= logp1) & (self.logp < logp3)
        layer3 =  self.logp >= logp3

        self.temperature[layer1] = \
            T0 + ((self.logp[layer1]-self.logp0) / (a1*self.loge))**2
        self.temperature[layer2] = \
            T2 + ((self.logp[layer2]-logp2) / (a2*self.loge))**2
        self.temperature[layer3] = np.tile(T3, np.sum(layer3))

        # Smooth PT profile:
        self.temperature = gaussian_filter1d(
            self.temperature, sigma=self.fsmooth, mode='nearest',
        )
        return np.copy(self.temperature)

    def __str__(self):
        with np.printoptions(formatter={'float':'{:.3e}'.format}):
            str_pressure = str(self.pressure)
        return (
            f'Model name: {self.name}\n'
            f'Number of parameters (npars): {self.npars}\n'
            f'Parameter names (pnames): {self.pnames}\n'
            f'Parameter Latex names (texnames): {self.texnames}\n'
            f'Pressure array (pressure, bar):\n {str_pressure}\n'
            f'Last evaluated profile (temperature, K):\n {self.temperature}\n'
        )



def get_model(name, *args, **kwargs):
    """Get a temperature-profile model by its name."""
    if name == 'isothermal':
        return Isothermal(*args, kwargs['pressure'])
    if name in ['tcea', 'guillot']:
        return Guillot(*args, kwargs['pressure'])
    if name == 'madhu':
        return Madhu(*args, kwargs['pressure'])

