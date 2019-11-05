# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = [
    'Isothermal',
    'TCEA',
    'Madhu',
    ]

import sys
import numpy as np
from scipy.ndimage import gaussian_filter1d
from collections import Iterable

from ... import tools     as pt
from ... import constants as pc

sys.path.append(pc.ROOT + 'pyratbay/lib/')
import pt as PT


class Isothermal(object):
    """Isothermal temperature profile model."""
    def __init__(self, nlayers):
        """
        Parameters
        ----------
        nlayers: Integer
            Number of layers in temperature profile.
        """
        self.temp = np.zeros(nlayers, np.double)

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
    model from Line et al. (2013).
    """
    def __init__(self, pressure, rstar, tstar, tint, gplanet, smaxis,
                 runits='cm'):
        """
        Parameters
        ----------
        pressure: 1D float ndarray
            Atmospheric pressure profile in barye units.
        rstar: String or float
            Stellar radius. If string, may specify units (default in cm units).
        tstar: String or float
            Stellar temperature in Kelvin degrees.
        tint: String or float
            Planetary internal temperature in Kelvin degrees.
        gplanet: String or float
            Planetary surface gravity in cm s-2.
        smaxis: String or float
            Orbital semi-major axis.
            If string, may specify units (default in cm units).
        runits: String
            Default units for rstar and smaxis.
        """
        self.pressure = pressure
        self.rstar   = pt.get_param('rstar',   rstar,   runits,   gt=0.0)
        self.tstar   = pt.get_param('tstar',   tstar,   'kelvin', gt=0.0)
        self.tint    = pt.get_param('tint',    tint,    'kelvin', ge=0.0)
        self.gplanet = pt.get_param('gplanet', gplanet, 'none',   gt=0.0)
        self.smaxis  = pt.get_param('smaxis',  smaxis,  runits,   gt=0.0)
        self.temp = np.zeros_like(pressure)

    def __call__(self, params):
        """
        Parameters
        ----------
        params: 1D iterable
            TCEA model parameters:
            log10(kappa):  Planck thermal IR opacity in units cm2 g-1
            log10(gamma1): Visible-to-thermal stream Planck mean opacity ratio
            log10(gamma2): Visible-to-thermal stream Planck mean opacity ratio
            alpha: Visible-stream partition (0.0--1.0)
            beta:  'catch-all' for albedo, emissivity, and day--night
                   redistribution (on the order of unity)

        Returns
        -------
        temp: 1D float ndarray
            Temperature profile in K.

        Examples
        --------
        >>> import pyratbay.atmosphere as pa
        >>> pressure = pa.pressure(1e-8, 1e2, 10, 'bar')
        >>> rstar   = '0.756 rsun'
        >>> tstar   = 5040.0       # K
        >>> tint    =  100.0       # K
        >>> gplanet = 2200.0       # cm s-2
        >>> smaxis  = '0.031 au'
        >>> tcea = pa.tmodels.TCEA(pressure, rstar, tstar, tint, gplanet, smaxis)
        >>> params = [-1.5, -0.8, -0.8, 0.5, 1.0]
        >>> print(tcea(params))
        [1047.04157312 1047.04189805 1047.04531644 1047.08118784 1047.45648563
         1051.34469989 1088.69956369 1311.86379107 1640.12857767 1660.02396061
         1665.30121021]
        """
        # Ensure Numpy array:
        if isinstance(params, (list, tuple)):
            params = np.array(params, np.double)
        self.temp[:] = PT.TCEA(params, self.pressure, self.rstar,
            self.tstar, self.tint, self.smaxis, self.gplanet)
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
        self.logp = np.log10(pressure)
        self.temp = np.zeros_like(pressure)
        self.logp0 = np.amin(self.logp)
        # Standard deviation of smoothing kernel (~0.3 dex in pressure):
        self.fsmooth = 0.33/np.ediff1d(self.logp)[0]
        self.loge = np.log10(np.e)

    def __call__(self, params):
        """
        Parameters
        ----------
        params: 1D float ndarray
            Temperature model parameters: log10(p1), log10(p2), log10(p3),
            a1, a2, and T0; as defined in MS (2009), with pressure values
            in barye.

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
        >>> params = 2.39,  5.23,  7.45, 0.85, 0.67, 870.0
        >>> inv = madhu(params)
        >>> # Non thermally-inverted profile (p1 > p2):
        >>> params = 5.23,  2.39,  7.45, 0.85, 0.67, 870.0
        >>> non_inv = madhu(params)
        >>> # Plot the results:
        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.semilogy(inv,     pressure/pc.bar, color='orange')
        >>> plt.semilogy(non_inv, pressure/pc.bar, color='b')
        >>> plt.ylim(1e3, 1e-6)
        """
        logp1, logp2, logp3, a1, a2, T0 = params
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
