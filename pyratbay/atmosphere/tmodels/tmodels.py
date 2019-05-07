# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = [
    'tcea',
    'isothermal',
    'madhu_inv',
    'madhu_noinv',
    ]

import sys
import numpy as np
from scipy.ndimage import gaussian_filter1d

from ... import tools     as pt
from ... import constants as pc

sys.path.append(pc.ROOT + 'pyratbay/lib/')
import pt as PT


def tcea(tparams, pressure, rstar, tstar, tint, gplanet, smaxis,
         runits='cm'):
    """
    Compute Three-channel Eddington Approximation (TCEA) temperature
    profile model.

    tparams: 1D iterable
        TCEA model parameters:
        log10(kappa):  Planck thermal IR opacity in units cm^2/gr
        log10(gamma1): Visible-to-thermal stream Planck mean opacity ratio
        log10(gamma2): Visible-to-thermal stream Planck mean opacity ratio
        alpha: Visible-stream partition (0.0--1.0)
        beta:  'catch-all' for albedo, emissivity, and day--night
               redistribution (on the order of unity)
    pressure: 1D float ndarray
        Atmospheric pressure profile in barye units.
    rstar: String or float
        Stellar radius (default in cm). If string, may specify units.
    tstar: String or float
        Stellar temperature in Kelvin degrees.
    tint: String or float
        Planetary internal temperature in Kelvin degrees.
    gplanet: String or float
        Planetary atmospheric temperature in cm s-2.
    smaxis: String or float
        Orbital semi-major axis (default in cm). If string, may specify units.
    runits: String
        Default units for rstar and smaxis.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants  as pc
    >>> nlayers = 11
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, 'bar')
    >>> rstar   = 0.756 * pc.rsun
    >>> tstar   = 5040.0  # K
    >>> tint    = 100.0   # K
    >>> gplanet = 2200.0   # cm s-2
    >>> smaxis  = 0.031 * pc.au
    >>> tparams = [-1.5, -0.8, -0.8, 0.5, 1.0]
    >>> temp = pa.tmodels.tcea(tparams, pressure, rstar, tstar, tint,
                               gplanet, smaxis)
    >>> print(temp)
    [1047.04157312 1047.04189805 1047.04531644 1047.08118784 1047.45648563
     1051.34469989 1088.69956369 1311.86379107 1640.12857767 1660.02396061
     1665.30121021]
    """
    # Ensure Numpy array:
    if isinstance(tparams, (list, tuple)):
        tparams = np.array(tparams, np.double)
    # Parse inputs:
    rstar   = pt.get_param('rstar',   rstar,   runits,   gt=0.0)
    tstar   = pt.get_param('tstar',   tstar,   'kelvin', gt=0.0)
    tint    = pt.get_param('tint',    tint,    'kelvin', ge=0.0)
    gplanet = pt.get_param('gplanet', gplanet, 'none',   gt=0.0)
    smaxis  = pt.get_param('smaxis',  smaxis,  runits,   gt=0.0)
    # Define model and arguments:
    targs  = [pressure, rstar, tstar, tint, smaxis, gplanet]
    return PT.TCEA(tparams, *targs)


def isothermal(tparams, nlayers):
    """
    Compute isothermal temperature profile.

    Parameters
    ----------
    tparams: scalar or iterable
        Temperature of the isothermal profile (Kelvin degree).
    nlayers: Integer
        Number of layers in temperature profile.

    Returns
    -------
    temp: 1D float ndarray
        Temperature profile.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> nlayers = 8
    >>> temp = pa.tmodels.isothermal(1500.0, nlayers)
    >>> print(temp)
    [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.]
    >>> temp = pa.tmodels.isothermal(np.array([1500.0]), nlayers)
    >>> print(temp)
    [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.]
    """
    # Cast to Numpy double array:
    if isinstance(tparams, (list, tuple)):
        tparams = np.array(tparams, np.double)
    else:
        tparams = np.array([tparams], np.double)

    # Actually, this is a wrapper of a CPython routine:
    return PT.isothermal(tparams, nlayers)


def madhu_inv(params, pressure):
    """
    Calculates PT profile for inversion case based on Equation (2) from
    Madhusudhan & Seager 2009.

    Parameters
    ----------
    pressure: 1D float ndarray
        Pressure array, needs to be equally spaced in log space from bottom
        to top of the atmosphere. Must be given in bars.
    params: 1D float ndarray
        Temperature model parameters:
            a1 - float, exponential factor in Layer 1,
            a2 - float, exponential factor in Layer 2,
            p1 - floa, pressure boundary between Layer 1 and 2 (in bars).
            p2 - float, pressure in the middle of tLayer 2
            p3 - float, pressure boundary between Layers 2 and 3 (in bars).
            T3 - float, temperature in the Layer 3.

    Returns
    -------
    temp: 1D float ndarray
        Gaussian smoothed temperatures.

    Example
    -------
    >>> import pyratbay.atmosphere as pa
    >>> # array of pressures, equally spaced in log space
    >>> press = pa.pressure(1e-5, 1e2, 100, 'bar')
    >>> # Params [a1,  a2,   p1,     p2,   p3,  T3]
    >>> params = 0.51, 0.25, 4.5e-3, 0.01, 1.0, 1500.0
    >>> temp = pa.tmodels.madhu_inv(params, press/pc.bar)
    """
    # Unpack params:
    a1, a2, p1, p2, p3, T3 = params
    # Top of the atmosphere:
    p0 = np.amin(pressure)

    # Calculate temperatures at layer boundaries:
    T2 = T3 - (np.log(p3/p2) / a2)**2
    T1 = T2 + (np.log(p1/p2) / a2)**2
    T0 = T1 - (np.log(p1/p0) / a1)**2

    layer1 =  pressure <  p1
    layer2 = (pressure >= p1) & (pressure < p3)
    layer3 =  pressure >= p3

    temp = np.concatenate([
        T0 + (np.log(pressure[layer1]/p0) / a1)**2,
        T2 + (np.log(pressure[layer2]/p2) / a2)**2,
        np.repeat(T3, np.sum(layer3))])

    # Smoothed PT profile:
    return gaussian_filter1d(temp, sigma=4.0, mode='nearest')


def madhu_noinv(params, pressure):
    """
    Calculates PT profile for inversion case based on Equation (2) from
    Madhusudhan & Seager 2009.

    Parameters
    ----------
    pressure: 1D float ndarray
        Pressure array, needs to be equally spaced in log space from bottom
        to top of the atmosphere. Must be given in bars.
    params: 1D float ndarray
        Temperature model parameters:
            a1 - float, exponential factor in Layer 1,
            a2 - float, exponential factor in Layer 2,
            p1 - floa, pressure boundary between Layer 1 and 2 (in bars).
            p3 - float, pressure boundary between Layers 2 and 3 (in bars).
            T3 - float, temperature in the Layer 3.

    Returns
    -------
    T_smooth:  1D array of floats, Gaussian smoothed temperatures,
               no kinks on Layer boundaries

    Example
    -------
    >>> import pyratbay.atmosphere as pa
    >>> # array of pressures, equally spaced in log space
    >>> press = pa.pressure(1e-5, 1e2, 100, 'bar')
    >>> # Params [a1,  a2,   p1,     p3,  T3]
    >>> params = 0.51, 0.25, 4.5e-3, 1.0, 1500.0
    >>> temp = pa.tmodels.madhu_noinv(params, press/pc.bar)
    """
    # Unpack params:
    a1, a2, p1, p3, T3 = params
    # Top of the atmosphere:
    p0 = np.amin(pressure)

    # Calculate temperatures at layer boundaries:
    T1 = T3 - (np.log(p3/p1) / a2)**2.0
    T0 = T1 - (np.log(p1/p0) / a1)**2.0

    layer1 =  pressure < p1
    layer2 = (pressure >= p1) & (pressure < p3)
    layer3 =  pressure >= p3

    temp = np.concatenate([
        T0 + (np.log(pressure[layer1]/p0) / a1)**2,
        T1 + (np.log(pressure[layer2]/p1) / a2)**2,
        np.repeat(T3, np.sum(layer3))])

    # Smoothed PT profile:
    return gaussian_filter1d(temp, sigma=4.0, mode='nearest')
