# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'convective_flux',
    ]

import numpy as np

from .. import constants as pc


def convective_flux(
    pressure, temperature, cp, gravity, mu, rho, alpha=1.0, beta=1.0):
    """
    Estimate the convective flux for an atmosphere following mixing-
    length theory as described in 'Modern Astrophysics (Carrol & Ostlie).

    Parameters
    ----------
    pressure: 1D float array
        Atmospheric pressure profile (barye).
    temperature: 1D float array
        Atmospheric temperature profile (kelvin degree).
    cp: 1D float array
        Atmospheric specific heat capacity at constant pressure
        (erg K-1 mol-1).
    gravity: 1D float array
        Atmospheric gravity profile (cm s-2).
    mu: 1D float array
        Atmospheric mean molecular mass profile (g mol-1).
    rho: 1D float array
        Atmospheric mass-density profile (g cm-3).
    alpha: Float
        Mixing-length scaling parameter: alpha = l/H, with l the
        mixing length and H the pressure scale height.
    beta: Float
        Free parameter from the average kinetic energy velocity
        estimation. Should take a value between 0 and 1.

    Returns
    -------
    F_conv: 1D float array
        Estimated convective flux (erg s-1 cm-2).
        This is not zero only where the actual temperature gradient
        is larger than the adiabatic gradient.
    """
    dpress = np.ediff1d(np.log(pressure), to_begin=1.0)
    # Actual (radiative) gradient:
    grad_t = np.ediff1d(np.log(temperature), to_begin=0.0) / dpress
    # Adiabatic gradient:
    cv = cp - pc.k/pc.amu
    gamma = cp / cv
    grad_ad = 1.0 - 1.0/gamma
    delta_grad = np.clip(grad_t - grad_ad, 0, np.inf)

    H = pc.k * temperature / (mu*pc.amu * gravity)

    F_conv = (
        alpha**2 * np.sqrt(beta)
        * cp/mu * rho * temperature
        * np.sqrt(gravity*H)
        * delta_grad**1.5
    )

    return F_conv

