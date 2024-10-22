# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'radiative_equilibrium',
]

import time

import numpy as np
from scipy.ndimage import gaussian_filter1d as gaussf

from .. import constants as pc
from .. import tools as pt
from . import convection as ps


def radiative_equilibrium(
        pressure,
        radeq_temps,
        nsamples,
        chem_model,
        two_stream_rt,
        wavenumber,
        spec,
        atm,
        convection=False,
        tmin=0.0, tmax=6000.0,
    ):
    """
    Compute radiative-thermochemical equilibrium atmosphere.
    Currently there is no convergence criteria implemented,
    some 100--300 iterations are typically sufficient to converge
    to a stable temperature-profile solution.

    Parameters
    ----------
    pressure: 1D float array
        Pressure profile in barye units.
    nsamples: Integer
        Number of radiative-equilibrium iterations to run.
    convection: Bool
        If True, include convective-flux transport in the radiative
        equilibrium calculation.

    Returns
    -------
    radeq_temps: 2D float array
        Temperature profiles at each iteration.
    """
    nlayers = len(pressure)
    n_prev = len(radeq_temps)

    # Total temperature samples:
    temp = radeq_temps = np.vstack((
        radeq_temps,
        np.zeros((nsamples, nlayers)),
    ))

    maxf = 1.0e08  # Maximum temperature scale factor
    # Scale factor for temperature jumps:
    dt_scale = atm._dt_scale
    dt_scale_tmp = np.copy(dt_scale)

    dpress = np.ediff1d(np.log(pressure), to_begin=1.0)
    dpress[0] = dpress[1]
    df_sign = np.zeros((nsamples, nlayers))

    t0 = time.time()
    for i in range(nsamples):
        k = n_prev + i - 1
        timeleft = pt.eta(time.time()-t0, i+1, nsamples, fmt='.2f')
        eta_text = (
            f'{i+1}/{nsamples} iterations, {100*(i+1)/nsamples:.2f} % done, '
            f'ETA: {timeleft}'
        )
        print(f'{eta_text:80s}', end='\r', flush=True)
        # Compute RT for input atmosphere:
        vmr = chem_model.thermochemical_equilibrium(temp[k])
        two_stream_rt(temp=temp[k], vmr=vmr)

        flux_up = spec.flux_up
        flux_down = spec.flux_down
        # Bolometric net fluxes through each layer:
        Qup = np.trapezoid(flux_up, wavenumber, axis=1)
        Qdown = np.trapezoid(flux_down, wavenumber, axis=1)
        Q_net = Qup - Qdown
        dF = np.ediff1d(Q_net, to_begin=0)

        # Update scaling factor:
        df_sign[k] = np.sign(dF)
        lo = np.amax([k-4, 0])
        wobble = np.any(df_sign[lo:k] - df_sign[k], axis=0)
        dt_scale_tmp = np.copy(dt_scale)
        dt_scale_tmp[wobble] *= 0.5
        dt_scale_tmp[~wobble] *= 1.15
        dt_scale_tmp = gaussf(np.clip(dt_scale_tmp, 1.0, maxf), 1.5)
        # Adaptive temperature update:
        dT = (
            dt_scale_tmp * np.sign(dF) * np.abs(dF)**0.1
            / (pc.sigma * temp[k]**3 * dpress)
        )
        temp[k+1] = temp[k] + dT
        temp[k+1,0] = temp[k+1,1]  # Isothermal top

        # Smooth out kinks and wiggles:
        avg_dT = np.mean(np.abs(dT))
        sigma = np.clip(avg_dT/10.0, 0.75, 2.0)
        temp[k+1,:-1] = gaussf(temp[k+1], sigma)[:-1]
        temp[k+1] = np.clip(temp[k+1], tmin, tmax)

        if not convection:
            dt_scale[:] = np.copy(dt_scale_tmp)
            continue

        # Radiative flux balance with convection:
        heat_capacity = chem_model.heat_capacity(temp[k+1])
        cp = np.sum(heat_capacity * vmr, axis=1) * pc.k/pc.amu
        gplanet = pc.G * atm.mplanet / atm.radius**2.0
        rho = np.sum(atm.d * atm.mol_mass, axis=1) * pc.amu

        conv_flux = ps.convective_flux(
            pressure, temp[k+1], cp,
            gplanet, atm.mm, rho,
        )
        if np.all(conv_flux == 0.0):
            dt_scale[:] = np.copy(dt_scale_tmp)
            continue

        # Update scaling factor:
        dF = np.ediff1d(Q_net + conv_flux, to_begin=0)
        df_sign[k] = np.sign(dF)
        wobble = np.any(df_sign[lo:k] - df_sign[k], axis=0)
        dt_scale[wobble] *= 0.5
        dt_scale[~wobble] *= 1.15
        dt_scale = gaussf(np.clip(dt_scale, 1.0, maxf), 1.5)
        # Adaptive temperature update:
        dT = (
            dt_scale * np.sign(dF) * np.abs(dF)**0.1
            / (pc.sigma * temp[k]**3 * dpress)
        )
        temp[k+1] = temp[k] + dT
        temp[k+1,0] = temp[k+1,1]

        # Smooth out kinks and wiggles:
        avg_dT = np.mean(np.abs(dT))
        sigma = np.clip(avg_dT/10.0, 0.75, 2.0)
        temp[k+1,:-1] = gaussf(temp[k+1], sigma)[:-1]
        temp[k+1] = np.clip(temp[k+1], tmin, tmax)

    return radeq_temps
