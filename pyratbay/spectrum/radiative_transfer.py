# Copyright (c) 2021-2022 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'radiative_equilibrium',
]

import sys

import numpy  as np
from scipy.ndimage import gaussian_filter1d as gaussf

from .. import constants as pc
from . import convection as ps


def radiative_equilibrium(
        pressure, temperature, nsamples,
        chem_model,
        two_stream_rt,
        wavenumber,
        spec, atm,
        radeq_temps=None,
        convection=False,
        tmin=0.0, tmax=6000.0, f_scale=None,
        mplanet=None,
        mol_mass=None,
    ):
    """
    Compute radiative-thermochemical equilibrium atmosphere.
    Currently there is no convergence criteria implemented,
    some 100--300 iterations are typically sufficient to converge
    to a stable temperature-profile solution.

    Parameters
    ----------
    nsamples: Integer
        Number of radiative-equilibrium iterations to run.
    continue_run: Bool
        If True, continue from a previous radiative-equil. run.
    no_convection: Bool
        If True, skip convective flux calculation in the radiative
        equilibrium calculation.

    Returns
    -------
    There are no returned values, but this method updates the
    temperature profile (self.atm.temp) and abundances (self.atm.vmr)
    with the values from the last radiative-equilibrium iteration.

    This method also defines pyrat.atm.radeq_temps, a 2D array
    containing all temperature-profile iterations.
    """
    nlayers = len(pressure)

    # Pre-existing iteration samples:
    if radeq_temps is None:
        radeq_temps = np.copy(np.atleast_2d(temperature))
    k = len(radeq_temps) - 1

    # Total number of samples:
    temp = radeq_temps = np.vstack((
        radeq_temps,
        np.zeros((nsamples, nlayers)),
    ))
    nsamples = len(temp)

    maxf = 1.0e08  # Maximum temperature scale factor
    # Initial delta-temperature scale factor
    if f_scale is None:
        f_scale = np.tile(1.0e5, nlayers)
    f_scale_tmp = np.copy(f_scale)

    dpress = np.ediff1d(np.log(pressure), to_begin=1.0)
    df_sign = np.zeros((nsamples, nlayers))

    for k in range(k, nsamples-1):
        sys.stdout.write(f"\rIteration {k+1:3d}/{nsamples-1}.")
        sys.stdout.flush()
        # Update atmosphere:
        temperature = temp[k]
        vmr = chem_model.thermochemical_equilibrium(temperature)
        two_stream_rt(temp=temperature, abund=vmr)

        flux_up = spec.flux_up
        flux_down = spec.flux_down
        # Bolometric net fluxes through each layer:
        Qup = np.trapz(flux_up, wavenumber, axis=1)
        Qdown = np.trapz(flux_down, wavenumber, axis=1)
        Q_net = Qup - Qdown
        dF = np.ediff1d(Q_net, to_begin=0)

        # Update scaling factor:
        f_scale_tmp = np.copy(f_scale)
        df_sign[k] = np.sign(dF)
        lo = np.amax([k-4, 0])
        wobble = np.any(df_sign[lo:k] - df_sign[k], axis=0)
        f_scale_tmp[wobble] *= 0.5
        f_scale_tmp[~wobble] *= 1.15
        f_scale_tmp = gaussf(np.clip(f_scale_tmp, 1.0, maxf), 1.5)
        # Adaptive temperature update:
        dT = (
            f_scale_tmp * np.sign(dF) * np.abs(dF)**0.1
            / (pc.sigma * temperature**3 * dpress)
        )
        temp[k+1] = temperature + dT
        temp[k+1,0] = temp[k+1,1]  # Isothermal top

        # Smooth out kinks and wiggles:
        avg_dT = np.mean(np.abs(dT))
        sigma = np.clip(avg_dT/10.0, 0.75, 2.0)
        temp[k+1,:-1] = gaussf(temp[k+1], sigma)[:-1]
        temp[k+1] = np.clip(temp[k+1], tmin, tmax)

        if not convection:
            f_scale = np.copy(f_scale_tmp)
            continue

        # Radiative flux balance with convection:
        radius = atm.radius
        number_density = atm.d
        mean_molecular_mass = atm.mm

        heat_capacity = chem_model.heat_capacity(temp[k+1])
        cp = np.sum(heat_capacity * vmr, axis=1) * pc.k/pc.amu
        gplanet = pc.G * mplanet / radius**2
        rho = np.sum(number_density * mol_mass, axis=1) * pc.amu

        conv_flux = ps.convective_flux(
            pressure, temp[k+1], cp,
            gplanet, mean_molecular_mass, rho,
        )
        if np.all(conv_flux == 0.0):
            f_scale = np.copy(f_scale_tmp)
            continue

        # Update scaling factor:
        dF = np.ediff1d(Q_net + conv_flux, to_begin=0)
        df_sign[k] = np.sign(dF)
        wobble = np.any(df_sign[lo:k] - df_sign[k], axis=0)
        f_scale[wobble] *= 0.5
        f_scale[~wobble] *= 1.15
        f_scale = gaussf(np.clip(f_scale, 1.0, maxf), 1.5)
        # Adaptive temperature update:
        dT = (
            f_scale * np.sign(dF) * np.abs(dF)**0.1
            / (pc.sigma * temperature**3 * dpress)
        )
        temp[k+1] = temperature + dT
        temp[k+1,0] = temp[k+1,1]

        # Smooth out kinks and wiggles:
        dtm = np.mean(np.abs(np.diff(temp[:k+1], n=1, axis=0)), axis=1)
        if np.size(dtm) > 0:
            sigma = np.clip(0.5*np.mean(dtm[-5:]), 0.5, 1.0)
        else:
            sigma = 1.0
        temp[k+1,:-1] = gaussf(temp[k+1], sigma)[:-1]
        temp[k+1] = np.clip(temp[k+1], tmin, tmax)

    return radeq_temps, f_scale
