# Copyright (c) 2021-2025 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'transmission',
    'plane_parallel_rt',
    'radiative_equilibrium',
]

import time

import numpy as np
from scipy.ndimage import gaussian_filter1d as gaussf
import scipy.interpolate as si

from .. import constants as pc
from .. import tools as pt
from . import convection as ps
from ..lib import _trapz as t
from .blackbody import blackbody_wn


def transmission(
        depth, radius, rstar, ideep=None, atm_itop=0,
        deck_rsurf=None, deck_itop=None,
    ):
    """
    Compute a transmission spectrum for transit geometry

    Parameters
    ----------
    depth: 2D float array
        Optical depth at each layer and wavelength channel [nlayers,nwave].
        Atmospheric layers are sorted from top to bottom (i.e.,
        depth[0] is the top-most layer)
    radius: 1D float array
        Radius profile of atmospheric layers (cm)
    rstar: Float
        Stellar radius (cm).
    ideep: 1D integer array
        Index of the 'bottom' of the atmosphere at each wavelength,
        from which to start the integration (e.g., at optically thick regime)
    atm_itop: Integer
        Index of the top of the atmosphere (to integrate only up to the layer)
    cloud_rsurf: Float
        If not None, the radius (cm) of an opaque cloud deck, which
        becomes the bottom integration boundary.
    cloud_itop: Integer
        If not None, index of the atmosphere layer right below the
        opaque cloud deck.

    Returns
    -------
    spectrum: 1D float array
        The transmission spectrum.
    """
    nlayers = ideep - atm_itop + 1

    # Get Delta radius (and integration variables):
    h = np.ediff1d(radius[atm_itop:])
    integ = np.exp(-depth[atm_itop:]) * np.expand_dims(radius[atm_itop:],1)

    # Replace (by interpolating) last layer with cloud top:
    if deck_rsurf is not None and deck_itop > atm_itop:
        h[deck_itop-atm_itop-1] = deck_rsurf - radius[deck_itop-1]
        f_interp = si.interp1d(radius[atm_itop:], integ, axis=0)
        integ[deck_itop-atm_itop] = f_interp(deck_rsurf)

    # Number of layers for integration at each wavelength:
    spectrum = t.trapz2D(integ, h, nlayers-1)
    spectrum = (radius[atm_itop]**2 + 2*spectrum) / rstar**2

    return spectrum


def plane_parallel_rt(
        depth, blackbody, wn, quadrature_mu, quadrature_weights=None,
        ideep=None, atm_itop=0,
        cloud_tsurf=None, cloud_itop=None,
    ):
    """
    Compute the intensity (and optionally flux) spectra under
    plane-parallel geometry for a range of slant angles

    Parameters
    ----------
    depth: 2D float array
        Optical depth at each layer and wavelength channel [nlayers,nwave].
        Atmospheric layers are sorted from top to bottom (i.e.,
        depth[0] is the top-most layer)
    blackbody: 2D float array
        Plank fuction evaluated at each layer and wavelength channel
        (same shape as depth).
    wn: 1D float array
        Wavenumber array (cm-1).
    quadrature_mu: 1D float array
        Cosine of the slant-path angles repect to the normal.
    quadrature_weights: 1D float array
        Weights for each quadrature_mu.  If given, compute and return the
        Gaussian-quadrature integral of the intensity over mu (i.e., flux)
    ideep: 1D integer array
        Index of the 'bottom' of the atmosphere at each wavelength,
        from which to start the integration.
    atm_itop: Integer
        Index of the top of the atmosphere (to integrate only up to the layer)
    cloud_tsurf: Float
        If not None, the temperature at cloud_itop, an opaque cloud
        deck that becomes the bottom integration boundary.
    cloud_itop: Integer
        If not None, index of the atmosphere layer right below an
        opaque cloud deck.

    Returns
    -------
    intensity: 2D float array
        The intensity spectra at each angle mu (erg s-1 cm-2 cm sr-1).
    flux: 1D float array [optional]
        The flux spectrum of the integrated intensities (erg s-1 cm-2 cm).
        (typically, a day-side integrated via the Gaussian-quadrature)
    """
    if ideep is None:
        nlayers, nwave = depth.shape
        ideep = np.tile(nlayers-1, nwave)

    if cloud_tsurf is not None:
        blackbody[cloud_itop] = blackbody_wn(wn, cloud_tsurf)
        ideep = np.clip(ideep, 0, cloud_itop)

    # Plane-parallel intensity integration
    intensity = t.intensity(
        depth, ideep, blackbody, quadrature_mu, atm_itop,
    )

    # Flux integral using Gaussian quadrature
    if quadrature_weights is not None:
        flux = np.sum(intensity * quadrature_weights, axis=0)
        return intensity, flux
    return intensity


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
