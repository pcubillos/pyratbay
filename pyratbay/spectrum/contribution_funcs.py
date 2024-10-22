# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'contribution_function',
    'transmittance',
    'band_cf',
]

import numpy as np


def contribution_function(optdepth, pressure, B):
    """
    Evaluate the contribution function equation as in Knutson et al. (2009)
    ApJ, 690, 822; Equation (2).

    Parameters
    ----------
    optdepth: 2D float ndarray
        Optical depth at each layer and wavenumber [nlayers, nwave].
    pressure: 1D float ndarray
        Atmospheric pressure array [nlayers].
    B: 2D float ndarray
        Plank emission at each layer and wavenumber [nlayers, nwave].

    Returns
    -------
    cf: 2D float ndarray
        The contribution function at each layer and wavenumber
        of shape [nlayers, nwave].
    """
    # Differential along layers:
    detau = np.diff(np.exp(-optdepth), axis=0)
    # Fix discontinuity:
    detau[np.where(detau > 0.1)] = 0.0

    # Log-pressure differential:
    dlogp = np.ediff1d(np.log(pressure))
    # Contribution functions at each layer and wavelength:
    cf = B[:-1] * detau / np.expand_dims(dlogp, axis=1)

    # Append row to preserve number of layers:
    cf = np.vstack([cf, np.zeros(np.shape(B)[1])])

    # Normalize to sum(cf) = 1 at each wavelength:
    cf = cf / np.sum(cf, axis=0)

    return cf


def transmittance(optdepth, ideep):
    """
    Compute the transmittance spectrum for the impact-parameter
    raypaths of a transmission model.

    Parameters
    ----------
    optdepth: 2D float ndarray
        Optical depth at each layer and wavenumber [nlayers, nwave].
    ideep: 1D integer ndarray
        Impact-parameter indices of deepest-calculated optical depth
        at each wavenumber.
    """
    # Transmittance:
    transmit = np.exp(-optdepth)
    # Fill-in values beyond ideep (completely opaque):
    for i in range(len(ideep)):
        transmit[ideep[i]:,i] = 0.0

    return transmit


def band_cf(cf, bands_response, wn, bands_idx):
    """
    Compute band-averaged contribution functions or transmittances.

    Parameters
    ----------
    cf: 2D float ndarray
        The contribution function or transmittance of
        shape [nlayers, nwave].
    bands_response: List of 1D ndarrays
        List of band transmission response curves.
    wn: 1D float ndarray
        The wavenumber sampling (in cm-1).
    bands_idx: List of 1D ndarrays
        List of wavenumber-indices in wn sampled by bands_response.

    Returns
    -------
    bands_cf: 2D float ndarray
        The band-integrated contribution functions of
        shape [nlayers, nbands].
    """
    nfilters = len(bands_response)
    nlayers = np.shape(cf)[0]

    # Allocate arrays for filter cf:
    bands_cf = np.zeros((nlayers, nfilters))

    # Number of filters
    for i in range(nfilters):
        response = bands_response[i]
        wn_idx = bands_idx[i]
        # Weighted CF (by filter response function):
        wcf = cf[:,wn_idx] * response
        # Integrated CF across bandpass at each layer:
        bands_cf[:,i] = np.trapezoid(wcf, wn[wn_idx], axis=1)

    return bands_cf
