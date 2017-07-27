# Copyright (c) 2016-2017 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import numpy as np

__all__ = ["cf", "transmittance", "bandcf"]

def cf(optdepth, pressure, B):
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
     The contribution function at each layer and wavenumber [nlayers, nwave].
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
  # Fill-in values beyond ideep:
  for i in np.arange(len(ideep)):
    transmit[ideep[i]:,i] = transmit[ideep[i],i]

  return transmit


def bandcf(cf, bandtrans, wn, bandidx):
  """
  Compute band-averaged contribution functions or transmittances.

  Parameters
  ----------
  cf: 2D float ndarray
     The contribution function or transmittance [nlayers, nwave]
  bandtrans: List of 1D ndarrays
     List of band transmission curves.
  wn: 1D float ndarray
     The wavenumber sampling (in cm-1).
  bandidx: List of 1D ndarrays
     List of wavenumber-index arrays for each band transmission curve.

  Returns
  -------
  bandcf: 2D float ndarray
     The band-integrated contribution functions.
  """
  nfilters = len(bandtrans)
  nlayers  = np.shape(cf)[0]

  # Allocate arrays for filter cf:
  bandcf = np.zeros((nfilters, nlayers))

  # Number of filters
  for i in np.arange(nfilters):
    # Weighted CF (by filter response function):
    wcf = cf[:,bandidx[i]] * bandtrans[i]
    # Integrated CF across bandpass at each layer:
    bandcf[i] = np.trapz(wcf, wn[bandidx[i]], axis=1)

    # Normalize to 1:
    #filt_cf[i] /= np.sum(filt_cf[i])

  return bandcf
