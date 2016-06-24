# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

"""
Pyrat plotting utilities.
"""

__all__ = ["spectrum", "cf"]

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
from scipy.ndimage.filters import gaussian_filter1d as gaussf

from .. import constants as pc


def spectrum(wlength=None, spectrum=None, data=None, uncert=None,
    bandflux=None, bandtrans=None, bandidx=None, bandwl=None,
    starflux=None, rprs=None, path=None,
    pyrat=None, gaussbin=2, filename=None):
  """
  Plot a transmission or emission model spectrum with (optional) data
  points with error bars and band-integrated model.

  Parameters
  ----------
  wlength: 1D flaot ndarray
     The wavelength of the model in microns.
  spectrum: 1D float ndarray
     Planetary spectrum evaluated at wlength.
  data: 1D float ndarray
     Observing data points at each bandwl.
  uncert: 1D float ndarray
     Uncertainty of the data points.
  bandflux: 1D float ndarray
     Band-integrated model spectrum at each bandwl.
  bandtrans: List of 1D float ndarrays
     Transmission curve for each band.
  bandidx: List of 1D float ndarrays.
     The indices in wlength for each bandtrans.
  bandwl: 1D float ndarray
     The mean wavelengths for each band.
  starflux: 1D float ndarray
     Stellar spectrum evaluated at wlength.
  rprs: Float
     Planet-to-star radius ratio.
  path: String
     Observing-geometry path: transit or eclipse.
  pyrat: Pyrat instance
  gaussbin: Integer
     Standard deviation for Gaussian-kernel smoothing (in number of samples).
  filename: String
     Filename of the output figure.
  """
  # Unpack variables from Pyrat object:
  if pyrat is not None:
    wlength   = 1.0/(pyrat.spec.wn*pc.um)
    spectrum  = pyrat.spec.spectrum
    starflux  = pyrat.spec.starflux
    data      = pyrat.obs.data
    uncert    = pyrat.obs.uncert
    bandflux  = pyrat.obs.bandflux
    bandtrans = pyrat.obs.bandtrans
    bandidx   = pyrat.obs.bandidx
    if pyrat.obs.bandwn is not None:
      bandwl    = 1/(pyrat.obs.bandwn*pc.um)
    rprs      = pyrat.phy.rprs
    path      = pyrat.od.path

  if bandtrans is None:
    nfilters = 0
  else:
    nfilters = len(bandtrans)

  # Plotting setup:
  fs  = 14
  ms  =  7
  lw  = 1.5
  mew = 1.0

  plt.figure(-20, (8.5, 5))
  plt.clf()
  ax = plt.subplot(111)

  # Setup according to geometry:
  if   path == "eclipse":
    if starflux is not None:
      fscale = 1e3
      gmodel = gaussf(spectrum/starflux * rprs**2.0, gaussbin)
      plt.ylabel(r"$F_{\rm p}/F_{\rm s}\ (10^{-3})$", fontsize=fs)
    else:
      fscale = 1.0
      gmodel = gaussf(spectrum, gaussbin)
      plt.ylabel(r"$F_{\rm p}\ ({\rm erg\, s^{-1}cm^{-2}cm})$", fontsize=fs)
  elif path == "transit":
    fscale = 1.0
    gmodel = gaussf(spectrum, gaussbin)
    plt.ylabel(r"${\rm Modulation}\ \ (R_p/R_s)^2$", fontsize=fs)

  # Plot model:
  plt.plot(wlength, gmodel*fscale, lw=lw, label="Model", color="orange")
  # Plot band-integrated model:
  if bandwl is not None:
    plt.plot(bandwl, bandflux*fscale, "o", ms=ms, color="orange", mew=mew)
  # Plot data:
  if data is not None:
    plt.errorbar(bandwl, data*fscale, uncert*fscale, fmt="ob", label="Data",
                 ms=ms, elinewidth=lw, capthick=lw, zorder=3)

  # Transmission filters:
  ylim = ax.get_ybound()
  bandh = 0.06*(ylim[1] - ylim[0])
  for i in np.arange(nfilters):
    bandtr = bandh * bandtrans[i]/np.amax(bandtrans[i])
    plt.plot(wlength[bandidx[i]], ylim[0]+bandtr, "k")

  #logxtics = [0.7, 1.0, 2.0, 3.0, 4.0, 5.0]
  #if logxtics:
  #  ax.set_xscale('log')
  #  ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
  #  #ax.set_xticks(np.arange(round(min(wlength)),max(wlength),1))
  #  ax.set_xticks(logxtics)

  plt.xlabel(r"${\rm Wavelength\ \ (um)}$", fontsize=fs)
  plt.xlim(np.amin(wlength), np.amax(wlength))
  leg = plt.legend(loc="best", numpoints=1)
  if filename is not None:
    plt.savefig(filename)




def cf(bandcf, bandwl, pressure, filename=None, filters=None):
  """
  Plot the band-integrated contribution functions.

  Parameters
  ----------
  bandcf: 2D float ndarray
     Band-integrated contribution functions [nfilters, nlayers].
  bandwl: 1D float ndarray
     Mean wavelength of the bands in microns.
  pressure: 1D float ndarray
     Layer's pressure array in barye.
  filters: 1D string ndarray
     Name of the filter bands (optional).
  filename: String
     Filename of the output figure.
  """
  nfilters = len(bandwl)
  xran   = 0, np.amax(bandcf)
  press  = pressure/pc.bar
  wlsort = np.argsort(bandwl)

  fs  = 14
  lw  = 1.5
  plt.figure(-21)
  plt.clf()
  for i in np.arange(nfilters):
    idx = wlsort[i]
    ax = plt.subplot(1, nfilters, i+1)
    fname = " {:5.2f} um .".format(bandwl[idx])
    # Strip root and file extension:
    if filters is not None:
      fname = os.path.split(os.path.splitext(filters[idx])[0])[1] + " @" + fname
    c = int(10 + i / (nfilters-1.0) * 240)
    ax.semilogy(bandcf[idx], press, '-', lw=lw, color=plt.cm.rainbow(c))
    ax.set_ylim(np.amax(press), np.amin(press))
    plt.text(0.9*xran[1], np.amin(press), fname, rotation=90,
             ha="right", va="top")
    ax.set_xlim(xran)
    ax.set_xticklabels([])
    if i == 0:
      ax.set_ylabel(r'${\rm Pressure\ \ (bar)}$' , fontsize=fs)
    else:
      ax.set_yticklabels([])

  plt.subplots_adjust(0.1, 0.11, 0.95, 0.95, 0, 0)
  plt.suptitle(r'${\rm Contribution\ functions}$', fontsize=fs, y=0.09, x=0.52)

  if filename is not None:
    plt.savefig(filename)
