# Copyright (c) 2016-2017 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

"""
Pyrat plotting utilities.
"""

__all__ = ["spectrum", "cf", "TCEA"]

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
from scipy.ndimage.filters import gaussian_filter1d as gaussf

from .. import constants  as pc
from .. import atmosphere as atm
from .. import wine       as w

rootdir = os.path.realpath(os.path.dirname(__file__) + "/../../")
sys.path.append(rootdir + "/pyratbay/lib/")
import pt as PT

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

  Returns
  -------
  ax: AxesSubplot instance
    The matplotlib Axes of the figure.
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
    if bandflux is None or np.all(bandflux==0):
      bandflux = w.bandintegrate(pyrat=pyrat)

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

  return ax


def cf(bandcf, bandwl, path, layers, rtop=0, filename=None, filters=None):
  """
  Plot the band-integrated contribution functions (emission) or
  transmittance (transmission).

  Parameters
  ----------
  bandcf: 2D float ndarray
     Band-integrated contribution functions [nfilters, nlayers].
  bandwl: 1D float ndarray
     Mean wavelength of the bands in microns.
  path: String
     Observing geometry (transit or eclipse).
  layers: 1D float ndarray
     Layer's pressure (eclipse) or impact parameter (transit) array (CGS units).
  rtop: Integer
     Index of topmost valid layer.
  filters: 1D string ndarray
     Name of the filter bands (optional).
  filename: String
     Filename of the output figure.
  """
  nfilters = len(bandwl)
  wlsort   = np.argsort(bandwl)

  if   path == "eclipse":
    press = layers/pc.bar
    xran = -0.03*np.amax(bandcf), 1.03*np.amax(bandcf)
    yran = np.amax(press), np.amin(press)
    xlabel = r'${\rm Contribution\ functions}$'
    ylabel = r'${\rm Pressure\ \ (bar)}$'
  elif path == "transit":
    rad  = layers[rtop:]/pc.km
    xran = -0.03, 1.03
    yran = np.amin(rad), np.amax(rad)
    xlabel = r'${\rm Band-averaged\ transmittance}$'
    ylabel = r'${\rm Impact\ parameter\ \ (km)}$'
  else:
    print("Invalid geometry.  Select from: 'eclipse' or 'transit'.")
    return

  fs  = 14
  lw  = 2.0
  colors = np.asarray(np.linspace(10, 240, nfilters), np.int)

  plt.figure(-21)
  plt.clf()
  for i in np.arange(nfilters):
    idx = wlsort[i]
    ax = plt.subplot(1, nfilters, i+1)
    fname = " {:5.2f} um .".format(bandwl[idx])
    # Strip root and file extension:
    if filters is not None:
      fname = os.path.split(os.path.splitext(filters[idx])[0])[1] + " @" + fname
    c = colors[i]
    if    path == "eclipse":
      ax.semilogy(bandcf[idx], press, '-', lw=lw, color=plt.cm.rainbow(c))
    elif  path == "transit":
      ax.plot(bandcf[idx], rad,       '-', lw=lw, color=plt.cm.rainbow(c))

    ax.set_ylim(yran)
    ax.set_xlim(xran)
    plt.text(0.9*xran[1], yran[1], fname, rotation=90, ha="right", va="top")
    ax.set_xticklabels([])
    if i == 0:
      ax.set_ylabel(ylabel, fontsize=fs)
    else:
      ax.set_yticklabels([])

  plt.subplots_adjust(0.12, 0.11, 0.97, 0.95, 0, 0)
  plt.suptitle(xlabel, fontsize=fs, y=0.09, x=0.52)

  if filename is not None:
    plt.savefig(filename)


def TCEA(posterior, pressure=None, tparams=None, tstepsize=None,
         besttpars=None, rstar=None, tstar=None, tint=None, smaxis=None,
         gplanet=None, pyrat=None):
  """
  Plot the posterior TCEA PT profile.

  Parameters
  ----------
  posterior: 2D float ndarray
     The MCMC posterior array of shape [nsamples, nparams]
  pressure: 1D float ndarray
     The atmospheric pressure profile in barye.
  tparams: 1D float ndarray
     The list of temperature-profile parameters.
  tstepsize: 1D float ndarray
     Stepsize of the temperature-profile parameters.
  besttpars: 1D float ndarray
     Array of best-fitting temperature-profile parameters.
  rstar: Float
     Stellar radius in cm.
  tstar: Float
     Stellar effective temperature in K.
  tint: Float
     Planetary internal temperature in K.
  smaxis: Float
     Semi-major axis in cm.
  gplanet: Float
     Planetary surface gravity in cm s-2.
  pyrat: Pyrat instance
  """
  if pyrat is not None:
    pressure  = pyrat.atm.press
    targs     = pyrat.ret.targs
    tparams   = pyrat.ret.params[pyrat.ret.itemp]
    tstepsize = pyrat.ret.stepsize[pyrat.ret.itemp]
    tmodel    = pyrat.ret.tmodel
  elif (pressure is None  or  tparams is None  or  tstepsize is None or
        rstar    is None  or  tstar   is None  or  tint      is None or
        smaxis   is None  or  gplanet is None):
    print("One or more input parameters is missing (pressure, tparams, "
          "tstepsize, rstar, tstar, tint, smaxis, gplanet).")

  if pyrat is None:
    targs = [pressure, rstar, tstar, tint, smaxis, gplanet]
    tmodel = PT.TCEA

  ifree = tstepsize > 0
  nfree = np.sum(ifree)
  ipost = np.arange(nfree)

  nlayers = len(pressure)
  nsamples, npars = np.shape(posterior)

  # Evaluate posterior PT profiles:
  PTprofiles = np.zeros((nsamples, nlayers), np.double)
  for i in np.arange(nsamples):
    tparams[ifree] = posterior[i, ipost]
    PTprofiles[i] = tmodel(tparams, *targs)

  # Get percentiles (for 1,2-sigma boundaries):
  low1 = np.percentile(PTprofiles, 16.0, axis=0)
  hi1  = np.percentile(PTprofiles, 84.0, axis=0)
  low2 = np.percentile(PTprofiles,  2.5, axis=0)
  hi2  = np.percentile(PTprofiles, 97.5, axis=0)
  median = np.median(PTprofiles, axis=0)

  # Plot figure:
  plt.figure(2)
  plt.clf()
  ax=plt.subplot(111)
  ax.fill_betweenx(pressure/pc.bar, low2, hi2,
                   facecolor="#62B1FF", edgecolor="0.5")
  ax.fill_betweenx(pressure/pc.bar, low1, hi1,
                   facecolor="#1873CC", edgecolor="#1873CC")
  plt.semilogy(median, pressure/pc.bar, "k-", lw=2, label="Median")
  if besttpars is not None:
    bestpt = tmodel(besttpars, *targs)
    plt.semilogy(bestpt, pressure/pc.bar, "r-", lw=2, label="Best fit")
  plt.ylim(np.amax(pressure/pc.bar), np.amin(pressure/pc.bar))
  plt.legend(loc="best")
  plt.xlabel("Temperature  (K)", size=15)
  plt.ylabel("Pressure  (bar)",  size=15)

  # Save figure:
  plt.savefig("MCMC_PT-profiles.pdf")
