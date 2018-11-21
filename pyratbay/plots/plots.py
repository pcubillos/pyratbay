# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

"""
Pyrat plotting utilities.
"""

__all__ = ["spectrum", "cf", "PT"]

import os, sys
import matplotlib

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
from scipy.ndimage.filters import gaussian_filter1d as gaussf

from .. import constants  as pc

rootdir = os.path.realpath(os.path.dirname(__file__) + "/../../")
sys.path.append(rootdir + "/pyratbay/lib/")
import pt as PT

sys.path.append(rootdir + "/pyratbay/atmosphere/")
import MadhuTP


def spectrum(wlength=None, spectrum=None, data=None, uncert=None,
    bandflux=None, bandtrans=None, bandidx=None, bandwl=None,
    starflux=None, rprs=None, path=None, logxticks=None, yran=None,
    pyrat=None, gaussbin=2.0, filename=None, fignum=-11):
  """
  Plot a transmission or emission model spectrum with (optional) data
  points with error bars and band-integrated model.

  Parameters
  ----------
  wlength: 1D float ndarray
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
  logxticks: 1D float ndarray
     If not None, switch the X-axis scale from linear to log, and set
     the X-axis ticks at the locations given by logxticks.
  yran: 1D float ndarray
     If not None, set the spectrum's Y-axis boundaries.
  pyrat: Pyrat instance
  gaussbin: Integer
     Standard deviation for Gaussian-kernel smoothing (in number of samples).
  filename: String
     Filename of the output figure.
  fignum: Integer
     The Figure number.

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
  ms  =  5
  lw  = 1.5
  mew = 0.4

  plt.figure(fignum, (8, 5))
  plt.clf()
  ax = plt.subplot(111)
  plt.subplots_adjust(0.14, 0.12, 0.97, 0.95)

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
    plt.ylabel("$(R_p/R_s)^2$", fontsize=fs)

  # Plot model:
  plt.plot(wlength, gmodel*fscale, lw=lw, label="Model", color="orange")
  # Plot band-integrated model:
  if bandwl is not None:
    plt.plot(bandwl, bandflux*fscale, "o", ms=ms, color="orange",
             mec="k", mew=mew)
  # Plot data:
  if data is not None:
    plt.errorbar(bandwl, data*fscale, uncert*fscale, fmt="ob", label="Data",
                 ms=ms, elinewidth=lw, capthick=lw, zorder=3)

  # Set Y-axis limits:
  if yran is not None:
    ax.set_ylim(np.array(yran)*fscale)
  yran = ax.get_ylim()  # Note this yran may differ from input (fscale).

  # Transmission filters:
  bandh = 0.06*(yran[1] - yran[0])
  for i in np.arange(nfilters):
    bandtr = bandh * bandtrans[i]/np.amax(bandtrans[i])
    plt.plot(wlength[bandidx[i]], yran[0]+bandtr, "0.4", zorder=-100)
  ax.set_ylim(yran)
  if logxticks is not None:
    ax.set_xscale('log')
    plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticks(logxticks)

  ax.tick_params(labelsize=fs-2)
  plt.xlabel("Wavelength  (um)", fontsize=fs)
  leg = plt.legend(loc="best", numpoints=1, fontsize=fs-1)
  plt.xlim(np.amin(wlength), np.amax(wlength))

  if filename is not None:
    plt.savefig(filename)

  return ax


def cf(bandcf, bandwl, path, pressure, radius, rtop=0,
       filename=None, filters=None, fignum=-21):
  """
  Plot the band-integrated normalized contribution functions
  (emission) or transmittance (transmission).

  Parameters
  ----------
  bandcf: 2D float ndarray
     Band-integrated contribution functions [nfilters, nlayers].
  bandwl: 1D float ndarray
     Mean wavelength of the bands in microns.
  path: String
     Observing geometry (transit or eclipse).
  pressure: 1D float ndarray
     Layer's pressure array (barye units).
  radius: 1D float ndarray
     Layer's impact parameter array (cm units).
  rtop: Integer
     Index of topmost valid layer.
  filename: String
     Filename of the output figure.
  filters: 1D string ndarray
     Name of the filter bands (optional).
  fignum: Integer
     Figure number.

  Notes
  -----
  - The dashed lines denote the 0.16 and 0.84 percentiles of the
    cumulative contribution function or the transmittance (i.e.,
    the boundaries of the central 68% of the respective curves).
  - If there are more than 80 filters, this code will thin the
    displayed filter names.
  """
  nfilters = len(bandwl)
  nlayers  = len(pressure)

  wlsort = np.argsort(bandwl)
  press = pressure[rtop:]/pc.bar
  rad   = radius  [rtop:]/pc.km

  press = pressure[rtop:]/pc.bar
  rad   = radius[rtop:]/pc.km
  if   path == "eclipse":
    yran = np.amax(np.log10(press)), np.amin(np.log10(press))
    zz = bandcf/np.amax(bandcf)
    xlabel = 'contribution function'
    ylabel = ''
    yright = 0.9
    cbtop  = 0.5
  elif path == "transit":
    zz = bandcf/np.amax(bandcf)
    yran = np.amin(rad), np.amax(rad)
    xlabel = r'transmittance'
    ylabel = r'Impact parameter (km)'
    yright = 0.84
    cbtop  = 0.8
  else:
    print("Invalid geometry.  Select from: 'eclipse' or 'transit'.")
    return

  fs  = 12
  colors = np.asarray(np.linspace(0, 255, nfilters), np.int)
  # 68% percentile boundaries of the central cumulative function:
  lo = 0.5*(1-0.683)
  hi = 1.0 - lo
  # Filter fontsize and thinning:
  ffs = 8.0 + (nfilters<50) + (nfilters<65)
  thin = (nfilters>80) + (nfilters>125) + (nfilters<100) + nfilters//100

  # Colormap and percentile limits:
  z = np.empty((nlayers, nfilters, 4), dtype=float)
  plo = np.zeros(nfilters+1)
  phi = np.zeros(nfilters+1)
  for i in np.arange(nfilters):
    z[:,i, :] = plt.cm.rainbow(colors[i])
    z[:,i,-1] = zz[i]**(0.5+0.5*(path=='transit'))
    if path == "eclipse":
      cumul = np.cumsum(zz[i])/np.sum(zz[i])
      plo[i], phi[i] = press[cumul>lo][0], press[cumul>hi][0]
    elif path == "transit":
      plo[i], phi[i] = press[zz[i]<lo][0], press[zz[i]<hi][0]
  plo[-1] = plo[-2]
  phi[-1] = phi[-2]

  fig = plt.figure(fignum, (8.5, 5))
  plt.clf()
  plt.subplots_adjust(0.105, 0.10, yright, 0.95)
  ax = plt.subplot(111)
  pax = ax.twinx()
  if   path == "eclipse":
    ax.imshow(z[:,wlsort], aspect='auto', extent=[0,nfilters,yran[0],yran[1]],
              origin='upper', interpolation='nearest')
    ax.yaxis.set_visible(False)
    pax.spines["left"].set_visible(True)
    pax.yaxis.set_label_position('left')
    pax.yaxis.set_ticks_position('left')
  elif path == "transit":
    ax.imshow(z[:,wlsort], aspect='auto', extent=[0,nfilters,yran[0],yran[1]],
              origin='upper', interpolation='nearest')
    # Setting the right radius tick labels requires some sorcery:
    fig.canvas.draw()
    ylab = [l.get_text() for l in ax.get_yticklabels()]
    rint = si.interp1d(rad, press, bounds_error=False)
    pticks = rint(ax.get_yticks())
    bounds = np.isfinite(pticks)
    pint = si.interp1d(press, np.linspace(yran[1], yran[0], nlayers),
                       bounds_error=False)
    ax.set_yticks(pint(pticks[bounds]))
    ax.set_yticklabels(np.array(ylab)[bounds])

  pax.plot(plo, drawstyle="steps-post", color="0.25", lw=0.75, ls="--")
  pax.plot(phi, drawstyle="steps-post", color="0.25", lw=0.75, ls="--")
  pax.set_yscale('log')
  pax.set_ylim(np.amax(press), np.amin(press))
  pax.set_ylabel(r'Pressure (bar)', fontsize=fs)

  ax.set_xlim(0, nfilters)
  ax.set_ylim(yran)
  ax.set_xticklabels([])
  ax.set_ylabel(ylabel, fontsize=fs)
  ax.set_xlabel("Band-averaged {:s}".format(xlabel), fontsize=fs)

  # Print filter names/wavelengths:
  for i in np.arange(0, nfilters-thin//2, thin):
    idx = wlsort[i]
    fname = " {:5.2f} um ".format(bandwl[idx])
    # Strip root and file extension:
    if filters is not None:
      fname = os.path.split(os.path.splitext(filters[idx])[0])[1] + " @" + fname
    ax.text(i+0.1, yran[1], fname, rotation=90, ha="left", va="top",
            fontsize=ffs)

  # Color bar:
  cbar = plt.axes([0.925, 0.10, 0.015, 0.85])
  cz = np.zeros((100, 2, 4), dtype=float)
  cz[:,0,3] = np.linspace(0.0,cbtop,100)**(0.5+0.5*(path=='transit'))
  cz[:,1,3] = np.linspace(0.0,cbtop,100)**(0.5+0.5*(path=='transit'))
  cbar.imshow(cz, aspect='auto', extent=[0, 1, 0, 1],
              origin='lower', interpolation='nearest')
  if path == "transit":
    cbar.axhline(0.1585, color="k", lw=1.0, dashes=(2.5,1))
    cbar.axhline(0.8415, color="w", lw=1.0, dashes=(2.5,1))
  cbar.spines["right"].set_visible(True)
  cbar.yaxis.set_label_position('right')
  cbar.yaxis.set_ticks_position('right')
  cbar.set_ylabel(xlabel.capitalize(), fontsize=fs)
  cbar.xaxis.set_visible(False)

  fig.canvas.draw()
  if filename is not None:
    plt.savefig(filename)


def PT(posterior, pressure=None, tparams=None, tstepsize=None,
       besttpars=None, rstar=None, tstar=None, tint=None, smaxis=None,
       gplanet=None, filename=None, pyrat=None):
  """
  Plot the posterior PT profile.

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
  filename: String
     Filename of the output figure.
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

  if pyrat is None and pyrat.ret.tmodel=="MadhuInv":
    targs = [pyrat.atm.press*1e-6]
    tmodel = MadhuTP.inversion

  if pyrat is None and pyrat.ret.tmodel=="MadhuNoInv":
    targs = [pyrat.atm.press*1e-6]
    tmodel = MadhuTP.no_inversion

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
  if filename is not None:
    plt.savefig(filename)
