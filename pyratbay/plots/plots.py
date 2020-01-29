# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = [
    'spectrum',
    'cf',
    'posterior_pt',
]

import os

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
from scipy.ndimage.filters import gaussian_filter1d as gaussf

from .. import constants as pc
from .. import tools as pt


def spectrum(spectrum, wavelength, path,
             data=None, uncert=None, bandwl=None, bandflux=None,
             bandtrans=None, bandidx=None,
             starflux=None, rprs=None, label='model', bounds=None,
             logxticks=None,
             gaussbin=2.0, yran=None, filename=None, fignum=501,
             axis=None):
  """
  Plot a transmission or emission model spectrum with (optional) data
  points with error bars and band-integrated model.

  Parameters
  ----------
  spectrum: 1D float ndarray
      Planetary spectrum evaluated at wavelength.
  wavelength: 1D float ndarray
      The wavelength of the model in microns.
  path: String
      Observing-geometry path: transit or eclipse.
  data: 1D float ndarray
      Observing data points at each bandwl.
  uncert: 1D float ndarray
      Uncertainties of the data points.
  bandwl: 1D float ndarray
      The mean wavelength for each band/data point.
  bandflux: 1D float ndarray
      Band-integrated model spectrum at each bandwl.
  bandtrans: List of 1D float ndarrays
      Transmission curve for each band.
  bandidx: List of 1D float ndarrays.
      The indices in wavelength for each bandtrans.
  starflux: 1D float ndarray
      Stellar spectrum evaluated at wavelength.
  rprs: Float
      Planet-to-star radius ratio.
  label: String
      Label for spectrum curve.
  bounds: Tuple
      Tuple with -2, -1, +1, and, +2 sigma boundaries of spectrum.
      If not None, plot shaded area between +/-1sigma and +/-2sigma
      boundaries.
  logxticks: 1D float ndarray
      If not None, switch the X-axis scale from linear to log, and set
      the X-axis ticks at the locations given by logxticks.
  gaussbin: Integer
      Standard deviation for Gaussian-kernel smoothing (in number of samples).
  yran: 1D float ndarray
      Figure's Y-axis boundaries.
  filename: String
      If not None, save figure to filename.
  fignum: Integer
      Figure number.
  axis: TBD
      TBD

  Returns
  -------
  ax: AxesSubplot instance
      The matplotlib Axes of the figure.
  """
  # Plotting setup:
  fs = 14.0
  ms =  6.0
  lw =  1.25

  if axis is None:
      plt.figure(fignum, (8, 5))
      plt.clf()
      ax = plt.subplot(111)
  else:
      ax = axis

  #fscale = {'':1.0, '%':100.0, 'ppt':1e3, 'ppm':1e6}

  spec_kw = {'label':label}
  if bounds is None:
      spec_kw['color'] = 'orange'
  else:
      spec_kw['color'] = 'orangered'


  # Setup according to geometry:
  if path == 'eclipse':
      if starflux is not None and rprs is not None:
          spectrum = spectrum/starflux * rprs**2.0
          if bounds is not None:
              bounds = [bound/starflux * rprs**2.0 for bound in bounds]
          fscale = 1e3
          plt.ylabel(r'$F_{\rm p}/F_{\rm s}\ (10^{-3})$', fontsize=fs)
      else:
          fscale = 1.0
          plt.ylabel(r'$F_{\rm p}$ (erg s$^{-1}$ cm$^{-2}$ cm)', fontsize=fs)
  elif path == 'transit':
      fscale = 100.0
      plt.ylabel(r'$(R_{\rm p}/R_{\rm s})^2$  (%)', fontsize=fs)

  gmodel = gaussf(spectrum, gaussbin)
  if bounds is not None:
      gbounds = [gaussf(bound, gaussbin) for bound in bounds]
      ax.fill_between(wavelength, fscale*gbounds[0], fscale*gbounds[3],
          facecolor='gold', edgecolor='none')
      ax.fill_between(wavelength, fscale*gbounds[1], fscale*gbounds[2],
          facecolor='orange', edgecolor='none')

  # Plot model:
  plt.plot(wavelength, gmodel*fscale, lw=lw, **spec_kw)
  # Plot band-integrated model:
  if bandflux is not None and bandwl is not None:
      plt.plot(bandwl, bandflux*fscale, 'o', ms=ms, color='tomato',
               mec='maroon', mew=lw)
  # Plot data:
  if data is not None and uncert is not None and bandwl is not None:
      plt.errorbar(bandwl, data*fscale, uncert*fscale, fmt='o', label='data',
                   color='blue', ms=ms, elinewidth=lw, capthick=lw, zorder=3)

  if yran is not None:
      ax.set_ylim(np.array(yran))
  yran = ax.get_ylim()

  # Transmission filters:
  if bandtrans is not None and bandidx is not None:
      bandh = 0.06*(yran[1] - yran[0])
      for btrans, bidx in zip(bandtrans, bandidx):
          btrans = bandh * btrans/np.amax(btrans)
          plt.plot(wavelength[bidx], yran[0]+btrans, '0.4', zorder=-10)
      ax.set_ylim(yran)

  if logxticks is not None:
      ax.set_xscale('log')
      plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
      ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
      ax.set_xticks(logxticks)

  ax.tick_params(labelsize=fs-2)
  plt.xlabel('Wavelength  (um)', fontsize=fs)
  plt.legend(loc='best', numpoints=1, fontsize=fs-1)
  plt.xlim(np.amin(wavelength), np.amax(wavelength))

  ax2 = ax.twinx()
  ax2.tick_params(right=True, direction='in')
  ax2.set_yticks(ax.get_yticks())
  ax2.set_yticklabels([])
  ax2.set_ylim(ax.get_ylim())
  plt.tight_layout()

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

  Returns
  -------
  ax: AxesSubplot instance
      The matplotlib Axes of the figure.

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
  if path == 'eclipse':
      yran = np.amax(np.log10(press)), np.amin(np.log10(press))
      zz = bandcf/np.amax(bandcf)
      xlabel = 'contribution function'
      ylabel = ''
      yright = 0.9
      cbtop  = 0.5
  elif path == 'transit':
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
      if path == 'eclipse':
          cumul = np.cumsum(zz[i])/np.sum(zz[i])
          plo[i], phi[i] = press[cumul>lo][0], press[cumul>hi][0]
      elif path == 'transit':
          plo[i], phi[i] = press[zz[i]<lo][0], press[zz[i]<hi][0]
  plo[-1] = plo[-2]
  phi[-1] = phi[-2]

  fig = plt.figure(fignum, (8.5, 5))
  plt.clf()
  plt.subplots_adjust(0.105, 0.10, yright, 0.95)
  ax = plt.subplot(111)
  pax = ax.twinx()
  if path == 'eclipse':
      ax.imshow(z[:,wlsort], aspect='auto',
                extent=[0, nfilters, yran[0], yran[1]],
                origin='upper', interpolation='nearest')
      ax.yaxis.set_visible(False)
      pax.spines['left'].set_visible(True)
      pax.yaxis.set_label_position('left')
      pax.yaxis.set_ticks_position('left')
  elif path == 'transit':
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

  pax.plot(plo, drawstyle='steps-post', color='0.25', lw=0.75, ls='--')
  pax.plot(phi, drawstyle='steps-post', color='0.25', lw=0.75, ls='--')
  pax.set_ylim(np.amax(press), np.amin(press))
  pax.set_yscale('log')
  pax.set_ylabel(r'Pressure (bar)', fontsize=fs)

  ax.set_xlim(0, nfilters)
  ax.set_ylim(yran)
  ax.set_xticklabels([])
  ax.set_ylabel(ylabel, fontsize=fs)
  ax.set_xlabel('Band-averaged {:s}'.format(xlabel), fontsize=fs)

  # Print filter names/wavelengths:
  for i in np.arange(0, nfilters-thin//2, thin):
      idx = wlsort[i]
      fname = ' {:5.2f} um '.format(bandwl[idx])
      # Strip root and file extension:
      if filters is not None:
          fname = (os.path.split(os.path.splitext(filters[idx])[0])[1]
                   + ' @' + fname)
      ax.text(i+0.1, yran[1], fname, rotation=90, ha='left', va='top',
              fontsize=ffs)

  # Color bar:
  cbar = plt.axes([0.925, 0.10, 0.015, 0.85])
  cz = np.zeros((100, 2, 4), dtype=float)
  cz[:,0,3] = np.linspace(0.0,cbtop,100)**(0.5+0.5*(path=='transit'))
  cz[:,1,3] = np.linspace(0.0,cbtop,100)**(0.5+0.5*(path=='transit'))
  cbar.imshow(cz, aspect='auto', extent=[0, 1, 0, 1],
              origin='lower', interpolation='nearest')
  if path == 'transit':
      cbar.axhline(0.1585, color='k', lw=1.0, dashes=(2.5,1))
      cbar.axhline(0.8415, color='w', lw=1.0, dashes=(2.5,1))
  cbar.spines['right'].set_visible(True)
  cbar.yaxis.set_label_position('right')
  cbar.yaxis.set_ticks_position('right')
  cbar.set_ylabel(xlabel.capitalize(), fontsize=fs)
  cbar.xaxis.set_visible(False)

  fig.canvas.draw()
  if filename is not None:
      plt.savefig(filename)
  return ax


def posterior_pt(posterior, tmodel, tpars, ifree, pressure,
                 bestpars=None, filename=None):
  """
  Plot the posterior PT profile.

  Parameters
  ----------
  posterior: 2D float ndarray
      MCMC posterior distribution for tmodel (of shape [nparams, nfree]).
  tmodel: Callable
      Temperature-profile model.
  tpars: 1D float ndarray
      Temperature-profile parameters (including fixed parameters).
  ifree: 1D bool ndarray
      Mask of free (True) and fixed (False) parameters in tpars.
      The number of free parameters must match nfree in posterior.
  pressure: 1D float ndarray
      The atmospheric pressure profile in barye.
  bestpars: 1D float ndarray
      Best-fitting temperature-profile parameters.
  filename: String
      If not None, save figure to filename.

  Returns
  -------
  ax: AxesSubplot instance
      The matplotlib Axes of the figure.
  """
  nlayers = len(pressure)

  u, uind, uinv = np.unique(posterior[:,0], return_index=True,
      return_inverse=True)
  nsamples = len(u)

  # Evaluate posterior PT profiles:
  profiles = np.zeros((nsamples, nlayers), np.double)
  for i in range(nsamples):
      tpars[ifree] = posterior[uind[i]]
      profiles[i] = tmodel(tpars)

  # Get percentiles (for 1,2-sigma boundaries and median):
  low1   = np.zeros(nlayers, np.double)
  low2   = np.zeros(nlayers, np.double)
  median = np.zeros(nlayers, np.double)
  high1  = np.zeros(nlayers, np.double)
  high2  = np.zeros(nlayers, np.double)
  for i in range(nlayers):
      tpost = profiles[uinv,i]
      low2[i]   = np.percentile(tpost,  2.275)
      low1[i]   = np.percentile(tpost, 15.865)
      median[i] = np.percentile(tpost, 50.000)
      high1[i]  = np.percentile(tpost, 84.135)
      high2[i]  = np.percentile(tpost, 97.725)

  # alpha != 0 does not work for ps/eps figures:
  if filename is not None and filename.endswith('ps'):
      fc1, fc2 = '#3366ff', '#77aaff'
      alpha1, alpha2 = 1.0, 1.0
  else:
      fc1, fc2 = 'royalblue', 'royalblue'
      alpha1, alpha2 = 0.8, 0.6

  # Plot figure:
  plt.figure(500)
  plt.clf()
  ax = plt.subplot(111)
  ax.fill_betweenx(pressure/pc.bar, low2, high2, facecolor=fc2,
      edgecolor='none', alpha=alpha2)
  ax.fill_betweenx(pressure/pc.bar, low1, high1, facecolor=fc1,
      edgecolor='none', alpha=alpha1)
  plt.plot(median, pressure/pc.bar, 'navy', lw=2, label='Median')
  if bestpars is not None:
      bestpt = tmodel(bestpars)
      plt.plot(bestpt, pressure/pc.bar, 'r-', lw=2, label='Best fit')
  ax.set_ylim(np.amax(pressure/pc.bar), np.amin(pressure/pc.bar))
  ax.set_yscale('log')
  plt.legend(loc='best')
  plt.xlabel('Temperature  (K)', size=15)
  plt.ylabel('Pressure  (bar)',  size=15)
  ax.tick_params(labelsize=12)

  plt.tight_layout()
  if filename is not None:
      plt.savefig(filename)
  return ax
