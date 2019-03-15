# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["resample", "bandintegrate"]

import numpy as np
import scipy.constants   as sc
import scipy.interpolate as si

from .. import constants as pc


def resample(specwn, filterwn, filtertr, starwn=None, starfl=None):
  """
  Resample the filtertr curve from the filterwn sampling into specwn

  Parameters
  ----------
  specwn: 1D ndarray
     A wavenumber sampling array (in cm-1).
  filterwn: 1D ndarray
     Filter wavenumber sampling array (in cm-1).
  filtertr: 1D ndarray
     Filter transmission curve sampled as filterwn.
  starwn: 1D ndarray
     Stellar model wavenumber sampling array (in cm-1).
  starfl: 1D ndarray
     Stellar flux.

  Returns
  -------
  nifilter: 1D ndarray
     The normalized interpolated filter transmission curve.
  wnidx: 1D ndarray
     The indices of specwn covered by the filter.
  istarfl: 1D ndarray
     The interpolated stellar flux.
  """
  # Indices in the spectrum wavenumber array included in the band
  # wavenumber range:
  wnidx = np.where((specwn < np.amax(filterwn))
                 & (np.amin(filterwn) < specwn))[0]

  # Make function to spline-interpolate the filter and stellar flux:
  finterp = si.interp1d(filterwn, filtertr)
  # Evaluate over the spectrum wavenumber array:
  ifilter = finterp(specwn[wnidx])
  # Normalize to integrate to 1.0:
  nifilter = ifilter/np.trapz(ifilter, specwn[wnidx])

  if starfl is not None and starwn is not None:
      sinterp = si.interp1d(starwn, starfl)
      istarfl = sinterp(specwn[wnidx])
  else:
      istarfl = None

  # Return the normalized interpolated filter and the indices:
  return nifilter, wnidx, istarfl


def bandintegrate(spectrum=None, specwn=None, wnidx=None, bandtrans=None,
                  starflux=None, rprs=None, path='transit', pyrat=None):
  """
  Integrate a spectrum over the band transmission.

  Parameters
  ----------
  spectrum: 1D float ndarray
     Spectral signal to be integrated.
  specwn: 1D float ndarray
     Wavenumber of spectrum in cm-1.
  wnidx: List of 1D integer ndarrays
     List of indices of specwn covered by each filter.
  bandtrans: List of 1D float ndarray
     List  of normalized interpolated band transmission values in each filter.
  starflux: List of 1D float ndarray
     List of interpolated stellar flux values in each filter (eclipse geometry).
  rprs: Float
     Planet-to-star radius ratio (used for eclipse geometry).
  path: String
     Observing geometry.  Select from: transit or eclipse.
  pyrat: A Pyrat object.
     If not None, extract input values from pyrat object.

  Returns
  -------
  bflux: 1D ndarray
     Array of band-integrated values.

  Example
  -------
  >>> import sys
  >>> pbpath = "../pyratbay/"
  >>> sys.path.append(pbpath)
  >>> import pyratbay as pb
  >>> import pyratbay.io as io
  
  >>> # Get a stellar spectrum:
  >>> kmodel = "fp00k2odfnew.pck"
  >>> sflux, swn, tm, gm = pb.starspec.readkurucz(kmodel, 5800, 4.43)
  
  >>> # Load Spitzer IRAC filters:
  >>> wn1, irac1 = io.read_filter(pbpath +
                       "inputs/filters/spitzer_irac1_sa.dat")
  >>> wn2, irac2 = io.read_filter(pbpath +
                       "inputs/filters/spitzer_irac2_sa.dat")
  
  >>> # Resample the filters into the stellar wavenumber array:
  >>> nifilter1, wnidx1 = pb.wine.resample(swn, wn1, irac1)
  >>> nifilter2, wnidx2 = pb.wine.resample(swn, wn2, irac2)
  >>> # Integrate the spectrum over the filter band:
  >>> bandflux = pb.wine.bandintegrate(spectrum=sflux, specwn=swn,
       wnidx=[wnidx1,wnidx2], bandtrans=[nifilter1, nifilter2],
       path='transit')

  >>> # Plot the results:
  >>> plt.figure(1, (8,5))
  >>> plt.clf()
  >>> plt.semilogy(1e4/swn, sflux, "b")
  >>> plt.plot(np.mean(1e4/wn1), bandflux[0], "o", color="red")
  >>> plt.plot(np.mean(1e4/wn2), bandflux[1], "o", color="limegreen")
  >>> plt.plot(1e4/wn1, (irac1+1)*2e5, "red")
  >>> plt.plot(1e4/wn2, (irac2+1)*2e5, "limegreen")
  >>> plt.xlim(0.3, 6.0)
  >>> plt.ylim(2e5, 6e6)
  >>> plt.xlabel("Wavelength  (um)")
  >>> plt.ylabel(r"Flux  (erg s$^{-1}$ cm$^{-2}$ cm)")
  """

  if pyrat is not None:
    # Unpack variables from pyrat object:
    spectrum  = pyrat.spec.spectrum
    specwn    = pyrat.spec.wn
    bflux     = pyrat.obs.bandflux
    wnidx     = pyrat.obs.bandidx
    nfilters  = pyrat.obs.nfilters
    bandtrans = pyrat.obs.bandtrans
    starflux  = pyrat.obs.starflux
    rprs      = pyrat.phy.rprs
    path      = pyrat.od.path
  else:
    if np.isscalar(wnidx[0]):  # A single filter
      nfilters  = 1
      wnidx     = [wnidx]
      bandtrans = [bandtrans]
      starflux  = [starflux]
    else:  # Multiple filters
      nfilters = len(wnidx)
    bflux = np.zeros(nfilters)

  if nfilters == 0:
    return None

  # Band-integrate spectrum:
  for i in np.arange(nfilters):
    # Integrate the spectrum over the filter band:
    if   path == "transit":
      bflux[i] = np.trapz(spectrum[wnidx[i]]*bandtrans[i], specwn[wnidx[i]])
    elif path == "eclipse":
      fluxrat = spectrum[wnidx[i]]/starflux[i] * rprs**2.0
      bflux[i] = np.trapz(fluxrat*bandtrans[i], specwn[wnidx[i]])

  return bflux
