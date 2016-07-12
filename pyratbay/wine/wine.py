# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["readfilter", "resample", "bandintegrate"]

import numpy as np
import scipy.constants   as sc
import scipy.interpolate as si

from .. import constants as pc


def readfilter(filt):
  """
  Load a filter bandpass from file.

  Parameters
  ----------
  filt: String
      Filter file name.

  Return
  ------
  wavenumber: 1D ndarray
     The filter pass band wavenumber in cm-1.
  transmission: 1D ndarray 
     The filter spectral response. No specific units.

  Notes
  -----
  - The file can contains empty lines and comments (with '#' character)
    before the data.  No comments or empty lines after the data.
  - The data must come in two columns.  The first column must contain
    the wavelength in microns, the second the filter response, other
    columns will be ignored.
  """
  # Open and read the filter file:
  data = open(filt, "r")
  lines = data.readlines()
  data.close()

  # Remove header comments and empty lines:
  while lines[0].startswith("#") or not lines[0].strip():
    comment = lines.pop(0)

  # Allocate arrays for the wavelength and response:
  nlines = len(lines)
  wavel  = np.zeros(nlines, np.double) # filter's wavelengths  (in microns)
  transm = np.zeros(nlines, np.double) # filter's pass bands

  # Read the data and store in reverse order:
  for i in np.arange(nlines):
    wavel[nlines-1-i], transm[nlines-1-i] = lines[i].strip().split()[0:2]

  m2cm  = 1e-4 # Microns to cm conversion factor
  # Get wavenumber in cm-1:
  waven = 1.0 / (wavel*m2cm)
  # Return statement:
  return waven, transm


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
  wnidx = np.where((specwn < filterwn[-1]) & (filterwn[0] < specwn))[0]

  # Make function to spline-interpolate the filter and stellar flux:
  finterp = si.interp1d(filterwn, filtertr)
  # Evaluate over the spectrum wavenumber array:
  ifilter = finterp(specwn[wnidx])
  # Normalize to integrate to 1.0:
  nifilter = ifilter/np.trapz(ifilter, specwn[wnidx])

  if starfl is not None and starwn is not None:
    sinterp = si.interp1d(starwn,   starfl)
    istarfl = sinterp(specwn[wnidx])
  else:
    istarfl = None

  # Return the normalized interpolated filter and the indices:
  return nifilter, wnidx, istarfl


def bandintegrate(spectrum=None, specwn=None, wnidx=None, nfilters=None,
                  bandtrans=None, starflux=None, rprs=None, path=None,
                  pyrat=None):
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
  nfilters: Integer
     Number of filters.
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
  
  >>> # Get a stellar spectrum:
  >>> kmodel = "fp00k2odfnew.pck"
  >>> sflux, swn, tm, gm = pb.starspec.readkurucz(kmodel, 5800, 4.43)
  
  >>> # Load Spitzer IRAC filters:
  >>> wn1, irac1 = pb.wine.readfilter(pbpath +
                       "inputs/filters/spitzer_irac1_sa.dat")
  >>> wn2, irac2 = pb.wine.readfilter(pbpath +
                       "inputs/filters/spitzer_irac2_sa.dat")
  
  >>> # Resample the filters into the stellar wavenumber array:
  >>> nifilter1, wnidx1 = pb.wine.resample(swn, wn1, irac1)
  >>> nifilter2, wnidx2 = pb.wine.resample(swn, wn2, irac2)
  >>> # Integrate the spectrum over the filter band:
  >>> bandflux = pb.wine.bandintegrate(spectrum=sflux, specwn=swn,
       wnidx=[wnidx1,wnidx2], nfilters=2, bandtrans=[nifilter1, nifilter2],
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
    wn        = pyrat.spec.wn
    bflux     = pyrat.obs.bandflux
    wnidx     = pyrat.obs.bandidx
    nfilters  = pyrat.obs.nfilters
    bandtrans = pyrat.obs.bandtrans
    starflux  = pyrat.obs.starflux
    rprs      = pyrat.phy.rprs
    path      = pyrat.od.path

  if bflux is None:
    bflux = np.zeros(nfilters)

  # Band-integrate spectrum:
  for i in np.arange(nfilters):
    # Integrate the spectrum over the filter band:
    if   path == "transit":
      bflux[i] = np.trapz(spectrum[wnidx[i]]*bandtrans[i], wn[wnidx[i]])
    elif path == "eclipse":
      fluxrat = spectrum[wnidx[i]]/starflux[i] * rprs**2.0
      bflux[i] = np.trapz(fluxrat*bandtrans[i], wn[wnidx[i]])

  return bflux
