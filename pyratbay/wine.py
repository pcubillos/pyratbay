__all__ = ["readfilter", "resample", "bandintegrate"]

import numpy as np
import scipy.constants   as sc
import scipy.interpolate as si

from . import constants as pc

"""
WINE: Waveband INtegrated Emission module

This set of routines read waveband filters and compute band-integrated
fluxes over the filter transmission curve.
"""


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
  wnindices: 1D ndarray
     The indices of specwn covered by the filter.
  istarfl: 1D ndarray
     The interpolated stellar flux.
  """
  # Indices in the spectrum wavenumber array included in the band
  # wavenumber range:
  wnindices = np.where((specwn < filterwn[-1]) & (filterwn[0] < specwn))

  # Make function to spline-interpolate the filter and stellar flux:
  finterp = si.interp1d(filterwn, filtertr)
  # Evaluate over the spectrum wavenumber array:
  ifilter = finterp(specwn[wnindices])
  # Normalize to integrate to 1.0:
  nifilter = ifilter/np.trapz(ifilter, specwn[wnindices])

  if starfl is not None and starwn is not None:
    sinterp = si.interp1d(starwn,   starfl)
    istarfl = sinterp(specwn[wnindices])
  else:
    istarfl = None

  # Return the normalized interpolated filter and the indices:
  return nifilter, wnindices, istarfl


def bandintegrate(spectrum, specwn, nifilter):
  """
  Integrate a spectrum over the band transmission.

  Parameters
  ----------
  spectrum: 1D ndarray
     Spectral signal to be integrated
  specwn: 1D ndarray
     Wavenumber of spectrum in cm-1
  nifilter: 1D ndarray
     The normalized interpolated filter transmission curve.

  Example
  -------
  >>> import sys
  >>> pbpath = "../Pyrat-Bay/"
  >>> sys.path.append(pbpath)
  >>> import pyratbay.pyratbay as pb
  
  >>> # Get a stellar spectrum:
  >>> kmodel = "fp00k2odfnew.pck"
  >>> sflux, swn, tm, gm = pb.starspec.readkurucz(kmodel, 5800, 4.43)
  
  >>> # Load Spitzer IRAC filters:
  >>> wn1, irac1 = pb.w.readfilter(pbpath+"inputs/filters/spitzer_irac1_sa.dat")
  >>> wn2, irac2 = pb.w.readfilter(pbpath+"inputs/filters/spitzer_irac2_sa.dat")
  
  >>> # Resample the filters into the stellar wavenumber array:
  >>> nifilter, wnindices = pb.w.resample(swn, wn1, irac1)
  >>> # Integrate the spectrum over the filter band:
  >>> bandflux1 = pb.w.bandintegrate(sflux[wnindices], swn[wnindices], nifilter)
  >>> nifilter, wnindices = pb.w.resample(swn, wn2, irac2)
  >>> bandflux2 = pb.w.bandintegrate(sflux[wnindices], swn[wnindices], nifilter)
  
  >>> # Plot the results:
  >>> plt.figure(1, (8,5))
  >>> plt.clf()
  >>> plt.semilogy(1e4/swn, sflux, "b")
  >>> plt.plot(np.mean(1e4/wn1), bandflux1, "o", color="red")
  >>> plt.plot(np.mean(1e4/wn2), bandflux2, "o", color="limegreen")
  >>> plt.plot(1e4/wn1, (irac1+1)*2e5, "red")
  >>> plt.plot(1e4/wn2, (irac2+1)*2e5, "limegreen")
  >>> plt.xlim(0.3, 6.0)
  >>> plt.ylim(2e5, 6e6)
  >>> plt.xlabel("Wavelength  (um)")
  >>> plt.ylabel(r"Flux  (erg s$^{-1}$ cm$^{-2}$ cm)")
  """

  return np.trapz(spectrum*nifilter, specwn)


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this module
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__ ):
        del locals()[varname]
del(varname)
