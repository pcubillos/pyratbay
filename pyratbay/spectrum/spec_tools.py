# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'tophat',
    'resample',
    'band_integrate',
    ]

from collections import Iterable

import numpy as np
import scipy.interpolate as si

from .. import constants as pc
from .. import io as io


def tophat(wl0, width, margin=None, dlambda=None, resolution=None, ffile=None):
    r"""
    Generate a top-hat filter function, with transmission = 1.0 from
    wl0-width/2 to wl0+width/2, and an extra margin with transmission
    = 0.0 at each end.

    Parameters
    ----------
    ffile: String
        Name of the output file.
    wl0:  Float
        Filter central wavelength in microns.
    width: Float
        Filter width in microns.
    margin: Float
        Margin (in microns) with zero-valued transmission, to append
        at each end of the filter.
    dlambda: Float
        Spectral sampling rate in microns.
    resolution: Float
        Spectral sampling resolution (used if dlambda is None).
    ffile: String
        If not None, save filter to file.

    Examples
    --------
    >>> import pyratbay.spectrum as ps
    >>> wl0     = 1.50
    >>> width   = 0.50
    >>> margin  = 0.10
    >>> dlambda = 0.05
    >>> wl, trans = ps.tophat(wl0, width, margin, dlambda)
    >>> print(wl, trans, sep='\n')
    [1.15 1.2  1.25 1.3  1.35 1.4  1.45 1.5  1.55 1.6  1.65 1.7  1.75 1.8
     1.85]
    [0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0.]
    """
    if margin is None:
        margin = 0.1 * width

    if dlambda is None and resolution is None:
        raise ValueError('Either dlambda or resolution must be defined.')

    # Wavelength array:
    wllow  = wl0 - 0.5*width - margin
    wlhigh = wl0 + 0.5*width + margin
    if dlambda is not None:
        wl = np.arange(wllow, wlhigh, dlambda)
    else:
        f = 0.5 / resolution
        g = (1.0-f) / (1.0+f)
        imax = int(np.ceil(np.log(wllow/wlhigh) / np.log(g))) + 1
        dwl = wlhigh * g**np.arange(imax)
        wl = 0.5 * np.flip(dwl[1:] + dwl[:-1], axis=0)

    transmission = np.array(np.abs(wl-wl0) < 0.5*width, np.double)

    if ffile is not None:
        io.write_spectrum(wl*pc.um, transmission, ffile, type='filter')

    return wl, transmission


def resample(signal, wn, specwn, normalize=False):
    r"""
    Resample signal from wn to specwn wavenumber sampling using a linear
    interpolation.

    Parameters
    ----------
    signal: 1D ndarray
        A spectral signal sampled at wn.
    wn: 1D ndarray
        Signal's wavenumber sampling, in cm-1 units.
    specwn: 1D ndarray
        Wavenumber sampling to resample into, in cm-1 units.
    normalize: Bool
        If True, normalized the output resampled signal to integrate to
        1.0 (note that a normalized curve when specwn is a decreasing
        function results in negative values for resampled).

    Returns
    -------
    resampled: 1D ndarray
        The interpolated signal.
    wnidx: 1D ndarray
        The indices of specwn covered by the input signal.

    Examples
    --------
    >>> import pyratbay.spectrum as ps
    >>> import numpy as np
    >>> wn     = np.linspace(1.3, 1.7, 11)
    >>> signal = np.array(np.abs(wn-1.5)<0.1, np.double)
    >>> specwn = np.linspace(1, 2, 51)
    >>> resampled, wnidx = ps.resample(signal, wn, specwn)
    >>> print(wnidx, specwn[wnidx], resampled, sep='\n')
    [16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34]
    [1.32 1.34 1.36 1.38 1.4  1.42 1.44 1.46 1.48 1.5  1.52 1.54 1.56 1.58
     1.6  1.62 1.64 1.66 1.68]
    [0.  0.  0.  0.  0.5 1.  1.  1.  1.  1.  1.  1.  1.  1.  0.5 0.  0.  0.
     0. ]
    """
    if np.amin(wn) < np.amin(specwn) or np.amax(wn) > np.amax(specwn):
        raise ValueError("Resampling signal's wavenumber is not contained "
                         "in specwn.")

    # Indices in the spectrum wavenumber array included in the band
    # wavenumber range:
    wnidx = np.where((specwn < np.amax(wn)) & (np.amin(wn) < specwn))[0]

    # Spline-interpolate:
    resampled = si.interp1d(wn,signal)(specwn[wnidx])

    if normalize:
        resampled /= np.trapz(resampled, specwn[wnidx])

    # Return the normalized interpolated filter and the indices:
    return resampled, wnidx


def band_integrate(spectrum, specwn, bandtrans, bandwn):
    """
    Integrate a spectrum over the band transmission.

    Parameters
    ----------
    spectrum: 1D float iterable
        Spectral signal to be integrated.
    specwn: 1D float iterable
        Wavenumber of spectrum in cm-1.
    bandtrans: 1D float iterable
        List of normalized interpolated band transmission values in each filter.
    bandwn: 1D float iterable

    Returns
    -------
    bflux: 1D float list
        Band-integrated values.

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pyratbay.io as io
    >>> import pyratbay.spectrum as ps
    >>> import pyratbay.constants as pc
    >>> # Load Spitzer IRAC filters:
    >>> wn1, irac1 = io.read_spectrum(
    >>>     pc.ROOT+'pyratbay/data/filters/spitzer_irac1_sa.dat')
    >>> wn2, irac2 = io.read_spectrum(
    >>>     pc.ROOT+'pyratbay/data/filters/spitzer_irac2_sa.dat')
    >>> # Spectrum to integrate:
    >>> wn = np.arange(1500, 5000.1, 1.0)
    >>> sflux = ps.bbflux(wn, 1800.0)
    >>> # Integrate over single filter:
    >>> bandflux = ps.band_integrate(sflux, wn, irac1, wn1)
    >>> # Integrate over multiple:
    >>> bandfluxes = ps.band_integrate(sflux, wn, [irac1,irac2], [wn1, wn2])
    >>> # Plot the results:
    >>> meanwn = [np.mean(wn1), np.mean(wn2)]
    >>> width = 0.5*(np.amax(wn1)-np.amin(wn1)), 0.5*(np.amax(wn2)-np.amin(wn2))
    >>> plt.figure(1)
    >>> plt.clf()
    >>> plt.semilogy(wn, sflux, 'k')
    >>> plt.plot(wn1, (irac1+1)*4e4, 'red')
    >>> plt.plot(wn2, (irac2+1)*4e4, 'blue')
    >>> plt.errorbar(meanwn[0], bandflux, xerr=width[0], fmt='o', color='red')
    >>> plt.errorbar(meanwn, bandfluxes, xerr=width, fmt='o', color='none',
    >>>              mec='k', ecolor='k')
    >>> plt.xlim(np.amin(wn), np.amax(wn))
    >>> plt.ylim(4e4, 1.2e5)
    >>> plt.xlabel('Wavenumber  (cm$^{-1}$)')
    >>> plt.ylabel(r'Flux  (erg s$^{-1}$ cm$^{-2}$ cm)')
    """
    if not isinstance(bandtrans[0], Iterable):
        bandtrans = [bandtrans]
        bandwn    = [bandwn]

    bflux = []
    for btrans, wn in zip(bandtrans, bandwn):
        # Resample bandpasses into spectrum wavenumber sampling:
        resampled, wnidx = resample(btrans, wn, specwn, normalize=True)
        # Band-integrate spectrum:
        bflux.append(np.trapz(spectrum[wnidx]*resampled, specwn[wnidx]))

    return bflux

