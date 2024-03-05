# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'PassBand',
    'Tophat',
    'constant_resolution_spectrum',
    'bin_spectrum',
    'tophat',
    'resample',
    'band_integrate',
]

from collections.abc import Iterable
import operator
import os

import numpy as np
import scipy.interpolate as si

from .. import constants as pc
from .. import io as io


class PassBand():
    """
    A Filter passband object.
    """
    def __init__(self, filter_file, wn=None):
        """
        Parameters
        ----------
        filter_file: String
            Path to filter file containing wavelength (um) and passband
            response in two columns.
            Comment and blank lines are ignored.

        Examples
        --------
        >>> import pyratbay.spectrum as ps
        >>> import pyratbay.constants as pc
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np

        >>> filter_file = f'{pc.ROOT}pyratbay/data/filters/spitzer_irac2_sa.dat'
        >>> band = ps.PassBand(filter_file)

        >>> # Evaluate over a wavelength array (um):
        >>> wl = np.arange(3.5, 5.5, 0.001)
        >>> out_wl, out_response = band(wl)

        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.plot(out_wl, out_response)
        >>> plt.plot(band.wl, band.response)  # Same variables
        >>> # Note wl differs from band.wl, but original array can be used as:
        >>> plt.plot(wl[band.idx], band.response)

        >>> # Evaluate over a wavenumber array:
        >>> wn = 1e4 / wl
        >>> band(wn=wn)
        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.plot(band.wn, band.response, dashes=(5,3))
        >>> plt.plot(wn[band.idx], out_response)
        """
        # Read filter wavenumber and transmission curves:
        filter_file = filter_file.replace('{ROOT}', pc.ROOT)
        self.filter_file = os.path.realpath(filter_file)
        input_wl, input_response = io.read_spectrum(self.filter_file, wn=False)

        self.wl0 = np.sum(input_wl*input_response) / np.sum(input_response)
        input_wn = 1.0 / (input_wl * pc.um)
        self.wn0 = 1.0 / (self.wl0 * pc.um)

        # Sort it in increasing wavenumber order, store it:
        wn_sort = np.argsort(input_wn)
        self.input_response = input_response[wn_sort]
        self.input_wn = input_wn[wn_sort]
        self.input_wl = input_wl[wn_sort]

        self.response = np.copy(input_response[wn_sort])
        self.wn = np.copy(input_wn[wn_sort])
        self.wl = 1.0 / (self.wn * pc.um)

        self.name = os.path.splitext(os.path.basename(filter_file))[0]

        # Resample the filters into the planet wavenumber array:
        if wn is not None:
            self.__eval__(wn=wn)

    def __repr__(self):
        return f"pyratbay.spectrum.PassBand('{self.filter_file}')"

    def __str__(self):
        return f'{self.name}'


    def __call__(self, wl=None, wn=None):
        """
        Interpolate filter response function at specified spectral array.
        The response funciton is normalized such that the integral over
        wavenumber equals one.

        Parameters
        ----------
        wl: 1D float array
            Wavelength array at which evaluate the passband response's
            in micron units.
            (only one of wl or wn should be provided on call)
        wn: 1D float array
            Wavenumber array at which evaluate the passband response's
            in cm-1 units.
            (only one of wl or wn should be provided on call)

        Defines
        -------
        self.response  Normalized interpolated response function
        self.idx       IndicesWavenumber indices
        self.wn        Passband's wavenumber array
        self.wl        Passband's wavelength array

        Returns
        -------
        out_wave: 1D float array
            Same as self.wl or self.wn depending on the input argument.
        out_response: 1D float array
            Same as self.response

        Examples
        --------
        >>> # See examples in help(ps.PassBand.__init__)
        """
        if wl is None and wn is None:
            raise ValueError(
                'Neither of wavelength (wl) nor wavenumber (wn) were provided'
            )
        if wl is not None and wn is not None:
            raise ValueError(
                'Either provide wavelength or wavenumber array, not both'
            )
        input_is_wl = wn is None
        if input_is_wl:
            wn = 1.0 / (wl*pc.um)

        sign = np.sign(np.ediff1d(wn))
        if not (np.all(sign == 1) or np.all(sign == -1)):
            raise ValueError(
                'Input wavelength/wavenumber array must be strictly '
                'increasing or decreasing'
            )
        sign = sign[0]

        response, wn_idx = resample(
            self.input_response, self.input_wn, wn, normalize=True,
        )

        # Internally, wavenumber is always monotonically increasing:
        wn_sort = np.argsort(wn[wn_idx])

        self.response = response[wn_sort] * sign
        self.wn = wn[wn_idx][wn_sort]
        self.idx = wn_idx[wn_sort]

        self.wl = 1.0 / (self.wn * pc.um)
        if input_is_wl:
            out_wave = self.wl
        else:
            out_wave = self.wn

        return out_wave, self.response

    def save_filter(self, save_file):
        """
        Write filter response function data to file, into two columns
        the wavelength (um) and the response function.

        Parameters
        ----------
        save_file: String
            File where to save the filter data.
        """
        io.write_spectrum(
            self.wl*pc.um, self.response, save_file, type='filter',
        )


class Tophat(PassBand):
    """
    A Filter passband object with a tophat-shaped passband.
    """
    def __init__(
        self, wl0, half_width,
        name='tophat',
    ):
        """
        Parameters
        ----------
        wl0: Float
            The passband's central wavelength (um units).
        half_width: Float
            The passband's half-width (um units).
        name: Str
            A user-defined name for the filter when calling str(self),
            e.g., to identify the instrument provenance of this filter.

        Examples
        --------
        >>> import pyratbay.spectrum as ps
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np

        >>> hat = ps.Tophat(4.5, 0.5)

        >>> # Evaluate over a wavelength array (um units):
        >>> wl = np.arange(3.5, 5.5, 0.001)
        >>> out_wl, out_response = hat(wl)

        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.plot(out_wl, out_response)
        >>> plt.plot(hat.wl, hat.response)  # Same variables
        >>> # Note wl differs from hat.wl, but original array can be used as:
        >>> plt.plot(wl[hat.idx], hat.response)

        >>> # Evaluate over a wavenumber array:
        >>> wn = 1e4 / wl
        >>> hat(wn=wn)
        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.plot(hat.wn, hat.response, dashes=(5,3))
        >>> plt.plot(wn[hat.idx], out_response)
        """
        self.wl0 = wl0
        self.half_width = half_width
        # Read filter wavenumber and transmission curves:
        self.wn0 = 1.0 / (self.wl0 * pc.um)

        self.name = name

    def __repr__(self):
        return f'pyratbay.spectrum.Tophat({self.wl0}, {self.half_width})'

    def __str__(self):
        return f'{self.name}_{self.wl0}um'


    def __call__(self, wl=None, wn=None):
        """
        Interpolate filter response function at specified spectral array.
        The response funciton is normalized such that the integral over
        wavenumber equals one.

        Parameters
        ----------
        wl: 1D float array
            Wavelength array at which evaluate the passband response's
            in micron units.
            (only one of wl or wn should be provided on call)
        wn: 1D float array
            Wavenumber array at which evaluate the passband response's
            in cm-1 units.
            (only one of wl or wn should be provided on call)

        Defines
        -------
        self.response  Normalized interpolated response function
        self.idx       IndicesWavenumber indices
        self.wn        Passband's wavenumber array
        self.wl        Passband's wavelength array

        Returns
        -------
        out_wave: 1D float array
            Same as self.wl or self.wn depending on the input argument.
        out_response: 1D float array
            Same as self.response

        Examples
        --------
        >>> # See examples in help(ps.Tophat.__init__)
        """
        if not operator.xor(wl is None, wn is None):
            raise ValueError(
                'Either provide wavelength or wavenumber array, not both'
            )
        input_is_wl = wn is None
        if input_is_wl:
            wn = 1.0 / (wl*pc.um)

        sign = np.sign(np.ediff1d(wn))
        if not (np.all(sign == 1) or np.all(sign == -1)):
            raise ValueError(
                'Input wavelength/wavenumber array must be strictly '
                'increasing or decreasing'
            )
        sign = sign[0]
        nwave = len(wn)

        wl_low  = self.wl0 - self.half_width
        wl_high = self.wl0 + self.half_width

        wn_low = 1.0 / (wl_high*pc.um)
        wn_high = 1.0 / (wl_low*pc.um)
        idx = (wn >= wn_low) & (wn <= wn_high)
        indices = np.where(idx)[0]

        # One spectral point as margin:
        idx_first = indices[0]
        idx_last = indices[-1] + 1

        if idx_first > 0:
            idx_first -= 1
        if idx_last < nwave:
            idx_last += 1

        if sign < 0.0:
            self.idx = np.flip(np.arange(idx_first, idx_last))
        else:
            self.idx = np.arange(idx_first, idx_last)

        self.wn = wn[self.idx]
        self.response = np.array(idx[self.idx], np.double)
        self.response /= np.trapz(self.response, self.wn)

        self.wl = 1.0 / (self.wn * pc.um)
        if input_is_wl:
            out_wave = self.wl
        else:
            out_wave = self.wn

        return out_wave, self.response


def constant_resolution_spectrum(wave_min, wave_max, resolution):
    """
    Compute a constant resolving-power sampling array.

    Parameters
    ----------
    wave_min: Float
        Lower spectral boundary.  This could be either a wavelength
        or a wavenumber. This is agnositc of units.
    wave_max: Float
        Upper spectral boundary.  This could be either a wavelength
        or a wavenumber. This is agnositc of units.
    resolution: Float
        The sampling resolving power: R = wave / delta_wave.

    Returns
    -------
    wave: 1D float array
        A spectrum array with the given resolving power.

    Examples
    --------
    >>> import numpy as np
    >>> import pyratbay.spectrum as ps

    >>> # A low-resolution wavelength sampling:
    >>> wl_min = 0.5
    >>> wl_max = 4.0
    >>> resolution = 5.5
    >>> wl = ps.constant_resolution_spectrum(wl_min, wl_max, resolution)
    >>> print(wl)
    [0.5        0.6        0.72       0.864      1.0368     1.24416
     1.492992   1.7915904  2.14990848 2.57989018 3.09586821 3.71504185]
    >>> # The actual resolution matches the input:
    >>> wl_mean = 0.5*(wl[1:]+wl[:-1])
    >>> print(wl_mean/np.ediff1d(wl))
    [5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5]
    """
    f = 0.5 / resolution
    g = (1.0+f) / (1.0-f)

    nwave = int(np.ceil(-np.log(wave_min/wave_max) / np.log(g)))
    wave = wave_min * g**np.arange(nwave)
    return wave


def bin_spectrum(bin_wl, wl, spectrum, half_widths=None):
    """
    Bin down a spectrum.

    Parameters
    ----------
    bin_wl: 1D float array
        Central wavelength (um) of the desired binned down spectra.
    wl: 1D float array
        Wavelength samples of the original spectrum.
    spectrum: 1D float array
        Spectral values to be binned down.
    half_widths: 1D float array
        The bin half widths (um).
        If None, assume that the bin edges are at the mid-points
        of the bin_wl array.

    Returns
    -------
    bin_spectrum: 1D float array
        The binned spectrum.

    Notes
    -----
    Probably bad things will happen if bin_wl has a similar
    or coarser resolution than wl.

    Examples
    --------
    >>> import pyratbay.spectrum as ps
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt

    >>> # Make a noisy high-resolution signal
    >>> wl = ps.constant_resolution_spectrum(1.0, 3.0, resolution=5000)
    >>> spectrum = np.sin(3.14*wl) + np.random.normal(1.5, 0.1, len(wl))
    >>> # Bin it down:
    >>> bin_wl = ps.constant_resolution_spectrum(1.0, 3.0, resolution=125)
    >>> bin_spectrum = ps.bin_spectrum(bin_wl, wl, spectrum)

    >>> # Compare original and binned signals
    >>> plt.figure(0)
    >>> plt.clf()
    >>> plt.plot(wl, spectrum, '.', ms=2, color='gray')
    >>> plt.plot(bin_wl, bin_spectrum, color='red')
    """
    if half_widths is None:
        half_widths = np.ediff1d(bin_wl, 0, 0)
        half_widths[0] = half_widths[1]
        half_widths[-1] = half_widths[-2]
        half_widths /= 2.0
    bands = [
        Tophat(wl0, half_width)
        for wl0, half_width in zip(bin_wl, half_widths)
    ]
    nbands = len(bands)
    band_flux = np.zeros(nbands)
    for i,band in enumerate(bands):
        band_wl, response = band(wl)
        band_flux[i] = np.trapz(spectrum[band.idx]*response, band.wn)
    return band_flux


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
        raise ValueError(
            "Resampling signal's wavenumber is not contained in specwn."
        )

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
        bandwn = [bandwn]

    bflux = []
    for btrans, wn in zip(bandtrans, bandwn):
        # Resample bandpasses into spectrum wavenumber sampling:
        resampled, wnidx = resample(btrans, wn, specwn, normalize=True)
        # Band-integrate spectrum:
        bflux.append(np.trapz(spectrum[wnidx]*resampled, specwn[wnidx]))

    return bflux

