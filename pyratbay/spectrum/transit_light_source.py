# Copyright (c) 2021-2026 Cubillos & Blecic
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'TransitLightSource',
]

import os

import numpy as np
import scipy.interpolate as si

from .phoenix import read_phoenix
from .. import spectrum as ps


def _tophat_binning(bands, bin_wl, spectrum):
    """
    Bin spectra over a same set of tophat bands (this is called multiple
    times in TransitLightSource().__init__, avoiding constructing the
    same band objects every time.
    """
    nbands = len(bands)
    band_flux = np.zeros(nbands)
    for i,band in enumerate(bands):
        if band.idx is None:
            band_flux[i] = np.nan
        else:
            band_flux[i] = band(spectrum)

    # Patch gaps if needed
    mask = np.isnan(band_flux)
    band_flux[mask] = np.interp(bin_wl[mask], bin_wl[~mask], band_flux[~mask])

    return band_flux


class TransitLightSource():
    """
    A Transit-light-source model.
    """
    def __init__(self, sed_folder, teff, wl=None, wl_range=None):
        r"""
        Parameters
        ----------
        sed_folder: String
            Folder containing a list of PHOENIX New-Era models.
            Note that all files should have a same logg and metallicity,
            TransitLightSource() will sort according to their Teff.
        teff: Float
            Effective temperature (K) of the star.
        wl: 1D float array
            Wavelength array where to evaluate the TLS effect (microns).
            Strongly suggested to use a regular grid (as in example below).
        wl_range: Two-element float pair
            If provided, trim output wavelengh range to be between
            wl_range[0] <= wl <= wl_range[-1].

        Example
        -------
        >>> import pyratbay.spectrum as ps
        >>> import numpy as np
        >>> import matplotlib
        >>> import matplotlib.pyplot as plt
        >>>
        >>> # A folder containing a list of phoenix new-era models
        >>> # (see ps.fetch_phoenix() function)
        >>> sed_folder = 'phoenix/'
        >>> # Initialize TLS model
        >>> teff = 4800.0
        >>> wl = ps.constant_resolution_spectrum(0.3, 12.0, resolution=300.0)
        >>> tls = ps.TransitLightSource(sed_folder, teff, wl)

        >>> # Evaluate TLS effect for a range of star spot/faculae temperature
        >>> f_spot = 0.05
        >>> t_spots = teff + np.linspace(-1500, 1500, 11)
        >>> epsilon = [tls(t_spot, f_spot) for t_spot in t_spots]
        >>>
        >>> fig = plt.figure(0)
        >>> plt.clf()
        >>> fig.set_size_inches(8,4)
        >>> ax = plt.axes([0.1, 0.12, 0.89, 0.87])
        >>> for i,t in enumerate(t_spots):
        >>>     col = 'red' if t==teff else plt.cm.viridis(i/10)
        >>>     label = f'Tspot = {t_spots[i]:.0f} K'
        >>>     plt.plot(tls.wl, epsilon[i], color=col, label=label)
        >>> plt.legend(loc='lower left', fontsize=9, framealpha=0.75)
        >>> ax.set_xlim(0.45, 12)
        >>> ax.set_ylim(0.95, 1.055)
        >>> ax.set_xlabel(r"Wavelength ($\mathrm{\mu}$m)", fontsize=12)
        >>> ax.set_ylabel(r"TLS $\epsilon$", fontsize=12)
        >>> ax.set_xscale('log')
        >>> ax.tick_params(which='both', direction='in', labelsize=11)
        >>> ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        >>> ax.set_xticks([0.5, 0.7, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0])
        """
        self.teff = teff
        if not os.path.exists(sed_folder):
            msg = f'SED folder for TLS model not found: {repr(sed_folder)}'
            raise ValueError(msg)
        models = sorted([
            pfile
            for pfile in os.listdir(sed_folder)
            if pfile.startswith('lte')
        ])
        if len(models) == 0:
            msg = f'No PHOENIX models found in TLS folder: {repr(sed_folder)}'
            raise ValueError(msg)

        if wl_range is None:
            wl_range = -np.inf, np.inf

        # Setup wavelength array:
        file = f'{sed_folder}/{models[0]}'
        sed_wl, sed_flux = read_phoenix(file)

        has_custom_wl = wl is not None
        if has_custom_wl:
            wl_mask = (wl>=wl_range[0]) & (wl<=wl_range[1])
            wl = wl[wl_mask]
            # Make top-hat bands
            half_widths = 0.5 * np.ediff1d(wl, 0, 0)
            half_widths[0] = half_widths[1]
            half_widths[-1] = half_widths[-2]
            bands = [
                ps.Tophat(wl0, half_width, wl=sed_wl, ignore_gaps=True)
                for wl0, half_width in zip(wl, half_widths)
            ]
            self.wl = wl
        else:
            wl_mask = (sed_wl>=wl_range[0]) & (sed_wl<=wl_range[1])
            self.wl = sed_wl[wl_mask]

        nwave = len(self.wl)
        ntemps = len(models)
        self.temps = np.zeros(ntemps)
        fluxes = np.zeros((ntemps, nwave))
        for i in range(ntemps):
            file = f'{sed_folder}/{models[i]}'
            sed_wl, sed_flux = read_phoenix(file)
            if has_custom_wl:
                fluxes[i] = _tophat_binning(bands, wl, sed_flux)
            else:
                fluxes[i] = sed_flux[wl_mask]
            self.temps[i] = float(models[i][3:8])

        # Temperature interpolation
        i = np.searchsorted(self.temps, self.teff)
        star_flux = (
            (self.teff-self.temps[i-1]) * fluxes[i] +
            (self.temps[i] - self.teff) * fluxes[i-1]
        ) / (self.temps[i] - self.temps[i-1])

        # A temperature interpolator of flux ratios
        self.flux_ratios = si.interp1d(
            self.temps, fluxes/star_flux,
            axis=0,
            bounds_error=False,
            fill_value=1e100,
        )

    def __call__(self, t_spot, f_spot):
        """
        Evaluate transit-light-source contamination factor into
        a transit-depth spectrum (Equation 2 of Rackham+2018, ApJ, 853)
        """
        epsilon = 1.0 / (1.0 - f_spot * (1.0 - self.flux_ratios(t_spot)))
        return epsilon

    def __repr__(self):
        return "pyratbay.spectrum.TransitLightSource()"

    def __str__(self):
        with np.printoptions(threshold=100):
            return (
                f"Transit Light Source model\n"
                f"T_eff (K) = {self.teff}\n"
                f"temps (K) = {self.temps}\n"
                f"wl (micron) = {self.wl}\n"
            )

