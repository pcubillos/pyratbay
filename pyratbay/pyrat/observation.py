# Copyright (c) 2021-2026 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Observation',
]

import numpy as np

from .. import constants as pc
from .. import io as io
from .. import tools as pt


class Observation():
    def __init__(self, inputs, wn, log):
        # Transit or eclipse data point
        self.data = []
        self.uncert = []
        self.bands = []
        self.band_wl = []
        self.half_widths = []
        self.units = inputs.dunits
        self._dunits = pt.u(self.units)
        self._wn = wn

        # Low-resolution data (pass bands)
        if inputs.obsfile is not None:
            for obs_file in inputs.obsfile:
                obs_data = io.read_observations(obs_file)
                self.bands += obs_data[0]
                self.band_wl += list(obs_data[1])
                self.half_widths += list(obs_data[2])
                if len(obs_data) == 5:
                    self.data += list(obs_data[3])
                    self.uncert += list(obs_data[4])

        n_low_res = len(self.bands)

        # High-resolution data (sampled at nyquist frequency)
        self.inst_resolution = inputs.inst_resolution
        if inputs.obsfile_hires is not None:
            # TBD: Check spec.wn is at constant resolution
            if self.inst_resolution is None:
                raise ValueError(
                    'Undefined instrumental resolution is required for '
                    'convolution of high-resolution data'
                )
            for obs_file in inputs.obsfile_hires:
                obs_data = io.read_observations(obs_file)
                bands = obs_data[0]
                for band in bands:
                    band.half_width = band.wl0 / self.inst_resolution / 2.0
                self.bands += bands
                self.band_wl += list(obs_data[1])
                self.half_widths += [band.half_width for band in bands]
                if len(obs_data) == 5:
                    self.data += list(obs_data[3])
                    self.uncert += list(obs_data[4])

        # Low- and high-resolution flags
        self.is_low_res = np.arange(len(self.bands)) < n_low_res
        self.has_low_res = np.any(self.is_low_res)
        self.has_hi_res = np.any(~self.is_low_res)

        self.data = np.array(self.data)
        self.uncert = np.array(self.uncert)

        # Number of datapoints and filters
        self.ndata = len(self.data)
        self.nbands = len(self.bands)
        if self.ndata > 0 and self.ndata != self.nbands:
            log.error(
                f'Number of filter bands ({self.nbands}) does not '
                f'match the number of data points ({self.ndata})'
            )

        # Per-band variables
        self.band_wl = np.array(self.band_wl)
        self.half_widths = np.array(self.half_widths)
        self.bandwn = np.array([band.wn0 for band in self.bands])
        self.bandwn[~self.is_low_res] = 1.0/(self.band_wl[~self.is_low_res]*pc.um)
        self.bandflux = np.zeros(self.nbands, np.double)
        # Resample the filters into the planet wavenumber array:
        for band in self.bands:
            band.set_sampling(wn=wn)

        # Instrumental offsets and error-scaling parameters
        band_names = [band.name for band in self.bands]
        self.offset_inst = inputs.offset_inst
        self.offset_pars = inputs.offset_pars
        self.uncert_scaling = inputs.uncert_scaling
        self.uncert_pars = inputs.uncert_pars

        # This object contains the original transit/eclipse depth data
        # self.data and self.uncert can be modified and computed with
        # the methods of self.depth
        self.depth = pt.Data(
            self.data, self.uncert, band_names,
            self.offset_inst, self.uncert_scaling,
            self.units,
        )
        if len(self.offset_pars) > 0:
            self.data = self.depth.offset_data(self.offset_pars, self.units)
        if len(self.uncert_pars) > 0:
            # default values (zero scaling):
            for i in range(self.depth.n_epars):
                scaling = self.depth.scaling_modes[i]
                if scaling == 'scale' and self.uncert_pars[i] is None:
                    self.uncert_pars[i] = 0.0
                elif scaling == 'quadrature' and self.uncert_pars[i] is None:
                    self.uncert_pars[i] = -100.0
            self.uncert = self.depth.scale_errors(self.uncert_pars, self.units)
        # RV parameters will be initialized from pyrat.Retrieval()


    def parse_rv_parameters(self, pnames, log):
        """
        Setup RV model parameters.
        """
        band_names = [band.name for band in self.bands]
        rv_inst = []
        rv_texnames = []
        rv_mask = []
        for var in pnames:
            if var == 'rv_shift':
                rv_texnames.append('RV (km/s)')
                rv_mask.append(np.ones(self.nbands, bool))
                rv_inst.append(var)
            elif var.startswith('rv_'):
                inst = var[3:].replace('_', ' ')
                rv_texnames.append(fr'RV$_{{\rm {inst}}}$ (km/s)')
                mask = np.array([inst in name for name in band_names])
                rv_mask.append(mask)
                rv_inst.append(var)
                if np.sum(mask) == 0:
                    log.error(
                        f"Invalid RV parameter '{var}'. There is "
                        f"no instrument matching the name '{inst}'"
                    )
        rv_mask_count = np.sum(rv_mask, axis=0)
        if np.any(rv_mask_count>1):
            log.error('Multiple RV parameters apply to a same data point')
        # Add no-RV to mask
        rv_mask.append(rv_mask_count == 0)

        self.rv_inst = rv_inst
        self.rv_texnames = rv_texnames
        self.rv_mask = rv_mask
        self.rv_wn = np.zeros((len(self.rv_mask), len(self._wn)))
        self.rv_wn[-1] = self._wn


    def __str__(self):
        units = pt.u(self.units)
        fw = pt.Formatted_Write()
        fw.write('Observing information:')
        if self.ndata>0 or self.bands is not None:
            fw.write('Data/bandflux display units (units): {}', self.units)
            fw.write('Data/bandflux internal units: none')
        fw.write('Number of data points (ndata): {}', self.ndata)
        if self.ndata > 0:
            fw.write('        Data  Uncertainty   Wavenumber  Wavelength\n'
                     '     {:>7s}      {:>7s}         cm-1          um\n'
                     '      (data)     (uncert)     (bandwn)',
                     self.units, self.units)
            for data, uncert, wn in zip(self.data, self.uncert, self.bandwn):
                fw.write('  {:10.5f}   {:10.5f}    {:9.2f}  {:10.3f}',
                data/units, uncert/units, wn, 1.0/(wn*pc.um))

        fw.write('\nNumber of filter pass bands (nbands): {}', self.nbands)
        if self.nbands == 0:
            return fw.text
        fw.write(
            'Wavelength    Bandflux  Filter name\n'
            '        um     {:>7s}\n'
            ' (band_wl)  (bandflux)  (bands)',
            self.units,
        )
        for i,band in enumerate(self.bands):
            band_flux = self.bandflux[i] / units
            fw.write(
                '{:10.3f}  {:10.5f}  {:s}',
                band.wl0, band_flux, band.name,
            )
        return fw.text

