# Copyright (c) 2021-2026 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Observation',
]

import numpy as np

from .. import constants as pc
from .. import io as io
from .. import spectrum as ps
from .. import tools as pt


class Observation():
    def __init__(self, inputs, wn, log):
        # Transit or eclipse data point
        self.data = inputs.data
        self.uncert = inputs.uncert
        self.data_hires = None
        self.uncert_hires = None
        self.band_wl = None
        self.units = inputs.dunits
        self._dunits = pt.u(self.units)

        if inputs.filters is None:
            self.bands = []
        else:
            self.bands = [
                ps.PassBand(filter_file)
                for filter_file in inputs.filters
            ]

        # Low-resolution data (pass bands)
        if inputs.obsfile is not None:
            # TBD: Throw error if filters already exist
            obs_data = io.read_observations(inputs.obsfile)
            if len(obs_data) == 5:
                # TBD: Throw error if data or uncert already exist
                self.data, self.uncert = obs_data[3:]
            self.bands, self.band_wl, self.half_widths = obs_data[0:3]

        # Number of datapoints and filters:
        self.ndata = 0
        if self.data is None and self.uncert is not None:
            log.error("Undefined transit/eclipse data (data)")

        if self.data is not None:
            self.ndata = len(self.data)

        if self.uncert is not None:
            n_uncert = len(self.uncert)
            if self.ndata != n_uncert:
                log.error(
                    f'Number of data uncertainty values ({n_uncert}) does not '
                    f'match the number of data points ({self.ndata})'
                )

        self.nbands = len(self.bands)
        if self.nbands > 0 and self.ndata > 0 and self.ndata != self.nbands:
            log.error(
                f'Number of filter bands ({self.nbands}) does not '
                f'match the number of data points ({self.ndata})'
            )

        # Resample the filters into the planet wavenumber array:
        for band in self.bands:
            band.set_sampling(wn=wn)
        # Per-band variables:
        self.bandwn = np.array([band.wn0 for band in self.bands])
        self.bandflux = np.zeros(self.nbands, np.double)


        # High-resolution data (sampled at nyquist frequency)
        self.inst_resolution = inputs.inst_resolution
        self.bands_hires = []
        if inputs.obsfile_hires is not None:
            # TBD: Check spec.wn is at constant resolution
            if self.inst_resolution is None:
                raise ValueError(
                    'Undefined instrumental resolution is required for '
                    'convolution of high-resolution data'
                )
            obs_data = io.read_observations(inputs.obsfile_hires)
            if len(obs_data) == 5:
                self.data_hires, self.uncert_hires = obs_data[3:]
            self.bands_hires, self.band_wl_hires, hw = obs_data[0:3]
            self.wn_hires = 1.0 / (self.band_wl_hires*pc.um)
            for band in self.bands_hires:
                band.half_width = band.wl0 / self.inst_resolution / 2.0
                band.set_sampling(wn=wn)

        if self.data_hires is not None:
            self.ndata_hires = len(self.data_hires)
        else:
            self.ndata_hires = 0

        if self.uncert_hires is not None:
            self.n_uncert_hires = len(self.uncert_hires)

        self.nbands_hires = len(self.bands_hires)
        self.bandflux_hires = np.zeros(self.nbands_hires, np.double)

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


    def __str__(self):
        units = pt.u(self.units)
        fw = pt.Formatted_Write()
        fw.write('Observing information:')
        if self.data is not None or self.bands is not None:
            fw.write('Data/bandflux display units (units): {}', self.units)
            fw.write('Data/bandflux internal units: none')
        fw.write('Number of data points (ndata): {}', self.ndata)
        if self.data is not None:
            fw.write('        Data  Uncertainty   Wavenumber  Wavelength\n'
                     '     {:>7s}      {:>7s}         cm-1          um\n'
                     '      (data)     (uncert)     (bandwn)',
                     self.units, self.units)
            for data, uncert, wn in zip(self.data, self.uncert, self.bandwn):
                fw.write('  {:10.5f}   {:10.5f}    {:9.2f}  {:10.3f}',
                data/units, uncert/units, wn, 1.0/(wn*pc.um))
        # TBD: add hires data

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

