# Copyright (c) 2021-2023 Patricio Cubillos
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
        self.units = inputs.dunits
        self._dunits = pt.u(self.units)

        if inputs.filters is None:
            self.filters = []
        else:
            self.filters = [
                ps.PassBand(filter_file)
                for filter_file in inputs.filters
            ]

        if inputs.obsfile is not None:
            # TBD: Throw error if filters already exist
            obs_data = io.read_observations(inputs.obsfile)
            if np.ndim(obs_data) == 2:
                # TBD: Throw error if data or uncert already exist
                self.filters, self.data, self.uncert = obs_data
            elif np.ndim(obs_data) == 1:
                self.filters = obs_data
        self.offset_instruments = inputs.offset_instruments

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

        self.nfilters = len(self.filters)
        if self.nfilters > 0 and self.ndata > 0 and self.ndata != self.nfilters:
            log.error(
                f'Number of filter bands ({self.nfilters}) does not '
                f'match the number of data points ({self.ndata})'
            )

        # Resample the filters into the planet wavenumber array:
        for band in self.filters:
            band(wn=wn)
        # Per-band variables:
        self.bandwn = np.array([band.wn0 for band in self.filters])
        self.bandflux = np.zeros(self.nfilters, np.double)

        # Parse instrumental offsets if any
        if self.offset_instruments is None:
            self.offset_instruments = []

        n_offsets = len(self.offset_instruments)
        self.offset_indices = []
        self.offset_pars = np.zeros(n_offsets)

        band_names = [band.name for band in self.filters]
        for var in self.offset_instruments:
            inst = var.replace('offset_', '', 1)
            flags = [inst in name for name in band_names]
            self.offset_indices.append(flags)
            if np.sum(flags) == 0:
                log.error(
                    f"Invalid retrieval parameter '{var}'. "
                    f"There is no instrument matching the name '{inst}'"
                )
        offsets = np.sum(self.offset_indices, axis=0)
        if np.any(offsets > 1):
            log.error('Multiple instrumental offsets apply to a same bandpass')


    def __str__(self):
        units = pt.u(self.units)
        fw = pt.Formatted_Write()
        fw.write('Observing information:')
        if self.data is not None or self.filters is not None:
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

        fw.write('\nNumber of filter pass bands (nfilters): {}', self.nfilters)
        if self.nfilters == 0:
            return fw.text
        fw.write(
            'Wavenumber  Wavelength    Bandflux  Filter name\n'
            '      cm-1          um     {:>7s}\n'
            '  (bandwn)              (bandflux)  (filters)',
            self.units,
        )
        for i,band in enumerate(self.filters):
            band_flux = self.bandflux[i] / units
            fw.write(
                ' {:9.2f}  {:10.3f}  {:10.5f}  {:s}',
                band.wn0, 1.0/(band.wn0*pc.um), band_flux, band.name,
            )
        return fw.text

