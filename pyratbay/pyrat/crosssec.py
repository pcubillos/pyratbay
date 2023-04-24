# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import numpy as np

from .. import opacity as op
from .. import tools as pt


class CIA(object):
    def __init__(self, cia_files, wn, atm_species, log):
        """
        Read collision induced absorption (CIA) files.
        """
        self.files = None  # CS file names
        self.nfiles = 0
        self.models = []
        self.mol_indices = []
        self.tmin = 0.0   # Minimum temperature sampled by all CS files
        self.tmax = 1e6   # Maximum temperature sampled by all CS files
        self.ec = None  # extinction coefficient in cm-1 [nlayer, nwave]

        log.head('\nReading collision induced absorption files.')

        if cia_files is None:
            log.head('No CIA files to read.', indent=2)
            return

        self.nfiles = len(cia_files)
        self.files = cia_files
        self.nwave = wn

        for cia_file in cia_files:
            log.head(f"Read CIA file: '{cia_file}'.", indent=2)
            cia = op.Collision_Induced(cia_file, wn)
            self.models.append(cia)

            # Check that CIA species are in the atmospheric file:
            absent = np.setdiff1d(cia.species, atm_species)
            if len(absent) > 0:
                log.error(
                    f'These cia species {absent} are not listed in '
                    'the atmospheric file\n'
                )
            # Indices for CIA species in the atmosphere:
            mol_indices = [list(atm_species).index(mol) for mol in cia.species]
            self.mol_indices.append(mol_indices)

            # Screen output:
            species = '-'.join(cia.species)
            wn_min = cia.wn[cia._wn_lo_idx]
            wn_max = cia.wn[cia._wn_hi_idx-1]
            log.msg(
                f'CIA opacity for {species}:\n'
                f'Read {cia.nwave} wavenumber and {cia.ntemp} temperature samples.\n'
                f'Temperature sample limits: {cia.tmin:.1f}--{cia.tmax:.1f} K\n'
                f'Wavenumber sample limits: {wn_min:.1f}--{wn_max:.1f} cm-1',
                indent=4,
            )

        # Update temperature boundaries:
        self.tmin = np.amax([cia.tmin for cia in self.models])
        self.tmax = np.amin([cia.tmax for cia in self.models])

        log.head('CIA read done.')

    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write('CIA extinction information:')
        fw.write('Number of CIA files (nfiles): {:d}', self.nfiles)
        for i,cia in enumerate(self.models):
            fw.write("\nCIA file name (files[{}]): '{:s}'", i, self.files[i])
            fw.write('CIA species: {:s}', '-'.join(cia.species),
            )
            fw.write('Number of temperature samples: {}', cia.ntemp)
            fw.write('Number of wavenumber samples: {}', cia.nwave)
            fw.write(
                'Temperature array (temp, K):\n  {}',
                cia.temps,
                prec=1, lw=800, edge=50,
            )
            fw.write(
                'Wavenumber array (wavenumber, cm-1):\n  {}',
                cia.wn,
                fmt={'float':'{:.1f}'.format}, lw=80, edge=3,
            )
        fw.write(
            '\nMinimum and maximum temperatures (tmin, tmax) in K: '
            '[{:.1f}, {:.1f}]',
            self.tmin, self.tmax,
        )
        if self.ec is not None:
            fw.write(
                'CIA extinction coefficient (ec, cm-1):\n{}',
                self.ec,
                fmt={'float':'{:.2e}'.format},
            )
        return fw.text


    def calc_extinction_coefficient(self, temperature, densities):
        """
        Interpolate the CS absorption into the planetary model temperature.
        """
        ec = 0.0
        for i,cia in enumerate(self.models):
            density = densities[:,self.mol_indices[i]]
            ec += cia.calc_extinction_coefficient(temperature, density)

        if not np.isscalar(ec):
            self.ec = ec


    def get_ec(self, temperature, densities, layer):
        """
        Interpolate the CS absorption into the planetary model temperature.
        """
        ec = []
        label = []
        for i,cia in enumerate(self.models):
            imol = self.mol_indices[i]
            ext_coeff = cia.calc_extinction_coefficient(
                temperature[layer], densities[layer,imol],
            )
            ec.append(ext_coeff)
            label.append(cia.name)

        return ec, label
