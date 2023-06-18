# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Line_Sample',
]

import importlib

import numpy as np
import mc3

from .. import io as io
from .. import tools as pt
from ..lib import _extcoeff as ec


class Line_Sample():
    """Line-by-line sampled opacities"""
    def __init__(self, cs_files, min_wn=-np.inf, max_wn=np.inf, log=None):
        """
        Read line-sampled cross-section table(s), with units of cm2 molec-1.

        Parameters
        ----------
        cs_files: String
            Line-sampled cross section files.
        min_wn: 1D float ndarray
            Minimum wavenumber value to extract from line-sample files (cm-1)
        max_wn: 1D float ndarray
            Maximum wavenumber value to extract from line-sample files (cm-1)

        Examples
        --------
        >>> import pyratbay.opacity as op
        >>> # TDB: how to get these files
        >>> cs_files = [
        >>>     'cs_table_H2O_150-3000K_0.5-5.0um.npz',
        >>>     'cs_table_CO_150-3000K_0.5-5.0um.npz',
        >>> ]
        >>> # Line-sample opacities at wavenumbers > 1.0 cm-1:
        >>> ls = op.Line_Sample(cs_files, max_wn=1e4/1.0)

        >>> # Total or per-mol cross sections (cm2 molec-1):
        >>> temp = np.tile(1400.0, ls.nlayers)
        >>> cs = ls.calc_cross_section(temp)
        >>> cs_per_mol = ls.calc_cross_section(temp, per_mol=True)

        >>> # Cross section per species (cm2 molec-1):
        >>> ec = ls.calc_cross_section(temp, per_mol=True)

        >>> # Take a look
        >>> wl = 1e4/ls.wn
        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.plot(wl, ec_per_mol[0,50], color='blue', lw=1.0)
        >>> plt.plot(wl, cs_per_mol[1,50], alpha=0.6, color='orange', lw=1.0)
        >>> plt.yscale('log')
        """
        self.name = 'line sampling'
        self.nspec = 0
        if log is None:
            log = mc3.utils.Log(width=80)

        if isinstance(cs_files, str):
            cs_files = [cs_files]
        self.cs_files = cs_files

        missing_files = [
            cs_file
            for cs_file in self.cs_files
            if pt.isfile(cs_file) == 0
        ]
        if len(missing_files) > 0:
            log.error(f'Missing opacity files: {missing_files}')

        # Load opacities into shared memory if and only if possible and needed:
        use_shared_memory = False
        mpi_exists = importlib.util.find_spec('mpi4py') is not None
        if mpi_exists:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            use_shared_memory = comm.size > 1

        # Get dimensions first:
        # Species, temperature (K), pressure (barye), and wavenumber (cm-1)
        species, self.temp, self.press, wn = io.read_opacity(
            self.cs_files[0], extract='arrays',
        )
        self.ntemp = len(self.temp)
        self.nlayers = len(self.press)
        self.wn = wn[(wn >= min_wn) & (wn <= max_wn)]
        self.nwave = len(self.wn)

        self.species = []
        species_per_file = []
        wn_masks = []
        for cs_file in self.cs_files:
            #cs_file = os.path.basename(cs_file)
            species, temp, press, wn = io.read_opacity(cs_file,extract='arrays')
            wn_mask = (wn >= min_wn) & (wn <= max_wn)
            wn = wn[wn_mask]
            wn_masks.append(wn_mask)
            species_per_file.append(list(species))
            ntemp = len(temp)
            nlayers = len(press)
            nwave = len(wn)

            # TBD: need to reconsider if applying p-resampling
            shape_mismatch = (
                ntemp != self.ntemp or
                nlayers != self.nlayers or
                nwave != self.nwave
            )
            if shape_mismatch:
                log.error(
                    f"Shape of the cross-section file '{cs_file}' "
                    "does not match with previous file shapes."
                )

            value_mismatch = [
                np.any(np.abs(1.0-temp/self.temp) > 0.01),
                np.any(np.abs(1.0-press/self.press) > 0.01),
                np.any(np.abs(1.0-wn/self.wn) > 0.01),
            ]
            if np.any(value_mismatch):
                vals = np.array(['temperature', 'pressure', 'wavenumber'])
                mismatch = ', '.join(vals[value_mismatch])
                log.error(
                    f"Tabulated {mismatch} values in file '{cs_file}' "
                    "do not match with previous arrays"
                )
            # Add new species:
            self.species += [
                spec
                for spec in species
                if spec not in self.species
            ]

        spec_indices = []
        for species in species_per_file:
            spec_indices.append([
                self.species.index(spec) for spec in species
            ])
        self.species = np.array(self.species)
        self.nspec = len(self.species)

        # Cross-sections table (cm2 molecule-1):
        cs_shape = (self.nspec, self.ntemp, self.nlayers, self.nwave)
        if not use_shared_memory:
            self.cs_table = np.zeros(cs_shape)
            for i,cs_file in enumerate(self.cs_files):
                idx = spec_indices[i]
                mask = wn_masks[i]
                self.cs_table[idx] += \
                    io.read_opacity(cs_file, extract='opacity')[:,:,:,mask]
        else:
            itemsize = MPI.DOUBLE.Get_size()
            if comm.rank == 0:
                nbytes = np.prod(cs_shape) * itemsize
            else:
                nbytes = 0
            # on rank 0, create the shared block
            # else get a handle to it
            win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)

            buf, itemsize = win.Shared_query(0)
            assert itemsize == MPI.DOUBLE.Get_size()
            self.cs_table = np.ndarray(buffer=buf, dtype='d', shape=cs_shape)
            if comm.rank == 0:
                for i,cs_file in enumerate(self.cs_files):
                    idx = spec_indices[i]
                    mask = wn_masks[i]
                    self.cs_table[idx] += \
                        io.read_opacity(cs_file, extract='opacity')[:,:,:,mask]
            comm.Barrier()

        # Set tabulated temperature extrema:
        self.tmin = np.amin(self.temp)
        self.tmax = np.amax(self.temp)


    def calc_cross_section(self, temperature, layer=None, per_mol=False):
        """
        Calculate cross-section spectra (cm2 molec-1) over temperature
        profiles by interpolating from tabulated values.

        Parameters
        ----------
        temperature: 1D float array
            Temperature array (Kelvin) at which to interpolate the
            cross section (must match nlayers size).
        layer: Integer
            If not None, compute cross-sections only at selected layer.
            In this case, the output array dimension is reduced by 1.
        per_mol: bool
            If True, compute cross sections individually per species.
            If False, co-add cross section contributions from all species.

        Returns
        -------
        cross_section: 1D, 2D, or 3D float array
            Cross section spectra (cm2 molec-1)
            Output array has shape [nspec, nlayers, nwave].
            If per_mol is False, nspec dimension is removed.
            If layer is not None, nlayers dimension is removed.
        """
        if np.amax(temperature) > self.tmax or np.amin(temperature) < self.tmin:
            raise ValueError('Temperatures are out of line-sample bounds')
        if layer is None:
            layer1 = 0
            layer2 = self.nlayers
        elif np.isscalar(layer):
            layer1 = layer
            layer2 = layer+1
        elif len(layer) == 2:
            layer1 = layer[0]
            layer2 = layer[1]
        else:
            raise ValueError('Invalid layer input')

        if per_mol:
            cross_section = np.zeros((self.nspec, self.nlayers, self.nwave))
            interp_ec = ec.interp_ec_per_mol
        else:
            cross_section = np.zeros((self.nlayers, self.nwave))
            interp_ec = ec.interp_ec

        densities = np.ones((self.nlayers, self.nspec))
        interp_ec(
            cross_section,
            self.cs_table, self.temp,
            temperature, densities,
            layer1, layer2,
        )

        if np.isscalar(layer):
            if per_mol:
                cross_section = cross_section[:,layer]
            else:
                cross_section = cross_section[layer]
        return cross_section


    def calc_extinction_coefficient(
        self, temperature, densities, layer=None, per_mol=False,
    ):
        """
        Calculate extinction-coefficient spectra (cm-1) for temperature
        and number-density profiles.

        Parameters
        ----------
        temperature: 1D float array
            Temperature array (Kelvin) at which to interpolate the
            cross section (must match nlayers size).
        densities: 1D float array
            Number-density profiles (molec cm-3) at each layer.
            Array has shape [nspec,nlayers].
        layer: Integer
            If not None, compute extinction coefficient only at selected layer.
            In this case, the output array dimension is reduced by 1.
        per_mol: bool
            If True, compute extinction coefficients individually per species.
            If False, co-add extinction contributions from all species.

        Returns
        -------
        extinction: 1D, 2D, or 3D float array
            Extinction-coefficient spectra (cm-1)
            Output array has shape [nspec, nlayers, nwave].
            If per_mol is False, nspec dimension is removed.
            If layer is not None, nlayers dimension is removed.
        """
        if np.amax(temperature) > self.tmax or np.amin(temperature) < self.tmin:
            raise ValueError('Temperatures are out of line-sample bounds')
        if layer is None:
            layer1 = 0
            layer2 = self.nlayers
        elif np.isscalar(layer):
            layer1 = layer
            layer2 = layer+1
        elif len(layer) == 2:
            layer1 = layer[0]
            layer2 = layer[1]
        else:
            raise ValueError('Invalid layer input')

        if per_mol:
            extinction = np.zeros((self.nspec, self.nlayers, self.nwave))
            interp_ec = ec.interp_ec_per_mol
        else:
            extinction = np.zeros((self.nlayers, self.nwave))
            interp_ec = ec.interp_ec

        interp_ec(
            extinction,
            self.cs_table, self.temp,
            temperature, densities,
            layer1, layer2,
        )

        if np.isscalar(layer):
            if per_mol:
                extinction = extinction[:,layer]
            else:
                extinction = extinction[layer]
        return extinction


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write(
            f"Line-sampling cross-section files (cs_files):\n{self.cs_files}"
        )
        fw.write(f'Number of species (nspec): {self.nspec}')
        fw.write(f'Number of temperature samples (ntemp): {self.ntemp}')
        fw.write(f'Number of pressure layers (nlayers): {self.nlayers}')
        fw.write(f'Number of wavenumber samples (nwave): {self.nwave}')
        fw.write(
            '\nMinimum and maximum temperatures (tmin, tmax) in K: '
            f'[{self.tmin:.1f}, {self.tmax:.1f}]'
        )
        fw.write(f'\nLine-sample species (species): {self.species}')
        fw.write(f'Temperature array (temps, K):\n{self.temp}')
        fw.write(f'Pressure layers (pressure, barye):\n{self.press}')
        fw.write(f'Wavenumber array (wn, cm-1):\n{self.wn}')

        fw.write(
            'Tabulated cross sections (cs_table, cm2 molecule-1):\n{}',
            self.cs_table,
            edge=1,
        )
        return fw.text

