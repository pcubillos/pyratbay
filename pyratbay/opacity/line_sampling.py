# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Line_Sample',
]

import numpy as np
import mc3

from .. import io as io
from .. import tools as pt
from .. import constants as pc
from ..lib import _extcoeff as ec


class Line_Sample():
    """Line-by-line sampled opacities"""
    def __init__(
        self, cs_files, *, pressure=None, temperature=None,
        min_wl=None, max_wl=None, min_wn=None, max_wn=None,
        wn_thinning=1,
        log=None,
    ):
        """
        Read line-sampled cross-section table(s), with units of cm2 molec-1.

        Parameters
        ----------
        cs_files: String or iterable of strings
            Line-sampled cross section file(s) to read.
        pressure: 1D floar array
            Pressure profile where to resample the opacities (bar).
            If None, use the tabulated pressure array from the opacities.
            If not None, it is allowed to extrapolate to lower pressures
            but not to higher pressures.
        min_wl: 1D float ndarray
            Minimum wavelength value to extract from line-sample files (um)
            (only one of min_wl or max_wn should be provided).
        max_wl: 1D float ndarray
            Maximum wavelength value to extract from line-sample files (um)
            (only one of min_wn or max_wl should be provided).
        min_wn: 1D float ndarray
            Minimum wavenumber value to extract from line-sample files (cm-1)
            (only one of min_wn or max_wl should be provided).
        max_wn: 1D float ndarray
            Maximum wavenumber value to extract from line-sample files (cm-1)
            (only one of min_wl or max_wn should be provided).
        wn_thinning: Integer
            Thinning factor to take every n-th value of the wavenumber array

        Examples
        --------
        >>> import pyratbay.opacity as op
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt

        >>> # Generate these files from
        >>> # pyratbay.rtfd.io/en/latest/cookbooks/line_list_hitran.html
        >>> cs_files = [
        >>>     'cross_section_R020K_0150-3000K_0.5-5.0um_hitran_H2O.npz',
        >>>     'cross_section_R020K_0150-3000K_0.5-5.0um_hitemp_CO.npz',
        >>> ]
        >>> # Line-sample opacities at wavenumbers > 1.0 cm-1:
        >>> ls = op.Line_Sample(cs_files, min_wl=1.0)

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
        >>> plt.plot(wl, ec[0,35], color='blue', lw=1.0)
        >>> plt.plot(wl, ec[1,35], alpha=0.6, color='orange', lw=1.0)
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
        use_shared_memory = pt.get_mpi_size() > 1
        if use_shared_memory:
            from mpi4py import MPI
            rank = pt.get_mpi_rank()

        # Get dimensions first:
        # Species, temperature (K), pressure (bar), and wavenumber (cm-1)
        species, temp, press, wn = io.read_opacity(
            self.cs_files[0], extract='arrays',
        )

        if temperature is None:
            self.temp = temp
        else:
            self.temp = temperature
        self.ntemp = len(self.temp)

        if pressure is None:
            self.press = press
        else:
            self.press = pressure
        self.nlayers = len(self.press)

        if min_wn is not None and max_wl is not None:
            log.error('Either define min_wn or max_wl, not both')
        if max_wn is not None and min_wl is not None:
            log.error('Either define min_wl or max_wn, not both')

        if min_wn is None:
            min_wn = 0.0 if max_wl is None else 1.0/(max_wl*pc.um)

        if max_wn is None:
            max_wn = np.inf if min_wl is None else 1.0/(min_wl*pc.um)

        self.wn = wn[(wn >= min_wn) & (wn <= max_wn)][::wn_thinning]
        self.nwave = len(self.wn)

        self.species = []
        species_per_file = []
        wn_masks = []
        for cs_file in self.cs_files:
            #cs_file = os.path.basename(cs_file)
            species, temp, press, wn = io.read_opacity(cs_file,extract='arrays')
            wn_mask = (wn >= min_wn) & (wn <= max_wn)
            wn = wn[wn_mask][::wn_thinning]
            wn_masks.append(wn_mask)
            species_per_file.append(list(species))
            ntemp = len(temp)
            nwave = len(wn)

            wave_mismatch = (
               nwave != self.nwave or
               np.any(np.abs(1.0-wn/self.wn) > 0.01)
            )
            if wave_mismatch:
                log.error(
                    f"Wavenumber array of cross-section file '{cs_file}' "
                    "does not match with previous arrays"
                )

            check_pressure_boundaries(self.press, press)
            # TBD: Implement or raise warning
            #check_temperature_boundaries(self.temp, temp)

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
                self.cs_table[idx] += pt.interpolate_opacity(
                    cs_file, self.temp, self.press, mask, wn_thinning,
                )
        else:
            itemsize = MPI.DOUBLE.Get_size()
            if rank == 0:
                nbytes = np.prod(cs_shape) * itemsize
            else:
                nbytes = 0
            # on rank 0, create the shared block
            # else get a handle to it
            win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=MPI.COMM_WORLD)

            buf, itemsize = win.Shared_query(0)
            assert itemsize == MPI.DOUBLE.Get_size()
            self.cs_table = np.ndarray(buffer=buf, dtype='d', shape=cs_shape)
            if rank == 0:
                for i,cs_file in enumerate(self.cs_files):
                    idx = spec_indices[i]
                    mask = wn_masks[i]
                    self.cs_table[idx] += pt.interpolate_opacity(
                        cs_file, self.temp, self.press, mask, wn_thinning,
                    )
            pt.mpi_barrier()

        # Set tabulated temperature extrema:
        self.tmin = np.amin(self.temp)
        self.tmax = np.amax(self.temp)

    def get_wl(self, units='um'):
        """
        Get the wavelength array in the desired units.

        Parameters
        ----------
        units: String
            Select one from: ['A', 'nm', 'um', 'mm', 'cm', 'm', 'km']

        Returns
        -------
        wl: 1D float ndarray
            The wavelength array from the line-sample grid.
        """
        valid_units = 'A nm um mm cm m km'.split()
        if units not in valid_units:
            raise ValueError(
                f"Invalid wavelength units '{units}', "
                f"select one from {valid_units}"
            )
        return 1.0/(self.wn * pt.u(units))


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

        density = np.ones((self.nlayers, self.nspec))
        interp_ec(
            cross_section,
            self.cs_table, self.temp,
            temperature, density,
            layer1, layer2,
        )

        if np.isscalar(layer):
            if per_mol:
                cross_section = cross_section[:,layer]
            else:
                cross_section = cross_section[layer]
        return cross_section


    def calc_extinction_coefficient(
        self, temperature, density, layer=None, per_mol=False,
    ):
        """
        Calculate extinction-coefficient spectra (cm-1) for temperature
        and number-density profiles.

        Parameters
        ----------
        temperature: 1D float array
            Temperature array (Kelvin) at which to interpolate the
            cross section (must match nlayers size).
        density: 1D float array
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
            temperature, density,
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
        fw.write(
            'Minimum and maximum pressures in bar: '
            f'[{np.amin(self.press):.3e}, {np.amax(self.press):.3e}]'
        )
        fw.write(
            'Minimum and maximum wavelengths in um: '
            f'[{np.amin(self.get_wl()):.3f}, {np.amax(self.get_wl()):.3f}]'
        )
        fw.write(f'\nLine-sample species (species): {self.species}')
        fw.write(f'Temperature array (temps, K):\n{self.temp}')
        with np.printoptions(precision=3):
            fw.write(f'Pressure layers (pressure, bar):\n{self.press}')
            fw.write(f'Wavenumber array (wn, cm-1):\n  {self.wn}')

        fw.write(
            'The tabulated cross sections (cs_table, cm2 molecule-1) '
            'are an array\nof dimensions [nspec,ntemp,nlayers,nwave] and '
            f'shape {self.cs_table.shape}'
        )
        return fw.text


def check_pressure_boundaries(press, tabulated_press):
    """
    Check that values in pressure profile are no larger than the maximum
    tabulated pressure (to relative precission of 1e-3).
    """
    pmax = np.amax(press)
    pmax_table = np.amax(tabulated_press)
    if pmax/pmax_table - 1 > 1e-3:
        raise ValueError(
            'Pressure profile extends beyond the maximum tabulated pressure'
        )
