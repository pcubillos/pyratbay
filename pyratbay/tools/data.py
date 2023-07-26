# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'parse_error_param',
    'Data',
]

import numpy as np
from . import tools as pt


def parse_error_param(var):
    r"""
    Parse error-scaling parameters. There are two options:
    - err_scale_name: Scale uncertainties as a multiplicative factor.
    - err_quad_name: Scale uncertainties by adding in quadrature.

    Parameters
    ----------
    var: String
        Parameter name. Must begin with either "err_scale_" or "err_quad_".

    Returns
    -------
    inst: String
        Instrument name, parsed from the end of var.
    texname: String
        Instrument name as a tex string, for plotting purposes.
    scaling: String
        What type of uncertainty scaling: "scale" or "quadrature".

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> # Valid options:
    >>> print(pt.parse_error_param('err_scale_WFC3'))
    ('WFC3', '$S_\\sigma^{\\rm WFC3}$', 'scale')
    >>> print(pt.parse_error_param('err_quad_NIRSpec_PRISM'))
    ('NIRSpec PRISM', '$\\sigma_{\\rm NIRSpec PRISM}$', 'quadrature')

    >>> # Not valid scaling throws error:
    >>> print(pt.parse_error_param('err_fudging_IRAC1'))
    ValueError: Invalid error scaling parameter 'fudging_IRAC1'. Valid options begin with: ['scale_', 'quad_']
    """
    error_scalings = {
        'err_scale_': r'$\log\ S^\sigma_{\rm INST}$',
        'err_quad_': r'$\log\ \sigma_{\rm INST}$',
    }

    inst = None
    for scaling, tex in error_scalings.items():
        if var.startswith(scaling):
            inst = var.replace(scaling, '', 1)
            inst = inst.replace('_', ' ')
            texname = tex.replace('INST', inst)
            mode = scaling.split('_')[1]
            break

    if inst is None:
        available_error_scalings = list(error_scalings.keys())
        raise ValueError(
            f"Invalid error scaling parameter '{var}'. "
            f"Valid options begin with: {available_error_scalings}"
        )

    mode = mode.replace('quad', 'quadrature')
    return inst, texname, mode


class Data():
    def __init__(
        self, data, uncert, band_names,
        offset_models=None, err_models=None,
    ):
        """
        Parameters
        ----------
        data: 1D float iterable
            The data values.
        uncert: 1D float iterable
            The uncertainty values.
        band_names: 1D string iterable
            Names for the uncertainties.
        offset_models: String or 1D iterable of strings
            List of data offset models, the strings must match
            a substring of at least one of the band_names, these specific
            data points will be affected by the respective offset model.
        err_models: String or 1D iterable of strings
            List of error inflation model names, must begin with either
            - "err_scale_" to scale uncertainties as a multiplicative factor
            - "err_quad_" to scale uncertainties by adding in quadrature
            followed by a string matching a substring of at least one of
            the band_names, these specific data points will be affected by
            the error scaling model.

        Examples
        --------
        >>> import pyratbay.tools as pt
        >>> import pyratbay.constants as pc
        >>> import pyratbay.io as io
        >>> import matplotlib.pyplot as plt

        >>> # Load a set of uncertainties obtained with JWST instruments:
        >>> obs_file = '/Users/pato/Documents/compendia/ERS/WASP39b/data/synthesis_v02/wasp39b_g395h_lrs.dat'
        >>> bands, depths, uncert = io.read_observations(obs_file)
        >>> band_names = [band.name for band in bands]
        >>> wl = [band.wl0 for band in bands]
        >>> print(set(band_names))
        {'nirspec_g395h_nrs1', 'nirspec_g395h_nrs2', 'miri_lrs'}

        >>> # Offsets for MIRI data only:
        >>> offsets = 'offset_miri'
        >>> obs = pt.Data(depths, uncert, band_names, offset_models=offsets)
        >>> offset_depths = obs.offset_data([400.0], 'ppm')

        >>> fig = plt.figure(0, (8,4))
        >>> plt.clf()
        >>> plt.semilogx(wl, obs.data/pc.ppm, 'o', color='darkorange')
        >>> plt.plot(wl, offset_depths/pc.ppm, '^' , color='blue', mfc='none')
        >>> plt.xlabel('Wavelength (um)')
        >>> plt.ylabel('Depths (ppm)')

        >>> # Multiplicative error scaling (increase by 1.5x):
        >>> # Target NRS1 uncertainties specifically of the NIRSpec instrument
        >>> obs = pt.Data(depths, uncert, band_names, err_models='err_scale_nrs1')
        >>> inflated_err = obs.scale_errors([0.3])

        >>> fig = plt.figure(1, (8,4))
        >>> plt.clf()
        >>> plt.semilogx(wl, obs.uncert/pc.ppm, 'o', color='xkcd:green')
        >>> plt.plot(wl, inflated_err/pc.ppm, '^' , color='blue', mfc='none')
        >>> plt.xlabel('Wavelength (um)')
        >>> plt.ylabel('Uncertainties (ppm)')

        >>> # Add-in-quadrature error scaling:
        >>> # Target all NIRSpec uncertainties (200ppm), for MIRI add 300ppm
        >>> err_models = ['err_quad_nirspec', 'err_quad_miri']
        >>> obs2 = pt.Data(depths, uncert, band_names, err_models=err_models)
        >>> inflated_err2 = obs2.scale_errors([2.3, 2.5], 'ppm')
        >>> plt.plot(wl, inflated_err2/pc.ppm, 's', color='orange', mfc='none')
        >>> plt.tight_layout()
        """
        # Input data
        self.data = np.copy(data)
        self.uncert = np.copy(uncert)
        if uncert is not None:
            self.ndata = len(self.data)

        # Data offset models
        if offset_models is None:
            offset_models = []
        elif isinstance(offset_models, str):
            offset_models = [offset_models]
        self.offset_models = offset_models
        self.n_offsets = len(self.offset_models)

        self.offset_indices = []
        self.offset_texnames = []
        for var in self.offset_models:
            inst = var.replace('offset_', '').replace('_',' ')
            texname = r'$\Delta$ {inst}'
            indices = np.array([inst in name for name in band_names])
            if np.sum(indices) == 0:
                raise ValueError(
                    f"Invalid instrumental offset parameter '{var}'. "
                    f"There is no instrument matching the name '{inst}'"
                )
            self.offset_indices.append(indices)
            self.offset_texnames.append(texname)

        offsets_per_data = np.sum(self.offset_indices, axis=0)
        if np.any(offsets_per_data > 1):
            raise ValueError(
                'Multiple instrumental offsets apply to a same data point'
            )

        # Error scaling models
        if err_models is None:
            err_models = []
        elif isinstance(err_models, str):
            err_models = [err_models]
        self.n_epars = len(err_models)

        self.err_inst = []
        self.err_texnames = []
        self.err_indices = []
        self.scaling_modes = []
        for var in err_models:
            inst, texname, scaling = parse_error_param(var)
            self.err_inst.append(inst)
            self.err_texnames.append(texname)
            self.scaling_modes.append(scaling)

            indices = np.array([inst in name for name in band_names])
            if np.sum(indices) == 0:
                raise ValueError(
                    f"Invalid retrieval parameter '{var}'. "
                    f"There is no instrument matching the name '{inst}'"
                )
            self.err_indices.append(indices)

        n_scales_per_data = np.sum(self.err_indices, axis=0)
        if np.any(n_scales_per_data > 1):
            raise ValueError(
                'Multiple uncertainty scaling apply to a same data point',
            )


    def offset_data(self, vals=None, units='none'):
        """
        Offset data values

        Parameters
        ----------
        vals: 1D float iterable
            Data offset values for each model.
        units: String
            The units of the scaling value. Options are:
            'none', 'percent', 'ppt', or 'ppm'.

        Returns
        -------
        data: 1D float array
            The offset data array.
        """
        data = np.copy(self.data)
        if vals is None or self.n_offsets==0:
            return data

        for j, val in enumerate(vals):
            indices = self.offset_indices[j]
            data[indices] += val * pt.u(units)

        return data


    def scale_errors(self, vals=None, units='none'):
        """
        Scale uncertainties according to input scaling values.

        Parameters
        ----------
        vals: 1D float iterable
            Uncertainty scaling value (in log10).
        units: String
            The units of the scaling value for quadrature errors.
            Options are: 'none', 'percent', 'ppt', or 'ppm'.

        Returns
        -------
        uncert: 1D float array
            The scaled uncertainty array.
        """
        uncert = np.copy(self.uncert)
        if vals is None or self.n_epars==0:
            return uncert

        for j, val in enumerate(vals):
            indices = self.err_indices[j]
            mode = self.scaling_modes[j]

            e_scaling = 10.0**val
            if mode == 'scale':
                uncert[indices] *= e_scaling
            elif mode == 'quadrature':
                e_scaling *= pt.u(units)
                uncert[indices] = np.sqrt(uncert[indices]**2.0 + e_scaling**2.0)

        return uncert

