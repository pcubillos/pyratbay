# Copyright (c) 2021-2025 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'optical_depth',
]

import numpy as np

from .. import atmosphere as pa
from .. import constants as pc
from ..lib import _trapezoid as t
from ..lib import cutils as cu


def optical_depth(
        rt_path, extinction, radius=None,
        itop=0, ibottom=None, maxdepth=np.inf,
        extinction_cloudy=None, raypath=None,
    ):
    """
    Calculate the optical depth.

    Parameters
    ----------
    rt_path: String
        Radiative-transfer path.
    extinction: 2D float array
        Extinction coefficient (cm-1) at each layer and wavelength channel.
    radius: 1D float array
        Radius (altitude) array of the atmospheric layers (cm), from top
        to bottom.  Used to compute the raypath.
    itop: Integer
        Index at top of the atmosphere, opacity at layer above this
        will be ignored.
    ibottom: Integer
        Index at the bottom of the atmosphere. The optical depth will
        be calculated down to the layer pointed by this index.
    maxdepth: Float
        Maximum optical depth to compute.  The optical depth calculation
        will stop when maxdepth is reached (checked independently in each
        wavelength channel).
    extinction_couldy: 2d float array
        If provided, signals that the optical depth should compute
        separately the depth of a cloudy and a clear atmosphere.
        This array contains the extinction coefficient of the 'cloud'
        absorbers.
    raypath:
        The ray paths, which can be directly provided instead of radius.

    Returns
    -------
    raypath:
        Path followed by the rays.
    depth: 2D float array
        The optical depth at each layer and wavelength channel.
    ideep: 1D integer array
        Indices of each wavelenght channel of the layer where the optical
        depth surpassed maxdepth.
    depth_clear: 2D float array
        If extinction_couldy is provided, the optical depth at each
        layer and wavelength channel for a cloud-less atmosphere.
    ideep_clear: 1D integer array
        Same as ideep but for depth_clear.
    """
    is_patchy = extinction_cloudy is not None
    nlayers, nwave = extinction.shape

    if rt_path in pc.transmission_rt:
        is_transit = True
    elif rt_path in pc.emission_rt + pc.eclipse_rt:
        is_transit = False
    else:
        raise ValueError('Invalid radiative-transfer path')

    if ibottom is None:
        ibottom = nlayers

    # Calculate the ray path:
    if radius is None and raypath is None:
        raise ValueError(
            'Need to provide either radius or raypath to compute the '
            'path for the optical-depth calculation'
        )

    if raypath is None and is_transit:
        raypath = pa.transit_path(radius, itop)
    if raypath is None and not is_transit:
        raypath = -cu.ediff(radius)

    # Evaluate the extinction coefficient at each layer:
    ec = np.copy(extinction)
    depth = np.zeros((nlayers, nwave))
    # If patchy, compute clear and cloudy spectra separately:
    if is_patchy:
        ec_clear = np.copy(extinction)
        ec[itop:] += extinction_cloudy[itop:]
        depth_clear = np.zeros((nlayers, nwave))
    else:
        depth_clear = None
        ideep_clear = None

    # Calculate the optical depth for each wavenumber:
    if is_transit:
        ideep = np.array(np.tile(-1, nwave), dtype=np.intc)
        r = itop
        # Optical depth at each level (tau = 2.0*integral e*ds):
        for r in range(itop, ibottom):
            depth[r] = t.optdepth(
                ec[itop:r+1], raypath[r], maxdepth, ideep, r,
            )
        ideep[ideep<0] = r

        if is_patchy:
            ideep_clear = np.array(np.tile(-1,nwave), dtype=np.intc)
            for r in range(itop, nlayers):
                depth_clear[r] = t.optdepth(
                    ec_clear[itop:r+1], raypath[r], maxdepth,
                    ideep_clear, r,
                )
            ideep_clear[ideep_clear<0] = r

    else:
        ideep = np.tile(nlayers-1, nwave)
        if 'two_stream' in rt_path:
            maxdepth = np.inf
        t.plane_parallel_optical_depth(
            depth, ideep,
            ec, raypath, maxdepth, itop, ibottom,
        )
        if is_patchy:
            ideep_clear = np.tile(nlayers-1, nwave)
            t.plane_parallel_optical_depth(
                depth_clear, ideep_clear,
                ec_clear, raypath, maxdepth, itop, nlayers,
            )

    return (
        raypath,
        depth,
        ideep,
        depth_clear,
        ideep_clear,
    )


