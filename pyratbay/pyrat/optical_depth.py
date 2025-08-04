# Copyright (c) 2021-2025 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import numpy as np

from .. import atmosphere as pa
from .. import constants as pc
from ..lib import _trapezoid as t
from ..lib import cutils as cu


def optical_depth(pyrat):
    """
    Calculate the optical depth.
    """
    od = pyrat.od
    nwave = pyrat.spec.nwave
    nlayers = pyrat.atm.nlayers
    rtop = pyrat.atm.rtop
    pyrat.log.head('\nBegin optical-depth calculation.')

    # Evaluate the extinction coefficient at each layer:
    od.ec = np.empty((nlayers, nwave))
    od.depth = np.zeros((nlayers, nwave))
    if pyrat.opacity.is_patchy:
        od.ec_clear = np.empty((nlayers, nwave))
        od.depth_clear = np.zeros((nlayers, nwave))

    # Calculate the ray path:
    if pyrat.od.rt_path in pc.emission_rt + pc.eclipse_rt:
        pyrat.od.raypath = -cu.ediff(pyrat.atm.radius)
    elif pyrat.od.rt_path in pc.transmission_rt:
        pyrat.od.raypath = pa.transit_path(pyrat.atm.radius, pyrat.atm.rtop)

    # Sum all contributions to the extinction:
    od.ec[rtop:] = pyrat.opacity.ec[rtop:]

    # If patchy, compute clear and cloudy spectra separately:
    if pyrat.opacity.is_patchy:
        od.ec_clear[rtop:] = np.copy(od.ec[rtop:])
        od.ec[rtop:] += pyrat.opacity.ec_cloud[rtop:]


    rbottom = nlayers
    for model in pyrat.opacity.models:
        if model.name == 'deck':
            rbottom = model.itop + 1
            break
    # Calculate the optical depth for each wavenumber:
    if od.rt_path in pc.emission_rt + pc.eclipse_rt:
        od.ideep = np.tile(nlayers-1, nwave)
        maxdepth = np.inf if 'two_stream' in od.rt_path else od.maxdepth
        t.plane_parallel_optical_depth(
            od.depth, od.ideep,
            od.ec, od.raypath, maxdepth, rtop, rbottom,
        )
        if pyrat.opacity.is_patchy:
            od.ideep_clear = np.tile(nlayers-1, nwave)
            t.plane_parallel_optical_depth(
                od.depth_clear, od.ideep_clear,
                od.ec_clear, od.raypath, maxdepth, rtop, nlayers,
            )

    elif od.rt_path in pc.transmission_rt:
        od.ideep = ideep = np.array(np.tile(-1, nwave), dtype=np.intc)
        r = rtop
        # Optical depth at each level (tau = 2.0*integral e*ds):
        for r in range(rtop, rbottom):
            od.depth[r] = t.optdepth(
                od.ec[rtop:r+1], od.raypath[r], od.maxdepth, ideep, r,
            )
        ideep[ideep<0] = r

        if pyrat.opacity.is_patchy:
            od.ideep_clear = ideep = np.array(np.tile(-1,nwave), dtype=np.intc)
            for r in range(rtop, nlayers):
                od.depth_clear[r] = t.optdepth(
                    od.ec_clear[rtop:r+1], od.raypath[r], od.maxdepth,
                    ideep, r,
                )
            ideep[ideep<0] = r

    pyrat.log.head('Optical depth done.')
