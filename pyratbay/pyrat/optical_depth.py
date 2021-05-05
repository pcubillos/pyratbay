# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import ctypes
import multiprocessing as mpr

import numpy as np

from . import extinction as ex
from .. import atmosphere as pa
from .. import constants as pc
from ..lib import _extcoeff as ec
from ..lib import _trapz as t
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
    pyrat.ex.ec = np.zeros((nlayers, nwave))
    od.ec       = np.empty((nlayers, nwave))
    od.depth    = np.zeros((nlayers, nwave))
    if pyrat.cloud.fpatchy is not None:
        od.ec_clear = np.empty((nlayers, nwave))
        od.depth_clear = np.zeros((nlayers, nwave))

    # Calculate the ray path:
    if pyrat.od.rt_path in pc.emission_rt:
        pyrat.od.raypath = -cu.ediff(pyrat.atm.radius)
    elif pyrat.od.rt_path in pc.transmission_rt:
        pyrat.od.raypath = pa.transit_path(pyrat.atm.radius, pyrat.atm.rtop)

    # Interpolate extinction coefficient from table:
    if pyrat.ex.extfile is not None:
        r = rtop
        imol = [list(pyrat.mol.name).index(mol) for mol in pyrat.ex.species]
        while r < nlayers:
            ec.interp_ec(
                pyrat.ex.ec[r], pyrat.ex.etable[:,:,r,:],
                pyrat.ex.temp, pyrat.atm.temp[r], pyrat.atm.d[r,imol])
            r += 1

    # Calculate the extinction coefficient on the spot:
    elif pyrat.lt.tlifile is not None:
        sm_ext = mpr.Array(ctypes.c_double,
            np.zeros(nlayers*nwave, np.double))
        pyrat.ex.ec = np.ctypeslib.as_array(
            sm_ext.get_obj()).reshape((nlayers, nwave))
        processes = []
        indices = np.arange(rtop, nlayers) % pyrat.ncpu
        for i in range(pyrat.ncpu):
            proc = mpr.Process(target=ex.extinction,   #      grid   add
                        args=(pyrat, np.where(indices==i)[0], False, True))
            processes.append(proc)
            proc.start()
        for proc in processes:
            proc.join()

    # Sum all contributions to the extinction:
    od.ec[rtop:] = (
        + pyrat.ex.ec[rtop:]
        + pyrat.cs.ec[rtop:]
        + pyrat.rayleigh.ec[rtop:]
        + pyrat.cloud.ec[rtop:]
        + pyrat.alkali.ec[rtop:]
    )
    # If fpatchy, compute a separate spectrum with clear skies:
    if pyrat.cloud.fpatchy is not None:
        od.ec_clear[rtop:] = np.copy(od.ec[rtop:]) - pyrat.cloud.ec[rtop:]

    rbottom = nlayers
    if 'deck' in (m.name for m in pyrat.cloud.models):
        deck = pyrat.cloud.models[pyrat.cloud.model_names.index('deck')]
        rbottom = deck.itop + 1
    # Calculate the optical depth for each wavenumber:
    if od.rt_path in pc.emission_rt:
        od.ideep = np.tile(nlayers-1, nwave)
        i = 0
        while i < nwave:
            od.ideep[i] = rtop - 1 + t.cumtrapz(
                od.depth[rtop:,i], od.ec[rtop:,i], od.raypath[rtop:rbottom],
                od.maxdepth)
            i += 1

    elif od.rt_path in pc.transmission_rt:
        od.ideep = ideep = np.array(np.tile(-1, nwave), dtype=np.intc)
        r = rtop
        # Optical depth at each level (tau = 2.0*integral e*ds):
        for r in range(rtop, rbottom):
            od.depth[r] = t.optdepth(
                od.ec[rtop:r+1], od.raypath[r], od.maxdepth, ideep, r)
        ideep[ideep<0] = r

        if pyrat.cloud.fpatchy is not None:
            rbottom = nlayers
            od.ideep_clear = ideep = np.array(np.tile(-1,nwave), dtype=np.intc)
            for r in range(rtop, rbottom):
                od.depth_clear[r] = t.optdepth(
                    od.ec_clear[rtop:r+1], od.raypath[r], od.maxdepth,
                    ideep, r)
            ideep[ideep<0] = r

    pyrat.log.head('Optical depth done.')
