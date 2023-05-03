# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import ctypes
import multiprocessing as mp

import numpy as np

from . import extinction as ex_mod
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
    ex = pyrat.ex
    lbl = pyrat.lbl
    nwave = pyrat.spec.nwave
    nlayers = pyrat.atm.nlayers
    rtop = pyrat.atm.rtop
    good_status = True

    pyrat.log.head('\nBegin optical-depth calculation.')

    # Evaluate the extinction coefficient at each layer:
    ex.ec = np.zeros((nlayers, nwave))
    od.ec = np.empty((nlayers, nwave))
    od.depth = np.zeros((nlayers, nwave))
    if pyrat.cloud.fpatchy is not None:
        od.ec_clear = np.empty((nlayers, nwave))
        od.depth_clear = np.zeros((nlayers, nwave))

    # Calculate the ray path:
    if pyrat.od.rt_path in pc.emission_rt:
        pyrat.od.raypath = -cu.ediff(pyrat.atm.radius)
    elif pyrat.od.rt_path in pc.transmission_rt:
        pyrat.od.raypath = pa.transit_path(pyrat.atm.radius, pyrat.atm.rtop)

    # Interpolate extinction coefficient from table:
    if ex.extfile is not None:
        if np.any(pyrat.atm.temp > ex.tmax) or np.any(pyrat.atm.temp < ex.tmin):
            pyrat.log.warning(
                "Atmospheric temperature values lie out of the cross-section "
                f"boundaries (K): [{ex.tmin:6.1f}, {ex.tmax:6.1f}]"
            )
            good_status = False
            return good_status

        r = rtop
        imol = [list(pyrat.atm.species).index(mol) for mol in ex.species]
        ec.interp_ec(
            ex.ec, ex.etable,
            ex.temp, pyrat.atm.temp, pyrat.atm.d[:,imol], rtop,
        )

    # Calculate the extinction coefficient on the spot:
    elif lbl.tlifile is not None:
        imol = lbl.mol_index[lbl.mol_index>=0]
        ex.species = pyrat.atm.species[np.unique(imol)]
        ex.nspec = len(ex.species)

        if np.any(pyrat.atm.temp>lbl.tmax) or np.any(pyrat.atm.temp<lbl.tmin):
            pyrat.log.warning(
                "Atmospheric temperature values lie out of the line-transition "
                f"boundaries (K): [{lbl.tmin:6.1f}, {lbl.tmax:6.1f}]"
            )
            good_status = False
            return good_status

        # Update partition functions:
        lbl.iso_pf = np.zeros((lbl.niso, nlayers))
        for i in range(lbl.niso):
            lbl.iso_pf[i] = lbl.iso_pf_interp[i](pyrat.atm.temp)

        sm_ext = mp.Array(
            ctypes.c_double, np.zeros(nlayers*nwave, np.double))
        ex.ec = np.ctypeslib.as_array(
            sm_ext.get_obj()).reshape((nlayers, nwave))
        processes = []
        indices = np.arange(rtop, nlayers) % pyrat.ncpu
        for i in range(pyrat.ncpu):
            grid = False
            add = True
            proc = mp.get_context('fork').Process(
                target=ex_mod.extinction,
                args=(pyrat, np.where(indices==i)[0], grid, add))
            processes.append(proc)
            proc.start()
        for proc in processes:
            proc.join()

    # Sum all contributions to the extinction:
    od.ec[rtop:] = (
        + ex.ec[rtop:]
        + pyrat.cloud.ec[rtop:]
    )
    if pyrat.rayleigh.ec is not None:
        od.ec[rtop:] += pyrat.rayleigh.ec[rtop:]
    if pyrat.alkali.ec is not None:
        od.ec[rtop:] += pyrat.alkali.ec[rtop:]
    if pyrat.h_ion.model is not None:
        od.ec[rtop:] += pyrat.h_ion.ec[rtop:]
    if pyrat.cs.ec is not None:
        od.ec[rtop:] += pyrat.cs.ec[rtop:]

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
        maxdepth = np.inf if od.rt_path=='emission_two_stream' else od.maxdepth
        i = 0
        while i < nwave:
            od.ideep[i] = rtop - 1 + t.cumtrapz(
                od.depth[rtop:,i],
                od.ec[rtop:,i],
                od.raypath[rtop:rbottom],
                maxdepth)
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
    return good_status
