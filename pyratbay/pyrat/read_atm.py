# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import numpy as np

from .. import atmosphere as pa
from .. import constants as pc
from .. import tools as pt


def update_atm(
        pyrat, temp=None, vmr=None, radius=None,
        # Deprecated parameters:
        abund=None,
    ):
    """
    Update temperature, abundances, and radius profiles of a pyrat object.

    Parameters
    ----------
    pyrat: A Pyrat instance
    temp: 1D float ndarray
        Layer's temperature profile (Kelvin) sorted from top to bottom.
    abund: 2D float ndarray
        Species mole mixing ratio profiles [nlayers, nmol].
    radius: 1D float ndarray
        Layer's altitude profile (in cm), same order as temp.
    """
    atm = pyrat.atm
    phy = pyrat.phy
    ex = pyrat.ex
    lt = pyrat.lt

    # Temperature profile:
    if temp is not None:
        # Need to null tpars since it does not represent temp anymore
        atm.tpars = None
    elif atm.tpars is not None:
        temp = atm.temp_model(atm.tpars)
    else:
        temp = atm.temp

    # Check that the dimensions match:
    if np.size(temp) != atm.nlayers:
        pyrat.log.error(
            f"The temperature array size ({np.size(temp)}) doesn't match "
            f"the Pyrat's temperature size ({np.size(atm.temp)})"
        )

    # Check temperature boundaries:
    msg = (
        "Atmospheric temperature values lie out of the {:s} "
        "boundaries (K): [{:6.1f}, {:6.1f}]."
    )
    # TBD: Check if retrieval runs need to pass a verbosity to mute warnings
    if ex.extfile is not None:
        if np.any(temp > ex.tmax) or np.any(temp < ex.tmin):
            msg = msg.format('extinction-coefficient', ex.tmin, ex.tmax)
            pyrat.log.warning(msg)
            return 0
    elif lt.ntransitions > 0:
        if np.any(temp > lt.tmax) or np.any(temp < lt.tmin):
            msg = msg.format('line-transition', lt.tmin, lt.tmax)
            pyrat.log.warning(msg)
            return 0
    if pyrat.cs.nfiles > 0:
        if np.any(temp > pyrat.cs.tmax) or np.any(temp < pyrat.cs.tmin):
            msg = msg.format('cross-section', pyrat.cs.tmin, pyrat.cs.tmax)
            pyrat.log.warning(msg)
            return 0
    atm.temp = temp

    # Volume mixing ratios:
    if vmr is not None:
        atm.molpars = []

    elif atm.chemistry == 'tea':
        net = atm.chem_model
        metallicity = None
        e_abundances = {}
        e_ratio = {}
        #e_scale = {}
        if np.any(atm._equil_var):
            equil_vars = np.array(atm.mol_pnames)[atm._equil_var]
            equil_pars = np.array(atm.molpars)[atm._equil_var]
        else:
            equil_vars, equil_pars = [], []
        for var,val in zip(equil_vars, equil_pars):
            if var == 'metal':
                metallicity = val
            elif var.startswith('[') and var.endswith('/H]'):
                element = var[1:-3]
                idx = list(net._base_composition).index(element)
                solar_abundance = net._base_dex_abundances[idx]
                e_abundances[var] = solar_abundance + val
            elif '/' in var:
                e_ratio[var.replace('/','_')] = val
        vmr = net.thermochemical_equilibrium(
            atm.temp,
            metallicity=metallicity,
            e_abundances=e_abundances,
            e_ratio=e_ratio,
            #e_scale=e_scale,
        )
    elif np.any(~atm._equil_var) and atm.molpars is not None:
        vmr_vars = np.array(atm.mol_pnames)[~atm._equil_var]
        vmr_pars = np.array(atm.molpars)[~atm._equil_var]
        vmr = pa.qscale(
            np.copy(atm.base_vmr),
            pyrat.atm.species,
            vmr_vars, vmr_pars, atm.bulk,
            iscale=atm.ifree, ibulk=atm.ibulk,
            bratio=atm.bulkratio, invsrat=atm.invsrat,
        )
    else:
        vmr = np.copy(atm.base_vmr)

    if np.shape(vmr) != np.shape(atm.vmr):
        pyrat.log.error(
            f"The shape of the abundances array {np.shape(vmr)} doesn't "
             "match the shape of the Pyrat's abundance size "
            f"{np.shape(atm.vmr)}"
        )
    atm.vmr = vmr

    # Mean molecular mass:
    atm.mm = np.sum(atm.vmr * atm.mol_mass, axis=1)

    # Radius profile (take input, re-compute, or keep current):
    if radius is not None:
        atm.radius = radius
    elif atm.rmodelname is not None:
        atm.radius = atm.rad_model(
            atm.press, atm.temp, atm.mm,
            atm.mplanet, atm.gplanet,
            atm.refpressure, atm.rplanet,
        )
    else:
        pass

    # Number density (molecules cm-3):
    atm.d = pa.ideal_gas_density(atm.vmr, atm.press, atm.temp)

    # Check radii lie within Hill radius:
    phy.rhill = pa.hill_radius(phy.smaxis, atm.mplanet, phy.mstar)
    atm.rtop = pt.ifirst(atm.radius<phy.rhill, default_ret=0)
    if atm.rtop > 0:
        rhill = phy.rhill / pt.u(atm.runits)
        if pyrat.runmode == 'mcmc':
            log_call = pyrat.log.msg
        else:
            log_call = pyrat.log.warning
        log_call(
            "The atmospheric pressure array extends beyond the Hill radius "
            f"({rhill:.5f} {atm.runits}) at pressure "
            f"{atm.press[atm.rtop]/pc.bar:.3e} bar (layer {atm.rtop}).  "
            "Extinction beyond this layer will be neglected.",
        )

    # Partition function:
    for i in range(pyrat.iso.niso):
        pyrat.iso.z[i] = pyrat.iso.zinterp[i](atm.temp)

