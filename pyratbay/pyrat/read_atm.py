# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import numpy as np
import scipy.interpolate as sip

from .. import atmosphere as pa
from .. import constants as pc
from .. import io as io
from .. import tools as pt


def read_atm(pyrat):
    """Read the atmospheric file, store variables into pyrat object."""
    # Check atmfile:
    pyrat.log.head(f"\nReading atmospheric file: '{pyrat.atm.atmfile}'.")

    # User-input atmospheric-data object:
    atm    = pyrat.atm
    atm_in = pyrat.inputs.atm

    with pt.log_error(pyrat.log):
        atm_inputs = io.read_atm(pyrat.atm.atmfile)

    # Atmospheric-file units, species, and profiles:
    punits, tunits, qunits, runits = atm_inputs[0]
    pyrat.mol.name = atm_inputs[1]
    pyrat.mol.nmol = len(pyrat.mol.name)

    # Last resort to set the pressure units:
    if atm.punits is None:
        atm.punits = punits

    # Planetary radius units (if not set by rplanet):
    if atm.runits is None:
        atm.runits = runits
    if atm.runits is None:
        atm.runits = 'km'

    # Read molecular constant values:
    get_constants(pyrat)

    # Store values in CGS system of units:
    atm_in.press = atm_inputs[2] * pt.u(punits)
    atm_in.temp = atm_inputs[3] * pt.u(tunits)
    atm_in.q = atm_inputs[4]
    if atm_inputs[5] is not None:
        atm_in.radius = atm_inputs[5] * pt.u(runits)
    atm_in.nlayers = len(atm_in.press)

    # Store the abundances as volume mixing ratio:
    if qunits == 'mass':
        atm_in.q /= pyrat.mol.mass * np.sum(atm_in.q/pyrat.mol.mass,axis=1)
    elif qunits not in ['volume', 'number']:
        pyrat.log.error(f"Invalid input abundance units '{qunits}'.")

    # Calculate the mean molecular mass per layer:
    atm_in.mm = np.sum(atm_in.q*pyrat.mol.mass, axis=1)

    # Calculate number density profiles for each molecule (in molecules cm-3):
    atm_in.d = pa.ideal_gas_density(atm_in.q, atm_in.press, atm_in.temp)

    pyrat.log.msg(f"Species list:\n  {pyrat.mol.name}", indent=2, si=4)

    pyrat.log.msg(
        f"Abundances are given by {qunits} mixing ratio.", indent=2)
    pyrat.log.msg(
        f"Unit factors: radius: {runits}, pressure: {punits}, "
        f"temperature: {tunits}", indent=2)

    pyrat.log.msg(
        f"Number of layers in the input atmospheric file: {atm_in.nlayers}",
        indent=2)
    pyrat.log.msg("Atmospheric file pressure limits: "
        f"{atm_in.press[ 0]/pt.u(atm.punits):.2e}--"
        f"{atm_in.press[-1]/pt.u(atm.punits):.2e} {atm.punits}.",
        indent=2)

    pyrat.log.msg(
        f"Median mean molecular mass: {np.median(atm_in.mm):.3f} g mol-1.",
        indent=2)

    pyrat.log.head("Read atmosphere done.")


def get_constants(pyrat):
    """Set molecules constant values (mass, radius)."""
    # Read file with molecular info:
    pyrat.log.msg(
        f"Taking species constant parameters from: '{pyrat.mol.molfile}'.",
        indent=2)
    symbol, mass, diam = io.read_molecs(pyrat.mol.molfile)

    # Check that all atmospheric species are listed in molfile:
    absent = np.setdiff1d(pyrat.mol.name, symbol)
    if len(absent) > 0:
        pyrat.log.error(
            f"These species: {absent} are not listed in the molecules "
            f"info file: {pyrat.mol.molfile}.")

    # Set molecule's values:
    pyrat.mol.symbol = np.zeros(pyrat.mol.nmol, 'U20')
    pyrat.mol.mass   = np.zeros(pyrat.mol.nmol)
    pyrat.mol.radius = np.zeros(pyrat.mol.nmol)

    pyrat.log.msg(
        'Molecule   Radius  Mass\n'
        '           (A)     (gr/mol)',
        indent=4)
    for i in range(pyrat.mol.nmol):
        # Find the molecule in the list:
        imol = np.where(symbol == pyrat.mol.name[i])[0]
        # Set molecule name, mass, and collision radius:
        pyrat.mol.symbol[i] = symbol[imol][0]
        pyrat.mol.mass[i]   = mass  [imol]
        pyrat.mol.radius[i] = 0.5*diam[imol] * pc.A
        pyrat.log.msg(
            f"{pyrat.mol.name[i]:>10s}:  "
            f"{pyrat.mol.radius[i]/pc.A:.3f}  "
            f"{pyrat.mol.mass[i]:8.4f}", indent=2)


def update_atm(pyrat, temp=None, abund=None, radius=None):
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
    # Recompute temperature profile:
    if temp is not None:
        # Need to null tpars since it does not represent temp anymore
        pyrat.atm.tpars = None
    elif pyrat.atm.tpars is not None:
        temp = pyrat.atm.tmodel(pyrat.atm.tpars)
    else:
        temp = pyrat.atm.temp

    # Check that the dimensions match:
    if np.size(temp) != np.size(pyrat.atm.temp):
        pyrat.log.error(
            f"The temperature array size ({np.size(temp)}) doesn't match "
            f"the Pyrat's temperature size ({np.size(pyrat.atm.temp)}).")

    # Check temperature boundaries:
    error = (
        "One or more input temperature values lies out of the {:s} "
        "temperature boundaries (K): [{:6.1f}, {:6.1f}].")
    if pyrat.ex.extfile is not None:
        if np.any(temp > pyrat.ex.tmax) or np.any(temp < pyrat.ex.tmin):
            pyrat.log.warning(error.format(
                'tabulated extinction-coefficient',
                pyrat.ex.tmin, pyrat.ex.tmax))
            return 0
    elif pyrat.lt.ntransitions > 0:
        if np.any(temp > pyrat.lt.tmax) or np.any(temp < pyrat.lt.tmin):
            pyrat.log.warning(error.format(
                'line-transition', pyrat.lt.tmin, pyrat.lt.tmax))
            return 0
    if pyrat.cs.nfiles > 0:
        if np.any(temp > pyrat.cs.tmax) or np.any(temp < pyrat.cs.tmin):
            pyrat.log.warning(error.format(
                'cross-section', pyrat.cs.tmin, pyrat.cs.tmax))
            return 0

    # Recompute abundance profiles:
    q0 = np.copy(pyrat.atm.qbase)
    if abund is not None:
        pyrat.atm.molpars = None
    elif pyrat.atm.molpars is not None:
        abund = pa.qscale(
            q0, pyrat.mol.name, pyrat.atm.molmodel,
            pyrat.atm.molfree, pyrat.atm.molpars, pyrat.atm.bulk,
            iscale=pyrat.atm.ifree, ibulk=pyrat.atm.ibulk,
            bratio=pyrat.atm.bulkratio, invsrat=pyrat.atm.invsrat)
    else:
        abund = q0
    if np.shape(abund) != np.shape(pyrat.atm.q):
        pyrat.log.error(
            f"The shape of the abundances array {np.shape(abund)} doesn't "
             "match the shape of the Pyrat's abundance size "
            f"{np.shape(pyrat.atm.q)}")

    # Update values:
    pyrat.atm.temp = temp
    pyrat.atm.q    = abund

    # Mean molecular mass:
    pyrat.atm.mm = np.sum(pyrat.atm.q*pyrat.mol.mass, axis=1)

    # Number density (molecules cm-3):
    pyrat.atm.d = pa.ideal_gas_density(
        pyrat.atm.q, pyrat.atm.press, pyrat.atm.temp)

    # Take radius if provided, else use hydrostatic-equilibrium equation:
    if radius is not None:
        pyrat.atm.radius = radius
    elif pyrat.atm.rmodelname is None:
        pyrat.atm.radius = pyrat.inputs.atm.radius
    else:
        pyrat.atm.radius = pyrat.hydro(
            pyrat.atm.press, pyrat.atm.temp, pyrat.atm.mm,
            pyrat.phy.gplanet, pyrat.phy.mplanet,
            pyrat.atm.refpressure, pyrat.phy.rplanet)

    # Check radii lie within Hill radius:
    pyrat.phy.rhill = pa.rhill(
        pyrat.phy.smaxis, pyrat.phy.mplanet, pyrat.phy.mstar)
    pyrat.atm.rtop = pt.ifirst(pyrat.atm.radius<pyrat.phy.rhill, default_ret=0)
    if pyrat.atm.rtop > 0:
        rhill = pyrat.phy.rhill/pt.u(pyrat.atm.runits)
        if pyrat.runmode == 'mcmc':
            logger = pyrat.log.msg
        else:
            logger = pyrat.log.warning
        logger(
            "The atmospheric pressure array extends beyond the Hill radius "
           f"({rhill:.5f} {pyrat.atm.runits}) at pressure "
           f"{pyrat.atm.press[pyrat.atm.rtop]/pc.bar:.3e} bar (layer "
           f"{pyrat.atm.rtop}).  "
            "Extinction beyond this layer will be neglected.")

    # Partition function:
    for db in pyrat.lt.db:            # For each Database
        for j in np.arange(db.niso):  # For each isotope in DB
            zinterp = sip.interp1d(db.temp, db.z[j], kind='slinear')
            pyrat.iso.z[db.iiso+j] = zinterp(pyrat.atm.temp)
