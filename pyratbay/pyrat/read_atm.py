# Copyright (c) 2021-2022 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import numpy as np
import scipy.interpolate as sip

from .. import atmosphere as pa
from .. import constants as pc
from .. import io as io
from .. import tools as pt


def make_atmosphere(pyrat):
    """
    Create an atmospheric model object for a given configuration file
    There are four main properties to compute (in this order):
    The pressure, the temperature, the volume mixing ratios (VMRs),
    and the radius profiles.

    Properties can be read from input profiles, computed by models,
    or skipped, depending on the configuration file.

    The rules are simple:
    - if there is a model in the config file, calculate the property
    - else, read from an input profile file
    - otherwise, skip the calculation
    - read VMR from profile only if already read T
    - if calculate p, any further read (T,VMR,r) will interpolate
    """
    log = pyrat.log
    atm = pyrat.atm
    inputs = pyrat.inputs
    log.head('\nGenerating atmospheric profile')

    # User-provided PT profile:
    if pt.isfile(atm.ptfile) == 1:
        input_atm_source = 'ptfile'
        input_atm_file = atm.ptfile
    # Existing atmospheric file:
    #elif pt.isfile(atm.atmfile) == 1 and pyrat.runmode != 'atmosphere':
    elif pt.isfile(atm.atmfile) == 1:
        input_atm_source = 'atmfile'
        input_atm_file = atm.atmfile
    else:
        input_atm_source = None

    # Input-atmosphere data:
    inputs.pressure = None
    inputs.temperature = None
    inputs.vmr = None
    inputs.radius = None
    inputs.atm_species = None
    if input_atm_source is not None:
        log.msg(f"\nReading atmospheric profile from: '{input_atm_file}'")
        input_atmosphere = check_input_atmosphere(input_atm_file, log)

        p_units, t_units, vmr_units, r_units = input_atmosphere[0]
        inputs.pressure = input_atmosphere[2]
        inputs.temperature = input_atmosphere[3]
        if input_atm_source == 'atmfile':
            inputs.atm_species = input_atmosphere[1]
            inputs.vmr = input_atmosphere[4]
            inputs.radius = input_atmosphere[5]
        # TBD: Check for reversed layers order from make_atmprofiles()

    # Figure out where to start:
    p_status = check_pressure(pyrat)
    t_status = check_temperature(pyrat)
    vmr_status = check_chemistry(pyrat, t_status)
    r_status = check_altitude(pyrat, vmr_status)
    #print(
    #    f"Status:\nP {p_status}\nT {t_status}\n"
    #    f"VMR {vmr_status}\nR {r_status}"
    #)

    # Pressure profile:
    if p_status == 'calculate':
        pressure = pa.pressure(
            atm.ptop, atm.pbottom, atm.nlayers, 'barye', log,
        )
    elif p_status == 'read':
        pressure = inputs.pressure
        atm.nlayers = len(pressure)
        if atm.punits is None:
            atm.punits = p_units
    atm.press = pressure

    if p_status == 'calculate' and 'read' in [t_status, vmr_status, r_status]:
        # Interpolate if needed:
        logp_input = np.log(inputs.pressure)
        if t_status == 'read':
            temp_input = np.copy(inputs.temperature)
            temp_extrap = inputs.temperature[0], inputs.temperature[-1]
            temp_interp = sip.interp1d(
                logp_input, temp_input,
                kind='slinear', bounds_error=False, fill_value=temp_extrap,
            )
            inputs.temperature = temp_interp(np.log(pressure))
        if vmr_status == 'read' and inputs.vmr is not None:
            log_vmr_input = np.log(inputs.vmr)
            vmr_extrap = np.log(inputs.vmr[0]), np.log(inputs.vmr[-1])
            vmr_interp = sip.interp1d(
                logp_input, log_vmr_input, axis=0,
                kind='slinear', bounds_error=False, fill_value=vmr_extrap,
            )
            inputs.vmr = np.exp(vmr_interp(np.log(pressure)))
        if r_status == 'read' and inputs.radius is not None:
            rad_interp = sip.interp1d(
                logp_input, inputs.radius, kind='slinear',
            )
            inputs.radius = rad_interp(np.log(pressure))


    # Temperature profile:
    if t_status == 'calculate':
        atm.tmodel = pa.tmodels.get_model(
            atm.tmodelname,
            pressure=atm.press,
            nlayers=atm.nlayers,
        )
        temperature = atm.tmodel(atm.tpars)
    elif t_status == 'read':
        temperature = inputs.temperature
    atm.temp = temperature

    # Composition (volume-mixing-ratio) profiles:
    species = None
    radius = None
    if vmr_status == 'calculate':
        chem_net = pa.chemistry(
            atm.chemistry,
            pressure, temperature, pyrat.inputs.species,
            metallicity=atm.metallicity,
            e_abundances=atm.e_abundances,
            e_scale=atm.e_scale,
            e_ratio=atm.e_ratio,
            solar_file=pyrat.inputs.solar,
            log=log,
            punits=atm.punits,
            q_uniform=pyrat.inputs.uniform,
        )
        atm.chem_model = chem_net
        vmr = chem_net.vmr
        species = chem_net.species
    elif vmr_status == 'read':
        species = inputs.atm_species
        vmr = inputs.vmr
    elif vmr_status == 'skip':
        vmr = None
    atm.vmr = vmr
    # TBD: Deprecate q

    # Set values of species properties:
    pyrat.atm.species = pyrat.mol.name = species
    if species is not None:
        pyrat.mol.nmol = len(pyrat.mol.name)
        get_constants(pyrat)


    # Radius profile:
    if r_status == 'calculate':
        # Mean molecular mass:
        mean_mass = pa.mean_weight(atm.vmr, species)
        # Altitude profile:
        radius = pyrat.hydro(
            pressure, temperature, mean_mass, pyrat.phy.gplanet,
            pyrat.phy.mplanet, atm.refpressure, pyrat.phy.rplanet,
        )
    elif r_status == 'read':
        radius = inputs.radius
    atm.radius = radius
    # TBD: Add extra bits/logs from makesample.make_atmosphere()
    # TBD: Same for read_atm() below?

    # Return atmospheric model if requested:
    if atm.atmfile is not None:
        # Guess radius units if not defined (by rplanet):
        #if radius is not None and atm.runits is None:
        #    atm.runits = 'rjup' if phy.rplanet > 0.5*pc.rjup else 'rearth'
        header = '# pyrat bay atmospheric model\n'
        io.write_atm(
            atm.atmfile, pressure, temperature, species,
            vmr, radius, atm.punits, atm.runits, header=header,
        )
        log.msg(f"Output atmospheric file: '{atm.atmfile}'.")

    return


# TBD: Delete, function no longer used
def read_atm(pyrat):
    """Read an atmospheric file, store variables into pyrat object."""
    log = pyrat.log

    # User-input atmospheric-data object:
    atm = pyrat.atm
    atm_in = pyrat.inputs.atm

    with pt.log_error(log):
        atm_inputs = io.read_atm(pyrat.atm.atmfile)

    punits, tunits, qunits, runits = atm_inputs[0]

    # Planetary radius units (if not set by rplanet):
    if atm.runits is None:
        atm.runits = runits
    if atm.runits is None:
        atm.runits = 'rearth'

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
        log.error(f"Invalid input abundance units '{qunits}'.")

    # Calculate the mean molecular mass per layer:
    atm_in.mm = np.sum(atm_in.q*pyrat.mol.mass, axis=1)

    # Calculate number density profiles for each molecule (in molecules cm-3):
    atm_in.d = pa.ideal_gas_density(atm_in.q, atm_in.press, atm_in.temp)

    log.msg(f"Species list:\n  {pyrat.mol.name}", indent=2, si=4)

    log.msg(
        f"Abundances are given by {qunits} mixing ratio.", indent=2)
    log.msg(
        f"Unit factors: radius: {runits}, pressure: {punits}, "
        f"temperature: {tunits}", indent=2)

    log.msg(
        f"Number of layers in the input atmospheric file: {atm_in.nlayers}",
        indent=2)
    log.msg("Atmospheric file pressure limits: "
        f"{atm_in.press[ 0]/pt.u(atm.punits):.2e}--"
        f"{atm_in.press[-1]/pt.u(atm.punits):.2e} {atm.punits}.",
        indent=2)

    log.msg(
        f"Median mean molecular mass: {np.median(atm_in.mm):.3f} g mol-1.",
        indent=2)


def get_constants(pyrat):
    """Set molecules constant values (mass, radius)."""
    # Read file with molecular info:
    pyrat.log.msg(
        f"Taking species constant parameters from: '{pyrat.mol.molfile}'.",
        indent=2)
    symbol, mass, radius = io.read_molecs(pyrat.mol.molfile)

    # Check that all atmospheric species are listed in molfile:
    absent = np.setdiff1d(pyrat.mol.name, symbol)
    if len(absent) > 0:
        pyrat.log.error(
            f"These species: {absent} are not listed in the molecules "
            f"info file: {pyrat.mol.molfile}.")

    # Set molecule's values:
    pyrat.mol.symbol = np.zeros(pyrat.mol.nmol, 'U20')
    pyrat.mol.mass = np.zeros(pyrat.mol.nmol)
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
        pyrat.mol.mass[i] = mass[imol]
        pyrat.mol.radius[i] = radius[imol] * pc.A
        pyrat.log.msg(
            f"{pyrat.mol.name[i]:>10s}:  "
            f"{pyrat.mol.radius[i]/pc.A:.3f}  "
            f"{pyrat.mol.mass[i]:8.4f}",
            indent=2,
        )


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
    # Temperature profile:
    if temp is not None:
        # Need to null tpars since it does not represent temp anymore
        atm.tpars = None
    elif atm.tpars is not None:
        temp = atm.tmodel(atm.tpars)
    else:
        temp = atm.temp

    if temp is not None:
        # Check that the dimensions match:
        if np.size(temp) != atm.nlayers:
            pyrat.log.error(
                f"The temperature array size ({np.size(temp)}) doesn't match "
                f"the Pyrat's temperature size ({np.size(atm.temp)}).")

        # Check temperature boundaries:
        msg = (
            "One or more input temperature values lies out of the {:s} "
            "temperature boundaries (K): [{:6.1f}, {:6.1f}].")
        if pyrat.ex.extfile is not None:
            if np.any(temp > pyrat.ex.tmax) or np.any(temp < pyrat.ex.tmin):
                msg = msg.format(
                    'tabulated extinction-coefficient',
                    pyrat.ex.tmin, pyrat.ex.tmax)
                pyrat.log.warning(msg)
                return 0
        elif pyrat.lt.ntransitions > 0:
            if np.any(temp > pyrat.lt.tmax) or np.any(temp < pyrat.lt.tmin):
                msg = msg.format('line-transition', pyrat.lt.tmin, pyrat.lt.tmax)
                pyrat.log.warning(msg)
                return 0
        if pyrat.cs.nfiles > 0:
            if np.any(temp > pyrat.cs.tmax) or np.any(temp < pyrat.cs.tmin):
                msg = msg.format('cross-section', pyrat.cs.tmin, pyrat.cs.tmax)
                pyrat.log.warning(msg)
                return 0
        atm.temp = temp


    # Volume mixing ratios:
    base_vmr = np.copy(atm.qbase)
    if vmr is not None:
        atm.molpars = None
    elif atm.molpars is not None:
        vmr = pa.qscale(
            base_vmr, pyrat.mol.name, atm.molmodel,
            atm.molfree, atm.molpars, atm.bulk,
            iscale=atm.ifree, ibulk=atm.ibulk,
            bratio=atm.bulkratio, invsrat=atm.invsrat,
        )
    else:
        vmr = base_vmr

    if vmr is not None:
        if np.shape(vmr) != np.shape(atm.q):
            pyrat.log.error(
                f"The shape of the abundances array {np.shape(vmr)} doesn't "
                 "match the shape of the Pyrat's abundance size "
                f"{np.shape(atm.q)}")
        atm.q = vmr

    # Mean molecular mass:
    atm.mm = np.sum(atm.q*pyrat.mol.mass, axis=1)

    # Radius profile:
    if radius is not None:
        atm.radius = radius
    elif atm.rmodelname is None:
        atm.radius = pyrat.inputs.atm.radius
    else:
        atm.radius = pyrat.hydro(
            atm.press, atm.temp, atm.mm,
            pyrat.phy.gplanet, pyrat.phy.mplanet,
            atm.refpressure, pyrat.phy.rplanet,
        )


    # Number density (molecules cm-3):
    atm.d = pa.ideal_gas_density(atm.q, atm.press, atm.temp)

    # Check radii lie within Hill radius:
    pyrat.phy.rhill = pa.rhill(
        pyrat.phy.smaxis, pyrat.phy.mplanet, pyrat.phy.mstar,
    )
    atm.rtop = pt.ifirst(atm.radius<pyrat.phy.rhill, default_ret=0)
    if atm.rtop > 0:
        rhill = pyrat.phy.rhill/pt.u(atm.runits)
        if pyrat.runmode == 'mcmc':
            logger = pyrat.log.msg
        else:
            logger = pyrat.log.warning
        logger(
            "The atmospheric pressure array extends beyond the Hill radius "
            f"({rhill:.5f} {atm.runits}) at pressure "
            f"{atm.press[atm.rtop]/pc.bar:.3e} bar (layer {atm.rtop}).  "
            "Extinction beyond this layer will be neglected.",
        )

    # Partition function:
    for database in pyrat.lt.db:
        for i_isotope in range(database.niso):
            zinterp = sip.interp1d(
                database.temp, database.z[i_isotope], kind='slinear')
            pyrat.iso.z[database.iiso+i_isotope] = zinterp(atm.temp)



def check_input_atmosphere(atm_file, log):
    """
    Make sure that the input atmospheric profile is sorted correctly
    (from low pressure to high pressure).
    """
    input_atm = io.read_atm(atm_file)
    p_units, t_units, vmr_units, r_units = units = input_atm[0]
    species = input_atm[1]
    pressure = input_atm[2] * pt.u(p_units)
    temperature = input_atm[3]
    vmr = input_atm[4]
    radius = input_atm[5]
    if radius is not None:
        radius *= pt.u(r_units)

    # Check that the layers are sorted from top to bottom:
    sort    = np.all(np.ediff1d(pressure) > 0)
    reverse = np.all(np.ediff1d(pressure) < 0)
    if radius is not None:
        sort    &= np.all(np.ediff1d(radius) < 0)
        reverse &= np.all(np.ediff1d(radius) > 0)

    if sort:  # Layers are in the correct order
        pass
    elif reverse:  # Layers in reverse order
        pressure = np.flip(pressure)
        temperature = np.flip(temperature)
        if vmr is not None:
            vmr = np.flip(vmr)
        if radius is not None:
            radius = np.flip(radius)
    else:
        log.error(
            'The layers of input atmosphere are neither sorted '
            'from top to bottom, nor from bottom to top'
        )
    return units, species, pressure, temperature, vmr, radius


def check_pressure(pyrat):
    """
    Determine whether to calculate or read the atmospheric pressure
    profile.
    """
    all_parameters_defined = (
        pyrat.atm.nlayers is not None and
        pyrat.atm.ptop is not None and
        pyrat.atm.pbottom is not None
    )
    if all_parameters_defined:
        return 'calculate'

    # User-provided PT profile:
    if pyrat.inputs.pressure is not None:
        return 'read'

    pyrat.log.error(
        'Cannot compute pressure profile, either set {ptop, pbottom, nlayers} '
        'parameters, or provide an input PT profile (ptfile) or atmospheric '
        'file (atmfile)'
    )


def check_temperature(pyrat):
    """
    Determine whether to calculate or read the atmospheric temperature
    profile.
    """
    # A model takes precedence:
    if pyrat.atm.tmodelname is not None:
        if pyrat.atm.tpars is None:
            pyrat.log.error('Undefined temperature-model parameters (tpars)')
        return 'calculate'

    # User-provided PT profile:
    if pyrat.inputs.temperature is not None:
        return 'read'

    pyrat.log.error(
        'Cannot compute temperature profile, either set a temperature model '
        '(tmodelname) and parameters (tpars), or provide an input PT '
        'profile (ptfile) or atmospheric file (atmfile)'
    )


def check_chemistry(pyrat, t_status):
    """
    Determine whether to calculate, read, or skip the atmospheric compostion
    profiles.
    """
    atm = pyrat.atm
    log = pyrat.log

    if atm.chemistry is None:
        if t_status == 'read' and pyrat.inputs.vmr is not None:
            if pyrat.inputs.species is not None:
                pyrat.log.warning(
                    "Composition will be taken from input atmospheric file, "
                    "user input 'species' will be ingnored"
                )
            return 'read'
        return 'skip'

    if pyrat.inputs.species is None:
        log.error('Undefined atmospheric species list (species).')

    # Uniform-abundances profile:
    if atm.chemistry == 'uniform':
        if pyrat.inputs.uniform is None:
            log.error(
                'Undefined list of uniform volume mixing ratios '
                f'(uniform) for {atm.chemistry} chemistry model.'
            )
        nuniform = len(pyrat.inputs.uniform)
        nspecies = len(pyrat.inputs.species)
        if nuniform != nspecies:
            pyrat.log.error(
                f'Number of uniform abundances ({nuniform}) does '
                f'not match the number of species ({nspecies}).'
            )

    return 'calculate'


def check_altitude(pyrat, vmr_status):
    """
    Determine whether to calculate, read, or skip the atmospheric
    altitude profile.
    """
    phy = pyrat.phy

    if pyrat.atm.rmodelname is None:
        if vmr_status == 'read' and pyrat.inputs.radius is not None:
            return 'read'
        return 'skip'

    if vmr_status == 'skip':
        pyrat.log.error(
            'Cannot compute hydrostatic-equilibrium radius profile.\n'
            'radius model needs to know the composition'
        )

    err = []
    if phy.rplanet is None:
        err += ['Undefined planet radius (rplanet).']
    if phy.mplanet is None and phy.gplanet is None:
        err += ['Undefined planet mass (mplanet) or surface gravity (gplanet).']
    if pyrat.atm.refpressure is None:
        err += ['Undefined reference pressure level (refpressure).']

    if len(err) != 0:
        error_message = '\n'.join(err)
        pyrat.log.error(
            'Cannot compute hydrostatic-equilibrium radius profile.\n'
            f'{error_message}'
        )

    return 'calculate'
