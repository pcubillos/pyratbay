# Copyright (c) 2021-2023 Patricio Cubillos
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
    - else if there is an input_atmfile or ptfile, read properties from file
    - else, skip the calculation
    - if calculate p, any further reads (T,VMR,r) will interpolate
    """
    log = pyrat.log
    atm = pyrat.atm
    mol = pyrat.mol
    inputs = pyrat.inputs
    log.head('\nGenerating atmospheric profile')

    # User-provided PT profile:
    if pt.isfile(atm.ptfile) == 1:
        input_atm_source = 'ptfile'
        input_atm_file = atm.ptfile
    # Existing atmospheric file:
    elif pt.isfile(atm.input_atmfile) == 1:
        input_atm_source = 'atmfile'
        input_atm_file = atm.input_atmfile
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
        # Always store the abundances as volume mixing ratios:
        if inputs.vmr is not None:
            if vmr_units == 'mass':
                inputs.vmr /= mol.mass * np.sum(inputs.vmr/mol.mass, axis=1)
            elif vmr_units not in ['volume', 'number']:
                log.error(f"Invalid input abundance units '{vmr_units}'.")

    # Figure out where to start:
    p_status = check_pressure(pyrat)
    t_status = check_temperature(pyrat)
    vmr_status = check_chemistry(pyrat, t_status)
    r_status = check_altitude(pyrat, vmr_status)

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
        # TBD: Do I want to update with input ptop/pbottom?
        atm.ptop = np.amin(pressure)
        atm.pbottom = np.amax(pressure)
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
            pressure, temperature, inputs.species,
            metallicity=atm.metallicity,
            e_abundances=atm.e_abundances,
            e_scale=atm.e_scale,
            e_ratio=atm.e_ratio,
            solar_file=inputs.solar,
            log=log,
            punits=atm.punits,
            q_uniform=inputs.uniform,
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
        if atm.runits is None:
            atm.runits = r_units
    atm.radius = radius
    # Planetary radius units (if not set by rplanet nor atmfile):
    if atm.runits is None:
        if pyrat.phy.rplanet is not None:
            atm.runits = 'rjup' if pyrat.phy.rplanet > 0.5*pc.rjup else 'rearth'
        else:
            atm.runits = 'rearth'

    # Mean molecular mass and number densities:
    mmm_text = ''
    if atm.vmr is not None:
        atm.mm = np.sum(atm.vmr*mol.mass, axis=1)
        # Number density profiles for each molecule (in molecules cm-3):
        atm.d = pa.ideal_gas_density(atm.vmr, atm.press, atm.temp)
        mmm_text += (
            f"\nMedian mean molecular mass: {np.median(atm.mm):.3f} g mol-1."
        )
        # Base abundance profiles:
        atm.base_vmr = np.copy(atm.vmr)

    # Print radius array:
    if atm.radius is not None:
        radius_arr = atm.radius / pt.u(atm.runits)
        log.msg(
            f'Radius array ({atm.runits}) = \n{radius_arr}', indent=2, si=4,
        )
        log.msg(
            'Upper/lower radius boundaries:    '
            f'{radius_arr[atm.rtop]:.5f} - {radius_arr[-1]:.5f} {atm.runits}.',
            indent=2,
        )

    #log.msg(
    #    'Lower/higher pressure boundaries: '
    #   f'{atm.press[atm.rtop]/punits:.2e} - {atm.pbottom/punits:.2e} '
    #   f'{atm.punits}.', indent=2)
    #log.msg(f'Number of model layers: {atm.nlayers-atm.rtop}.', indent=2)

    # Provide a summary of what happened here:
    log.msg(
        f'Provenance status for main atmospheric properties:\n'
        f'Pressure profile: {p_status}\n'
        f'Temperature profile: {t_status}\n'
        f'VMR profiles: {vmr_status}\n'
        f'Radius profile: {r_status}',
        indent=2,
    )
    log.msg(f"Species list:\n  {mol.name}", indent=2, si=4)
    min_p = atm.press[ 0] / pt.u(atm.punits)
    max_p = atm.press[-1] / pt.u(atm.punits)
    log.msg(
        f"Abundances are given by volume mixing ratio.\n"
        f"Unit factors: radius: {atm.runits}, pressure: {atm.punits}, "
        f"temperature: {atm.tunits}\n"
        f"Number of layers in atmospheric profile: {atm.nlayers}\n"
        f"Atmospheric pressure limits: {min_p:.2e}--{max_p:.2e} {atm.punits}."
        f"{mmm_text}",
        indent=2,
    )
    # TBD: Add extra bits/logs from makesample.make_atmosphere()

    # Return atmospheric model if requested:
    if atm.atmfile is not None:
        header = '# pyrat bay atmospheric model\n'
        io.write_atm(
            atm.atmfile, pressure, temperature, species,
            vmr, radius, atm.punits, atm.runits, header=header,
        )
        log.msg(f"Output atmospheric file: '{atm.atmfile}'.")

    return


def get_constants(pyrat):
    """Set molecules constant values (mass, radius)."""
    # Read file with molecular info:
    pyrat.log.msg(
        f"Taking species constant parameters from: '{pyrat.mol.molfile}'.",
        indent=2,
    )
    symbol, mass, radius = io.read_molecs(pyrat.mol.molfile)

    # Check that all atmospheric species are listed in molfile:
    absent = np.setdiff1d(pyrat.mol.name, symbol)
    if len(absent) > 0:
        pyrat.log.error(
            f"These species: {absent} are not listed in the molecules "
            f"info file: {pyrat.mol.molfile}."
        )

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
    phy = pyrat.phy
    ex = pyrat.ex
    lt = pyrat.lt

    # Temperature profile:
    if temp is not None:
        # Need to null tpars since it does not represent temp anymore
        atm.tpars = None
    elif atm.tpars is not None:
        temp = atm.tmodel(atm.tpars)
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
        atm.molpars = None

    elif atm.chemistry == 'tea':
        net = atm.chem_model
        metallicity = None
        e_abundances = {}
        e_ratio = {}
        e_scale = {}
        for model,var,val in zip(atm.molmodel, atm.molfree, atm.molpars):
            if model != 'equil':
                continue
            if var == 'metal':
                metallicity = val
            elif var.endswith('_metal'):
                var = var.rstrip('_metal')
                idx = list(net._base_composition).index(var)
                solar_abundance = net._base_dex_abundances[idx]
                e_abundances[var] = solar_abundance + val
            elif '_' in var:
                e_ratio[var] = val
            else:
                e_abundances[var] = val
        vmr = net.thermochemical_equilibrium(
            atm.temp,
            metallicity=metallicity,
            e_abundances=e_abundances,
            e_ratio=e_ratio,
            e_scale=e_scale,
        )
    # TBD: Check this is the right (best) criterion
    elif atm.molpars is not None: #atm.chemistry == 'uniform':
        vmr = pa.qscale(
            np.copy(atm.base_vmr),
            pyrat.mol.name, atm.molmodel,
            atm.molfree, atm.molpars, atm.bulk,
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
    atm.mm = np.sum(atm.vmr * pyrat.mol.mass, axis=1)

    # Radius profile (take input, re-compute, or keep current):
    if radius is not None:
        atm.radius = radius
    elif atm.rmodelname is not None:
        atm.radius = pyrat.hydro(
            atm.press, atm.temp, atm.mm,
            phy.gplanet, phy.mplanet,
            atm.refpressure, phy.rplanet,
        )
    else:
        pass

    # Number density (molecules cm-3):
    atm.d = pa.ideal_gas_density(atm.vmr, atm.press, atm.temp)

    # Check radii lie within Hill radius:
    phy.rhill = pa.hill_radius(phy.smaxis, phy.mplanet, phy.mstar)
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
    if pyrat.atm.ptop is not None and pyrat.atm.pbottom is not None:
        if pyrat.atm.pbottom <= pyrat.atm.ptop:
            pbottom = pyrat.atm.pbottom / pt.u(pyrat.atm.punits)
            ptop = pyrat.atm.ptop / pt.u(pyrat.atm.punits)
            pyrat.log.error(
               f'Bottom-layer pressure ({pbottom:.2e} {pyrat.atm.punits}) '
                'must be higher than the top-layer pressure '
               f'({ptop:.2e} {pyrat.atm.punits})'
            )

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
        'file (input_atmfile)'
    )


def check_temperature(pyrat):
    """
    Determine whether to calculate or read the atmospheric temperature
    profile.
    """
    # A model takes precedence:
    if pyrat.atm.tmodelname is not None:
        if pyrat.atm.tpars is None and pyrat.inputs.temperature is not None:
            return 'read'
        if pyrat.atm.tpars is None:
            pyrat.log.error('Undefined temperature-model parameters (tpars)')
        return 'calculate'

    # User-provided PT profile:
    if pyrat.inputs.temperature is not None:
        return 'read'

    pyrat.log.error(
        'Cannot compute temperature profile, either set a temperature model '
        '(tmodelname) and parameters (tpars), or provide an input PT '
        'profile (ptfile) or atmospheric file (input_atmfile)'
    )


def check_chemistry(pyrat, t_status):
    """
    Determine whether to calculate, read, or skip the atmospheric compostion
    profiles.
    """
    atm = pyrat.atm
    log = pyrat.log
    inputs = pyrat.inputs

    # No model:
    if atm.chemistry is None:
        if inputs.vmr is None:
            return 'skip'
        if inputs.species is not None:
            log.warning(
                "Composition will be taken from input atmospheric file. "
                "The input variable 'species' will be ingnored"
            )
        return 'read'

    if inputs.species is None:
        inputs.species = inputs.atm_species
    if inputs.species is None:
        log.error(
            'Cannot compute VMRs. Undefined atmospheric species list (species)'
        )

    # Uniform-abundances profile:
    if atm.chemistry == 'uniform':
        if inputs.uniform is None:
            log.error(
                'Undefined list of uniform volume mixing ratios '
                f'(uniform) for {atm.chemistry} chemistry model'
            )
        nuniform = len(inputs.uniform)
        nspecies = len(inputs.species)
        if nuniform != nspecies:
            log.error(
                f'Number of uniform abundances ({nuniform}) does '
                f'not match the number of species ({nspecies})'
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
