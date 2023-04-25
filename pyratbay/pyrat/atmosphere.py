# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Atmosphere',
]

import numpy as np
import scipy.interpolate as sip

from .. import atmosphere as pa
from .. import constants as pc
from .. import io as io
from .. import tools as pt


class MassGravity():
    """
    Descriptor object that keeps the planet mass and gravity
    consistent with each other by ensuring that
        gplanet = G * mplanet / rplanet**2
    whenever one of these two variables are modified.

    To understand this sorcery see:
    https://docs.python.org/3/howto/descriptor.html
    """
    def __set_name__(self, obj, name):
        self.private_name = '_' + name

    def __get__(self, obj, objtype=None):
        value = getattr(obj, self.private_name)
        return value

    def __set__(self, obj, value):
        priv_name = self.private_name
        var_name = self.private_name[1:]
        if hasattr(obj, priv_name) and value == getattr(obj, priv_name):
            return
        setattr(obj, priv_name, value)

        if obj.rplanet is not None and value is not None:
            if var_name == 'mplanet':
                obj.gplanet = pc.G * obj.mplanet / obj.rplanet**2
            elif var_name == 'gplanet':
                obj.mplanet = obj.gplanet * obj.rplanet**2 / pc.G


class Atmosphere():
    mplanet = MassGravity()  # Planetary mass
    gplanet = MassGravity()  # Planetary surface gravity (at rplanet)

    def __init__(self, inputs, log):
        """
        Initialize atmosphere.
        """
        self.rtop = 0  # Index of topmost layer (within Hill radius)

        # Check that input files exist:
        if inputs.molfile is None:
            self.molfile = pc.ROOT + 'pyratbay/data/molecules.dat'
        else:
            self.molfile = inputs.molfile
        with pt.log_error(log):
            pt.file_exists('molfile', 'Molecular-data', self.molfile)

        self.atmfile = inputs.atmfile
        # Display units (internally, these variables are in CGS units):
        self.qunits = 'vmr'
        self.tunits = 'kelvin'
        self.runits = inputs.runits
        self.punits = inputs.punits

        self.refpressure = inputs.refpressure
        self.tmodelname = inputs.tmodelname
        self.tpars = inputs.tpars

        self.chemistry = inputs.chemistry
        self.metallicity = inputs.metallicity
        self.e_scale = inputs.e_scale
        self.e_ratio = inputs.e_ratio

        #self.bulk = []
        #if inputs.bulk is not None:
        #    self.bulk = inputs.bulk
        self.bulk = inputs.bulk

        self.rmodelname = inputs.rmodelname

        self.rplanet = inputs.rplanet
        self.gplanet = inputs.gplanet
        self.mplanet = inputs.mplanet
        self.mass_units = inputs.mass_units

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
        log.head('\nGenerating atmospheric model')

        # User-provided PT profile:
        if pt.isfile(inputs.ptfile) == 1:
            input_atm_source = 'ptfile'
            self.input_atmfile = inputs.ptfile
        # Existing atmospheric file:
        elif pt.isfile(inputs.input_atmfile) == 1:
            input_atm_source = 'atmfile'
            self.input_atmfile = inputs.input_atmfile
        else:
            input_atm_source = None

        self.input_atm_source = input_atm_source

        # Input-atmosphere data:
        inputs.pressure = None
        inputs.temperature = None
        inputs.vmr = None
        inputs.radius = None
        inputs.atm_species = None
        if input_atm_source is not None:
            log.msg(f"\nReading atmospheric profile from: '{self.input_atmfile}'")
            input_atmosphere = check_input_atmosphere(self.input_atmfile, log)

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
        p_status = check_pressure(inputs, log)
        t_status = check_temperature(inputs, log)
        vmr_status = check_chemistry(inputs, log, t_status)
        r_status = check_altitude(inputs, log, vmr_status)

        log.msg(
            f'Provenance status for main atmospheric properties:\n'
            f'  Pressure profile: {p_status}\n'
            f'  Temperature profile: {t_status}\n'
            f'  VMR profiles: {vmr_status}\n'
            f'  Radius profile: {r_status}',
        )

        # Pressure profile:
        if p_status == 'calculate':
            self.ptop = inputs.ptop
            self.pbottom = inputs.pbottom
            self.nlayers = inputs.nlayers
            pressure = pa.pressure(
                self.ptop, self.pbottom, self.nlayers, 'barye', log,
            )
        elif p_status == 'read':
            pressure = inputs.pressure
            if self.punits is None:
                self.punits = p_units
            # TBD: Do I want to update with input ptop/pbottom?
            self.ptop = np.amin(pressure)
            self.pbottom = np.amax(pressure)
            self.nlayers = len(pressure)
        self.press = pressure


        if p_status == 'calculate' and 'read' in [t_status,vmr_status,r_status]:
            # Interpolate if needed:
            logp_input = np.log(inputs.pressure)
            if t_status == 'read':
                temp_input = np.copy(inputs.temperature)
                temp_extrap = inputs.temperature[0], inputs.temperature[-1]
                temp_interp = sip.interp1d(
                    logp_input, temp_input,
                    kind='slinear',
                    bounds_error=False, fill_value=temp_extrap,
                )
                inputs.temperature = temp_interp(np.log(pressure))
            if vmr_status == 'read' and inputs.vmr is not None:
                log_vmr_input = np.log(inputs.vmr)
                vmr_extrap = np.log(inputs.vmr[0]), np.log(inputs.vmr[-1])
                vmr_interp = sip.interp1d(
                    logp_input, log_vmr_input, axis=0,
                    kind='slinear',
                    bounds_error=False, fill_value=vmr_extrap,
                )
                inputs.vmr = np.exp(vmr_interp(np.log(pressure)))
            if r_status == 'read' and inputs.radius is not None:
                rad_interp = sip.interp1d(
                    logp_input, inputs.radius, kind='slinear',
                )
                inputs.radius = rad_interp(np.log(pressure))

        # Temperature profile:
        if self.tmodelname in pc.tmodels:
            self.temp_model = pa.tmodels.get_model(
                self.tmodelname,
                pressure=self.press,
                nlayers=self.nlayers,
            )
        else:
            self.temp_model = None

        if t_status == 'calculate':
            temperature = self.temp_model(self.tpars)
        elif t_status == 'read':
            temperature = inputs.temperature
        self.temp = temperature

        # Composition (volume-mixing-ratio) profiles:
        species = None
        radius = None
        if vmr_status == 'calculate':
            chem_net = pa.chemistry(
                self.chemistry,
                pressure, temperature, inputs.species,
                metallicity=self.metallicity,
                #e_abundances=self.e_abundances,
                e_scale=self.e_scale,
                e_ratio=self.e_ratio,
                solar_file=inputs.solar,
                log=log,
                punits=self.punits,
                q_uniform=inputs.uniform,
            )
            self.chem_model = chem_net
            vmr = chem_net.vmr
            species = chem_net.species
        elif vmr_status == 'read':
            species = inputs.atm_species
            vmr = inputs.vmr
        elif vmr_status == 'skip':
            vmr = None
        self.vmr = vmr

        # Set values of species properties:
        self.species = species
        if species is not None:
            self.nmol = len(species)
            log.msg(
                f"Read species physical properties from: '{self.molfile}'",
                indent=2,
            )
            mol_names, mol_mass, mol_radius = io.read_molecs(self.molfile)

            # Check that all atmospheric species are listed in molfile:
            absent = np.setdiff1d(species, mol_names)
            if len(absent) > 0:
                log.error(
                    f"These species: {absent} are not listed in the molecules "
                    f"info file: {self.molfile}"
                )

            # Set molecule's values:
            self.mol_mass = np.zeros(self.nmol)
            self.mol_radius = np.zeros(self.nmol)

            log.msg(
                'Molecule   Radius  Mass\n'
                '           (A)     (gr/mol)',
                indent=4,
            )
            for i in range(self.nmol):
                # Find the molecule in the list:
                imol = list(mol_names).index(self.species[i])
                # Set molecule name, mass, and collision radius:
                self.mol_mass[i] = mol_mass[imol]
                self.mol_radius[i] = mol_radius[imol] * pc.A
                log.msg(
                    f"{self.species[i]:>10s}:  "
                    f"{self.mol_radius[i]/pc.A:.3f}  "
                    f"{self.mol_mass[i]:8.4f}",
                    indent=2,
                )


        # Radius profile:
        if r_status == 'calculate':
            # Mean molecular mass:
            mean_mass = pa.mean_weight(self.vmr, species)
            # Altitude profile:
            radius = self.rad_model(
                pressure, temperature, mean_mass, self.mplanet,
                self.gplanet, self.refpressure, self.rplanet,
            )
        elif r_status == 'read':
            radius = inputs.radius
            if self.runits is None:
                self.runits = r_units
        self.radius = radius
        # Planetary radius units (if not set by rplanet nor atmfile):
        if self.runits is None:
            if self.rplanet is not None:
                self.runits = 'rjup' if self.rplanet > 0.5*pc.rjup else 'rearth'
            else:
                self.runits = 'rearth'

        # Mean molecular mass and number densities:
        mmm_text = ''
        if self.vmr is not None:
            self.mm = np.sum(self.vmr*self.mol_mass, axis=1)
            # Number density profiles for each molecule (in molecules cm-3):
            self.d = pa.ideal_gas_density(self.vmr, self.press, self.temp)
            median_mmm = np.median(self.mm)
            mmm_text += (
                f"\nMedian mean molecular mass: {median_mmm:.3f} g mol-1."
            )
            # Base abundance profiles:
            self.base_vmr = np.copy(self.vmr)

        # Print radius array:
        if self.radius is not None:
            radius_arr = self.radius / pt.u(self.runits)
            log.msg(
                f'Radius array ({self.runits}) = \n{radius_arr}',
                indent=2, si=4,
            )
            log.msg(
                f'Upper/lower radius boundaries: {radius_arr[self.rtop]:.5f}'
                f'--{radius_arr[-1]:.5f} {self.runits}.',
                indent=2,
            )

        # Provide a summary of what happened here:
        log.msg(f"Species list:\n  {self.species}", indent=2, si=4)
        min_p = self.press[ 0] / pt.u(self.punits)
        max_p = self.press[-1] / pt.u(self.punits)
        log.msg(
            f"Abundances are given by volume mixing ratio.\n"
            f"Unit factors: radius: {self.runits}, pressure: {self.punits}, "
            f"temperature: {self.tunits}\n"
            f"Number of layers in atmospheric profile: {self.nlayers}\n"
            f"Atmospheric pressure limits: {min_p:.2e}--{max_p:.2e} "
            f"{self.punits}."
            f"{mmm_text}",
            indent=2,
        )
        # TBD: Add extra bits/logs from makesample.make_atmosphere()

        # Return atmospheric model if requested:
        if self.atmfile is not None:
            header = '# pyrat bay atmospheric model\n'
            io.write_atm(
                self.atmfile, pressure, temperature, species,
                vmr, radius, self.punits, self.runits, header=header,
            )
            log.msg(f"Output atmospheric file: '{self.atmfile}'.")


    def rad_model(self, pressure, temperature, mu, mplanet, gplanet, p0, r0):
        """
        Calculate radius profile in hydrostatic equilibrium.

        Depending on self.rmodelname, select between a g(r)=GM/r**2
        (hydro_m) or constant-g (hydro_g) formula to compute
        the hydrostatic-equilibrium radii at each layer.

        Parameters
        ----------
        pressure: 1D float ndarray
            Atmospheric pressure for each layer (in barye).
        temperature: 1D float ndarray
            Atmospheric temperature for each layer (in K).
        mu: 1D float ndarray
            Mean molecular mass for each layer (in g mol-1).
        mplanet: Float
            Planetary mass in g (ignored for hydro_g).
        gplanet: Float
            Atmospheric gravity in cm s-2 (ignored for hydro_m).
        p0: Float
            Reference pressure level (in barye) where radius(p0) = r0.
        r0: Float
            Reference radius level (in cm) corresponding to p0.
        """
        if self.rmodelname is None:
            return None
        # H.E. with  g(r) = GM/r**2:
        elif self.rmodelname == 'hydro_m':
            return pa.hydro_m(pressure, temperature, mu, mplanet, p0, r0)
        # H.E. with constant g:
        elif self.rmodelname == 'hydro_g':
            return pa.hydro_g(pressure, temperature, mu, gplanet, p0, r0)


    def validate_species(self, var, log, molec=None, elements=[]):
        # Validate composition variables:
        if molec is not None:
            if molec not in self.species:
                log.error(
                    f"Invalid molvars variable '{var}', "
                    f"species {molec} is not in the atmosphere"
                )
            return

        in_equillibrium = self.chemistry == 'tea'
        if not in_equillibrium:
            log.error(f"molvars variable '{var}' requires chemistry=tea")
        for element in elements:
            if element not in self.chem_model.elements:
                log.error(
                    f"Invalid molvars variable '{var}', "
                    f"element '{element}' is not in the atmosphere"
                )


    def parse_abundance_parameters(self, molvars, log=None):
        # Sort out abundance free-parameters:
        self.mol_pnames = []
        self.mol_texnames = []
        self.mol_npars = len(molvars)

        free_vmr = []
        self._equil_var = np.zeros(self.mol_npars, bool)
        for i,var in enumerate(molvars):
            # VMR variables
            if var.startswith('log_'):
                molec = var[4:]
                free_vmr.append(molec)
                self.validate_species(var, log, molec=molec)
                self.mol_pnames.append(var)
                self.mol_texnames.append(fr'$\log\ X_{{\rm {molec}}}$')
            elif var.startswith('scale_'):
                molec = var[6:]
                free_vmr.append(molec)
                self.validate_species(var, log, molec=molec)
                self.mol_pnames.append(var)
                self.mol_texnames.append(fr'$\log\ X_{{\rm {molec}}}$')
            # Equillibrium variables
            elif var == 'metal':
                self._equil_var[i] = True
                self.validate_species(var, log)
                self.mol_pnames.append(var)
                self.mol_texnames.append('[M/H]')
            elif var.startswith('[') and var.endswith('/H]'):
                self._equil_var[i] = True
                self.validate_species(var, log, elements=[var[1:-3]])
                self.mol_pnames.append(var)
                self.mol_texnames.append(var)
            elif '/' in var:
                self._equil_var[i] = True
                idx = var.index('/')
                elements = [var[0:idx], var[idx+1:]]
                self.validate_species(var, log, elements=elements)
                self.mol_pnames.append(var)
                self.mol_texnames.append(var)
            else:
                log.error(f"Unrecognized molvars variable name: '{var}'")


        self.ifree = [
            list(self.species).index(mol)
            for mol in free_vmr
        ]

        if len(molvars) == 0:
            self.ibulk = None

    def __str__(self):
        fmt = {'float': '{: .3e}'.format}
        fw = pt.Formatted_Write()
        press  = self.press/pt.u(self.punits)
        radius = self.radius/pt.u(self.runits)
        fw.write('Atmospheric model information:')
        fw.write(
            f"Input atmospheric file name (input_atmfile): '{self.input_atmfile}'"
        )
        fw.write("Output atmospheric file name (atmfile): '{}'", self.atmfile)
        fw.write('Number of layers (nlayers): {:d}', self.nlayers)

        fw.write('\nPressure display units (punits): {}', self.punits)
        fw.write('Pressure internal units: barye')
        fw.write('Pressure at top of atmosphere (ptop):        {:.2e} {}',
            self.ptop/pt.u(self.punits), self.punits)
        fw.write('Pressure at bottom of atmosphere (pbottom):  {:.2e} {}',
            self.pbottom/pt.u(self.punits), self.punits)
        fw.write('Reference pressure at rplanet (refpressure): {:.2e} {}',
            self.refpressure/pt.u(self.punits), self.punits)
        fw.write(
            'Pressure profile (press, {}):\n    {}',
            self.punits, press, fmt=fmt, edge=3,
        )

        fw.write('Atmospheric species information:')
        fw.write('Number of species (nmol): {:d}\n', self.nmol)
        fw.write(
            '\nMolecule    Mass       Radius\n'
            '            g/mol      Angstrom\n'
            '(name)      (mass)     (radius)  '
        )
        for i in range(self.nmol):
            fw.write(
                '  {:8s}  {:8.4f}  {:10.3f}',
                self.species[i], self.mol_mass[i], self.mol_radius[i]/pc.A,
            )
        fw.write("Molecular data taken from (molfile): '{}'", self.molfile)

        rplanet = None if self.rplanet is None else self.rplanet/pc.rjup
        mplanet = None if self.mplanet is None else self.mplanet/pc.mjup
        gplanet = self.gplanet
        fw.write('\nPlanetary radius (rplanet, Rjup): {:.3f}', rplanet)
        fw.write('Planetary mass (mplanet, Mjup): {:.3f}', mplanet)
        fw.write('Planetary surface gravity (gplanet, cm s-2): {:.1f}', gplanet)

        fw.write('\nRadius display units (runits): {}', self.runits)
        fw.write('Radius internal units: cm', self.runits)
        fw.write('Radius model name (rmodelname): {}', self.rmodelname)
        fw.write('Radius profile (radius, {}):\n    {}', self.runits, radius,
            prec=4, edge=3, lw=800)

        fw.write('\nTemperature units (tunits): {}', self.tunits)
        fw.write('Temperature model name (tmodelname): {}', self.tmodelname)
        if self.temp_model is not None:
            fw.write('  tmodel parameters (tpars): {}', self.tpars)
        fw.write('Temperature profile (temp, K):\n    {}', self.temp,
            fmt={'float': '{:9.3f}'.format}, edge=3)


        fw.write('\nMean molecular mass (mm, amu):\n    {}', self.mm,
            fmt={'float': '{:8.4f}'.format}, edge=3)
        fw.write('\nAbundance units (qunits): {}', self.qunits)
        fw.write('Abundance internal units: mole mixing fraction')
        fw.write('Number of atmospheric species: {:d}', len(self.vmr[0]))
        if len(self.mol_pnames) > 0:
            molpars = self.molpars
            if self.molpars == []:
                molpars = [None for _ in self.mol_pnames]

            fw.write('Abundance models:\n  molvars    molpars  ifree')
            for var, val in zip(self.mol_pnames, molpars):
                fw.write(f'  {var:15s}  {val:10s}')
            fw.write('Bulk species:\n  ibulk  bulk')
            for ibulk, bulk in zip(self.ibulk, self.bulk):
                fw.write('     {:2d}  {:10s}', ibulk, bulk)

        fw.write('Abundance profiles (vmr, mole mixing fraction):')
        for i, q in enumerate(self.vmr.T):
            fw.write('    species [{:2d}]:   {}', i, q,    fmt=fmt, edge=2)
        fw.write('Density profiles (d, molecules cm-3):')
        for i, dens in enumerate(self.d.T):
            fw.write('    species [{:2d}]:   {}', i, dens, fmt=fmt, edge=2)

        return fw.text



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


def check_pressure(inputs, log):
    """
    Determine whether to calculate or read the atmospheric pressure
    profile.
    """
    if inputs.ptop is not None and inputs.pbottom is not None:
        if inputs.pbottom <= inputs.ptop:
            pbottom = inputs.pbottom / pt.u(inputs.punits)
            ptop = inputs.ptop / pt.u(inputs.punits)
            log.error(
               f'Bottom-layer pressure ({pbottom:.2e} {inputs.punits}) '
                'must be higher than the top-layer pressure '
               f'({ptop:.2e} {inputs.punits})'
            )

    all_parameters_defined = (
        inputs.nlayers is not None and
        inputs.ptop is not None and
        inputs.pbottom is not None
    )
    if all_parameters_defined:
        return 'calculate'

    # User-provided PT profile:
    if inputs.pressure is not None:
        return 'read'

    log.error(
        'Cannot compute pressure profile, either set {ptop, pbottom, nlayers} '
        'parameters, or provide an input PT profile (ptfile) or atmospheric '
        'file (input_atmfile)'
    )


def check_temperature(inputs, log):
    """
    Determine whether to calculate or read the atmospheric temperature
    profile.
    """
    # A model takes precedence:
    if inputs.tmodelname is not None:
        if inputs.tpars is None and inputs.temperature is not None:
            return 'read'
        if inputs.tpars is None:
            log.error('Undefined temperature-model parameters (tpars)')
        return 'calculate'

    # User-provided PT profile:
    if inputs.temperature is not None:
        return 'read'

    log.error(
        'Cannot compute temperature profile, either set a temperature model '
        '(tmodelname) and parameters (tpars), or provide an input PT '
        'profile (ptfile) or atmospheric file (input_atmfile)'
    )


def check_chemistry(inputs, log, t_status):
    """
    Determine whether to calculate, read, or skip the atmospheric compostion
    profiles.
    """
    # No model:
    if inputs.chemistry is None:
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
    if inputs.chemistry == 'uniform':
        if inputs.uniform is None:
            log.error(
                'Undefined list of uniform volume mixing ratios '
                f'(uniform) for {inputs.chemistry} chemistry model'
            )
        nuniform = len(inputs.uniform)
        nspecies = len(inputs.species)
        if nuniform != nspecies:
            log.error(
                f'Number of uniform abundances ({nuniform}) does '
                f'not match the number of species ({nspecies})'
            )

    return 'calculate'


def check_altitude(inputs, log, vmr_status):
    """
    Determine whether to calculate, read, or skip the atmospheric
    altitude profile.
    """
    if inputs.rmodelname is None:
        if vmr_status == 'read' and inputs.radius is not None:
            return 'read'
        return 'skip'

    if vmr_status == 'skip':
        log.error(
            'Cannot compute hydrostatic-equilibrium radius profile.\n'
            'radius model needs to know the composition'
        )

    err = []
    if inputs.rplanet is None:
        err += ['Undefined planet radius (rplanet).']
    if inputs.mplanet is None and inputs.gplanet is None:
        err += ['Undefined planet mass (mplanet) or surface gravity (gplanet).']
    if inputs.refpressure is None:
        err += ['Undefined reference pressure level (refpressure).']

    if len(err) != 0:
        error_message = '\n'.join(err)
        log.error(
            'Cannot compute hydrostatic-equilibrium radius profile.\n'
            f'{error_message}'
        )

    return 'calculate'

