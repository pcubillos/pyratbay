# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Atmosphere',
]

import numpy as np
import scipy.interpolate as sip

from .. import atmosphere as pa
from ..atmosphere import vmr_models
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

    def __init__(self, inputs, log, mstar=None):
        """
        Initialize an atmospheric model object
        There are four main properties to compute (in this order):
        The pressure, the temperature, the volume mixing ratios (VMRs),
        and the radius profiles.

        Properties can be read from input profiles, computed by models,
        or skipped, depending on the configuration file.

        The rules are simple:
        - if there is a model in the config file, calculate the property
        - else if there is an input atmfile or ptfile, read properties from file
        - else, skip the calculation
        - if calculate p, any further reads (T,VMR,r) will interpolate
        """
        log.head('\nGenerating atmospheric model')

        self.rtop = 0  # Index of topmost layer (within Hill radius)

        # Check that input files exist:
        if inputs.molfile is None:
            self.molfile = pc.ROOT + 'pyratbay/data/molecules.dat'
        else:
            self.molfile = inputs.molfile
        with pt.log_error(log):
            pt.file_exists('molfile', 'Molecular-data', self.molfile)

        # Display units (internally, these variables are in CGS units):
        self.qunits = 'vmr'
        self.tunits = 'kelvin'
        self.runits = inputs.runits
        self.punits = inputs.punits

        self.tint = inputs.tint
        self.beta_irr = inputs.beta_irr
        self.rhill = np.inf
        self.smaxis = inputs.smaxis

        self.refpressure = inputs.refpressure
        self.tmodelname = inputs.tmodelname
        self.tpars = inputs.tpars

        self.chemistry = inputs.chemistry
        vmr_vars = inputs.vmr_vars
        if vmr_vars is None:
            vmr_vars = ''
        vmr_vars = [
            par for par in vmr_vars.splitlines()
            if par != ''
        ]
        # if any item is a number, then assume {vars,pars} pairs, one per line
        has_pars = np.any([
            pt.is_number(val)
            for vars in vmr_vars
            for val in vars.split()
        ])
        self.vmr_vars = []
        self.vmr_pars = []
        if has_pars:
            for vmr_var in vmr_vars:
                vmr_var = vmr_var.split()
                vmr_model = vmr_var[0]
                self.vmr_vars.append(vmr_model)
                if len(vmr_var) == 1:
                    error = f'Unspecified parameter value for {vmr_model}'
                    raise ValueError(error)
                pars = np.array(vmr_var[1:], float)
                self.vmr_pars.append(pars)
        else:
            vmr_vars = ' '.join(vmr_vars)
            self.vmr_vars = vmr_vars.split()
            self.vmr_pars = None

        #self.bulk = []
        #if inputs.bulk is not None:
        #    self.bulk = inputs.bulk
        self.bulk = inputs.bulk

        self.rmodelname = inputs.rmodelname

        self.rplanet = inputs.rplanet
        self.gplanet = inputs.gplanet
        self.mplanet = inputs.mplanet
        self.mass_units = inputs.mass_units

        self.output_atmfile = inputs.output_atmfile
        # User-provided PT profile:
        if pt.isfile(inputs.ptfile) == 1:
            input_atm_source = 'ptfile'
            self.atmfile = inputs.ptfile
        # Existing atmospheric file:
        elif pt.isfile(inputs.atmfile) == 1:
            input_atm_source = 'atmfile'
            self.atmfile = inputs.atmfile
        else:
            input_atm_source = None
            self.atmfile = None

        self.input_atm_source = input_atm_source

        # Input-atmosphere data:
        inputs.pressure = None
        inputs.temperature = None
        inputs.vmr = None
        inputs.radius = None
        inputs.atm_species = None
        if input_atm_source is not None:
            log.msg(f"\nReading atmospheric profile from: '{self.atmfile}'")
            input_atmosphere = check_input_atmosphere(self.atmfile, log)

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
                    # TBD: get masses for these species
                    #inputs.vmr /= mol.mass * np.sum(inputs.vmr/mol.mass, axis=1)
                    pass
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
                self.ptop, self.pbottom, self.nlayers, 'bar', log,
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


        # Composition / volume mixing ratio profiles:
        if vmr_status == 'skip':
            self.species = None
            self.vmr = None
        elif vmr_status == 'read':
            self.species = inputs.atm_species
            self.vmr = inputs.vmr
        elif vmr_status == 'calculate':
            chem_net = pa.chemistry(
                self.chemistry,
                pressure, temperature, inputs.species,
                solar_file=inputs.solar,
                log=log,
                punits=self.punits,
                q_uniform=inputs.uniform,
            )
            self.chem_model = chem_net
            self.species = chem_net.species
            self.vmr = np.copy(chem_net.vmr)

        self.base_vmr = None
        if self.vmr is not None:
            self.base_vmr = np.copy(self.vmr)

        self.parse_abundance_parameters(log)
        # Set species' mass and collision radius:
        if self.species is not None:
            self.nmol = len(self.species)
            log.msg(
                f"Read species physical properties from: '{self.molfile}'",
                indent=2,
            )
            mol_names, mol_mass, mol_radius = io.read_molecs(self.molfile)
            absent = np.setdiff1d(self.species, mol_names)
            if len(absent) > 0:
                log.error(
                    f"These species: {absent} are not listed in the molecules "
                    f"info file: {self.molfile}"
                )

            log.msg(
                'Molecule   Radius  Mass\n'
                '           (A)     (gr/mol)',
                indent=4,
            )
            self.mol_mass = np.zeros(self.nmol)
            self.mol_radius = np.zeros(self.nmol)
            for i in range(self.nmol):
                imol = list(mol_names).index(self.species[i])
                self.mol_mass[i] = mol_mass[imol]
                self.mol_radius[i] = mol_radius[imol] * pc.A
                log.msg(
                    f"{self.species[i]:>10s}:  "
                    f"{self.mol_radius[i]/pc.A:.3f}  "
                    f"{self.mol_mass[i]:8.4f}",
                    indent=2,
                )


        if r_status == 'skip':
            self.radius = None
        elif r_status == 'read':
            self.radius = inputs.radius
            if self.runits is None:
                self.runits = r_units

        # Planetary radius units (if not set by rplanet nor atmfile):
        if self.runits is None:
            if self.rplanet is not None:
                self.runits = 'rjup' if self.rplanet > 0.5*pc.rjup else 'rearth'
            else:
                self.runits = 'rearth'

        # Compute VMR and radius profiles (when needed)
        # and other properties (mean molecular mass, number density, Hill radius)
        self.calc_profiles(mstar=mstar, log=log)


        # Screen outputs:
        mmm_text = ''
        if self.vmr is not None:
            median_mmm = np.median(self.mm)
            mmm_text += (
                f"\nMedian mean molecular mass: {median_mmm:.3f} g mol-1."
            )

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

        log.msg(f"Species list:\n  {self.species}", indent=2, si=4)
        min_p = self.press[ 0] * pc.bar/pt.u(self.punits)
        max_p = self.press[-1] * pc.bar/pt.u(self.punits)
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

        # Save atmospheric model to file if requested:
        if self.output_atmfile is not None:
            header = '# pyrat bay atmospheric model\n'
            io.write_atm(
                self.output_atmfile, pressure, temperature, self.species,
                self.vmr, self.radius, self.punits, self.runits, header=header,
            )
            log.msg(f"Output atmospheric file: '{self.output_atmfile}'.")


    def calc_profiles(
            self, temp=None, vmr=None, radius=None, mstar=None, log=None,
            # Deprecated parameters:
            abund=None,
        ):
        """
        Update temperature, abundances, and radius profiles

        Parameters
        ----------
        temp: 1D float ndarray
            Layer's temperature profile (Kelvin) sorted from top to bottom.
        vmr: 2D float ndarray
            Species mole mixing ratio profiles [nlayers, nmol].
        radius: 1D float ndarray
            Layer's altitude profile (in cm), same order as temp.
        """
        # Temperature profile:
        if temp is not None:
            # Need to null tpars since it does not represent temp anymore
            self.tpars = None
        elif self.temp_model is not None and self.tpars is not None:
            temp = self.temp_model(self.tpars)
        else:
            temp = self.temp

        # Check that the dimensions match:
        if np.size(temp) != self.nlayers:
            log.error(
                f"The temperature array size ({np.size(temp)}) doesn't match "
                f"the Pyrat's temperature size ({np.size(self.temp)})"
            )
        self.temp = temp

        # Volume mixing ratios:
        if vmr is not None:
            self.molpars = []
            if np.shape(vmr) != np.shape(self.vmr):
                log.error(
                    f"The shape of the input VMR array {np.shape(vmr)} "
                     "doesn't match the shape of the Pyrat VMR "
                    f"{np.shape(self.vmr)}"
                )
        elif self.chemistry == 'tea':
            metallicity = None
            e_ratio = {}
            e_scale = {}
            for i,equil_model in enumerate(self.vmr_models):
                if self.vmr_pars is None or not self._is_equil_model[i]:
                    continue
                val = self.vmr_pars[i][0]
                if equil_model.name == 'metal_equil':
                    metallicity = val
                elif equil_model.name == 'scale_equil':
                    e_scale[equil_model.element] = val
                elif equil_model.name == 'ratio_equil':
                    e_ratio[equil_model.element_ratio] = val
            vmr = self.chem_model.thermochemical_equilibrium(
                self.temp,
                metallicity=metallicity,
                e_ratio=e_ratio,
                e_scale=e_scale,
            )
        elif np.any(~self._is_equil_model) and self.vmr_pars is not None:
            vmr_pars = [
                self.vmr_pars[i]
                for i,is_equil in enumerate(self._is_equil_model)
                if not is_equil
            ]
            vmr = pa.vmr_scale(
                self.base_vmr,
                self.species,
                self._free_models, vmr_pars, self.bulk,
                iscale=self.ifree, ibulk=self.ibulk,
                bratio=self.bulkratio, invsrat=self.invsrat,
            )
        elif self.base_vmr is None:
            vmr = self.base_vmr
        else:
            vmr = np.copy(self.base_vmr)
        self.vmr = vmr

        if self.vmr is None:
            return

        # Number density (molecules cm-3):
        self.d = pa.ideal_gas_density(self.vmr, self.press, self.temp)
        # Mean molecular mass:
        self.mm = pa.mean_weight(self.vmr, mass=self.mol_mass)

        # Radius profile:
        if radius is not None:
            self.radius = radius
        elif self.rmodelname is not None:
            self.radius = self.rad_model(
                self.press, self.temp, self.mm,
                self.mplanet, self.gplanet,
                self.refpressure, self.rplanet,
            )
        else:
            pass

        # Check radii lie within Hill radius:
        self.rhill = pa.hill_radius(self.smaxis, self.mplanet, mstar)
        if self.radius is not None:
            self.rtop = pt.ifirst(self.radius<self.rhill, default_ret=0)
            if self.rtop > 0:
                rhill = self.rhill / pt.u(self.runits)
                log.warning(
                    "The atmospheric pressure array extends beyond the Hill "
                    f"radius ({rhill:.5f} {self.runits}) at pressure "
                    f"{self.press[self.rtop]:.3e} bar (layer {self.rtop})."
                    "  Extinction beyond this layer will be neglected."
                )


    def rad_model(self, pressure, temperature, mu, mplanet, gplanet, p0, r0):
        """
        Calculate radius profile in hydrostatic equilibrium.

        Depending on self.rmodelname, select between a g(r)=GM/r**2
        (hydro_m) or constant-g (hydro_g) formula to compute
        the hydrostatic-equilibrium radii at each layer.

        Parameters
        ----------
        pressure: 1D float ndarray
            Atmospheric pressure for each layer (in bar).
        temperature: 1D float ndarray
            Atmospheric temperature for each layer (in K).
        mu: 1D float ndarray
            Mean molecular mass for each layer (in g mol-1).
        mplanet: Float
            Planetary mass in g (ignored for hydro_g).
        gplanet: Float
            Atmospheric gravity in cm s-2 (ignored for hydro_m).
        p0: Float
            Reference pressure level (in bar) where radius(p0) = r0.
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
        """
        Validate composition variables.
        """
        if molec is not None:
            if molec not in self.species:
                log.error(
                    f"Invalid vmr_vars variable '{var}', "
                    f"species {molec} is not in the atmosphere"
                )
            return

        in_equillibrium = self.chemistry == 'tea'
        if not in_equillibrium:
            log.error(f"vmr_vars variable '{var}' requires chemistry=tea")
        for element in elements:
            if element not in self.chem_model.elements:
                log.error(
                    f"Invalid vmr_vars variable '{var}', "
                    f"element '{element}' is not in the atmosphere"
                )


    def parse_abundance_parameters(self, log=None):
        """
        Sort out variables related to the VMR modeling
        """
        species = [] if self.species is None else list(self.species)

        # Sort out abundance free-parameters:
        self.vmr_models = []
        self.mol_pnames = []  # --> self.vmr_pnames
        self.mol_texnames = []  # --> self.vmr_texnames

        for i,var in enumerate(self.vmr_vars):
            # VMR variables
            if var.startswith('log_'):
                molec = var[4:]
                self.validate_species(var, log, molec=molec)
                vmr_model = vmr_models.IsoVMR(molec, self.press)
            elif var.startswith('scale_'):
                molec = var[6:]
                self.validate_species(var, log, molec=molec)
                imol = species.index(molec)
                vmr_model = vmr_models.ScaleVMR(
                    molec, self.press, self.vmr[:,imol],
                )
            elif var.startswith('slant_'):
                molec = var[6:]
                self.validate_species(var, log, molec=molec)
                vmr_model = vmr_models.SlantVMR(molec, self.press)
            # Equillibrium variables
            elif var == '[M/H]':
                vmr_model = vmr_models.MetalEquil()
                self.validate_species(var, log)
            elif var.startswith('[') and var.endswith('/H]'):
                vmr_model = vmr_models.ScaleEquil(var)
                self.validate_species(var, log, elements=[vmr_model.element])
            elif '/' in var:
                vmr_model = vmr_models.RatioEquil(var)
                self.validate_species(var, log, elements=vmr_model.elements)
            else:
                log.error(f"Unrecognized VMR model (vmr_vars): '{var}'")
            self.mol_pnames += vmr_model.pnames
            self.mol_texnames += vmr_model.texnames
            self.vmr_models.append(vmr_model)

        self._equil_models = [
            model for model in self.vmr_models if model.type == 'equil'
        ]
        self._free_models = [
            model for model in self.vmr_models if model.type == 'free'
        ]

        is_equil_model = [
            model.type == 'equil'
            for model in self.vmr_models
        ]
        self._is_equil_model = np.array(is_equil_model, dtype=bool)

        free_vmr = [
            model.species
            for model in self.vmr_models
            if model.type=='free'
        ]
        # TBD: trigger this error in tests
        vmr_uniques, vmr_counts = np.unique(free_vmr, return_counts=True)
        if np.any(vmr_counts>1):
            duplicates = vmr_uniques[vmr_counts>1]
            log.error(f'There are repeated species {duplicates} in vmr_vars')

        self.ifree = [species.index(mol) for mol in free_vmr]

        self.mol_npars = len(self.mol_pnames)
        # Allow null molpars for now, check later after retrieval_params
        if self.vmr_pars is not None:
            for i, vmr_pars in enumerate(self.vmr_pars):
                vmr_model = self.vmr_models[i]
                vmr_model_name = self.vmr_vars[i]
                npars = np.size(vmr_pars)
                if vmr_model.npars != npars:
                    log.error(
                        f"The parameters for model '{vmr_model_name}' ({npars}) "
                        "does not match the expected number of values "
                        f"({vmr_model.npars})"
                    )

        # Obtain abundance ratios between the bulk species:
        self.ibulk = None
        if self.bulk is not None:
            # Ckeck bulk species:
            missing = np.setdiff1d(self.bulk, species)
            if len(missing) > 0:
                log.error(
                    'These bulk species are not present in the '
                    f'atmosphere: {missing}'
                )

            free_vmr = self.species[self.ifree]
            bulk_free_species = np.intersect1d(self.bulk, free_vmr)
            if len(bulk_free_species) > 0:
                log.error(
                    'These species were marked as both bulk and '
                    f'variable-abundance: {bulk_free_species}'
                )

            self.ibulk = [species.index(mol) for mol in self.bulk]
            self.bulkratio, self.invsrat = pa.ratio(self.vmr, self.ibulk)



    def __str__(self):
        fmt = {'float': '{:.3e}'.format}
        fw = pt.Formatted_Write()
        press = self.press*pc.bar/pt.u(self.punits)
        fw.write('Atmospheric model information:')
        fw.write(
            f"Input atmospheric file name (atmfile): '{self.atmfile}'"
        )
        fw.write(
            f"Output atmospheric file name (output_atmfile): '{self.output_atmfile}'"
        )
        fw.write('Number of layers (nlayers): {:d}', self.nlayers)

        rplanet = None if self.rplanet is None else self.rplanet/pc.rjup
        mplanet = None if self.mplanet is None else self.mplanet/pc.mjup
        gplanet = self.gplanet
        smaxis = pt.none_div(self.smaxis, pc.au)
        rhill = pt.none_div(self.rhill, pc.rjup)
        fw.write('\nPlanetary radius (rplanet, Rjup): {:.3f}', rplanet)
        fw.write('Planetary mass (mplanet, Mjup): {:.3f}', mplanet)
        fw.write('Planetary surface gravity (gplanet, cm s-2): {:.1f}', gplanet)
        fw.write(
            'Planetary internal temperature (tint, K):  {:.1f}',
            self.tint,
        )
        fw.write('Planetary Hill radius (rhill, Rjup):  {:.3f}', rhill)
        fw.write('Orbital semi-major axis (smaxis, AU): {:.4f}', smaxis)

        fw.write('\nPressure display units (punits): {}', self.punits)
        fw.write('Pressure internal units: bar')
        fw.write('Pressure at top of atmosphere (ptop):        {:.2e} {}',
            self.ptop*pc.bar/pt.u(self.punits), self.punits)
        fw.write('Pressure at bottom of atmosphere (pbottom):  {:.2e} {}',
            self.pbottom*pc.bar/pt.u(self.punits), self.punits)
        if self.refpressure is None:
            ref_pressure = None
        else:
            ref_pressure = self.refpressure*pc.bar/pt.u(self.punits)
        fw.write(
            'Reference pressure at rplanet (refpressure): {:.2e} {}',
            ref_pressure, self.punits,
        )
        fw.write(
            'Pressure profile (press, {}):\n    {}',
            self.punits, press, fmt=fmt, edge=3,
        )

        fw.write('\nTemperature units (tunits): {}', self.tunits)
        fw.write('Temperature model name (tmodelname): {}', self.tmodelname)
        if self.temp_model is not None:
            fw.write('  tmodel parameters (tpars): {}', self.tpars)
        fw.write('Temperature profile (temp, K):\n    {}', self.temp,
            fmt={'float': '{:9.3f}'.format}, edge=3)

        fw.write('\nAbundance units (qunits): {}', self.qunits)
        fw.write('Abundance internal units: VMR')
        fw.write('Abundance model (chemistry): {}', self.chemistry)
        fw.write('Number of species (nmol): {:d}\n', self.nmol)
        fw.write(
          '\nIndex   Molecule  Mass (g/mol)    Radius (A)\n'
            '       (species)    (mol_mass)  (mol_radius)'
        )
        for i in range(self.nmol):
            fw.write(
                '{:>5d}  {:>9s}  {:12.5f}  {:12.3f}',
                i, self.species[i], self.mol_mass[i], self.mol_radius[i]/pc.A,
            )
        fw.write("Molecular data taken from (molfile): '{}'", self.molfile)

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

        fw.write('Abundance profiles (vmr):')
        for i, q in enumerate(self.vmr.T):
            fw.write('{:>16s}:   {}', self.species[i], q, fmt=fmt, edge=2)
        fw.write('Density profiles (d, molecules cm-3):')
        for i, dens in enumerate(self.d.T):
            fw.write('{:>16s}:   {}', self.species[i], dens, fmt=fmt, edge=2)

        fw.write('\nRadius display units (runits): {}', self.runits)
        fw.write('Radius internal units: cm', self.runits)
        fw.write('Radius model name (rmodelname): {}', self.rmodelname)
        if self.radius is not None:
            fw.write(
                'Radius profile (radius, {}):\n    {}',
                self.runits, self.radius/pt.u(self.runits),
                prec=4, edge=3, lw=800,
            )


        fw.write('\nMean molecular mass (mm, amu):\n    {}', self.mm,
            fmt={'float': '{:8.4f}'.format}, edge=3)

        return fw.text



def check_input_atmosphere(atm_file, log):
    """
    Make sure that the input atmospheric profile is sorted correctly
    (from low pressure to high pressure).
    """
    input_atm = io.read_atm(atm_file)
    p_units, t_units, vmr_units, r_units = units = input_atm[0]
    species = input_atm[1]
    pressure = input_atm[2]
    if p_units != 'bar':
        pressure *= pt.u(p_units) / pc.bar
    temperature = input_atm[3]
    vmr = input_atm[4]
    radius = input_atm[5]
    if radius is not None:
        radius *= pt.u(r_units)

    # Check that the layers are sorted from top to bottom:
    sort = np.all(np.ediff1d(pressure) > 0)
    reverse = np.all(np.ediff1d(pressure) < 0)
    if radius is not None:
        sort &= np.all(np.ediff1d(radius) < 0)
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
            pbottom = inputs.pbottom * pc.bar / pt.u(inputs.punits)
            ptop = inputs.ptop * pc.bar / pt.u(inputs.punits)
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
        'file (atmfile)'
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
        'profile (ptfile) or atmospheric file (atmfile)'
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

