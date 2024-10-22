# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import numpy as np

from .. import opacity as op
from .. import tools as pt
from .line_by_line import Line_By_Line


def check_species_exists(species, atm_species, opa_model, log):
    """
    Check when a required species is not found in list of atmospheric species
    """
    absent = np.setdiff1d(species, atm_species)
    if len(absent) > 0:
        log.error(
            f'Species {absent}, required for opacity model {opa_model}, '
            'are not present in the atmosphere'
        )


class Opacity():
    """Interface between opacity models and pyrat object"""
    def __init__(self, inputs, wn, species, pressure, log, pyrat):
        """
        Read opacity sources.
        """
        self.models = []
        self.models_type = []
        self.mol_indices = []
        self.pnames = []
        # Min/max temperatures that can be sampled
        self.tmin = {}
        self.tmax = {}
        # extinction coefficient in cm-1
        self.nwave = len(wn)
        self.nlayers = len(pressure)
        self.ec = np.zeros((self.nlayers, self.nwave))
        self.ec_cloud = np.zeros((self.nlayers, self.nwave))

        self.fpatchy = inputs.fpatchy
        # The flag self.is_patchy will be set in pyrat.Retrieval() because
        # a patchy model may be set via the retrieval_params input

        min_wn = np.amin(wn)
        max_wn = np.amax(wn)
        species = list(species)

        self.nspec = []
        # TBD: without runmode?
        if inputs.extfile is not None and inputs.runmode != 'opacity':
            log.head("\nReading cross-section table file(s):")
            #for cs_file in self.extfile:
            #    log.head(f"  '{cs_file}'.")

            make_temp_array = (
                inputs.tmin is not None and
                inputs.tmax is not None and
                inputs.tstep is not None
            )
            if make_temp_array:
                # Instead of np.arange() to include inputs.tmax
                ntemp = int((inputs.tmax-inputs.tmin)/inputs.tstep) + 1
                tmax = inputs.tmin + (ntemp-1)*inputs.tstep
                temp_array = np.linspace(inputs.tmin, tmax, ntemp)
            else:
                temp_array = None

            # TBD: self.ls_files?
            ls = op.Line_Sample(
                inputs.extfile, pressure=pressure, temperature=temp_array,
                min_wn=min_wn, max_wn=max_wn, wn_thinning=inputs.wn_thinning,
                log=log,
            )
            self.models.append(ls)
            self.models_type.append('line_sample')
            check_species_exists(ls.species, species, ls.name, log)
            imol = [species.index(mol) for mol in ls.species]
            self.mol_indices.append(imol)
            self.tmin['line_sample'] = ls.tmin
            self.tmax['line_sample'] = ls.tmax

            log.msg(
                f"Line-sample opacity files have {ls.nspec} species, "
                f"{ls.ntemp} temperature samples, {ls.nlayers} layers, "
                f"and {ls.nwave} wavenumber samples.",
            )
            with np.printoptions(precision=1):
                str_temp = str(ls.temp)
            with np.printoptions(formatter={'float':'{:.2f}'.format}):
                str_wn = str(ls.wn)
            with np.printoptions(formatter={'float':'{:.3e}'.format}):
                str_press = str(ls.press)
            log.msg(
                f"Species names: {ls.species}\n"
                f"Temperatures (K):\n   {str_temp}\n"
                f"Pressure layers (bar):\n{str_press}\n"
                f"Wavenumber array (cm-1):\n   {str_wn}",
            )
            self.nspec.append(ls.nspec)
            self.pnames.append([])

        if inputs.tlifile is not None:
            lbl = Line_By_Line(inputs, species, min_wn, max_wn, log, pyrat)
            self.models.append(lbl)
            self.models_type.append('lbl')
            imol = [species.index(mol) for mol in lbl.species]
            self.mol_indices.append(imol)
            self.tmin['lbl'] = lbl.tmin
            self.tmax['lbl'] = lbl.tmax
            self.nspec.append(lbl.nspec)
            self.pnames.append([])

        if inputs.alkali_models is not None:
            for name in inputs.alkali_models:
                cutoff = inputs.alkali_cutoff
                model = op.alkali.get_model(name, pressure, wn=wn, cutoff=cutoff)
                self.models.append(model)
                self.models_type.append('alkali')
                check_species_exists(model.species, species, model.name, log)
                imol = species.index(model.species)
                self.mol_indices.append(imol)
                self.nspec.append(1)
                self.pnames.append([])

        if inputs.cia_files is not None:
            tmin = []
            tmax = []
            for cia_file in inputs.cia_files:
                log.head(f"Read CIA file: '{cia_file}'.", indent=2)
                cia = op.Collision_Induced(cia_file, wn=wn)
                log.msg(
                    f'{cia.name} opacity:\n'
                    f'Read {cia.nwave} wave and {cia.ntemp} temperature samples.\n'
                    f'Temperature ranges: {cia.tmin:.1f}--{cia.tmax:.1f} K\n'
                    f'Wavenumber ranges: {cia.wn[cia._wn_lo_idx]:.1f}--'
                    f'{cia.wn[cia._wn_hi_idx-1]:.1f} cm-1',
                    indent=4,
                )

                self.models.append(cia)
                self.models_type.append('cia')
                tmin.append(cia.tmin)
                tmax.append(cia.tmax)
                check_species_exists(cia.species, species, cia.name, log)
                imol = [species.index(mol) for mol in cia.species]
                self.mol_indices.append(imol)
                self.nspec.append(1)
                self.pnames.append([])
            self.tmin['cia'] = np.amax(tmin)
            self.tmax['cia'] = np.amin(tmax)

        if inputs.rayleigh is not None:
            npars = 0
            for name in inputs.rayleigh:
                if name.startswith('dalgarno_'):
                    mol = name.split('_')[1]
                    model = op.rayleigh.Dalgarno(wn, mol)
                    self.models_type.append('rayleigh')
                if name == 'lecavelier':
                    model = op.rayleigh.Lecavelier(wn)
                    self.models_type.append('cloud')
                self.models.append(model)
                self.nspec.append(1)
                self.pnames.append(model.pnames)

                check_species_exists(model.species, species, model.name, log)
                self.mol_indices.append(species.index(model.species))
                # Parse parameters:
                if inputs.rpars is None:
                    model.pars = np.tile(np.nan, model.npars)
                elif len(inputs.rpars) >= npars+model.npars:
                    model.pars = inputs.rpars[npars:npars+model.npars]
                npars += model.npars

            n_input = npars if inputs.rpars is None else len(inputs.rpars)
            if n_input != npars:
                log.error(
                    f'Number of input Rayleigh parameters ({n_input}) does not '
                    f'match the number of required model parameters ({npars})'
                )

        if inputs.clouds is not None:
            npars = 0
            for name in inputs.clouds:
                if name == 'ccsgray':
                    model = op.clouds.CCSgray(pressure, wn)
                if name == 'deck':
                    model = op.clouds.Deck(pressure, wn)
                self.models.append(model)
                self.models_type.append('cloud')
                self.mol_indices.append(None)
                self.nspec.append(1)
                self.pnames.append(model.pnames)

                # Parse parameters:
                if inputs.cpars is None:
                    model.pars = np.tile(np.nan, model.npars)
                elif len(inputs.cpars) >= npars+model.npars:
                    model.pars = inputs.cpars[npars:npars+model.npars]
                npars += model.npars

            n_input = npars if inputs.cpars is None else len(inputs.cpars)
            if n_input != npars:
                log.error(
                    f'Number of input cloud parameters ({n_input}) does not '
                    f'match the number of required model parameters ({npars})'
                )

        if inputs.h_ion is not None:
            model = op.Hydrogen_Ion(wn)
            self.models.append(model)
            self.models_type.append(model.name)

            check_species_exists(model.species, species, model.name, log)
            # For calculations only H and e- are necessary:
            imol = [species.index(mol) for mol in model.species]
            self.mol_indices.append(imol)
            self.nspec.append(1)
            self.pnames.append([])


    def calc_extinction_coefficient(
        self, temperature, radius, densities, skip=[],
    ):
        """
        Compute extinction coefficient over temperature and
        number density profiles.

        It is possible to neglect specific opacities using the skip
        argument, skipped sources can be tagged by type, name, or
        species (the latter for line_sample and lbl sources)
        """
        self.ec[:] = 0.0
        self.ec_cloud[:] = 0.0
        for i,model in enumerate(self.models):
            imol = self.mol_indices[i]
            model_type = self.models_type[i]
            density = densities[:,imol]
            args = {}

            # Skip opacities
            if model.name in skip or model_type in skip:
                if model.name == 'deck':
                    model.itop = self.nlayers - 1
                    model.tsurf = temperature[-1]
                    model.rsurf = radius[-1]
                continue
            if model_type == 'line_sample':
                intersect = np.isin(model.species, skip)
                if np.any(intersect) > 0:
                    density[:,intersect] = 0.0
            if model_type == 'lbl':
                args['skip_mol'] = skip

            if model.name == 'deck':
                args['temperature'] = temperature
                args['radius'] = radius
            elif model.name == 'ccsgray':
                args['temperature'] = temperature
            elif model_type in ['rayleigh', 'cloud']:
                args['density'] = density
            else:
                args['temperature'] = temperature
                args['density'] = density

            if hasattr(model, 'pars') and model.npars > 0:
                args['pars'] = model.pars
            if model_type == 'cloud' and self.is_patchy:
                self.ec_cloud += model.calc_extinction_coefficient(**args)
            else:
                self.ec += model.calc_extinction_coefficient(**args)

        return self.ec


    def get_ec(self, temperature, densities, layer):
        """
        Interpolate the CS absorption into the planetary model temperature.
        """
        ec = np.zeros((np.sum(self.nspec),self.nwave))
        label = []
        j = 0
        for i,model in enumerate(self.models):
            imol = self.mol_indices[i]
            model_type = self.models_type[i]
            density = densities[:,imol]
            args = {
                'layer': layer,
            }
            if model.name == 'deck':
                args['temperature'] = temperature
                args['radius'] = None
            elif model_type in ['rayleigh', 'cloud']:
                args['density'] = density
            elif model_type == 'cia':
                args['temperature'] = temperature[layer]
                args['density'] = density[layer]
                args.pop('layer')
            else:
                args['temperature'] = temperature
                args['density'] = density

            if self.models_type[i] == 'line_sample':
                args['per_mol'] = True
            if hasattr(model, 'pars'):
            #if hasattr(model, 'pars') and model.npars > 0:
                args['pars'] = model.pars

            extinction = model.calc_extinction_coefficient(**args)
            ec[j:j+self.nspec[i]] = extinction
            j += self.nspec[i]

            if model_type in ['line_sample', 'lbl']:
                label += list(model.species)
            elif model_type == 'alkali':
                label.append(model.species)
            else:
                label.append(model.name)

        return ec, label


    def check_temp_bounds(self, temperatures):
        """
        Check if any temperature lies out of bounds for any opacity
        model.

        Parameters
        ----------
        temperatures: 1D float iterable
            Temperatures to test.

        Returns
        -------
        oob: List
            List of models where temperaures are out of bounds.
        """
        min_temp = np.amin(temperatures)
        max_temp = np.amax(temperatures)

        oob = []
        for model,temp in self.tmin.items():
            if min_temp < temp:
                oob.append(model)

        for model,temp in self.tmax.items():
            if max_temp > temp:
                oob.append(model)

        return list(np.unique(oob))


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write('Opacity extinction information:')
        fw.write('Model           type           T_min   T_max')
        for i,model in enumerate(self.models):
            mtype = self.models_type[i]
            if mtype in self.tmin:
                t_lims = f' {self.tmin[mtype]:7.1f} {self.tmax[mtype]:7.1f}'
            else:
                t_lims = ''
            if mtype in ['line_sample', 'lbl']:
                for j in range(model.nspec):
                    fw.write(f'{model.species[j]:15} {mtype:11} {t_lims}')
            else:
                fw.write(f'{model.name:15} {mtype:11} {t_lims}')
        return fw.text

