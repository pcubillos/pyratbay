# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import numpy as np

from .. import constants as pc
from .. import opacity as op
from .. import tools as pt


class Opacity():
    """Interface between opacity models and pyrat object"""
    def __init__(self, inputs, wn, species, pressure, log):
        """
        Read collision induced absorption (CIA) files.
        """
        self.models = []
        self.models_type = []
        self.mol_indices = []
        # Min/max temperatures that can be sampled
        self.tmin = {}
        self.tmax = {}
        # extinction coefficient in cm-1
        nwave = len(wn)
        nlayers = len(pressure)
        self.ec = np.zeros((nlayers, nwave))

        min_wn = np.amin(wn)
        max_wn = np.amax(wn)
        species = list(species)

        #if len(inputs.extfile) > 0:
        # TBD: without runmode?
        if inputs.extfile is not None and inputs.runmode != 'opacity':
            log.head("\nReading cross-section table file(s):")
            #for cs_file in self.cs_files:
            #    log.head(f"  '{cs_file}'.")

            ls = op.Line_Sample(
                inputs.extfile, # TBD: self.ls_files?
                min_wn,
                max_wn,
                log,
            )
            self.models.append(ls)
            self.models_type.append('line_sample')
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
                str_press = str(ls.press/pc.bar)
            log.msg(
                f"Species names: {ls.species}\n"
                f"Temperatures (K):\n   {str_temp}\n"
                f"Pressure layers (bar):\n{str_press}\n"
                f"Wavenumber array (cm-1):\n   {str_wn}",
            )


    def calc_extinction_coefficient(self, temperature, densities):
        """
        Compute extinction coefficient over temperature and
        number density profiles.
        """
        self.ec[:] = 0.0
        for i,model in enumerate(self.models):
            density = densities[:,self.mol_indices[i]]
            self.ec += model.calc_extinction_coefficient(temperature, density)
        return self.ec


    def get_ec(self, temperature, densities, layer):
        """
        Interpolate the CS absorption into the planetary model temperature.
        """
        ec = []
        label = []
        for i,model in enumerate(self.models):
            imol = self.mol_indices[i]
            density = densities[:,imol]

            extinction = model.calc_extinction_coefficient(
                temperature, density, layer, per_mol=True,
            )
            ec += list(extinction)
            label += list(model.species)

        ## CIA
        ##for i,cia in enumerate(self.models):
        #    ext_coeff = model.calc_extinction_coefficient(
        #        temperature[layer], densities[layer,imol],
        #    )
        #    label.append(model.name)
        ## Alkali
        ##for imol,model in zip(self.imol, self.models):
        #    ext_coeff = model.calc_extinction_coefficient(
        #        temperature, density, layer,
        #    )
        #    label.append(model.mol)
        ## H-
        ##for imol,model in zip(self.imol, self.models):
        #    ext_coeff = model.calc_extinction_coefficient(
        #        temperature, density, layer,
        #    )
        #    label.append(model.mol)

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
            if mtype == 'line_sample':
                for j in range(model.nspec):
                    info = f'{model.species[j]:15} {mtype:10} {t_lims}'
                    fw.write(info)
        return fw.text

