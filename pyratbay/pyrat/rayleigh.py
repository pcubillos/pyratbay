# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Rayleigh',
]

import numpy as np
from .. import opacity as op
from .. import tools as pt


class Rayleigh():
    """Interface between Rayleigh opacity models and pyrat object"""
    def __init__(self, model_names, pars, species, wn, log, cloud_obj):
        self.model_names = model_names
        self.models = []
        self.pnames = []
        self.texnames = []
        self.mol_indices = []
        self.npars = 0
        self.pars = None
        self.ec = None
        self._cloud = cloud_obj

        if model_names is None:
            return

        for name in model_names:
            if name.startswith('dalgarno_'):
                mol = name.split('_')[1]
                model = op.rayleigh.Dalgarno(wn, mol)
            if name == 'lecavelier':
                model = op.rayleigh.Lecavelier(wn)
            self.models.append(model)
            self.npars += model.npars
            self.pnames += model.pnames
            self.texnames += model.texnames

        # Set index of species for each Rayleigh model:
        for model in self.models:
            mol = model.mol
            if mol not in species:
                self.mol_indices.append(None)
            else:
                self.mol_indices.append(list(species).index(mol))

        # Parse parameters:
        if pars is None:
            return
        self.pars = pars
        input_npars = len(self.pars)
        if self.npars != input_npars:
            log.error(
                f'Number of input Rayleigh parameters ({input_npars}) '
                'does not match the number of required '
                f'model parameters ({self.npars})'
            )
        j = 0
        for model in self.models:
            model.pars = self.pars[j:j+model.npars]
            j += model.npars


    def calc_extinction_coefficient(self, densities):
        """
        Evaluate the total Rayleigh absorption (cm-1) in the atmosphere.

        Parameters
        ----------
        densities: 2D float array
            Number-density atmospheric profiles (molecules cm-3)
        """
        self.ec = 0.0
        j = 0
        for idx,model in zip(self.mol_indices, self.models):
            if idx is None:
                continue

            args = dict()
            args['density'] = densities[:,idx]
            if model.npars > 0:
                args['pars'] = self.pars[j:j+model.npars]
                j += model.npars
            ec = model.calc_extinction_coefficient(**args)

            # Put into cloud.ec instead to apply fpatchy factor
            if model.name == 'lecavelier' and self._cloud.fpatchy is not None:
                self._cloud.ec += ec
            else:
                self.ec += ec

        if np.isscalar(self.ec):
            self.ec = None


    def get_ec(self, densities, layer):
        """
        Extract per-model extinction coefficient at requested layer.
        """
        ec, label = [], []
        j = 0
        for idx,model in zip(self.mol_indices, self.models):
            args = dict()
            args['density'] = densities[:,idx]
            args['layer'] = layer
            if model.npars > 0:
                args['pars'] = self.pars[j:j+model.npars]
                j += model.npars
            ec.append(model.calc_extinction_coefficient(**args))
            label.append(model.name)
        return ec, label


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write('Rayleigh-opacity models (models):')
        for model in self.models:
            fw.write('\n' + str(model))
        fw.write('\nTotal atmospheric Rayleigh extinction-coefficient '
                 '(ec, cm-1):\n{}', self.ec, fmt={'float': '{: .3e}'.format})
        return fw.text

