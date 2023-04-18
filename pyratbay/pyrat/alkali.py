# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import numpy as np
from .. import opacity as op
from .. import tools as pt


class Alkali(object):
    """Interface between Alkali opacity models and pyrat object"""
    def __init__(self, model_names, pressure, wn, cutoff, species, log):
        self.models = []
        self.imol = []
        self.ec = None

        if model_names is None:
            return

        for name in model_names:
            model = op.alkali.get_model(name, pressure, wn, cutoff)
            self.models.append(model)

            if model.mol in species:
                imol = list(species).index(model.mol)
            else:
                imol = None
            self.imol.append(imol)

        log.head("\nSetup Alkali opacity models.")



    def calc_extinction_coefficient(self, temperature, densities):
        """
        Evaluate the total alkali absorption in the atmosphere.
        """
        ec = 0.0
        for imol,model in zip(self.imol, self.models):
            if imol is None:
                continue
            # Number density of alkali species:
            dens = densities[:,imol]
            # Calculate extinction coefficient (cm-1):
            ec += model.calc_extinction_coefficient(temperature, dens)

        if not np.isscalar(ec):
            self.ec = ec


    def get_ec(self, temperature, densities, layer):
        """
        Extract per-species extinction coefficient (cm-1) at requested layer.
        """
        self.calc_extinction_coefficient(temperature, densities)
        ec, label = [], []
        for imol,model in zip(self.imol, self.models):
            if imol is None:
                continue
            ec.append(model.cross_section[layer] * densities[layer,imol])
            label.append(model.mol)
        return ec, label


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write('Alkali-opacity models (models):')
        for model in self.models:
            fw.write('\n' + str(model))
        fw.write(
            '\nTotal atmospheric alkali extinction-coefficient (ec, cm-1):\n{}',
            self.ec,
            fmt={'float': '{:.3e}'.format},
        )
        return fw.text

