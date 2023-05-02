# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

from .. import opacity as op
from .. import tools as pt


class H_Ion():
    """Interface between H- opacity models and pyrat object"""
    def __init__(self, inputs, wn, species, log):
        self.model = None
        self.ec = None

        if inputs.h_ion is None:
            return

        self.model = op.Hydrogen_Ion(wn)

        has_all_species = (
            'H' in species and
            'e-' in species and
            'H-' in species
        )
        if not has_all_species:
            log.error(
                "'h_ion' opacity model requires the atmosphere to "
                "contain H, H-, and e- species"
        )

        # For calculations only H and e- are necessary:
        self.mol_indices = [
            list(species).index(mol)
            for mol in ['H', 'e-']
        ]


    def calc_extinction_coefficient(self, temperature, densities):
        """
        Evaluate the H- bound-free/free-free absorption in the atmosphere.
        """
        if self.model is not None:
            self.ec = self.model.calc_extinction_coefficient(
                temperature, densities[:,self.mol_indices],
            )

    def get_ec(self, temperature, densities, layer):
        """
        Extract per-species extinction coefficient (cm-1) at requested layer.
        """
        ec, label = [], []
        for imol,model in zip(self.imol, self.models):
            if imol is None:
                continue
            ext_coeff = model.calc_extinction_coefficient(
                temperature, densities[:,imol], layer,
            )
            ec.append(ext_coeff)
            label.append(model.mol)
        return ec, label


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write('H- opacity model:')
        for model in self.models:
            fw.write('\n' + str(model))
        fw.write(
            '\nTotal atmospheric extinction-coefficient (ec, cm-1):\n{}',
            self.ec,
            fmt={'float': '{:.3e}'.format},
        )
        return fw.text

