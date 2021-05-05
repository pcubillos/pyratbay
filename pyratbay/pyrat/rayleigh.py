# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import numpy as np


def absorption(pyrat):
    """
    Evaluate the total Rayleigh absorption (cm-1) in the atmosphere.
    """
    pyrat.rayleigh.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

    for rmodel in pyrat.rayleigh.models:
        # Opacity cross section (in cm2 molecule-1):
        rmodel.extinction(pyrat.spec.wn)
        # Get molecule index:
        if rmodel.mol not in pyrat.mol.name:
            continue
        imol = np.where(pyrat.mol.name == rmodel.mol)[0][0]

        # Densities in molecules cm-3:
        dens = pyrat.atm.d[:,imol]
        # Extinction coefficient (cm-1):
        if rmodel.name == 'lecavelier' and pyrat.cloud.fpatchy is not None:
            # Put it in cloud.ec because we want to apply the fpatchy factor
            # to this 'unknown' haze source.
            pyrat.cloud.ec += rmodel.ec * np.expand_dims(dens, axis=1)
        else:
            pyrat.rayleigh.ec += rmodel.ec * np.expand_dims(dens, axis=1)


def get_ec(pyrat, layer):
    """
    Extract per-model extinction coefficient at requested layer.
    """
    ec, label = [], []
    for rmodel in pyrat.rayleigh.models:
        rmodel.extinction(pyrat.spec.wn)
        if rmodel.mol in pyrat.mol.name:
            imol = np.where(pyrat.mol.name == rmodel.mol)[0][0]
            ec.append(rmodel.ec * pyrat.atm.d[layer,imol])
        else:
            ec.append(np.zeros_like(pyrat.spec.wn))
        label.append(rmodel.name)
    return ec, label
