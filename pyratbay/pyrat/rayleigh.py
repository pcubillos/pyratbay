# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import numpy as np


def absorption(pyrat):
    """
    Evaluate the total Rayleigh absorption (cm-1) in the atmosphere.
    """
    pyrat.rayleigh.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

    for rmodel in pyrat.rayleigh.models:
        # Extinction coefficient (in cm2 molecule-1):
        rmodel.extinction(pyrat.spec.wn)
        # Get molecule index:
        if rmodel.mol not in pyrat.mol.name:
            continue
        imol = np.where(pyrat.mol.name == rmodel.mol)[0][0]

        # Densities in molecules cm-3:
        dens = pyrat.atm.d[:,imol]
        # Absorption (cm-1):
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
