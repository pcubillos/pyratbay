# Copyright (c) 2021-2022 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import numpy as np


def init(pyrat):
    """Setup alkali models for pyrat's atmosphere."""
    if pyrat.alkali.models != []:
        pyrat.log.head("\nSetup Alkali opacity models.")
    for alkali in pyrat.alkali.models:
        # Spectral sampling rate at alkali wn0:
        dwave = []
        for i,wn0 in enumerate(alkali.wn):
            i_wn = np.argmin(np.abs(pyrat.spec.wn-wn0))
            is_last = i_wn == pyrat.spec.nwave-1
            is_over = i_wn > 0 and pyrat.spec.wn[i_wn] > wn0
            if is_over or is_last:
                i_wn -= 1
            dwave.append(pyrat.spec.wn[i_wn+1]-pyrat.spec.wn[i_wn])

        alkali.setup(pyrat.mol.name, pyrat.mol.mass, dwave)


def absorption(pyrat):
    """
    Evaluate the total alkali absorption in the atmosphere.
    """
    # Initialize extinction coefficient:
    pyrat.alkali.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

    for alkali in pyrat.alkali.models:
        # Number density of alkali species:
        dens = np.expand_dims(pyrat.atm.d[:,alkali.imol], axis=1)
        # Calculate extinction coefficient (cm2 molecule-1):
        alkali.c_absorption(pyrat.atm.press, pyrat.atm.temp, pyrat.spec.wn)
        # Sum alkali extinction coefficient (cm-1):
        pyrat.alkali.ec += alkali.ec * dens


def get_ec(pyrat, layer):
    """
    Extract per-species extinction coefficient (cm-1) at requested layer.
    """
    absorption(pyrat)
    ec, label = [], []
    for alkali in pyrat.alkali.models:
        if alkali.imol >= 0:
            ec.append(alkali.ec[layer] * pyrat.atm.d[layer,alkali.imol])
            label.append(alkali.mol)
    return ec, label

