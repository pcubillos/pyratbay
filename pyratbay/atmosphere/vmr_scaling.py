# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'qcapcheck',
    'balance',
    'ratio',
    'qscale',
]

import numpy as np


def qcapcheck(vmr, qcap, ibulk):
    """
    Check if the cummulative abundance of traces exceeds qcap.

    Parameters
    ----------
    vmr: 2D float ndarray
        Volume mixing ratio of atmospheric species [nlayers, nspecies].
    qcap: Float
        Cap threshold for cummulative trace abundances.
    ibulk: 1D integer ndarray
        Indices of the bulk species to calculate the mixing ratio.

    Returns
    -------
    vmr_cap_flag: Bool
        Flag indicating whether trace abundances sum more than qcap.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> # Make an atmosphere:
    >>> pressure = pa.pressure(ptop=1e-8, pbottom=1e2, nlayers=11, units='bar')
    >>> temperature = np.tile(1500.0, 11)
    >>> species = ["H2", "He", "H2O"]
    >>> abundances = [0.8495, 0.15, 5e-4]
    >>> qprofiles = pa.uniform(pressure, temperature, species, abundances)
    >>> ibulk = [0,1]
    >>> # Sum of all metals (H2O) does not exceed qcap:
    >>> qcap = 1e-3
    >>> print(pa.qcapcheck(qprofiles, qcap, ibulk))
    False
    >>> # Sum of all metals (H2O) exceedes qcap:
    >>> qcap = 1e-4
    >>> print(pa.qcapcheck(qprofiles, qcap, ibulk))
    True
    """
    if qcap is None:
        return False

    # The shape of things:
    nlayers, nspecies = np.shape(vmr)

    # Get the indices of the species not in ibulk (trace species):
    itrace = np.setdiff1d(np.arange(nspecies), ibulk)

    # Sum the abundances of everything exept the ibulk species (per layer):
    qtrace = np.sum(vmr[:,itrace], axis=1)

    # Do sum of trace abundances exceed Qcap?
    vmr_cap_flag = np.any(qtrace > qcap)
    return vmr_cap_flag


def balance(vmr, ibulk, ratio, invsrat):
    r"""
    Balance the volume mixing ratios of bulk species, vmr[ibulk],
    such that sum(vmr) = 1.0 at each level.

    Parameters
    ----------
    vmr: 2D float ndarray
        Volume mixing ratio of atmospheric  species [nlayers, nspecies].
    ibulk: 1D integer ndarray
        Indices of the bulk species to calculate the mixing ratio.
    ratio: 2D float ndarray
        Abundance ratio between species indexed by ibulk.
    invsrat: 1D float ndarray
        Inverse of the sum of the ratios (at each layer).

    Notes
    -----
    Let the bulk abundance species be the remainder of the sum of the trace
    species:
        vmr_bulk = sum vmr_j = 1.0 - sum vmr_trace.
    This code assumes that the abundance ratio among bulk species
    remains constant in each layer:
        {\rm ratio}_j = vmr_j/vmr_0.
    The balanced abundance of the bulk species is then:
        vmr_j = \frac{{\rm ratio}_j * vmr_{\rm bulk}} {\sum {\rm ratio}}.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>>
    >>> vmr = np.tile([0.8, 0.2, 0.1], (5,1))
    >>> ibulk = [0, 1]
    >>> bratio, invsrat = pa.ratio(vmr, ibulk)
    >>> pa.balance(vmr, ibulk, bratio, invsrat)
    >>> # Balanced VMRs:
    >>> print(vmr[0])
    [0.72 0.18 0.1 ]
    >>> # Sum of VMRs equals one at each layer:
    >>> print(np.sum(vmr, axis=1))
    [1. 1. 1. 1. 1.]
    >>> # Ratio of 'bulk' species remains constant:
    >>> print(vmr[:,1]/vmr[:,0])
    [0.25 0.25 0.25 0.25 0.25]
    """
    # The shape of things:
    nlayers, nspecies = np.shape(vmr)
    nratio = len(ibulk)

    # Get the indices of the species not in ibulk (trace species):
    itrace = np.setdiff1d(np.arange(nspecies), ibulk)

    # Sum the abundances of everything exept the ibulk species (per layer):
    sum_vmr_metals = 1.0 - np.sum(vmr[:,itrace], axis=1)

    # Calculate the balanced mole mixing ratios:
    for j in range(nratio):
        vmr[:,ibulk[j]] = ratio[:,j] * sum_vmr_metals * invsrat


def ratio(vmr, ibulk):
    """
    Calculate the abundance ratios of the species indexed by ibulk, relative
    to the first species in the list.

    Parameters
    ----------
    vmr: 2D float ndarray
        Volume mixing ratio of atmospheric species [nlayers, nspecies].
    ibulk: 1D integer ndarray
        Indices of the species to calculate the ratio.

    Returns
    -------
    bratio: 2D float ndarray
        Abundance ratio between species indexed by ibulk.
    invsrat: 1D float ndarray
        Inverse of the sum of the ratios (at each layer).

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> vmr = np.tile([0.8, 0.2], (5,1))
    >>> ibulk = [0, 1]
    >>> bratio, invsrat = pa.ratio(vmr, ibulk)
    >>> print(bratio)
    [[ 1.    0.25]
     [ 1.    0.25]
     [ 1.    0.25]
     [ 1.    0.25]
     [ 1.    0.25]]
    >>> print(invsrat)
    [ 0.8  0.8  0.8  0.8  0.8]
    """
    # The shape of things:
    nlayers, nspecies = np.shape(vmr)
    nratio = len(ibulk)
    bratio = np.ones((nlayers, nratio))

    # Calculate the abundance ratio WRT first indexed species in ibulk:
    for j in range(1, nratio):
        bratio[:,j] = vmr[:,ibulk[j]] / vmr[:,ibulk[0]]

    # Inverse sum of ratio:
    invsrat = 1.0 / np.sum(bratio, axis=1)

    return bratio, invsrat


def qscale(
        vmr, species, molvars, molpars, bulk,
        qsat=None, iscale=None, ibulk=None, bratio=None, invsrat=None,
    ):
    """
    Scale specified species abundances and balance bulk abundances to
    conserve sum(vmr)=1 in each layer.

    Parameters
    ----------
    vmr: 2D float ndarray
        Volume mixing ratio of atmospheric species [nlayers, nspecies].
    species: 1D string ndarray
        Names of the species in the atmosphere.
    molvars: 1D string ndarray
        Model to vary the species abundances.
    molpars: 1D float ndarray
        Scaling factor (dex) for each species in molvars.
    bulk: 1D string ndarray
        Names of the bulk (dominant) species.
    qsat: Float
        Maximum allowed combined abundance for trace species.
    iscale: 1D integer ndarray
        Indices of molvars species in vmr.
    ibulk: 1D integer ndarray
        Indices of bulk species in vmr.
    bratio: 2D float ndarray
        Abundance ratios between the bulk species (relative to bulk[0]).
    invsrat: 1D float ndarray
        Inverse of the sum of the ratios (at each layer).

    Returns
    -------
    scaled_vmr: 2D float ndarray
       The modified atmospheric VMR profiles.

    Notes
    -----
    iscale, ibulk, bratio, and invsrat are optional parameters to
    speed up the routine.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa

    >>> nlayers = 2
    >>> vmr = np.tile([0.85, 0.15, 1e-4, 1e-4], (nlayers,1))
    >>> species = ['H2', 'He' ,'H2O', 'CH4']
    >>> # Set the H2O abundance to 1.0e-3:
    >>> molvars = ['log_H2O']
    >>> molpars = [-3.0]
    >>> bulk = ['H2', 'He']
    >>> scaled_vmr = pa.qscale(vmr, species, molvars, molpars, bulk)
    >>> print(scaled_vmr)
    [[8.49065e-01 1.49835e-01 1.00000e-03 1.00000e-04]
     [8.49065e-01 1.49835e-01 1.00000e-03 1.00000e-04]]
    """
    if iscale is None:
        molecs = [var.split('_')[-1] for var in molvars]
        iscale = [list(species).index(mol) for mol in molecs]
    if ibulk is None:
        ibulk  = [list(species).index(mol) for mol in bulk]
    if bratio is None:
        bratio, invsrat = ratio(vmr, ibulk)

    # Scale abundance of requested species:
    scaled_vmr = np.copy(vmr)
    for idx,model,value in zip(iscale, molvars, molpars):
        if model.startswith('scale_'):
            scaled_vmr[:,idx] = vmr[:,idx] * 10.0**value
        elif model.startswith('log_'):
            scaled_vmr[:,idx] = 10.0**value

    # Enforce saturation limit:
    if qsat is not None:
        indices = np.arange(len(species))
        ifixed = np.setdiff1d(indices, np.union1d(ibulk, iscale))
        sum_fixed = np.sum(scaled_vmr[:,ifixed], axis=1, keepdims=True)
        sum_scaled = np.sum(scaled_vmr[:,iscale], axis=1, keepdims=True)

        q0 = (qsat - sum_fixed) / sum_scaled
        q0 = np.clip(q0, 0.0, 1.0)
        scaled_vmr[:,iscale] *= q0

    # Scale abundance of bulk species to balance sum(vmr):
    balance(scaled_vmr, ibulk, bratio, invsrat)

    return scaled_vmr
