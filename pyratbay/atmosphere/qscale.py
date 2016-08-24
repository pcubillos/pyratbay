# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["balance", "ratio", "qscale"]

import numpy as np

def balance(Q, ibulk, ratio, invsrat):
  """
  Balance the mole mixing ratios of the bulk species, Q[ibulk],
  such that Sum(Q) = 1.0 at each level.

  Parameters
  ----------
  Q: 2D float ndarray
     Mole mixing ratio of the species in the atmosphere [Nlayers, Nspecies].
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
     Q_{\rm bulk} = \sum Q_j = 1.0 - \sum Q_{\rm trace}.
  This code assumes that the abundance ratio among bulk species
  remains constant in each layer:
     {\rm ratio}_j = Q_j/Q_0.
  The balanced abundance of the bulk species is then:
     Q_j = \frac{{\rm ratio}_j * Q_{\rm bulk}} {\sum {\rm ratio}}.
  """
  # The shape of things:
  nlayers, nspecies = np.shape(Q)
  nratio = len(ibulk)

  # Get the indices of the species not in ibulk (trace species):
  itrace = np.setdiff1d(np.arange(nspecies), ibulk)

  # Sum the abundances of everything exept the ibulk species (per layer):
  q = 1.0 - np.sum(Q[:,itrace], axis=1)

  # Calculate the balanced mole mixing ratios:
  for j in np.arange(nratio):
    Q[:,ibulk[j]] = ratio[:,j] * q * invsrat


def ratio(Q, ibulk):
  """
  Calculate the abundance ratios of the species indexed by ibulk, relative
  to the first species in the list.

  Parameters
  ----------
  Q: 2D float ndarray
     Mole mixing ratio of the species in the atmosphere [Nlayers, Nspecies].
  ibulk: 1D integer ndarray
     Indices of the species to calculate the ratio.

  Returns
  -------
  bratio: 2D float ndarray
     Abundance ratio between species indexed by ibulk.
  invsrat: 1D float ndarray
     Inverse of the sum of the ratios (at each layer).
  """
  # The shape of things:
  nlayers, nspecies = np.shape(Q)
  nratio = len(ibulk)
  bratio = np.ones((nlayers, nratio))

  # Calculate the abundance ratio WRT first indexed species in ibulk: 
  for j in np.arange(1, nratio):
    bratio[:,j] = Q[:,ibulk[j]] / Q[:,ibulk[0]]

  # Inverse sum of ratio:
  invsrat = 1.0 / np.sum(bratio, axis=1)

  return bratio, invsrat


def qscale(Q, spec, qscale, molscale, bulk, qsat=None,
           iscale=None, ibulk=None, bratio=None, invsrat=None):
  """
  Scale specified species abundances and balance bulk abundances to
  conserve sum(Q)=1 in each layer.

  Parameters
  ----------
  Q:  2D float ndarray
     Mole mixing ratio of the species in the atmosphere [Nlayers, Nspecies].
  spec:  1D string ndarray
     Names of the species in the atmosphere.
  qscale:  1D float ndarray
     Scaling factor (dex) for each species in molscale.
  molscale:  1D string ndarray
     Names of the species to scale their abundance profiles.
  bulk:  1D string ndarray
     Names of the bulk (dominant) species.
  qsat:  Float
     Maximum allowed combined abundance for trace species.
  iscale:  1D integer ndarray
     Indices of molscale species in Q.
  ibulk:  1D integer ndarray
     Indices of the bulk species in Q.
  bratio:  2D float ndarray
     Abundance ratios between the bulk species (relative to bulk[0]).
  invsrat:  1D float ndarray
     Inverse of the sum of the ratios (at each layer).

  Returns
  -------
  q:  2D float ndarray
     The modified atmospheric abundance profiles.

  Notes
  -----
  iscale, ibulk, bratio, and invsrat are optional parameters to speed up
     the routine.
  I'm not completely happy with qsat yet.
  """
  q = np.copy(Q)
  if iscale is None:
    iscale = []
    for mol in molscale:
      iscale += list(np.where(spec==mol)[0])
  if ibulk is None:
    ibulk = []
    for mol in bulk:
      ibulk  += list(np.where(spec==mol)[0])
  if bratio is None:
    bratio, invsrat = ratio(Q, ibulk)

  # Scale abundance of requested species:
  for i in np.arange(len(iscale)):
    m = iscale[i]
    q[:,m] = Q[:,m] * 10.0**qscale[i]

  # Enforce saturation limit:
  if qsat is not None:
    ifix = np.setdiff1d(np.arange(len(spec)), np.union1d(ibulk, iscale))
    #print(ifix)
    q0 = ((qsat - np.sum(q[:,ifix],   axis=1, keepdims=True)) /
                  np.sum(q[:,iscale], axis=1, keepdims=True))
    #print(q0)
    q0 = np.clip(q0, 0.0, 1.0)
    q[:,iscale] *= q0
    #print(np.sum(q[:,ifix], axis=1)+np.sum(q[:,iscale], axis=1))

  # Scale abundance of bulk species to balance sum(Q):
  balance(q, ibulk, bratio, invsrat)

  return q
