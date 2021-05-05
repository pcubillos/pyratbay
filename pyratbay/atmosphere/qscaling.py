# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = ["qcapcheck", "balance", "ratio", "qscale"]

import numpy as np


def qcapcheck(Q, qcap, ibulk):
  """
  Check if the cummulative abundance of traces exceed qcap.

  Parameters
  ----------
  Q: 2D float ndarray
     Mole mixing ratio of the species in the atmosphere [Nlayers, Nspecies].
  qcap: Float
     Cap threshold for cummulative trace abundances.
  ibulk: 1D integer ndarray
     Indices of the bulk species to calculate the mixing ratio.

  Returns
  -------
  qcapcheck: Bool
     Flag indicating whether trace abundances sum more than qcap.

  Examples
  --------
  >>> import pyratbay.atmosphere as pa
  >>> # Make an atmosphere:
  >>> pressure    = pa.pressure(ptop=1e-8, pbottom=1e2, nlayers=11, units='bar')
  >>> temperature = np.tile(1500.0, 11)
  >>> species     = ["H2", "He", "H2O"]
  >>> abundances  = [0.8495, 0.15, 5e-4]
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
  # The shape of things:
  nlayers, nspecies = np.shape(Q)

  # Get the indices of the species not in ibulk (trace species):
  itrace = np.setdiff1d(np.arange(nspecies), ibulk)

  # Sum the abundances of everything exept the ibulk species (per layer):
  qtrace = np.sum(Q[:,itrace], axis=1)

  # Do sum of trace abundances exceed Qcap?
  if np.any(qtrace > qcap):
      return True
  return False


def balance(Q, ibulk, ratio, invsrat):
  r"""
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

  Examples
  --------
  >>> import pyratbay.atmosphere as pa
  >>> q = np.tile([0.8, 0.2, 0.5], (5,1))
  >>> q[4] = 0.5, 0.5, 0.5
  >>> ibulk = [0, 1]
  >>> bratio, invsrat = pa.ratio(q, ibulk)
  >>> pa.balance(q, ibulk, bratio, invsrat)
  >>> print(np.sum(q,axis=1))
  [ 1.  1.  1.  1.  1.]
  >>> print(q[:,1]/q[:,0])
  [ 0.25  0.25  0.25  0.25  1.  ]
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

  Examples
  --------
  >>> import pyratbay.atmosphere as pa
  >>> q = np.tile([0.8, 0.2], (5,1))
  >>> q[4] = 0.5, 0.5
  >>> ibulk = [0, 1]
  >>> bratio, invsrat = pa.ratio(q, ibulk)
  >>> print(bratio)
  [[ 1.    0.25]
   [ 1.    0.25]
   [ 1.    0.25]
   [ 1.    0.25]
   [ 1.    1.  ]]
  >>> print(invsrat)
  [ 0.8  0.8  0.8  0.8  0.5]
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


def qscale(Q, spec, molmodel, molfree, molpars, bulk, qsat=None,
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
  molmodel: 1D string ndarray
     Model to vary the species abundances.
  molfree:  1D string ndarray
     Names of the species to vary their abundances.
  molpars:  1D float ndarray
     Scaling factor (dex) for each species in molfree.
  bulk:  1D string ndarray
     Names of the bulk (dominant) species.
  qsat:  Float
     Maximum allowed combined abundance for trace species.
  iscale:  1D integer ndarray
     Indices of molfree species in Q.
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
  spec = list(spec)
  q = np.copy(Q)
  if iscale is None:
      iscale = [spec.index(mol) for mol in molfree]
  if ibulk is None:
      ibulk  = [spec.index(mol) for mol in bulk]
  if bratio is None:
      bratio, invsrat = ratio(Q, ibulk)

  # Scale abundance of requested species:
  for idx,value,model in zip(iscale, molpars, molmodel):
      if model == 'scale':
          q[:,idx] = Q[:,idx] * 10.0**value
      elif model == 'vert':
          q[:,idx] = 10.0**value

  # Enforce saturation limit:
  if qsat is not None:
      ifix = np.setdiff1d(np.arange(len(spec)), np.union1d(ibulk, iscale))
      q0 = ((qsat - np.sum(q[:,ifix],   axis=1, keepdims=True)) /
                    np.sum(q[:,iscale], axis=1, keepdims=True))
      q0 = np.clip(q0, 0.0, 1.0)
      q[:,iscale] *= q0

  # Scale abundance of bulk species to balance sum(Q):
  balance(q, ibulk, bratio, invsrat)

  return q
