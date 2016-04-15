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
  q = 1.0 - np.sum(Q[:, itrace], axis=1)

  # Calculate the balanced mole mixing ratios:
  for j in np.arange(nratio):
    Q[:,ibulk[j]] = ratio[:,j] * q * invsrat


def ratio(Q, ibulk):
  """
  Calculate the abundance ratios of the species relative to the first
  species indexed in ibulk.

  Parameters
  ----------
  Q: 2D float ndarray
     Mole mixing ratio of the species in the atmosphere [Nlayers, Nspecies].
  ibulk: 1D integer ndarray
     Indices of the species to calculate the ratio.

  Returns
  -------
  ratio: 2D float ndarray
     Abundance ratio between species indexed by ibulk.
  invsrat: 1D float ndarray
     Inverse of the sum of the ratios (at each layer).
  """
  # The shape of things:
  nlayers, nspecies = np.shape(Q)
  nratio = len(ibulk)
  ratio = np.ones((nlayers, nratio))

  # Calculate the abundance ratio WRT first indexed species in ibulk: 
  for j in np.arange(1, nratio):
    ratio[:,j] = Q[:,ibulk[j]] / Q[:,ibulk[0]]

  # Inverse sum of ratio:
  invsrat = 1.0 / np.sum(ratio, axis=1)

  return ratio, invsrat
