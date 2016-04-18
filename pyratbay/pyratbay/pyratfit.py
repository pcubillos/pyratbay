from .  import qscale as qs

def init(pyrat):
  """
  Initialize variables that will be used in fit().
  """
  # Input converter:
  rat, invsrat = b.ratio(q, ibulk)
  ibulk  = np.where(np.in1d(spec, bulk))[0]
  iscale = np.where(np.in1d(spec, molscale))[0]

  # Output converter:
  wnindices = None
  istarfl   = None
  rprs      = None
  specwn    = None
  nifilter  = None


def fit(params, pyrat):
  """
  Fitting routine for MCMC.
  """

  # Parse fitting parameters:

  # Update atmospheric profile:
  temp = tmodel(tparams)
  q2 = pb.pyratbay.b.qscale(q, spec, qscale, molscale, bulk,
         iscale=iscale, ibulk=ibulk, bratio=rat, invsrat=invsrat)
  # Calculate spectrum:
  pyrat = pyrat.run(pyrat, temp, q2)
  # Band-integrate spectrum:
  bandflux = None

  return bandflux
