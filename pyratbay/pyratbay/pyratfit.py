import sys
import numpy as np

from .. import tools     as pt
from .. import constants as pc
from .. import pyrat     as py
from .. import wine      as w
from .. import atmosphere as atm


def init(pyrat, args, log):
  """
  Initialize variables that will be used in the atmospheric retrieval.
  Loads the stellar spectrum.
  Loads the waveband transmission filters.
  """

  # Check stellar spectrum model:
  if pyrat.phy.starflux is None:
    pt.error("Unspecified stellar flux model.", log)

  # Check filter files and data:
  if pyrat.obs.filter is None:
    pt.error("Undefined waveband transmission filters (filter).", log)
  if pyrat.obs.data is None:
    pt.error("Undefined data.", log)
  if pyrat.obs.uncert is None:
    pt.error("Undefined data uncertainties.", log)

  # Check rprs:
  if pyrat.od.path == "eclipse" and pyrat.phy.rprs is None:
    pt.error("Undefined Rp/Rs.", log)

  # Boundaries for temperature profile:
  pyrat.ret.tlow  = args.tlow
  pyrat.ret.thigh = args.thigh


def fit(params, pyrat, freeze=False):
  """
  Fitting routine for MCMC.

  Parameters
  ----------
  params: 1D float ndarray
     Array of fitting parameters that define the atmosphere.
  pyrat: Pyrat object instance
     Pyrat object.
  freeze: Bool
     If True, (after the spectrum calculation) reset the atmospheric
     temperature, abundance, and radius profiles to the original values.
     Note that in this case the pyrat object will contain inconsistent
     values between the atmospheric profiles and the spectrum.

  Returns
  -------
  bandflux: 1D float ndarray
     The waveband-integrated spectrum values.
  """

  if freeze:
    t0, q0, r0 = pyrat.atm.temp, pyrat.atm.q, pyrat.atm.radius

  # Update temperature profile:
  temp = pyrat.ret.tmodel(params[pyrat.ret.itemp], *pyrat.ret.targs)
  if np.any(temp < pyrat.ret.tlow) or np.any(temp > pyrat.ret.thigh):
    pyrat.obs.bandflux[:] = -1e10  # FINDME: what if np.inf? or nan?
    return pyrat.obs.bandflux
  # Update abundance profiles:
  q2 = atm.qscale(pyrat.atm.q, pyrat.mol.name, params[pyrat.ret.iabund],
                  pyrat.ret.molscale, pyrat.ret.bulk,
                  iscale=pyrat.ret.iscale, ibulk=pyrat.ret.ibulk,
                  bratio=pyrat.ret.bulkratio, invsrat=pyrat.ret.invsrat)
  # Update radius profile:
  if len(pyrat.ret.irad) > 0:
    radius = atm.hydro_equilibrium(pyrat.atm.press, temp, pyrat.atm.m,
        pyrat.gplanet, pyrat.refpressure, params[pyrat.irad][0]*pc.km)
  else:
    radius = None

  # Calculate spectrum:
  pyrat = py.run(pyrat, [temp, q2, radius])

  # Unpack some variables:
  spectrum = pyrat.spec.spectrum
  wn       = pyrat.spec.wn
  bflux    = pyrat.obs.bandflux
  wnidx    = pyrat.obs.bandidx
  # Band-integrate spectrum:
  for i in np.arange(pyrat.obs.nfilters):
    # Integrate the spectrum over the filter band:
    if   pyrat.od.path == "transit":
      bflux[i] = w.bandintegrate(spectrum[wnidx[i]], wn[wnidx[i]],
                                 pyrat.obs.bandtrans[i])
    elif pyrat.od.path == "eclipse":
      fluxrat = spectrum[wnidx[i]]/pyrat.obs.starflux[i] * pyrat.phy.rprs**2.0
      bflux[i] = w.bandintegrate(fluxrat, wn[wnidx[i]], pyrat.obs.bandtrans[i])

  # Revert changes in the atmospheric profile:
  if freeze:
    pyrat.atm.temp, pyrat.atm.q, pyrat.atm.radius = t0, q0, r0
  return pyrat.obs.bandflux
