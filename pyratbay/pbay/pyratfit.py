# Copyright (c) 2016-2017 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

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
  """

  # Check stellar spectrum model:
  if pyrat.od.path == "eclipse" and pyrat.phy.starflux is None:
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
    q0 = pyrat.atm.q

  rejectflag = False
  # Update temperature profile if requested:
  if pyrat.ret.itemp is not None:
    temp = pyrat.ret.tmodel(params[pyrat.ret.itemp], *pyrat.ret.targs)
  else:
    temp = pyrat.atm.temp
  # Turn-on reject flag if out-of-bounds temperature:
  if np.any(temp < pyrat.ret.tlow) or np.any(temp > pyrat.ret.thigh):
    temp[:] = 0.5*(pyrat.ret.tlow + pyrat.ret.thigh)
    rejectflag = True

  # Update abundance profiles if requested:
  if pyrat.ret.iabund is not None:
    q2 = atm.qscale(pyrat.atm.q, pyrat.mol.name, params[pyrat.ret.iabund],
                    pyrat.ret.molscale, pyrat.ret.bulk,
                    iscale=pyrat.ret.iscale, ibulk=pyrat.ret.ibulk,
                    bratio=pyrat.ret.bulkratio, invsrat=pyrat.ret.invsrat)
  else:
    q2 = pyrat.atm.q

  # Update reference radius if requested:
  if pyrat.ret.irad is not None:
    pyrat.phy.rplanet = params[pyrat.ret.irad][0]*pc.km

  # Update Rayleigh parameters if requested:
  if pyrat.ret.iray is not None:
    j = 0
    rpars = params[pyrat.ret.iray]
    for i in np.arange(pyrat.rayleigh.nmodels):
      pyrat.rayleigh.model[i].pars = rpars[j:j+pyrat.rayleigh.model[i].npars]
      j += pyrat.rayleigh.model[i].npars

  # Update haze parameters if requested:
  if pyrat.ret.ihaze is not None:
    j = 0
    hpars = params[pyrat.ret.ihaze]
    for i in np.arange(pyrat.haze.nmodels):
      pyrat.haze.model[i].pars = hpars[j:j+pyrat.haze.model[i].npars]
      j += pyrat.haze.model[i].npars

  # Update patchy fraction if requested:
  if pyrat.ret.ipatchy is not None:
    pyrat.haze.fpatchy = params[pyrat.ret.ipatchy]

  # Calculate spectrum:
  pyrat = py.run(pyrat, [temp, q2, None])

  # Band-integrate spectrum:
  pyrat.obs.bandflux = w.bandintegrate(pyrat=pyrat)

  # Turn-on reject flag if atm doesn't cover the transit-depth values:
  if (pyrat.od.path == "transit" and
      np.any(pyrat.obs.data >
             (pyrat.atm.radius[pyrat.atm.rtop]/pyrat.phy.rstar)**2)):
    rejectflag = True

  # Reject this iteration if there are invalid temperatures or radii:
  if rejectflag:
    pyrat.obs.bandflux[:] = np.inf

  # Revert abundances in the atmospheric profile:
  if freeze:
    pyrat.atm.q = q0
  return pyrat.obs.bandflux
