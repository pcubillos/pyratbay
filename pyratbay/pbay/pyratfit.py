# Copyright (c) 2016 Patricio Cubillos and contributors.
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
    t0, q0, r0 = pyrat.atm.temp, pyrat.atm.q, pyrat.atm.radius

  # Update temperature profile:
  temp = pyrat.ret.tmodel(params[pyrat.ret.itemp], *pyrat.ret.targs)
  # Turn-on out-of-bounds temperature flag:
  if np.any(temp < pyrat.ret.tlow) or np.any(temp > pyrat.ret.thigh):
    temp[:] = 0.5*(pyrat.ret.tlow + pyrat.ret.thigh)
    tempflag = True

  # Update abundance profiles:
  q2 = atm.qscale(pyrat.atm.q, pyrat.mol.name, params[pyrat.ret.iabund],
                  pyrat.ret.molscale, pyrat.ret.bulk,
                  iscale=pyrat.ret.iscale, ibulk=pyrat.ret.ibulk,
                  bratio=pyrat.ret.bulkratio, invsrat=pyrat.ret.invsrat)
  # Update radius profile:
  if len(pyrat.ret.irad) > 0:
    radius = pyrat.hydro(pyrat.atm.press, temp, pyrat.atm.mm,
                         pyrat.phy.gplanet, pyrat.phy.mplanet,
                         pyrat.refpressure, params[pyrat.ret.irad][0]*pc.km)
  else:
    radius = None
  # Update haze parameters:
  if len(pyrat.ret.ihaze) > 0:
    j = 0
    hpars = params[pyrat.ret.ihaze]
    for i in np.arange(pyrat.haze.nmodels):
      pyrat.haze.model[i].pars = hpars[j:j+pyrat.haze.model[i].npars]
      j += pyrat.haze.model[i].npars

  # Calculate spectrum:
  pyrat = py.run(pyrat, [temp, q2, radius])

  # Band-integrate spectrum:
  pyrat.obs.bandflux = w.bandintegrate(pyrat=pyrat)

  # Reject this iteration if there are invalid temperatures or radii:
  if tempflag or (pyrat.od.path == "transit" and
      np.any(pyrat.obs.data > (pyrat.atm.radius[0]/pyrat.phy.rstar)**2)):
    pyrat.obs.bandflux[:] = np.inf

  # Revert changes in the atmospheric profile:
  if freeze:
    pyrat.atm.temp, pyrat.atm.q, pyrat.atm.radius = t0, q0, r0
  return pyrat.obs.bandflux
