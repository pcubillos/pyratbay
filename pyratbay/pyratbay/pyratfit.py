import sys
import numpy as np

from .. import tools     as pt
from .. import constants as pc
from .. import pyrat     as py
from .. import wine      as w

from .  import qscale  as qs
from .  import makeatm as ma
from .  import argum   as ar

def init(pyrat, args, log):
  """
  Initialize variables that will be used in the atmospheric retrieval.
  Checks bulk and variable-abundance species.
  Loads the stellar spectrum.
  Loads the waveband transmission filters.
  """
  species = pyrat.mol.name
  # Process bulk-abundance species:
  if args.bulk is None:
    pt.error("Undefined bulk species list (bulk).", log)
  if len(np.setdiff1d(args.bulk, species)) > 0:
    pt.error("These bulk species are not present in the atmosphere: {:s}".
      format(str(np.setdiff1d(args.bulk, species))), log)

  # Process variable-abundance species:
  if args.molscale is None:
    pt.warning(pyrat.verb-2, "There are no variable-abundance species "
                             "(molscale).", log)
  else:
    if len(np.setdiff1d(args.molscale, species)) > 0:
      pt.error("These variable-abundance species are not present "
               "in the atmosphere: {:s}".
                format(str(np.setdiff1d(args.molscale, species))), log)
    if len(np.intersect1d(args.bulk, args.molscale)) > 0:
      pt.error("These species were marked as both bulk and "
               "variable-abundance: {:s}".
                format(np.intersect1d(args.bulk, args.molscale)))

  # Obtain abundance ratios between the bulk species:
  pyrat.bulk, pyrat.molscale = args.bulk, args.molscale
  pyrat.ibulk  = np.where(np.in1d(species, args.bulk))[0]
  pyrat.iscale = np.where(np.in1d(species, args.molscale))[0]
  pyrat.bulkratio, pyrat.invsrat = qs.ratio(pyrat.atm.q, pyrat.ibulk)

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

  if pyrat.od.path == "eclipse":
    pyrat.rprs = pyrat.rplanet/pyrat.rstar

  # Temperature model:
  # FINDME: Need to check args.tstar, tint, smaxis
  pyrat.tmodel, pyrat.targs, ntemp = ma.temperature(args.tmodel,
     eval=False, pressure=pyrat.atm.press, rstar=args.rstar, tstar=args.tstar,
     tint=args.tint, gplanet=args.gplanet, smaxis=args.smaxis,
     radunits=pyrat.radunits, nlayers=pyrat.atm.nlayers, log=log)

  # Boundaries for temperature profile:
  pyrat.tlow  = args.tlow
  pyrat.thigh = args.thigh

  # Indices to parse the array of fitting parameters:
  nrad   = int(pyrat.od.path == "transit")
  nabund = len(pyrat.iscale)
  nhaze  = 0
  nalk   = 0

  pyrat.itemp  = np.arange(0,          ntemp)
  pyrat.irad   = np.arange(ntemp,      ntemp+nrad)
  pyrat.iabund = np.arange(ntemp+nrad, ntemp+nrad+nabund)
  pyrat.ihaze  = np.arange(ntemp+nrad+nabund, ntemp+nrad+nabund+nhaze)
  pyrat.ialk   = np.arange(ntemp+nrad+nabund+nhaze,
                           ntemp+nrad+nabund+nhaze+nalk)

  if len(args.params) != ntemp+nrad+nabund:
    pt.error("The input number of fitting parameters ({:d}) does not "
             "match the number of model parameters ({:d}).".
              format(len(args.params), ntemp+nrad+nabund))


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
  temp = pyrat.tmodel(params[pyrat.itemp], *pyrat.targs)
  if np.any(temp < pyrat.tlow) or np.any(temp > pyrat.thigh):
    pyrat.obs.bandflux[:] = -1e10  # FINDME: what if np.inf? or nan?
    return pyrat.obs.bandflux
  # Update abundance profiles:
  q2 = qs.qscale(pyrat.atm.q, pyrat.mol.name, params[pyrat.iabund],
                 pyrat.molscale, pyrat.bulk,
                 iscale=pyrat.iscale, ibulk=pyrat.ibulk,
                 bratio=pyrat.bulkratio, invsrat=pyrat.invsrat)
  # Update radius profile:
  if len(pyrat.irad) > 0:
    radius = ma.hydro_equilibrium(pyrat.atm.press, temp, pyrat.atm.m,
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
      bflux[i] = w.bandintegrate(spectrum[wnidx[i]],
                          wn[wnidx[i]], pyrat.nifilter[i])
    elif pyrat.od.path == "eclipse":
      fluxrat = (spectrum[wnidx[i]]/pyrat.istarfl[i]) * pyrat.rprs**2.0
      bflux[i] = w.bandintegrate(fluxrat, wn[wnidx[i]], pyrat.nifilter[i])

  # Revert changes in the atmospheric profile:
  if freeze:
    pyrat.atm.temp, pyrat.atm.q, pyrat.atm.radius = t0, q0, r0
  return pyrat.obs.bandflux
