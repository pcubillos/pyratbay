# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import numpy as np
import scipy.interpolate as sip

from .. import tools     as pt
from .. import constants as pc


def makewavenumber(pyrat):
  """
  Make the wavenumber sample from user inputs.
  Store all values in CGS units.
  """
  # Alias for Pyrat's Spectrum object:
  spec = pyrat.spec

  pyrat.log.msg("\nGenerating wavenumber array.")
  # Low wavenumber boundary:
  if spec.wnlow is None:
    if spec.wlhigh is None:
      pyrat.log.error("Low wavenumber boundary is undefined.  Either set "
                      "wnlow or wlhigh.")
    else:
      spec.wnlow = 1.0 / spec.wlhigh
  elif spec.wlhigh is not None:
    pyrat.log.warning("Both wnlow ({:.2e} cm-1) and wlhigh ({:.2e} cm) were "
                      "defined.  Pyrat will take wnlow and ignore wlhigh.".
                      format(spec.wnlow, spec.wlhigh))

  # High wavenumber boundary:
  if spec.wnhigh is None:
    if spec.wllow is None:
      pyrat.log.error("High wavenumber boundary is undefined.  Either set "
                      "wnhigh or wllow.")
    else:
      spec.wnhigh = 1.0 / spec.wllow
  elif spec.wllow is not None:
    pyrat.log.warning("Both wnhigh ({:.2e} cm-1) and wllow ({:.2e} cm) were "
                      "defined.  Pyrat will take wnhigh and ignore wllow.".
                       format(spec.wnhigh, spec.wllow))

  # Consistency check (wnlow < wnhigh):
  if spec.wnlow > spec.wnhigh:
    pyrat.log.error("Wavenumber low boundary ({:.2e} cm-1) must be larger "
                    "than the high boundary ({:.2e} cm-1)".
                    format(spec.wnlow, spec.wnhigh))

  # Set wavelength limits based on the wavenumber limits:
  spec.wlhigh = 1.0 / spec.wnlow
  spec.wllow  = 1.0 / spec.wnhigh

  # Make the wavenumber array:
  spec.wn = np.arange(spec.wnlow, spec.wnhigh, spec.wnstep)

  # Re-set final boundary (stay inside given boundaries):
  if spec.wn[-1] != spec.wnhigh:
    pyrat.log.warning(
       "Final wavenumber boundary modified from {:10.4f} cm-1 (input)\n"
       "                                     to {:10.4f} cm-1 (Pyrat).".
       format(spec.wnhigh, spec.wn[-1]))
  # Set the number of spectral samples:
  spec.nwave  = len(spec.wn)

  # Make the fine-sampled (oversampled) wavenumber array:
  spec.ownstep = spec.wnstep / spec.wnosamp
  spec.onwave  = (spec.nwave - 1) *  spec.wnosamp + 1
  spec.own = np.linspace(spec.wn[0], spec.wn[-1], spec.onwave)

  # Get list of divisors:
  spec.odivisors = pt.divisors(spec.wnosamp)

  # Screen output:
  pyrat.log.msg("Initial wavenumber boundary:  {:.5e} cm-1  ({:.3e} "
                "{:s})".format(spec.wnlow, spec.wlhigh/pt.u(spec.wlunits),
                               spec.wlunits), verb=2, indent=2)
  pyrat.log.msg("Final   wavenumber boundary:  {:.5e} cm-1  ({:.3e} "
                "{:s})".format(spec.wnhigh, spec.wllow/pt.u(spec.wlunits),
                               spec.wlunits), verb=2, indent=2)
  pyrat.log.msg("Wavenumber sampling stepsize: {:.2g} cm-1\n"
                "Wavenumber sample size:      {:8d}\n"
                "Wavenumber fine-sample size: {:8d}\n".format(
                spec.wnstep, spec.nwave, spec.onwave), verb=2, indent=2, si=2)
  pyrat.log.msg("Wavenumber sampling done.")


def make_atmprofiles(pyrat):
  """
  Define atmospheric-profile layers sampling.

  Notes
  -----
  There are multiple things going on in this functions, here's a summary:
  - Read input atmospheric file.
  - Compute hydrostatic-equilibrium radius profile if necessary.
  - Reset (pressure/radius) boundaries if requested.
  - Resample (pressure/radius/abundance/temperature) profiles if requested.
  - Check whether top of atmosphere crosses the Hill radius.
  - Compute partition-function at layers temperatures.
  """
  pyrat.log.msg("\nGenerating atmospheric profile sample.")

  # Pyrat and user-input atmospheric-data objects:
  atm    = pyrat.atm
  atm_in = pyrat.inputs.atm

  # Copy atmospheric-model units:
  atm.tunits = "kelvin"
  atm.qunits = atm_in.qunits
  atm.punits = atm_in.punits
  atm.runits = atm_in.runits

  # FINDME: move this to readatm
  # Pressure limits from the atmospheric file:
  pyrat.log.msg("Atmospheric file pressure limits: {:.2e}--{:.2e} {:s}.".format(
         atm_in.press[ 0]/pt.u(atm_in.punits),
         atm_in.press[-1]/pt.u(atm_in.punits), atm_in.punits), verb=2, indent=2)

  # Check that the layers are sorted from the top to the bottom of
  #  the atmosphere:
  sort    = np.all(np.ediff1d(atm_in.press) > 0)  # Top to bottom
  reverse = np.all(np.ediff1d(atm_in.press) < 0)  # Bottom to top
  if atm_in.radius is not None:
    sort    *= np.all(np.ediff1d(atm_in.radius) < 0)
    reverse *= np.all(np.ediff1d(atm_in.radius) > 0)

  if   sort:     # Layers are in the correct order
    pass
  elif reverse:  # Layers in reverse order
    pyrat.warning("The atmospheric layers are in reversed order "
                  "(bottom-top).  Resorting to be from the top down.")
    if atm_in.radius is not None:
      atm_in.radius = atm_in.radius[::-1]
    atm_in.press  = atm_in.press [::-1]
    atm_in.temp   = atm_in.temp  [::-1]
    atm_in.mm     = atm_in.mm    [::-1]
    atm_in.q      = np.flipud(atm_in.q)
    atm_in.d      = np.flipud(atm_in.d)
  else:
    pyrat.log.error("The atmospheric layers are neither sorted from the "
                    "bottom up, nor from the top down.")

  if atm_in.radius is None:
    # Check that gplanet exists:
    if pyrat.phy.rplanet is None and pyrat.runmode != "opacity":
      pyrat.log.error("Undefined reference planetary radius (rplanet). Either "
          "include the radius profile in the atmospheric file or set rplanet.")
    if pyrat.phy.gplanet is None and pyrat.runmode != "opacity":
      pyrat.log.error("Undefined atmospheric gravity (gplanet).  Either "
          "include the radius profile in the atmospheric file, set the surface "
          "gravity, or set the planetary mass (mplanet).")
    if pyrat.refpressure is None and pyrat.runmode != "opacity":
      pyrat.log.error("Undefined reference pressure level (refpressure). "
          "Either include the radius profile in the atmospheric file or set "
          "refpress.")
    if (pyrat.phy.rplanet is not None and
        pyrat.phy.gplanet is not None and
        pyrat.refpressure is not None):
      # Atmopsheric reference pressure-radius level:
      pyrat.log.msg("Reference pressure: {:.3e} {:s}.".
             format(pyrat.refpressure/pt.u(pyrat.punits), pyrat.punits),
             verb=2, indent=2)
      pyrat.log.msg("Reference radius: {:8.1f} {:s}.".
             format(pyrat.phy.rplanet/pt.u(pyrat.radunits), pyrat.radunits),
             verb=2, indent=2)
      if not np.isinf(pyrat.phy.rhill):
        pyrat.log.msg("Hill radius:      {:8.1f} {:s}.".
             format(pyrat.phy.rhill/pt.u(pyrat.radunits), pyrat.radunits),
             verb=2, indent=2)

      # Calculate the radius profile using the hydostatic-equilibrium equation:
      atm_in.radius = pyrat.hydro(atm_in.press, atm_in.temp, atm_in.mm,
                                  pyrat.phy.gplanet, pyrat.phy.mplanet,
                                  pyrat.refpressure, pyrat.phy.rplanet)

  # Check if Hydrostatic Eq. breaks down:
  ibreak = 0  # Index and flag at the same time
  if np.any(np.ediff1d(atm_in.radius) > 0.0):
    ibreak = np.where(np.ediff1d(atm_in.radius)>0)[0][0] + 1

  # Set the interpolating function (for use later):
  try:
    radinterp   = sip.interp1d(atm_in.press [ibreak:],
                             atm_in.radius[ibreak:], kind='slinear')
    pressinterp = sip.interp1d(np.flipud(atm_in.radius[ibreak:]),
                             np.flipud(atm_in.press [ibreak:]), kind='slinear')
  except:
    pass

  # Set radius/pressure boundaries if exist:
  if pyrat.plow is not None:
    if pyrat.plow >= atm_in.press[ibreak]:
      ibreak = 0   # Turn-off break flag
  elif pyrat.radhigh is not None:
    if pyrat.radhigh <= atm_in.radius[ibreak]:
      ibreak = 0   # Turn-off break flag
      pyrat.plow = pressinterp(pyrat.radhigh)[0]
    #else:
    #  out-of-bounds error
  else:
    pyrat.plow = np.amin(atm_in.press)

  if pyrat.radlow is not None:
    pyrat.phigh = pressinterp(pyrat.radlow)[0]
  elif pyrat.phigh is None:
    pyrat.phigh = np.amax(atm_in.press)

  if ibreak != 0 and np.isinf(pyrat.phy.rhill):
    pyrat.log.error("Unbounded atmosphere.  Hydrostatic-equilibrium radius "
             "solution diverges at pressure {:.3e} bar.  Set mstar and smaxis "
             "to define a Hill radius (top boundary) and avoid error.".
              format(atm_in.press[ibreak]/pc.bar))

  # Out of bounds errors:
  if pyrat.plow < np.amin(atm_in.press):
    pyrat.log.error("User-defined bottom layer (p={:.3e} {:s}) is lower than "
             "the atmospheric-file bottom layer (p={:.3e} {:s}).".
             format(pyrat.plow/pt.u(pyrat.punits), pyrat.punits,
             np.amin(atm_in.press)/pt.u(pyrat.punits), pyrat.punits))
  if pyrat.phigh > np.amax(atm_in.press):
    pyrat.log.error("User-defined top layer (p={:.3e} {:s}) is higher than "
             "the atmospheric-file top layer (p={:.3e} {:s}).".
             format(pyrat.phigh/pt.u(pyrat.punits), pyrat.punits,
             np.amax(atm_in.press)/pt.u(pyrat.punits), pyrat.punits))

  pyrat.log.msg("User pressure boundaries: {:.2e}--{:.2e} bar.".
                format(pyrat.plow/pc.bar, pyrat.phigh/pc.bar), verb=2, indent=2)

  # Resample to equispaced log-pressure array if requested:
  if atm.nlayers is not None:
    atm.press = np.logspace(np.log10(pyrat.plow), np.log10(pyrat.phigh),
                            atm.nlayers)
    atm.radius = radinterp(atm.press)
    resample = True

  # Resample to equispaced radius array if requested:
  elif pyrat.radstep is not None:
    # Make equispaced radius array:
    if pyrat.radlow is None:
      pyrat.radlow  = radinterp(pyrat.phigh)[0]
    if pyrat.radhigh is None:
      pyrat.radhigh = radinterp(pyrat.plow)[0]

    atm.radius = np.arange(pyrat.radlow, pyrat.radhigh, pyrat.radstep)
    atm.nlayers = len(atm.radius)
    # Interpolate to pressure array:
    atm.press = pressinterp(atm.radius)
    resample = True

  else:  # Take the atmospheric-file sampling:
    # Get top-bottom indices:
    ilow  = np.where(atm_in.press >= pyrat.plow) [0][-1]
    ihigh = np.where(atm_in.press <= pyrat.phigh)[0][ 0]
    # Take values within the boundaries:
    atm.press   = atm_in.press [ihigh:ilow+1]
    if atm_in.radius is not None:
      atm.radius = atm_in.radius[ihigh:ilow+1]
    atm.temp    = atm_in.temp  [ihigh:ilow+1]
    atm.mm      = atm_in.mm    [ihigh:ilow+1]
    atm.q       = atm_in.q     [ihigh:ilow+1]
    atm.d       = atm_in.d     [ihigh:ilow+1]
    atm.nlayers = len(atm.press)
    resample = False

  # Check the radii lie within Hill radius:
  rtop = np.where(atm.radius > pyrat.phy.rhill)[0]
  if np.size(rtop) > 0:
    pyrat.atm.rtop = rtop[-1] + 1
    pyrat.log.warning("The atmospheric pressure array extends "
       "beyond the Hill radius ({:.1f} km) at pressure {:.2e} bar (layer "
       "#{:d}).  Extinction beyond this layer will be neglected.".
        format(pyrat.phy.rhill/pc.km, atm_in.press[pyrat.atm.rtop]/pc.bar,
               pyrat.atm.rtop))

  # Print radius array:
  if atm.radius is not None:
    radstr = '['+', '.join('{:9.2f}'.format(k) for k in atm.radius/pc.km)+']'
    pyrat.log.msg("Radius array (km) =   {:s}".format(radstr),
                  verb=2, indent=2, si=4)
    pyrat.log.msg("Valid upper/lower radius boundaries:    {:8.1f} - {:8.1f} "
         "{:s}.".format(atm.radius[atm.rtop]/pt.u(pyrat.radunits),
         atm.radius[-1]/pt.u(pyrat.radunits), pyrat.radunits), verb=2, indent=2)

  pyrat.log.msg("Valid lower/higher pressure boundaries: {:.2e} - "
                "{:.2e} {:s}.".format(atm.press[atm.rtop]/pt.u(pyrat.punits),
                pyrat.phigh/pt.u(pyrat.punits), pyrat.punits), verb=2, indent=2)
  pyrat.log.msg("Number of valid model layers: {:d}.".
                format(atm.nlayers-atm.rtop), verb=2, indent=2)

   # Interpolate to new atm-layer sampling if necessary:
  if resample:
    tempinterp = sip.interp1d(atm_in.press, atm_in.temp, kind='slinear')
    mminterp   = sip.interp1d(atm_in.press, atm_in.mm,   kind='slinear')
    atm.temp = tempinterp(atm.press)
    atm.m    =   mminterp(atm.press)
    # Interpolate abundance profiles:
    atm.q = np.zeros((atm.nlayers, pyrat.mol.nmol))
    atm.d = np.zeros((atm.nlayers, pyrat.mol.nmol))
    for i in np.arange(pyrat.mol.nmol):
      qinterp = sip.interp1d(atm_in.press, atm_in.q[:, i], kind='slinear')
      dinterp = sip.interp1d(atm_in.press, atm_in.d[:, i], kind='slinear')
      atm.q[:,i] = qinterp(atm.press)
      atm.d[:,i] = dinterp(atm.press)

  # Temperature boundaries check:
  if np.any(atm.temp < pyrat.lt.tmin):
    icold = np.where(atm.temp < pyrat.lt.tmin)[0][0]
    pyrat.log.error("The layer {:d} in the atmospheric model has a lower "
             "temperature ({:.1f} K) than the lowest allowed TLI temperature "
             "({:.1f} K).".format(icold, atm.temp[icold], pyrat.lt.tmin))
  if np.any(atm.temp > pyrat.lt.tmax):
    ihot = np.where(atm.temp > pyrat.lt.tmax)[0][0]
    pyrat.log.error("The layer {:d} in the atmospheric model has a higher "
             "temperature ({:.1f} K) than the highest allowed TLI temperature "
             "({:.1f} K).".format(ihot, atm.temp[ihot], pyrat.lt.tmax))

  # Interpolate isotopes partition function:
  pyrat.log.msg("Number of isotopes: {:d}".format(pyrat.iso.niso),
                verb=2, indent=2)
  # Initialize the partition-function array for pyrat.iso:
  pyrat.iso.z = np.zeros((pyrat.iso.niso, atm.nlayers))
  for i in np.arange(pyrat.lt.ndb):           # For each Database
    for j in np.arange(pyrat.lt.db[i].niso):  # For each isotope in DB
      pyrat.log.msg("Interpolating (isotope ID {:2d}) partition "
         "function.".format(pyrat.lt.db[i].iiso+j), verb=3, indent=4)
      zinterp = sip.interp1d(pyrat.lt.db[i].temp, pyrat.lt.db[i].z[j],
                             kind='slinear')
      pyrat.iso.z[pyrat.lt.db[i].iiso+j] = zinterp(atm.temp)

  pyrat.log.msg("Make atmosphere done.")
