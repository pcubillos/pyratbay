# Copyright (c) 2016-2017 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import numpy as np
import time
import matplotlib.pyplot as plt
import scipy.interpolate as sip

from .  import readatm   as ra
from .. import tools     as pt
from .. import constants as pc

def makewavenumber(pyrat):
  """
  Make the wavenumber sample from user inputs.
  Store all values in CGS units.
  """
  # Alias for Pyrat's Spectrum object:
  spec = pyrat.spec

  pt.msg(pyrat.verb-3, "\nGenerating wavenumber array.", pyrat.log)
  # Low wavenumber boundary:
  if spec.wnlow is None:
    if spec.wlhigh is None:
      pt.error("Low wavenumber boundary is undefined.  Set either wnlow or "
               "wlhigh.", pyrat.log)
    else:
      spec.wnlow = 1.0 / spec.wlhigh
  elif spec.wlhigh is not None:
    pt.warning(pyrat.verb-2, "Both wnlow ({:.2e} cm-1) and wlhigh ({:.2e} cm) "
        "were defined.  Pyrat will take wnlow and ignore wlhigh.".
         format(spec.wnlow, spec.wlhigh), pyrat.log, pyrat.wlog)

  # High wavenumber boundary:
  if spec.wnhigh is None:
    if spec.wllow is None:
      pt.error("High wavenumber boundary is undefined.  Set either wnhigh or "
               "wllow.", pyrat.log)
    else:
      spec.wnhigh = 1.0 / spec.wllow
  elif spec.wllow is not None:
    pt.warning(pyrat.verb-2, "Both wnhigh ({:.2e} cm-1) and wllow ({:.2e} cm) "
        "were defined.  Pyrat will take wnhigh and ignore wllow.".
         format(spec.wnhigh, spec.wllow), pyrat.log, pyrat.wlog)

  # Consistency check (wnlow < wnhigh):
  if spec.wnlow > spec.wnhigh:
    pt.error("Wavenumber low boundary ({:.2e} cm-1) must be larger than the "
             "high boundary ({:.2e} cm-1)".
              format(spec.wnlow, spec.wnhigh), pyrat.log)

  # Set wavelength limits based on the wavenumber limits:
  spec.wlhigh = 1.0 / spec.wnlow
  spec.wllow  = 1.0 / spec.wnhigh

  # Make the wavenumber array:
  spec.wn = np.arange(spec.wnlow, spec.wnhigh, spec.wnstep)

  # Re-set final boundary (stay inside given boundaries):
  if spec.wn[-1] != spec.wnhigh:
    pt.warning(pyrat.verb-2, "Final wavenumber boundary modified from "
                                               "{:10.4f} cm-1 (input)\n"
       "                                     to {:10.4f} cm-1 (Pyrat).".
               format(spec.wnhigh, spec.wn[-1]), pyrat.log, pyrat.wlog)
  # Set the number of spectral samples:
  spec.nwave  = len(spec.wn)

  # Make the fine-sampled (oversampled) wavenumber array:
  spec.ownstep = spec.wnstep / spec.wnosamp
  spec.onwave  = (spec.nwave - 1) *  spec.wnosamp + 1
  spec.own = np.linspace(spec.wn[0], spec.wn[-1], spec.onwave)

  # Get list of divisors:
  spec.odivisors = pt.divisors(spec.wnosamp)

  # Screen output:
  pt.msg(pyrat.verb-4, "Initial wavenumber boundary:  {:.5e} cm-1  ({:.3e} "
                    "{:s})".format(spec.wnlow, spec.wlhigh/pt.u(spec.wlunits),
                                   spec.wlunits), pyrat.log, 2)
  pt.msg(pyrat.verb-4, "Final   wavenumber boundary:  {:.5e} cm-1  ({:.3e} "
                    "{:s})".format(spec.wnhigh, spec.wllow/pt.u(spec.wlunits),
                                   spec.wlunits), pyrat.log, 2)
  pt.msg(pyrat.verb-4, "Wavenumber sampling stepsize: {:.2g} cm-1".
                            format(spec.wnstep),  pyrat.log, 2)
  pt.msg(pyrat.verb-4, "Wavenumber sample size:      {:8d}".format(spec.nwave),
                                                  pyrat.log,  2)
  pt.msg(pyrat.verb-4, "Wavenumber fine-sample size: {:8d}".format(spec.onwave),
                                                  pyrat.log, 2)
  pt.msg(pyrat.verb-3, "Done.", pyrat.log)


def makeradius(pyrat):
  """
  Generate atmospheric-profile layers sampling.

  If the user sets radius command-line arguments, make a sampling; else,
  default to the atmfile radius array.
  """
  pt.msg(pyrat.verb-3, "\nGenerating atmospheric profile sample.", pyrat.log)

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
  pt.msg(pyrat.verb-4, "Atmospheric file pressure limits: {:.2e}--{:.2e} {:s}.".
     format(atm_in.press[ 0]/pt.u(atm_in.punits),
            atm_in.press[-1]/pt.u(atm_in.punits), atm_in.punits), pyrat.log, 2)

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
    pt.warning(pyrat.verb-2, "The atmospheric layers are in reversed order "
     "(bottom-top).  Resorting to be from the top down.", pyrat.log, pyrat.wlog)
    if atm_in.radius is not None:
      atm_in.radius = atm_in.radius[::-1]
    atm_in.press  = atm_in.press [::-1]
    atm_in.temp   = atm_in.temp  [::-1]
    atm_in.mm     = atm_in.mm    [::-1]
    atm_in.q      = np.flipud(atm_in.q)
    atm_in.d      = np.flipud(atm_in.d)
  else:
    pt.error("The atmospheric layers are neither sorted from the bottom up, "
             "nor from the top down.", pyrat.log)

  if atm_in.radius is None:
    # Check that gplanet exists:
    if pyrat.phy.rplanet is None:
      pt.error("Undefined reference planetary radius (rplanet). Either include "
        "the radius profile in the atmospheric file or set rplanet.", pyrat.log)
    if pyrat.phy.gplanet is None:
      pt.error("Undefined atmospheric gravity (gplanet).  Either include "
        "the radius profile in the atmospheric file, set the surface "
        "gravity, or set the planetary mass (mplanet).", pyrat.log)
    if pyrat.refpressure is None:
      pt.error("Undefined reference pressure level (refpressure). Either "
        "include the radius profile in the atmospheric file or set refpress.",
        pyrat.log)
    # Atmopsheric reference pressure-radius level:
    pt.msg(pyrat.verb-4, "Reference pressure: {:.3e} {:s}.".
           format(pyrat.refpressure/pt.u(pyrat.punits), pyrat.punits),
           pyrat.log, 2)
    pt.msg(pyrat.verb-4, "Reference radius: {:8.1f} {:s}.".
           format(pyrat.phy.rplanet/pt.u(pyrat.radunits), pyrat.radunits),
           pyrat.log, 2)
    if not np.isinf(pyrat.phy.rhill):
      pt.msg(pyrat.verb-4, "Hill radius:      {:8.1f} {:s}.".
           format(pyrat.phy.rhill/pt.u(pyrat.radunits), pyrat.radunits),
           pyrat.log, 2)

    # Calculate the radius profile using the hydostatic-equilibrium equation:
    atm_in.radius = pyrat.hydro(atm_in.press, atm_in.temp, atm_in.mm,
                                pyrat.phy.gplanet, pyrat.phy.mplanet,
                                pyrat.refpressure, pyrat.phy.rplanet)

  # Check if Hydrostatic Eq. breaks down:
  ibreak = 0  # Index and flag at the same time
  if np.any(np.ediff1d(atm_in.radius) > 0.0):
    ibreak = np.where(np.ediff1d(atm_in.radius)>0)[0][0] + 1

  # Set the interpolating function (for use later):
  radinterp   = sip.interp1d(atm_in.press [ibreak:],
                             atm_in.radius[ibreak:], kind='slinear')
  pressinterp = sip.interp1d(np.flipud(atm_in.radius[ibreak:]),
                             np.flipud(atm_in.press [ibreak:]), kind='slinear')

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
    pt.error("Unbounded atmosphere.  Hydrostatic-equilibrium radius solution "
             "diverges at pressure {:.3e} bar.  Set mstar and smaxis to "
             "define a Hill radius (top boundary) and avoid error.".
              format(atm_in.press[ibreak]/pc.bar))

  # Out of bounds errors:
  if pyrat.plow < np.amin(atm_in.press):
    pt.error("User-defined bottom layer (p={:.3e} {:s}) is lower than the "
             "atmospheric-file bottom layer (p={:.3e} {:s}).".
             format(pyrat.plow/pt.u(pyrat.punits), pyrat.punits,
             np.amin(atm_in.press)/pt.u(pyrat.punits), pyrat.punits), pyrat.log)
  if pyrat.phigh > np.amax(atm_in.press):
    pt.error("User-defined top layer (p={:.3e} {:s}) is higher than the "
             "atmospheric-file top layer (p={:.3e} {:s}).".
             format(pyrat.phigh/pt.u(pyrat.punits), pyrat.punits,
             np.amax(atm_in.press)/pt.u(pyrat.punits), pyrat.punits), pyrat.log)

  pt.msg(pyrat.verb-4, "User pressure boundaries: {:.2e}--{:.2e} bar.".
         format(pyrat.plow/pc.bar, pyrat.phigh/pc.bar), pyrat.log, 2)

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
    # Avoid radinterp going out-of-bounds:
    # Interpolate to pressure array:
    atm.press = pressinterp(atm.radius)
    resample = True

  # Else, take the atmospheric-file sampling:
  else:
    # Get top-bottom indices:
    ilow  = np.where(atm_in.press >= pyrat.plow) [0][-1]
    ihigh = np.where(atm_in.press <= pyrat.phigh)[0][ 0]
    # Take values within the boundaries:
    atm.press   = atm_in.press [ihigh:ilow+1]
    atm.radius  = atm_in.radius[ihigh:ilow+1]
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
    pt.warning(pyrat.verb-2, "The atmospheric pressure array extends "
       "beyond the Hill radius ({:.1f} km) at pressure {:.2e} bar (layer "
       "#{:d}).  Extinction beyond this layer will be neglected.".
        format(pyrat.phy.rhill/pc.km, atm_in.press[pyrat.atm.rtop]/pc.bar,
               pyrat.atm.rtop), pyrat.log, pyrat.wlog)

  # Radius-vs-pressure from Atm. file and resampled array:
  # plt.figure(2)
  # plt.clf()
  # plt.semilogx(atm_in.press /pt.u(pyrat.punits),
  #              atm_in.radius/pt.u(pyrat.radunits),
  #              "o-r", mec="r", mfc='r')
  # plt.semilogx(atm.press /pt.u(pyrat.punits),
  #              atm.radius/pt.u(pyrat.radunits), "o-b",
  #              mec="b", mew=1, mfc='None')
  # plt.xlabel("Pressure  ({:s})".format(pyrat.punits))
  # plt.ylabel("Radius  ({:s})".format(pyrat.radunits))
  # plt.savefig("radpress.png")

  # Print radius array:
  radstr = '['+', '.join('{:9.2f}'.format(k) for k in atm_in.radius/pc.km)+']'
  pt.msg(pyrat.verb-4, "Radius array (km) =   {:s}".format(radstr),
         pyrat.log, 2, si=4)

  pt.msg(pyrat.verb-4, "Number of valid model layers: {:d}.".
                        format(atm.nlayers-atm.rtop), pyrat.log, 2)
  pt.msg(pyrat.verb-4, "Valid lower/higher pressure boundaries: {:.2e} - "
         "{:.2e} {:s}.".format(atm.press[atm.rtop]/pt.u(pyrat.punits),
                 pyrat.phigh/pt.u(pyrat.punits), pyrat.punits), pyrat.log, 2)
  pt.msg(pyrat.verb-4, "Valid upper/lower radius boundaries:    {:8.1f} - {:8.1f} "
         "{:s}.".format(atm.radius[atm.rtop]/pt.u(pyrat.radunits),
          atm.radius[-1]/pt.u(pyrat.radunits), pyrat.radunits), pyrat.log, 2)

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
    pt.error("The layer {:d} in the atmospheric model has a lower temperature "
             "({:.1f} K) than the lowest allowed TLI temperature ({:.1f} K).".
              format(icold, atm.temp[icold], pyrat.lt.tmin), pyrat.log)
  if np.any(atm.temp > pyrat.lt.tmax):
    ihot = np.where(atm.temp > pyrat.lt.tmax)[0][0]
    pt.error("The layer {:d} in the atmospheric model has a higher temperature "
             "({:.1f} K) than the highest allowed TLI temperature ({:.1f} K).".
              format(ihot, atm.temp[ihot], pyrat.lt.tmax), pyrat.log)

  # Interpolate isotopes partition function:
  pt.msg(pyrat.verb-4, "Number of isotopes: {:d}".format(pyrat.iso.niso),
         pyrat.log, 2)
  # Initialize the partition-function array for pyrat.iso:
  pyrat.iso.z = np.zeros((pyrat.iso.niso, atm.nlayers))
  for i in np.arange(pyrat.lt.ndb):           # For each Database
    for j in np.arange(pyrat.lt.db[i].niso):  # For each isotope in DB
      pt.msg(pyrat.verb-5, "Interpolating (isotope ID {:2d}) partition "
                        "function.".format(pyrat.lt.db[i].iiso+j), pyrat.log, 4)
      zinterp = sip.interp1d(pyrat.lt.db[i].temp, pyrat.lt.db[i].z[j],
                             kind='slinear')
      pyrat.iso.z[pyrat.lt.db[i].iiso+j] = zinterp(atm.temp)

  # # Plot interpolation:
  # plt.figure(3)
  # plt.clf()
  # plt.plot(pyrat.lt.db[0].temp, pyrat.lt.db[0].z[0], "o-r", mec="r")
  # plt.plot(atm.temp, pyrat.iso.z[0], "o-b", mec="b", mfc='None')
  # plt.xlabel("Temperature  (K)")
  # plt.ylabel("Partition function")
  # plt.xlim(0, 3000)
  # plt.ylim(0, 15000)
  # plt.savefig("PartitionFunction.png")

  pt.msg(pyrat.verb-3, "Done.", pyrat.log)
