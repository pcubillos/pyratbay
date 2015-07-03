import numpy as np
import time
import matplotlib.pyplot as plt
import scipy.constants   as sc
import scipy.integrate   as si
import scipy.interpolate as sip

import readatm    as ra
import pconstants as pc
import ptools     as pt


def makewavenumber(pyrat):
  """
  Make the wavenumber sample from user inputs.
  Store all values in CGS units.
  """
  # Alias for Pyrat's Spectrum object:
  spec = pyrat.spec

  pt.msg(pyrat.verb, "\nGenerating wavenumber array:", 0)
  # Low wavenumber boundary:
  if spec.wnlow is None:
    if spec.wlhigh is None:
      pt.error("Low wavenumber boundary is undefined.  Set either wnlow or "
               "wlhigh.")
    else:
      spec.wnlow = 1.0 / spec.wlhigh
  elif spec.wlhigh is not None:
    pt.warning("Both wnlow ({:.2e} cm-1) and wlhigh ({:.2e} cm) were "
               "defined.  Pyrat will take wnlow and ignore wlhigh.".
                format(spec.wnlow, spec.wlhigh))

  # High wavenumber boundary:
  if spec.wnhigh is None:
    if spec.wllow is None:
      pt.error("High wavenumber boundary is undefined.  Set either wnhigh or "
               "wllow.")
    else:
      spec.wnhigh = 1.0 / spec.wllow
  elif spec.wllow is not None:
    pt.warning("Both wnhigh ({:.2e} cm-1) and wllow ({:.2e} cm) were "
               "defined.  Pyrat will take wnhigh and ignore wllow.".
                format(spec.wnhigh, spec.wllow))

  # Consistency check (wnlow < wnhigh):
  if spec.wnlow > spec.wnhigh:
    pt.error("Wavenumber low boundary ({:.2e} cm-1) must be larger than the "
             "high boundary ({:.2e} cm-1)".
              format(spec.wnlow, spec.wnhigh))

  # Set wavelength limits based on the wavenumber limits:
  spec.wlhigh = 1.0 / spec.wnlow
  spec.wllow  = 1.0 / spec.wnhigh

  # Make the wavenumber array:
  spec.wn = np.arange(spec.wnlow, spec.wnhigh, spec.wnstep)

  # Re-set final boundary (stay inside given boundaries):
  if spec.wn[-1] != spec.wnhigh:
    pt.warning("Final wavenumber boundary modified from {:.4f} cm-1 (input)\n"
               "                                     to {:.4f} cm-1 (Pyrat).".
               format(spec.wnhigh, spec.wn[-1]))
  # Set the number of spectral samples:
  spec.nspec  = len(spec.wn)

  # Make the fine-sampled (oversampled) wavenumber array:
  spec.ownstep = spec.wnstep / spec.wnosamp
  spec.onspec  = (spec.nspec - 1) *  spec.wnosamp + 1
  spec.own = np.linspace(spec.wn[0], spec.wn[-1], spec.onspec)

  # Get list of divisors:
  spec.odivisors = pt.divisors(spec.wnosamp)

  # Screen output:
  pt.msg(pyrat.verb,"Initial wavenumber boundary:  {:.5e} cm-1  ({:.3e} "
                    "{:s})".format(spec.wnlow,
                          spec.wlhigh/pc.units[spec.wlunits], spec.wlunits), 2)
  pt.msg(pyrat.verb,"Final   wavenumber boundary:  {:.5e} cm-1  ({:.3e} "
                    "{:s})".format(spec.wnhigh,
                          spec.wllow /pc.units[spec.wlunits], spec.wlunits), 2)
  pt.msg(pyrat.verb,"Wavenumber sampling stepsize: {:.2g} cm-1".
                            format(spec.wnstep), 2)
  pt.msg(pyrat.verb,"Wavenumber sample size:      {:8d}".format(spec.nspec),  2)
  pt.msg(pyrat.verb,"Wavenumber fine-sample size: {:8d}".format(spec.onspec), 2)
  pt.msg(pyrat.verb, "Done.")


def makeradius(pyrat):
  """
  Generate atmospheric radius layers sampling.

  If the user sets radius command-line arguments, make a sampling; else,
  default to the atmfile radius array.
  """
  pt.msg(pyrat.verb, "\nGenerating atmospheric radius sample:")

  # Atmopsheric reference pressure-radius level:
  pt.msg(pyrat.verb, "Reference pressure: {:.3e} {:s}".format(
         pyrat.pressurebase/pc.units[pyrat.punits], pyrat.punits),   2)
  pt.msg(pyrat.verb, "Reference radius: {:8g} {:s}".format(
         pyrat.radiusbase/pc.units[pyrat.radunits], pyrat.radunits), 2)

  # FINDME: move this to readatm
  # Pressure limits from the atmospheric file:
  pt.msg(pyrat.verb, "Pressure limits: {:.3e} -- {:.3e} {:s}".format(
        pyrat.atmf.press[ 0]/pc.units[pyrat.punits],
        pyrat.atmf.press[-1]/pc.units[pyrat.punits], pyrat.punits), 2)

  if pyrat.atmf.radius is None:
    # Calculate radius array for given atmospheric profile by using the 
    # hydostatic-equilibrium equation:
    pyrat.atmf.radius = si.cumtrapz((-pc.k * sc.N_A * pyrat.atmf.temp) /
               (pyrat.atmf.mm * pyrat.surfgravity), np.log(pyrat.atmf.press))
    pyrat.atmf.radius  = np.concatenate(([0.0], pyrat.atmf.radius))
    pyrat.atmf.radius -= pyrat.atmf.radius[-1]

    # Get radius at the reference pressure level:
    if pyrat.pressurebase is not None and pyrat.radiusbase is not None:
      radinterp = sip.interp1d(pyrat.atmf.press, pyrat.atmf.radius,
                               kind='slinear')
      r0 = radinterp(pyrat.pressurebase)
      # Make correction to force: radius(ref. pressure) = ref. radius
      pyrat.atmf.radius += pyrat.radiusbase - r0

  # Reset the interpolating function (for use later):
  radinterp   = sip.interp1d(pyrat.atmf.press, pyrat.atmf.radius,
                             kind='slinear')
  pressinterp = sip.interp1d(pyrat.atmf.radius[::-1],
                             pyrat.atmf.press [::-1], kind='slinear')
  pt.msg(pyrat.verb, "Radius array (km) = {:s}".
                      format(pt.pprint(pyrat.atmf.radius/pc.units["km"],2)), 2)

  # Set pressure boundaries:
  if pyrat.radhigh is not None:
    pyrat.plow = pressinterp(pyrat.radhigh)[0]
  elif pyrat.phigh is None:
    pyrat.phigh = np.amax(pyrat.atmf.press)

  if pyrat.radlow is not None:
    pyrat.phigh = pressinterp(pyrat.radlow)[0]
  elif pyrat.plow is None:
    pyrat.plow  = np.amin(pyrat.atmf.press)

  pt.msg(pyrat.verb, "Pressure user boundaries: {:.3e} -- {:.3e} bar".format(
              pyrat.plow/pc.units["bar"], pyrat.phigh/pc.units["bar"]), 2)

  # Out of bounds errors:
  if pyrat.phigh > np.amax(pyrat.atmf.press):
    pt.error("User-defined top layer (p={:.3e} {:s}) is higher than the "
             "atmospheric-file top layer (p={:.3e} {:s}).".format(
              pyrat.phigh/pc.units[pyrat.punits], pyrat.punits,
              np.amax(pyrat.atmf.press)/pc.units[pyrat.punits], pyrat.punits))

  # Out of bounds errors:
  if pyrat.plow < np.amin(pyrat.atmf.press):
    pt.error("User-defined bottom layer (p={:.3e} {:s}) is lower than the "
             "atmospheric-file bottom layer (p={:.3e} {:s}).".format(
              pyrat.plow/pc.units[pyrat.punits], pyrat.punits,
              np.amin(pyrat.atmf.press)/pc.units[pyrat.punits], pyrat.punits))

  # Resample to equispaced log-pressure array if requested:
  if pyrat.atm.nlayers is not None:
    pyrat.atm.press = np.logspace(np.log10(pyrat.phigh), np.log10(pyrat.plow),
                                  pyrat.atm.nlayers)[::-1]
    pyrat.atm.radius = radinterp(pyrat.atm.press)
    resample = True

  # Resample to equispaced radius array if requested:
  elif pyrat.radstep is not None:
    # Make equispaced radius array:
    if pyrat.radlow is None:
      pyrat.radlow  = radinterp(pyrat.phigh)[0]
    if pyrat.radhigh is None:
      pyrat.radhigh = radinterp(pyrat.plow)[0]

    pyrat.atm.radius = np.arange(pyrat.radlow, pyrat.radhigh, pyrat.radstep)
    pyrat.atm.nlayers = len(pyrat.atm.radius)
    # Avoid radinterp going out-of-bounds:
    # radlow  = np.amax([rad[ 0], radinterp([pyrat.phigh])[0]])
    # radhigh = np.amin([rad[-1], radinterp([pyrat.plow ])[0]])
    # Interpolate to pressure array:
    pyrat.atm.press = pressinterp(pyrat.atm.radius)
    resample = True

  # Else, take the atmospheric-file sampling:
  else:
    # Get top-bottom indices:
    ilow  = np.where(pyrat.atmf.press >= pyrat.plow) [0][-1]
    ihigh = np.where(pyrat.atmf.press <= pyrat.phigh)[0][ 0]
    # Take values within the boundaries:
    pyrat.atm.press   = pyrat.atmf.press [ihigh:ilow+1]
    pyrat.atm.radius  = pyrat.atmf.radius[ihigh:ilow+1]
    pyrat.atm.temp    = pyrat.atmf.temp  [ihigh:ilow+1]
    pyrat.atm.mm      = pyrat.atmf.mm    [ihigh:ilow+1]
    pyrat.atm.q       = pyrat.atmf.q     [ihigh:ilow+1]
    pyrat.atm.d       = pyrat.atmf.d     [ihigh:ilow+1]
    pyrat.atm.nlayers = len(pyrat.atm.press)
    resample = False

  # Radius-vs-pressure from Atm. file and resampled array:
  # plt.figure(2)
  # plt.clf()
  # plt.semilogx(pyrat.atmf.press /pc.units[pyrat.punits],
  #              pyrat.atmf.radius/pc.units[pyrat.radunits],
  #              "o-r", mec="r", mfc='r')
  # plt.semilogx(pyrat.atm.press /pc.units[pyrat.punits],
  #              pyrat.atm.radius/pc.units[pyrat.radunits], "o-b",
  #              mec="b", mew=1, mfc='None')
  # plt.xlabel("Pressure  ({:s})".format(pyrat.punits))
  # plt.ylabel("Radius  ({:s})".format(pyrat.radunits))
  # plt.savefig("radpress.png")

  pt.msg(pyrat.verb, "Number of model layers: {:d}".format(pyrat.atm.nlayers),2)
  pt.msg(pyrat.verb, "Pressure lower/higher boundaries: {:.2e} - {:.2e} "
                     "{:s}".format(pyrat.plow /pc.units[pyrat.punits],
                           pyrat.phigh/pc.units[pyrat.punits], pyrat.punits), 2)
  pt.msg(pyrat.verb, "Radius lower/higher boundaries:   {:.1f} - {:.1f} "
         "{:s}".format(np.amin(pyrat.atm.radius)/pc.units[pyrat.radunits],
         np.amax(pyrat.atm.radius)/pc.units[pyrat.radunits], pyrat.radunits), 2)

   # Interpolate to new atm-layer sampling if necessary:
  if resample:
    tempinterp = sip.interp1d(pyrat.atmf.press, pyrat.atmf.temp,
                              kind='slinear')
    mminterp   = sip.interp1d(pyrat.atmf.press, pyrat.atmf.mm,
                              kind='slinear')
    pyrat.atm.temp = tempinterp(pyrat.atm.press)
    pyrat.atm.m    =   mminterp(pyrat.atm.press)
    # Interpolate abundance profiles:
    pyrat.atm.q = np.zeros((pyrat.atm.nlayers, pyrat.mol.nmol))
    pyrat.atm.d = np.zeros((pyrat.atm.nlayers, pyrat.mol.nmol))
    for i in np.arange(pyrat.mol.nmol):
      qinterp = sip.interp1d(pyrat.atmf.press, pyrat.atmf.q[:, i],
                             kind='slinear')
      dinterp = sip.interp1d(pyrat.atmf.press, pyrat.atmf.d[:, i],
                             kind='slinear')
      pyrat.atm.q[:,i] = qinterp(pyrat.atm.press)
      pyrat.atm.d[:,i] = dinterp(pyrat.atm.press)

  # Temperature boundaries check:
  if np.any(pyrat.atm.temp < pyrat.lt.tmin):
    icold = np.where(pyrat.atm.temp < pyrat.lt.tmin)[0][0]
    pt.error("The layer {:d} in the atmospheric model has a lower "
             "temperature ({:.1f} K) than the lowest allowed TLI temperature "
             "({:.1f} K).".format(icold, pyrat.atm.temp[icold], pyrat.lt.tmin))
  if np.any(pyrat.atm.temp > pyrat.lt.tmax):
    ihot = np.where(pyrat.atm.temp > pyrat.lt.tmax)[0][0]
    pt.error("The layer {:d} in the atmospheric model has a higher "
             "temperature ({:.1f} K) than the highest allowed TLI temperature "
             "({:.1f} K).".format(ihot, pyrat.atm.temp[ihot], pyrat.lt.tmax))

  # Interpolate isotopes partition function:
  pt.msg(pyrat.verb,"Number of isotopes: {:d}".format(pyrat.iso.niso), 2)
  # Initialize the partition-function array for pyrat.iso:
  pyrat.iso.z = np.zeros((pyrat.iso.niso, pyrat.atm.nlayers))
  for i in np.arange(pyrat.lt.ndb):           # For each Database
    for j in np.arange(pyrat.lt.db[i].niso):  # For each isotope in DB
      pt.msg(pyrat.verb, "Interpolating (isotope ID {:2d}) partition "
                         "function.".format(pyrat.lt.db[i].iiso+j), 4)
      zinterp = sip.interp1d(pyrat.lt.db[i].temp, pyrat.lt.db[i].z[j],
                             kind='slinear')
      pyrat.iso.z[pyrat.lt.db[i].iiso+j] = zinterp(pyrat.atm.temp)

  # # Plot interpolation:
  # plt.figure(3)
  # plt.clf()
  # plt.plot(pyrat.lt.db[0].temp, pyrat.lt.db[0].z[0], "o-r", mec="r")
  # plt.plot(pyrat.atm.temp, pyrat.iso.z[0], "o-b", mec="b", mfc='None')
  # plt.xlabel("Temperature  (K)")
  # plt.ylabel("Partition function")
  # plt.xlim(0, 3000)
  # plt.ylim(0, 15000)
  # plt.savefig("PartitionFunction.png")

  pt.msg(pyrat.verb, "Done.")


def reloadatm(pyrat, temp, abund):
  """
  Parameters:
  -----------
  pyrat: A Pyrat instance

  temp: 1D float ndarray
     Array with a temperature profile in Kelvin (from top to bottom layer).
  abund: 2D float ndarray
     Array with the species mole mixing ratio profiles [nlayers, nmol].

  Notes:
  ------
  This code assumes that the input temperature and abundances correspond
  to the final sampling of the atmospheric layers (after makeradius).
  """
  # Check that the dimensions match:
  if np.size(temp) != np.size(pyrat.atm.temp):
    pt.error("The temperature array size ({:d}) doesn't match the Pyrat's "
             "temperature size ({:d}).".format(np.size(temp),
                                               np.size(pyrat.atm.temp)))
  if np.shape(abund) != np.shape(pyrat.atm.q):
    pt.error("The shape of the abundances array {:s} doesn't match the "
             "shape of the Pyrat's abundance size {:s}".format(
              str(np.shape(abund)), str(np.shape(pyrat.atm.q))))

  # Put temperature and abundance data into the Pyrat object:
  pyrat.atm.temp = temp
  pyrat.atm.q    = abund

  # Mean molecular mass:
  pyrat.atm.mm = np.sum(pyrat.atm.q*pyrat.mol.mass, axis=1)

  # Density:
  for i in np.arange(pyrat.mol.nmol):
    pyrat.atm.d[:,i] = ra.IGLdensity(pyrat.atm.q[:,i], pyrat.mol.mass[i],
                                     pyrat.atm.press, pyrat.atm.temp)

  # Radius:
  pyrat.atm.radius[1:] = si.cumtrapz((-pc.k * sc.N_A * pyrat.atm.temp) /
             (pyrat.atm.mm * pyrat.surfgravity), np.log(pyrat.atm.press))
  pyrat.atm.radius[0] = 0.0
  radinterp = sip.interp1d(pyrat.atmf.press, pyrat.atmf.radius, kind='slinear')
  r0 = radinterp(pyrat.pressurebase)
  pyrat.atm.radius += pyrat.radiusbase - r0

  # Partition function:
  for i in np.arange(pyrat.lt.ndb):           # For each Database
    for j in np.arange(pyrat.lt.db[i].niso):  # For each isotope in DB
      zinterp = sip.interp1d(pyrat.lt.db[i].temp, pyrat.lt.db[i].z[j],
                             kind='slinear')
      pyrat.iso.z[pyrat.lt.db[i].iiso+j] = zinterp(pyrat.atm.temp)
