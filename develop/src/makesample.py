import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.integrate as si
import scipy.interpolate as sip

import pconstants as pc
import ptools as pt


def makewavenumber(pyrat):
  """
  Make the wavenumber sample from user inputs.
  Store all values in CGS units.

  Modification History:
  ---------------------
  2014-04-28  patricio  Initial python implementation.
  """
  pt.msg(pyrat.verb, "\nGenerating wavenumber array:", 0)
  # User-inputs object:
  inputs = pyrat.inputs

  # Accept wavenumber and wavelength units:
  pyrat.wnunits = inputs.wnunits
  pyrat.wlunits = inputs.wlunits

  # Initial wavenumber limit:
  if inputs.wnlow is not None:
    if inputs.wnlow < 0.0:
      pt.error(message="low wavenumber boundary (%.2e %s-1) must be >= 0.0"%(
              inputs.wnlow, inputs.wnunits))
    pyrat.wnlow = inputs.wnlow / pc.units[pyrat.wnunits]
  else:
    if inputs.wlhigh < 0.0:
      pt.error(message="high wavelength boundary (%.2e %s) must be >= 0.0"%(
                inputs.wlhigh, inputs.wlunits))
    pyrat.wnlow = 1.0 / (inputs.wlhigh * pc.units[pyrat.wlunits])

  # Final wavenumber limit:
  if inputs.wnhigh is not None:
    if inputs.wnhigh < 0.0:
      pt.error(message="high wavenumber boundary (%.2e %s-1) must be >= 0.0"%(
              inputs.wnhigh, inputs.wnunits))
    pyrat.wnhigh = inputs.wnhigh / pc.units[pyrat.wnunits]
  else:
    if inputs.wllow < 0.0:
      pt.error(message="low wavelength boundary (%.2e %s) must be >= 0.0"%(
                inputs.wllow, inputs.wlunits))
    pyrat.wnhigh = 1.0 / (inputs.wllow * pc.units[pyrat.wlunits])

  # Consistency check (wnlow < wnhigh):
  if pyrat.wnlow > pyrat.wnhigh:
    pt.error(message="Wavenumber low boundary (%.2e cm-1) must be larger "
             "than the high boundary (%.2e cm-1)"%(pyrat.wnlow, pyrat.wnhigh))

  # Set wavelength limits:
  pyrat.wlhigh = 1.0/pyrat.wnlow
  pyrat.wllow  = 1.0/pyrat.wnhigh

  # Set wavenumber step:
  if inputs.wnstep is not None:
    if inputs.wnstep <= 0.0:
      pt.error(message="Wavenumber step (%.2e %s-1) must be > 0.0"%
               (inputs.wnstep, pyrat.wnunits))
    pyrat.wnstep = inputs.wnstep / pc.units[pyrat.wnunits]
  else:
    if inputs.wlstep <= 0.0:
      pt.error(message="Wavelength step (%.2e %s-1) must be > 0.0"%
               (inputs.wlstep, pyrat.wlunits))
    wlstep = inputs.wlstep * pc.units[pyrat.wlunits]
    pyrat.wnsize = int((1/pyrat.wnlow - 1/pyrat.wnhigh)/wlstep)
    pyrat.wnstep = (pyrat.wnhigh - pyrat.wnlow)/pyrat.wnsize

  # Make the wavenumber array:
  pyrat.wn = np.arange(pyrat.wnlow, pyrat.wnhigh, pyrat.wnstep)  

  # Re-set final boundary (stay inside given boundaries):
  if pyrat.wn[-1] != pyrat.wnhigh:
    pt.warning("Final wavenumber boundary modified from %.4f cm^-1 (user)\n"
               "                                     to %.4f cm^-1 (Pyrat)."%
               (pyrat.wnhigh, pyrat.wn[-1]))
  #pyrat.wnhigh = pyrat.wn[-1]
  pyrat.nspec  = len(pyrat.wn)

  # Screen output:
  pt.msg(pyrat.verb-10, "Wavelength units: %s"%pyrat.wlunits, 2)
  pt.msg(pyrat.verb, "Initial wavenumber boundary:  %.5e cm-1  (%.3e %s)"%
     (pyrat.wnlow, pyrat.wlhigh/pc.units[pyrat.wlunits], pyrat.wlunits), 2)
  pt.msg(pyrat.verb, "Final   wavenumber boundary:  %.5e cm-1  (%.3e %s)"%
     (pyrat.wnhigh, pyrat.wllow/pc.units[pyrat.wlunits], pyrat.wlunits), 2)
  pt.msg(pyrat.verb, "Wavenumber sampling stepsize: %.2e cm-1"%pyrat.wnstep, 2)
  pt.msg(pyrat.verb, "Wavenumber sample size:  %d"%pyrat.nspec, 2)
  pt.msg(pyrat.verb, "Done.", 0)


def makeradius(pyrat):
  """
  Generate atmospheric radius layers sampling.

  Notes:
  ------
  If there are radius command-line arguments set, make a sampling; else,
  use atmfile radius array.

  Modification History:
  ---------------------
  2014-06-29  patricio  Initial implementation.
  2014-08-10  patricio  Added warning print for modified boundary.
                        Implemented the isotopic interpolation.
  2014-08-31  patricio  Added molecular abundances and density interpolation.
  """
  pt.msg(pyrat.verb, "\nGenerating atmospheric radius sample:")

  # User-inputs object:
  inputs = pyrat.inputs
  # Set user-defined units:
  pyrat.radunits = inputs.radunits
  pyrat.punits   = inputs.punits

  # Store atmopsheric base level arguments:
  if inputs.radiusbase <= 0:
    pt.error("Planetary radius base must be > 0.")
  pyrat.radiusbase = inputs.radiusbase * pc.units[pyrat.radunits]

  if inputs.pressurebase <= 0:
    pt.error("Planetary pressure base must be > 0.")
  pyrat.pressurebase = inputs.pressurebase * pc.units[pyrat.punits]

  if inputs.surfgravity <= 0:
    pt.error("Planetary surface gravity must be > 0.")
  pyrat.surfgravity = inputs.surfgravity

  print("Base pressure: {:.3e} {:s} ".format(
                     pyrat.pressurebase/pc.units[pyrat.punits], pyrat.punits))
  print("Base radius: {:8g} {:s} ".format(
                     pyrat.radiusbase/pc.units[pyrat.radunits], pyrat.radunits))
  # FINDME: move this to readatm
  # Pressure limits from the atmospheric file:
  print("Pressure limits: {:.3e} -- {:.3e} {:s}".format(
        pyrat.atmf.press[ 0]/pc.units[pyrat.punits],
        pyrat.atmf.press[-1]/pc.units[pyrat.punits], pyrat.punits))

  # Calculate radius array for given atmospheric profile by using the 
  # hydostatic-equilibrium equation:
  rad = si.cumtrapz((-pc.k*sc.N_A*pyrat.atmf.temp) /
          (pyrat.atmf.mm * pyrat.surfgravity), np.log(pyrat.atmf.press))
  rad = np.concatenate(([0.0], rad))
 
  # Get radius at the baseline pressure:
  radinterp = sip.interp1d(pyrat.atmf.press[::-1], rad[::-1], kind='cubic')
  r0 = radinterp(pyrat.pressurebase)
  # Make correction to have: rad(basepress) = baserad
  rad += pyrat.radiusbase - r0

  # Reset the interpolating function (for use later):
  radinterp   = sip.interp1d(pyrat.atmf.press[::-1], rad[::-1], kind='slinear')
  pressinterp = sip.interp1d(rad, pyrat.atmf.press,             kind='cubic')
  pt.msg(pyrat.verb, "Radius array (km) = \n{:s}".format(str(rad/1e5)), 2)


  # Set values for modeling atm object:
  pyrat.atm.nmol = pyrat.mol.nmol

  # Set pressure higher boundary:
  if inputs.phigh is not None:     # From user-defined pressure
    if inputs.phigh < 0.0:
      pt.error("High atm pressure boundary ({:.2e} {:s}) must be "
               ">= 0.0".format(inputs.phigh, inputs.punits))
    pyrat.phigh = inputs.phigh * pc.units[pyrat.punits]
  elif inputs.radlow is not None:  # From user-defined radius
    if inputs.radlow < 0.0:
      pt.error("Low atm radius boundary ({:.2e} {:s}) must be "
               ">= 0.0".format(inputs.radlow, pyrat.radunits))
    pyrat.phigh = pressinterp(inputs.radlow*pc.units[pyrat.radunits])[0]
  else:                                # From atm-file
    pyrat.phigh = np.amax(pyrat.atmf.press)
  # Out of bounds errors:
  if pyrat.phigh > np.amax(pyrat.atmf.press):
    pt.error("User defined top layer (p={:.3e} {:s}) is higher than the "
             "atmospheric-file top layer (p={:.3e} {:s}).".format(
              pyrat.phigh/pc.units[pyrat.punits], pyrat.punits,
              np.amax(pyrat.atmf.press)/pc.units[pyrat.punits], pyrat.punits))

  # Set pressure lower boundary:
  if inputs.plow is not None:       # From user-defined pressure
    if inputs.plow < 0.0:
      pt.error("Low atm pressure boundary ({:.2e} {:s}) must be >= 0.0".format(
               inputs.plow, inputs.punits))
    pyrat.plow = inputs.plow * pc.units[pyrat.punits]
  elif inputs.radhigh is not None:  # From user-defined radius
    if inputs.radhigh < 0.0:
      pt.error("High atm radius boundary ({:.2e} {:s}) must be >= 0.0".format(
               inputs.radhigh, pyrat.radunits))
    pyrat.plow = pressinterp([inputs.radhigh*pc.units[pyrat.radunits]])[0]
  else:                             # From atm-file
    pyrat.plow = np.amin(pyrat.atmf.press)
  # Out of bounds errors:
  if pyrat.plow < np.amin(pyrat.atmf.press):
    pt.error("User defined bottom layer (p={:.3e} {:s}) is lower than the "
             "atmospheric-file bottom layer (p={:.3e} {:s}).".format(
              pyrat.plow/pc.units[pyrat.punits], pyrat.punits,
              np.amin(pyrat.atmf.press)/pc.units[pyrat.punits], pyrat.punits))

  # Atmospheric-layers resampling:
  resample = True
  # Resample to equispaced log-pressure array:
  if inputs.nlayers != -1:
    pyrat.atm.layers = inputs.nlayers
    pyrat.atm.press = np.logspace(np.log10(pyrat.phigh),
                                  np.log10(pyrat.plow ), pyrat.atm.layers)
    pyrat.atm.radius = radinterp(pyrat.atm.press)

  # Resample to equispaced radius:
  elif inputs.radstep is not None:
    pyrat.atm.layers = int((np.amax(rad)-np.amin(rad))/
                           (inputs.radstep*pc.units[pyrat.radunits])) + 1
    # Avoid radinterp going out-of-bounds:
    radlow  = np.amax([rad[ 0], radinterp([pyrat.phigh])[0]])
    radhigh = np.amin([rad[-1], radinterp([pyrat.plow ])[0]])
    # Make equispaced radius array:
    pyrat.atm.radius = np.arange(radlow, radhigh,
                                 inputs.radstep*pc.units[pyrat.radunits])
    # Interpolate to pressure array:
    pyrat.atm.press = pressinterp(pyrat.atm.radius)
  # Take atm file sampling:
  else:
    pyrat.atm.layers = pyrat.atmf.layers
    pyrat.atm.press  = pyrat.atmf.press
    pyrat.atm.radius = rad
    resample = False

  # Radius-vs.-pressure from Atm. file and resampled array:
  plt.figure(2)
  plt.clf()
  #plt.semilogx(pyrat.atmf.press, pyrat.atmf.rad, "o-g")
  plt.semilogx(pyrat.atmf.press/pc.units[pyrat.punits],
               rad/pc.units[pyrat.radunits], "o-r", mec="r", mfc='r')
  plt.semilogx(pyrat.atm.press/pc.units[pyrat.punits],
               pyrat.atm.radius/pc.units[pyrat.radunits], "o-b", mec="b", mew=1,
               mfc='None')
  plt.xlabel("Pressure  (%s)"%pyrat.punits)
  plt.ylabel("Radius  (%s)"%pyrat.radunits)
  plt.savefig("radpress.png")

  pt.msg(pyrat.verb, "Number of model layers: {:d}".format(pyrat.atm.layers), 2)
  pt.msg(pyrat.verb, "Pressure lower/higher boundaries: {:.2e} - {:.2e} "
                     "{:s}".format(pyrat.plow /pc.units[pyrat.punits],
                           pyrat.phigh/pc.units[pyrat.punits], pyrat.punits), 2)
  pt.msg(pyrat.verb, "Radius lower/higher boundaries:   {:.2e} - {:.2e} "
         "{:s}".format(np.amin(pyrat.atm.radius)/pc.units[pyrat.radunits],
         np.amax(pyrat.atm.radius)/pc.units[pyrat.radunits], pyrat.radunits), 2)

  # If defaulting to atm file values, don't interpolate:
  if (pyrat.plow  == np.amin(pyrat.atmf.press) and 
      pyrat.phigh == np.amax(pyrat.atmf.press) and not resample):
    pyrat.atm.temp   = pyrat.atmf.temp
    pyrat.atm.mm     = pyrat.atmf.mm
    pyrat.atm.q      = pyrat.atmf.q
    pyrat.atm.d      = pyrat.atmf.d

  else:  # Interpolate to new atm-layer sampling:
    tempinterp = sip.interp1d(pyrat.atmf.press[::-1], pyrat.atmf.temp[::-1],
                              kind='cubic')
    mminterp   = sip.interp1d(pyrat.atmf.press[::-1], pyrat.atmf.mm[::-1],
                              kind='cubic')
    pyrat.atm.temp = tempinterp(pyrat.atm.press)
    pyrat.atm.m    =   mminterp(pyrat.atm.press)
    # Interpolate abundance profiles:
    pyrat.atm.q = np.zeros((pyrat.atm.layers, pyrat.atm.nmol))
    pyrat.atm.d = np.zeros((pyrat.atm.layers, pyrat.atm.nmol))
    for i in np.arange(pyrat.atmf.nmol):
      qinterp = sip.interp1d(pyrat.atmf.press[::-1], pyrat.atmf.q[::-1, i],
                             kind='cubic')
      dinterp = sip.interp1d(pyrat.atmf.press[::-1], pyrat.atmf.d[::-1, i],
                             kind='cubic')
      pyrat.atm.q[:,i] = qinterp(pyrat.atm.press)
      pyrat.atm.d[:,i] = dinterp(pyrat.atm.press)

  # Interpolate isotopes partition function:
  pt.msg(pyrat.verb,"Number of isotopes: {:d}".format(pyrat.iso.niso), 2)
  
  # Initialize the partition-function array for pyrat.iso:
  pyrat.iso.z = np.zeros((pyrat.iso.niso, pyrat.atm.layers))
  for i in np.arange(pyrat.lt.ndb):
    for j in np.arange(pyrat.lt.db[i].niso):
      zinterp = sip.interp1d(pyrat.lt.db[i].temp, pyrat.lt.db[i].z[j],
                             kind='cubic')
      pt.msg(pyrat.verb, "Interpolating (isotope ID {:2d}) partition "
                         "function.".format(pyrat.lt.db[i].iiso+j), 4)
      pyrat.iso.z[pyrat.lt.db[i].iiso+j] = zinterp(pyrat.atm.temp)

  # Plot interpolation:
  plt.figure(3)
  plt.clf()
  plt.plot(pyrat.lt.db[0].temp, pyrat.lt.db[0].z[0], "o-r", mec="r")
  plt.plot(pyrat.atm.temp, pyrat.iso.z[0], "o-b", mec="b", mfc='None')
  plt.xlabel("Temperature  (K)")
  plt.ylabel("Partition function")
  plt.xlim(0, 3000)
  plt.ylim(0, 15000)
  plt.savefig("PartitionFunction.png")

  pt.msg(pyrat.verb, "Done.")

