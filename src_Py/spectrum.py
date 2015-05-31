import numpy as np

import pconstants as pc
import ptools     as pt
import simpson    as s
import blackbody  as bb
import cutils     as cu

def spectrum(pyrat):
  """
  Spectrum calculation driver.
  """
  pt.msg(1, "\nCalculate the planetary spectrum.")

  # Initialize the spectrum array:
  pyrat.spectrum = np.empty(pyrat.nspec, np.double)

  # Call respective function depending on the geometry:
  if   pyrat.path == "transit":
    modulation(pyrat)

  elif pyrat.path == "eclipse":
    intensity(pyrat)
    flux(pyrat)

  # Print spectrum to file:
  printspec(pyrat)
  pt.msg(1, "Done.")


def modulation(pyrat):
  """
  Calculate the modulation spectrum for transit geometry.
  """
  pt.msg(1, "Modulation spectrum.", 2)
  # Get the stellar radius:
  rstar = pyrat.rstar
  h = np.ediff1d(pyrat.atm.radius)

  for i in np.arange(pyrat.nspec):
    # Layer index where the optical depth reached maxdepth:
    last = pyrat.od.ideep[i]
    # The integrand:
    integ = np.exp(-pyrat.od.depth[0:last+1,i]) * pyrat.atm.radius[0:last+1]

    # Integrate with Simpson's rule:
    pyrat.spectrum[i] = s.simps(integ, h[0:last], *s.geth(h[0:last]))

  pyrat.spectrum = (pyrat.atm.radius[0]**2 + 2*pyrat.spectrum) / pyrat.rstar**2


def intensity(pyrat):
  """
  Calculate the intensity spectrum [units] for eclipse geometry.
  """
  pt.msg(1, "Intensity spectrum.", 2)
  # Allocate intensity array:
  pyrat.nangles = len(pyrat.raygrid)
  pyrat.intensity = np.empty((pyrat.nangles, pyrat.nspec), np.double)

  # Calculate the Blackbody function:
  pyrat.B = np.empty((pyrat.nspec, pyrat.atm.nlayers), np.double)
  bb.planck(pyrat.B, pyrat.wn, pyrat.atm.temp, pyrat.od.ideep)

  # Allocate dtau:
  dtau = np.empty(pyrat.atm.nlayers, np.double)

  # Calculate the intensity for each angle in raygrid:
  j = 0
  while (j <pyrat.nangles):
    i = 0
    while (i < pyrat.nspec):
      # Layer index where the optical depth reached maxdepth:
      last = pyrat.od.ideep[i]
      # Optical depth:
      tau  = pyrat.od.depth[:last+1,i]
      cu.ediff(tau, dtau, last+1)
      # The integrand:
      integ = pyrat.B[i,:last+1] * np.exp(-tau/np.cos(pyrat.raygrid[j]))
      # Integrate 
      pyrat.intensity[j,i] = (s.simps(integ, dtau, *s.geth(dtau[:last])) /
                              np.cos(pyrat.raygrid[j])          )
      i += 1
    j += 1

def flux(pyrat):
  """
  Calculate the hemisphere-integrated flux spectrum [units] for eclipse
  geometry.
  """
  pt.msg(1, "Flux spectrum: '{:s}'".format(pyrat.outspec), 2)
  # Calculate the projected area:
  boundaries = np.linspace(0, 0.5*np.pi, pyrat.nangles+1)
  boundaries[1:pyrat.nangles] = 0.5 * (pyrat.raygrid[:-1] + pyrat.raygrid[1:])
  area = np.pi * (np.sin(boundaries[1:])**2 - np.sin(boundaries[:-1])**2)

  # Weight-sum the intensities to get the flux:
  pyrat.spectrum[:] = np.sum(pyrat.intensity*np.expand_dims(area,1), axis=0)


def printspec(pyrat):
  """
  Print the planetary spectrum to file.
  """

  # Type of spectrum and units:
  if   pyrat.path == "transit":
    spectype  = "Modulation"
    specunits = "[unitless]"
  elif pyrat.path == "eclipse":
    spectype  = "Flux"
    specunits = "[erg/s/cm]"

  # Wavelength units in brackets:
  wlunits = "[{:s}]".format(pyrat.wlunits)

  # Open-write file:
  specfile = open(pyrat.outspec, "w")
  # Write header:
  specfile.write("# Wavenumber   {:s}\n# {:>10s}   {:s}\n".
                                        format(spectype, wlunits, specunits))

  # Write the spectrum values:
  for i in np.arange(pyrat.nspec):
    specfile.write("  {:>10.5f}   {:>.8e}\n".
           format(1.0/pyrat.wn[i]/pc.units[pyrat.wlunits], pyrat.spectrum[i]))

  specfile.close()

