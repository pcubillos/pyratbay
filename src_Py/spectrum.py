import numpy as np

import pconstants as pc
import ptools     as pt
import simpson    as s

def spectrum(pyrat):
  """
  Spectrum calculation driver.
  """
  pt.msg(1, "\nCalculate the planetary spectrum.")

  # Initialize the spectrum array:
  pyrat.spectrum = np.zeros(pyrat.nspec, np.double)

  # Call respective function depending on the geometry:
  if   pyrat.path == "transit":
    modulation(pyrat)

  elif pyrat.path == "eclipse":
    pass

  # Print spectrum to file:
  printspec(pyrat)

  pt.msg(1, "Done.")

def modulation(pyrat):
  """
  Calculate the modulation spectrum for transit geometry.
  """
  # Get the stellar radius:
  rstar = pyrat.rstar
  h = np.ediff1d(pyrat.atm.radius)
  print(h)

  for i in np.arange(pyrat.nspec):
    # Layer index where the optical depth reached maxdepth:
    last = pyrat.od.ideep[i]
    # The integrand:
    integ = np.exp(-pyrat.od.depth[0:last+1,i]) * pyrat.atm.radius[0:last+1]
    # Integrate with Simpson's rule:
    pyrat.spectrum[i] = s.simps(integ, h[0:last], *s.geth(h[0:last]))

  print(pyrat.atm.radius[0])
  pyrat.spectrum = (pyrat.atm.radius[0]**2 - 2*pyrat.spectrum) / pyrat.rstar**2


def intensity():
  """
  Calculate the intensity spectrum [units] for eclipse geometry.
  """
  pass


def flux():
  """
  Calculate the hemisphere-integrated flux spectrum [units] for eclipse
  geometry.
  """
  pass


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
  specfile.write("# Wavenumber   {:s}\n"
                 "# {:>10s}   {:s}\n".format(spectype, wlunits, specunits))

  # Write the spectrum values:
  for i in np.arange(pyrat.nspec):
    specfile.write("  {:>10.5f}   {:>.8e}\n".format(
                            1.0/pyrat.wn[i]/pc.units[pyrat.wlunits],
                                                 pyrat.spectrum[i]))

  specfile.close()

