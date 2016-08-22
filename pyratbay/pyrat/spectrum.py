# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys, os
import numpy as np

from .. import tools     as pt
from .. import constants as pc

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../lib')
import simpson    as s
import blackbody  as bb
import cutils     as cu

def spectrum(pyrat):
  """
  Spectrum calculation driver.
  """
  pt.msg(pyrat.verb-3, "\nCalculate the planetary spectrum.", pyrat.log)

  # Initialize the spectrum array:
  pyrat.spec.spectrum = np.empty(pyrat.spec.nwave, np.double)

  # Call respective function depending on the geometry:
  if   pyrat.od.path == "transit":
    modulation(pyrat)

  elif pyrat.od.path == "eclipse":
    intensity(pyrat)
    flux(pyrat)

  # Print spectrum to file:
  printspec(pyrat)
  pt.msg(pyrat.verb-3, "Done.", pyrat.log)


def modulation(pyrat):
  """
  Calculate the modulation spectrum for transit geometry.
  """
  rtop = pyrat.atm.rtop
  # Get the stellar radius:
  h = np.ediff1d(pyrat.atm.radius[rtop:])

  for i in np.arange(pyrat.spec.nwave):
    # Layer index where the optical depth reached maxdepth:
    last = pyrat.od.ideep[i]
    # The integrand:
    integ = (np.exp(-pyrat.od.depth[rtop:last+1,i]) *
                   pyrat.atm.radius[rtop:last+1])

    # Integrate with Simpson's rule:
    pyrat.spec.spectrum[i] = s.simps(integ, h[0:last], *s.geth(h[0:last]))

  pyrat.spec.spectrum = ((pyrat.atm.radius[rtop]**2 + 2*pyrat.spec.spectrum) /
                         pyrat.phy.rstar**2)
  pt.msg(pyrat.verb-3, "Computed transmission spectrum: '{:s}'.".
                        format(pyrat.outspec), pyrat.log, 2)


def intensity(pyrat):
  """
  Calculate the intensity spectrum [units] for eclipse geometry.
  """
  pt.msg(pyrat.verb-4, "Computing intensity spectrum.", pyrat.log, 2)
  if pyrat.quadrature is not None:
    pyrat.raygrid = np.arccos(np.sqrt(pyrat.qnodes))

  # Allocate intensity array:
  pyrat.nangles = len(pyrat.raygrid)
  pyrat.spec.intensity = np.empty((pyrat.nangles, pyrat.spec.nwave), np.double)

  # Calculate the Planck Emission:
  pyrat.od.B = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave), np.double)
  bb.planck(pyrat.od.B, pyrat.spec.wn, pyrat.atm.temp, pyrat.od.ideep)

  # Allocate dtau:
  dtau  = np.empty(pyrat.atm.nlayers, np.double)
  dltau = np.empty(pyrat.atm.nlayers, np.double)

  # Calculate the intensity for each angle in raygrid:
  rtop = pyrat.atm.rtop
  i = 0
  while (i < pyrat.spec.nwave):
    # Layer index where the optical depth reached maxdepth:
    last = pyrat.od.ideep[i]
    # Optical depth:
    tau  = pyrat.od.depth[rtop:last+1,i]
    cu.ediff(tau, dtau, last+1)
    hsum, hratio, hfactor = s.geth(dtau[rtop:last])
    j = 0
    while (j < pyrat.nangles):
      # The integrand:
      integ = (pyrat.od.B[rtop:last+1,i] *
               np.exp(-tau/np.cos(pyrat.raygrid[j])) / np.cos(pyrat.raygrid[j]))
      # Simpson integration:
      pyrat.spec.intensity[j,i] = s.simps(integ, dtau, hsum, hratio, hfactor)
      #ltau = np.log(tau[1:])
      #integ=(pyrat.od.B[i,1:last+1]*np.exp(-tau[1:]/np.cos(pyrat.raygrid[j]))/
      #         tau[1:]/np.cos(pyrat.raygrid[j]))
      #pyrat.spec.intensity[j,i] = si.simps(integ, ltau)
      j += 1
    i += 1


def flux(pyrat):
  """
  Calculate the hemisphere-integrated flux spectrum [units] for eclipse
  geometry.
  """
  # Calculate the projected area:
  boundaries = np.linspace(0, 0.5*np.pi, pyrat.nangles+1)
  boundaries[1:pyrat.nangles] = 0.5 * (pyrat.raygrid[:-1] + pyrat.raygrid[1:])
  area = np.pi * (np.sin(boundaries[1:])**2 - np.sin(boundaries[:-1])**2)

  if pyrat.quadrature is not None:
    area = pyrat.qweights * np.pi
  # Weight-sum the intensities to get the flux:
  pyrat.spec.spectrum[:] = np.sum(pyrat.spec.intensity *
                                  np.expand_dims(area,1), axis=0)
  pt.msg(pyrat.verb-3, "Computed flux spectrum: '{:s}'.".
                        format(pyrat.outspec), pyrat.log, 2)


def printspec(pyrat):
  """
  Print the planetary spectrum to file.
  """

  # Type of spectrum and units:
  if   pyrat.od.path == "transit":
    spectype  = "Modulation"
    specunits = "[unitless]"
  elif pyrat.od.path == "eclipse":
    spectype  = "Flux"
    specunits = "[erg/s/cm]"

  # Wavelength units in brackets:
  wlunits = "[{:s}]".format(pyrat.spec.wlunits)

  # Open-write file:
  specfile = open(pyrat.outspec, "w")
  # Write header:
  specfile.write("# Wavenumber   {:s}\n# {:>10s}   {:s}\n".
                                        format(spectype, wlunits, specunits))

  # Write the spectrum values:
  for i in np.arange(pyrat.spec.nwave):
    specfile.write("  {:>10.5f}   {:>.8e}\n".
                    format(1.0/pyrat.spec.wn[i]/pt.u(pyrat.spec.wlunits),
                           pyrat.spec.spectrum[i]))

  specfile.close()

