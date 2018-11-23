# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys
import os
import numpy as np

from .. import tools     as pt
from .. import blackbody as bb
from .. import constants as pc

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../lib')
import simpson    as s
import trapz      as t
import cutils     as cu


def spectrum(pyrat):
  """
  Spectrum calculation driver.
  """
  pt.msg(pyrat.verb-3, "\nCalculate the planetary spectrum.", pyrat.log)

  # Initialize the spectrum array:
  pyrat.spec.spectrum = np.empty(pyrat.spec.nwave, np.double)
  if pyrat.haze.fpatchy is not None:
    pyrat.spec.clear  = np.empty(pyrat.spec.nwave, np.double)
    pyrat.spec.cloudy = np.empty(pyrat.spec.nwave, np.double)

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
  # The integrand:
  integ = (np.exp(-pyrat.od.depth[rtop:,:]) *
           np.expand_dims(pyrat.atm.radius[rtop:],1))
  if pyrat.haze.fpatchy is not None:
    pinteg = (np.exp(-pyrat.od.pdepth[rtop:,:]) *
              np.expand_dims(pyrat.atm.radius[rtop:],1))
  # Get Delta radius (and simps' integration variables):
  h = np.ediff1d(pyrat.atm.radius[rtop:])
  hseven, hreven, hfeven = s.geth(h[:-1])
  hsodd,  hrodd,  hfodd  = s.geth(h)
  hsum = [hseven, hsodd]  # Re-packed for even/odd values of ideep
  hrat = [hreven, hrodd]
  hfac = [hfeven, hfodd]

  nlayers = pyrat.od.ideep-rtop         # Number of layers for integration
  nhalf   = np.asarray(nlayers/2, int)  # Half-size of used layers
  par     = pyrat.od.ideep % 2          # Parity of ideep

  for i in np.arange(pyrat.spec.nwave):
    nl = nlayers[i]
    nh = nhalf[i]
    p  = par[i]
    # Integrate with Simpson's rule:
    pyrat.spec.spectrum[i] = s.simps(integ[0:nl+1,i], h[0:nl],
                             hsum[p][0:nh], hrat[p][0:nh], hfac[p][0:nh])
    # Extra spectrum for patchy model:
    if pyrat.haze.fpatchy is not None:
      pyrat.spec.cloudy[i] = s.simps(pinteg[0:nl,i], h[0:nl],
                             hsum[p][0:nh], hrat[p][0:nh], hfac[p][0:nh])

  pyrat.spec.spectrum = ((pyrat.atm.radius[rtop]**2 + 2*pyrat.spec.spectrum) /
                         pyrat.phy.rstar**2)
  if pyrat.haze.fpatchy is not None:
    pyrat.spec.cloudy = ((pyrat.atm.radius[rtop]**2 + 2*pyrat.spec.cloudy) /
                         pyrat.phy.rstar**2)
    pyrat.spec.clear = pyrat.spec.spectrum
    pyrat.spec.spectrum = (   pyrat.haze.fpatchy  * pyrat.spec.cloudy +
                           (1-pyrat.haze.fpatchy) * pyrat.spec.clear  )

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
  bb.Bwn2D(pyrat.spec.wn, pyrat.atm.temp, pyrat.od.B, pyrat.od.ideep)

  # Plane-parallel radiative-transfer intensity integration:
  pyrat.spec.intensity = t.intensity(pyrat.od.depth, pyrat.od.ideep,
                           pyrat.od.B, np.cos(pyrat.raygrid), pyrat.atm.rtop)


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
  if pyrat.outspec is None:
    return

  # Type of spectrum and units:
  if   pyrat.od.path == "transit":
    spectype  = "Modulation"
    specunits = "[unitless]"
  elif pyrat.od.path == "eclipse":
    spectype  = "Flux"
    specunits = "[erg s-1 cm-2 cm]"

  # Wavelength units in brackets:
  wlunits = "[{:s}]".format(pyrat.spec.wlunits)
  wlength = 1.0/pyrat.spec.wn/pt.u(pyrat.spec.wlunits)
  # Precision of 5 decimal places (or better if needed):
  precision = -np.floor(np.amin(np.log10(np.abs(np.ediff1d(wlength)))))
  precision = int(np.clip(precision+1, 5, np.inf))

  # Open-write file:
  specfile = open(pyrat.outspec, "w")
  # Write header:
  specfile.write("# Wavenumber   {:s}\n# {:>10s}   {:s}\n".
                                        format(spectype, wlunits, specunits))

  # Write the spectrum values:
  for i in np.arange(pyrat.spec.nwave):
    specfile.write("  {:>15.{:d}f}   {:>.8e}\n".
                    format(wlength[i], precision, pyrat.spec.spectrum[i]))
  specfile.close()
