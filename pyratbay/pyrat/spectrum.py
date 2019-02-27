# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys
import numpy as np

from .. import blackbody as bb
from .. import constants as pc
from .. import io        as io

sys.path.append(pc.ROOT + 'lib')
import simpson    as s
import trapz      as t


def spectrum(pyrat):
  """
  Spectrum calculation driver.
  """
  pyrat.log.msg("\nCalculate the planetary spectrum.")

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
  io.write_spectrum(1.0/pyrat.spec.wn, pyrat.spec.spectrum, pyrat.outspec,
                    pyrat.od.path)
  pyrat.log.msg("Done.")


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
  if (len(h) % 2) == 1:
        hseven, hreven, hfeven = s.geth(h[0:-1])
        hsodd,  hrodd,  hfodd  = s.geth(h[0:])
  else:
        hseven, hreven, hfeven = s.geth(h[0:])
        hsodd,  hrodd,  hfodd  = s.geth(h[0:-1])
  hsum = [hsodd, hseven]
  hrat = [hrodd, hreven]
  hfac = [hfodd, hfeven]

  nlayers = pyrat.od.ideep - rtop + 1    # Number of layers for integration
  nhalf   = np.asarray(nlayers//2, int)  # Half-size of used layers
  # Parity of nlayers (if nlayers is odd, parity is 1, h is even):
  parity  = nlayers % 2

  for i in np.arange(pyrat.spec.nwave):
      nl = nlayers[i]
      nh = nhalf  [i]
      p  = parity [i]
      # Integrate with Simpson's rule:
      pyrat.spec.spectrum[i] = s.simps(integ[0:nl,i], h[0:nl-1],
                               hsum[p][0:nh], hrat[p][0:nh], hfac[p][0:nh])
      # Extra spectrum for patchy model:
      if pyrat.haze.fpatchy is not None:
          pyrat.spec.cloudy[i] = s.simps(pinteg[0:nl,i], h[0:nl-1],
                               hsum[p][0:nh], hrat[p][0:nh], hfac[p][0:nh])

  pyrat.spec.spectrum = ((pyrat.atm.radius[rtop]**2 + 2*pyrat.spec.spectrum)
                         / pyrat.phy.rstar**2)
  if pyrat.haze.fpatchy is not None:
      pyrat.spec.cloudy = ((pyrat.atm.radius[rtop]**2 + 2*pyrat.spec.cloudy)
                           / pyrat.phy.rstar**2)
      pyrat.spec.clear = pyrat.spec.spectrum
      pyrat.spec.spectrum = (   pyrat.haze.fpatchy  * pyrat.spec.cloudy +
                             (1-pyrat.haze.fpatchy) * pyrat.spec.clear  )

  pyrat.log.msg("Computed transmission spectrum: '{:s}'.".
                format(pyrat.outspec), indent=2)


def intensity(pyrat):
  """
  Calculate the intensity spectrum [units] for eclipse geometry.
  """
  pyrat.log.msg("Computing intensity spectrum.", verb=2, indent=2)
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
  pyrat.log.msg("Computed flux spectrum: '{:s}'.".
                format(pyrat.outspec), indent=2)
