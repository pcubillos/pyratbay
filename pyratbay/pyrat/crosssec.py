# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import numpy as np

from .. import constants as pc
from .. import io as io
from ..lib import _spline as sp


def read(pyrat):
  """
  Read a Cross-section (CS) file.
  """
  pyrat.log.head('\nReading cross-section files.')
  pyrat.cs = pyrat.cs.clone_new(pyrat)

  if pyrat.cs.files is None:
      pyrat.log.head('No CS files to read.', indent=2)
      return

  pyrat.cs.nfiles = len(pyrat.cs.files)

  for csfile in pyrat.cs.files:
      pyrat.log.head("Read CS file: '{:s}'.".format(csfile), indent=2)
      absorption, molecules, temp, wavenumber = io.read_cs(csfile)

      pyrat.cs.absorption.append(absorption)
      pyrat.cs.molecules.append(molecules)
      pyrat.cs.temp.append(temp)
      pyrat.cs.wavenumber.append(wavenumber)

      ntemp = len(temp)
      nwave = len(wavenumber)

      # Check that CS species are in the atmospheric file:
      absent = np.setdiff1d(molecules, pyrat.mol.name)
      if len(absent) > 0:
          pyrat.log.error('These cross-section species {} are not listed in '
                          'the atmospheric file.\n'.format(str(absent)))

      # Update temperature boundaries:
      pyrat.cs.tmin = np.amax((pyrat.cs.tmin, temp[ 0]))
      pyrat.cs.tmax = np.amin((pyrat.cs.tmax, temp[-1]))

      # Wavenumber range check:
      if (wavenumber[ 0] > pyrat.spec.wn[ 0] or
          wavenumber[-1] < pyrat.spec.wn[-1] ):
          pyrat.log.warning("The wavenumber range [{:.3f}, {:.3f}] cm-1 "
              "of the CS file:\n  '{:s}',\ndoes not cover the Pyrat's "
              "wavenumber range: [{:.3f}, {:.3f}] cm-1.".
              format(wavenumber[0], wavenumber[-1], csfile,
                     pyrat.spec.wn[0], pyrat.spec.wn[-1]))

      # Screen output:
      pyrat.log.msg('Cross-section opacity for {:s}:\n'
          'Read {:d} wavenumber and {:d} temperature samples.'.
          format('-'.join(molecules), nwave, ntemp), indent=4)
      pyrat.log.msg('Temperature sample limits: {:g}--{:g} K'.
          format(temp[0], temp[-1]), indent=4)
      pyrat.log.msg('Wavenumber sample limits: {:.1f}--{:.1f} cm-1'.
          format(wavenumber[0], wavenumber[-1]), indent=4)

      # Wavenumber-interpolated CS:
      iabsorp = np.zeros((ntemp, pyrat.spec.nwave), np.double)
      for j in range(ntemp):
          z = sp.second_deriv(absorption[j], wavenumber)
          iabsorp[j] = sp.splinterp_1D(absorption[j], wavenumber, z,
              pyrat.spec.wn, 0.0)
      pyrat.cs.iabsorp.append(iabsorp)
      # Array with second derivatives:
      iz   = np.zeros((pyrat.spec.nwave, ntemp), np.double)
      wnlo = np.flatnonzero(pyrat.spec.wn >= wavenumber[ 0])[ 0]
      wnhi = np.flatnonzero(pyrat.spec.wn <= wavenumber[-1])[-1] + 1
      for j in range(wnlo, wnhi):
          iz[j] = sp.second_deriv(iabsorp[:,j], temp)
      pyrat.cs.iz.append(iz.T)
      pyrat.cs.iwnlo.append(wnlo)
      pyrat.cs.iwnhi.append(wnhi)

  pyrat.log.head('Cross-section read done.')


def interpolate(pyrat, layer=None):
  """
  Interpolate the CS absorption into the planetary model temperature.
  """
  pyrat.log.head('\nBegin CS interpolation.')

  # Allocate output extinction-coefficient array:
  if layer is None:   # Take a single layer
      ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))
      li, lf = 0, pyrat.atm.nlayers
  else:               # Take whole atmosphere
      ec = np.zeros((pyrat.cs.nfiles, pyrat.spec.nwave))
      li, lf = layer, layer+1
      label = []

  for i in np.arange(pyrat.cs.nfiles):
      cs_absorption = np.zeros((lf-li, pyrat.spec.nwave))
      sp.splinterp_2D(pyrat.cs.iabsorp[i], pyrat.cs.temp[i], pyrat.cs.iz[i],
                      pyrat.atm.temp[li:lf], cs_absorption,
                      pyrat.cs.iwnlo[i], pyrat.cs.iwnhi[i])

      # Get density scale factor in amagat:
      dens = 1.0
      for mol in pyrat.cs.molecules[i]:
          imol = np.where(pyrat.mol.name == mol)[0][0]
          dens *= pyrat.atm.d[li:lf,imol] / pc.amagat

      # Compute CS absorption in cm-1 units:
      if layer is None:
          ec += cs_absorption * np.expand_dims(dens, axis=1)
      else:
          ec[i] = cs_absorption * dens
          if len(pyrat.cs.molecules[i]) == 2:
              label.append('CIA ' + '-'.join(pyrat.cs.molecules[i]))
          else:
              label.append(pyrat.cs.molecules[i][0])

  # Return per-database EC if single-layer run:
  if layer is not None:
      return ec, label
  # Else, store cumulative result into pyrat object:
  pyrat.cs.ec = ec
  pyrat.log.head('Cross-section interpolate done.')
