# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import numpy as np

from .. import constants as pc

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../lib')
import spline as sp


def read(pyrat):
  """
  Read a Cross-section (CS) file.
  """

  pyrat.log.msg("\nReading cross-section files.")
  # Number of CS files:
  if pyrat.cs.files is None:
    pyrat.cs.nfiles = 0
  else:
    pyrat.cs.nfiles = len(pyrat.cs.files)
    # Allocate molecules array:
    pyrat.cs.molecules = np.zeros((pyrat.cs.nfiles, 2), pc.strfmt)
    pyrat.cs.nmol      = np.zeros(pyrat.cs.nfiles, np.int)
    # Allocate the number of temperature and wavenumber samples per file:
    pyrat.cs.ntemp     = np.zeros(pyrat.cs.nfiles, np.int)
    pyrat.cs.nwave     = np.zeros(pyrat.cs.nfiles, np.int)

  if pyrat.cs.nfiles == 0:
    pyrat.log.msg("No CS files to read.", indent=2)
  else:
    for i in np.arange(pyrat.cs.nfiles):
      pyrat.log.msg("Read CS file: '{:s}'.".format(pyrat.cs.files[i]),
                    pyrat.log, indent=2)
      f = open(pyrat.cs.files[i], "r")
      lines = f.readlines()
      f.seek(0)

      nheader = 0  # Number of lines in header
      # Read the file, extract the keywords info:
      while True:
        line = f.readline().strip()
        nheader += 1

        if   line == "@DATA":
          break

        # Skip empty and comment lines:
        elif line == '' or line.startswith('#'):
          pass

        # Get the molecules involved:
        elif line == "@SPECIES":
          molecs = f.readline().strip().split()
          nmol = pyrat.cs.nmol[i] = len(molecs)
          pyrat.cs.molecules[i, 0:nmol] = molecs
          # Check that CS species are in the atmospheric file:
          absent = np.setdiff1d(pyrat.cs.molecules[i, 0:nmol], pyrat.mol.name)
          if len(absent) > 0:
            pyrat.log.error("These species: {:s} are not listed in the "
                            "atmospheric file.\n".format(str(absent)))
          nheader += 1

        # Get the sampled temperatures:
        elif line == "@TEMPERATURES":
          pyrat.cs.temp.append(np.asarray(f.readline().strip().split(),
                                          np.double))
          pyrat.cs.ntemp[i] = len(pyrat.cs.temp[i])
          # Update temperature boundaries:
          pyrat.cs.tmin = np.amax((pyrat.cs.tmin, pyrat.cs.temp[i][ 0]))
          pyrat.cs.tmax = np.amin((pyrat.cs.tmax, pyrat.cs.temp[i][-1]))
          nheader += 1

        else:
          pyrat.log.error("CS file has unexpected line: \n'{:s}'".format(line))

      # Read the data:
      # Get number of wavenumber samples:
      pyrat.cs.nwave[i] = len(lines) - nheader
      # Allocate the wavenumber and absorption arrays:
      wavenumber = np.zeros(pyrat.cs.nwave[i], np.double)
      # Allocate the absorption (in cm-1 amagat-2):
      absorption = np.zeros((pyrat.cs.ntemp[i], pyrat.cs.nwave[i]), np.double)
      for j in np.arange(pyrat.cs.nwave[i]):
        data = f.readline().split()
        wavenumber[j]   = data[0]
        absorption[:,j] = data[1:]
      f.close()

      # Wavenumber range check:
      if (wavenumber[ 0] > pyrat.spec.wn[ 0] or
          wavenumber[-1] < pyrat.spec.wn[-1] ):
        pyrat.log.warning("The wavenumber range [{:.3f}, {:.3f}] cm-1 "
            "of the CS file:\n  '{:s}',\ndoes not cover the Pyrat's "
            "wavenumber range: [{:.3f}, {:.3f}] cm-1.".format(
            wavenumber[0], wavenumber[-1], pyrat.cs.files[i],
            pyrat.spec.wn[ 0], pyrat.spec.wn[-1]))

      # Add arrays to pyrat:
      pyrat.cs.wavenumber.append(wavenumber)
      pyrat.cs.absorption.append(absorption)

      # Screen output:
      pyrat.log.msg("For {:s} CS,\nRead {:d} wavenumber and {:d} temperature "
          "samples.".format("-".join(pyrat.cs.molecules[i,0:nmol]),
          pyrat.cs.nwave[i], pyrat.cs.ntemp[i]), verb=2, indent=4)
      pyrat.log.msg("Temperature sample limits: {:g}--{:g} K".
          format(pyrat.cs.temp[i][0], pyrat.cs.temp[i][-1]),
          verb=2, indent=4)
      pyrat.log.msg("Wavenumber sample limits: {:.1f}--{:.1f} cm-1".
          format(pyrat.cs.wavenumber[i][0], pyrat.cs.wavenumber[i][-1]),
          verb=2, indent=4)

      # Wavenumber-interpolated CS:
      iabsorp = np.zeros((pyrat.cs.ntemp[i], pyrat.spec.nwave), np.double)
      for j in np.arange(pyrat.cs.ntemp[i]):
        iabsorp[j] = sp.splinterp(absorption[j], wavenumber, pyrat.spec.wn, 0.0)
      pyrat.cs.iabsorp.append(iabsorp)
      # Array with second derivatives:
      iz   = np.zeros((pyrat.spec.nwave, pyrat.cs.ntemp[i]), np.double)
      wnlo = np.flatnonzero(pyrat.spec.wn >= wavenumber[ 0])[ 0]
      wnhi = np.flatnonzero(pyrat.spec.wn <= wavenumber[-1])[-1]
      for j in np.arange(wnlo, wnhi):
        iz[j] = sp.spline_init(iabsorp[:,j], pyrat.cs.temp[i])
      pyrat.cs.iz.append(iz.T)
      pyrat.cs.iwnlo.append(wnlo)
      pyrat.cs.iwnhi.append(wnhi)

  pyrat.log.msg("Cross-section read done.")


def interpolate(pyrat, layer=None):
  """
  Interpolate the CS absorption into the planetary model temperature.
  """
  pyrat.log.msg("\nBegin CS interpolation.")

  # Check temperature boundaries:
  if np.any(pyrat.atm.temp < pyrat.cs.tmin):
    icold = np.where(pyrat.atm.temp < pyrat.cs.tmin)[0][0]
    pyrat.log.error("Layer {:d} in the atmospheric model has a lower "
             "temperature ({:.1f} K) than the lowest allowed CS temperature "
             "({:.1f} K).".format(icold, pyrat.atm.temp[icold], pyrat.cs.tmin))
  if np.any(pyrat.atm.temp > pyrat.cs.tmax):
    ihot = np.where(pyrat.atm.temp > pyrat.cs.tmax)[0][0]
    pyrat.log.error("Layer {:d} in the atmospheric model has a higher "
             "temperature ({:.1f} K) than the highest allowed CS temperature "
             "({:.1f} K).".format(ihot, pyrat.atm.temp[ihot], pyrat.cs.tmax))

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

    # Apply density scale factor:
    dens = 1.0
    for j in np.arange(pyrat.cs.nmol[i]):
      # Get index from the pyrat list of molecules:
      imol = np.where(pyrat.mol.name == pyrat.cs.molecules[i,j])[0][0]
      # Densities in amagat:
      dens *= pyrat.atm.d[li:lf,imol] / pc.amagat

    # Compute CS absorption in cm-1 units:
    if layer is None:
      ec += cs_absorption * np.expand_dims(dens, axis=1)
    else:
      ec[i] = cs_absorption * dens
      label.append("-".join(pyrat.cs.molecules[i]))

    pyrat.log.msg("CS extinction: {}".format(ec[:,0]), verb=4, indent=2)
  # Return per-database EC if single-layer run:
  if layer is not None:
    return ec, label
  # Else, store cumulative result into pyrat object:
  pyrat.cs.ec = ec
  pyrat.log.msg("Cross-section interpolate done.")
