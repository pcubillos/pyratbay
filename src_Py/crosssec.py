import numpy as np
import scipy.constants   as sc
import scipy.interpolate as sip

import pconstants as pc
import ptools as pt

def read(pyrat):
  """
  Read a Cross-section (CS) file.
  """

  pt.msg(pyrat.verb, "\nReading the cross-section files:", 0)
  # Number of CS files:
  pyrat.cs.nfiles = len(pyrat.cs.files)
  # Allocate molecules array:
  pyrat.cs.molecules = np.zeros((pyrat.cs.nfiles, 2), pc.strfmt)
  pyrat.cs.nmol      = np.zeros(pyrat.cs.nfiles, np.int)
  # Allocate the number of temperature and wavenumber samples per file:
  pyrat.cs.ntemp     = np.zeros(pyrat.cs.nfiles, np.int)
  pyrat.cs.nwave     = np.zeros(pyrat.cs.nfiles, np.int)

  if pyrat.cs.nfiles == 0:
    pt.msg(pyrat.verb, "No CS files to read.", 2)
  else:
    for i in np.arange(pyrat.cs.nfiles):
      pt.msg(pyrat.verb, "Read CS file: '{:s}'.".format(pyrat.cs.files[i]), 2)
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
          nmol = pyrat.cs.nmol = len(molecs)
          pyrat.cs.molecules[i, 0:nmol] = molecs
          # Check that CS species are in the atmospheric file:
          absent = np.setdiff1d(pyrat.cs.molecules[i, 0:nmol], pyrat.mol.name)
          if len(absent) > 0:
            pt.error("These species: {:s} are not listed in the atmospheric "
                     "file.\n".format(str(absent)))
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
          pt.error("CS file has an unexpected line: \n'{:s}'".format(line))

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

      # Wavenumber range check:
      if (wavenumber[ 0] > pyrat.spec.wn[ 0] or
          wavenumber[-1] < pyrat.spec.wn[-1] ):
        pt.error("The wavenumber range [{:.3f}, {:.3f}] cm-1 of the CS "
                 "file:\n  '{:s}',\ndoes not cover the Pyrat's wavenumber "
                 "range: [{:.3f}, {:.3f}] cm-1.".
                 format(wavenumber[0], wavenumber[-1], pyrat.cs.files[i],
                        pyrat.spec.wn[ 0], pyrat.spec.wn[-1]))

      # Add arrays to pyrat:
      pyrat.cs.wavenumber.append(wavenumber)
      pyrat.cs.absorption.append(absorption)

      # Screen output:
      pt.msg(pyrat.verb, "For {:s} CS,\nRead {:d} wavenumber and {:d} "
        "temperature samples.".format("-".join(pyrat.cs.molecules[i,0:nmol]),
                                      pyrat.cs.nwave[i], pyrat.cs.ntemp[i]), 4)
      pt.msg(pyrat.verb, "Temperature sample limits: {:g}--{:g} K".
                          format(pyrat.cs.temp[i][0], pyrat.cs.temp[i][-1]), 4)
      pt.msg(pyrat.verb, "Wavenumber sample limits: {:.1f}--{:.1f} cm-1".
              format(pyrat.cs.wavenumber[i][0], pyrat.cs.wavenumber[i][-1]), 4)
      f.close()

  pt.msg(pyrat.verb, "Done.", 0)


def interpolate(pyrat):
  """
  Interpolate the CS absorption to the planetary model temperature and
  wavenumber samples.
  """
  pt.msg(pyrat.verb, "\nBegin CS interpolation.")

  # Check temperature boundaries:
  if np.any(pyrat.atm.temp < pyrat.cs.tmin):
    icold = np.where(pyrat.atm.temp < pyrat.cs.tmin)[0][0]
    pt.error("The layer {:d} in the atmospheric model has a lower "
             "temperature ({:.1f} K) than the lowest allowed CS temperature "
             "({:.1f} K).".format(icold, pyrat.atm.temp[icold], pyrat.cs.tmin))
  if np.any(pyrat.atm.temp > pyrat.cs.tmax):
    ihot = np.where(pyrat.atm.temp > pyrat.cs.tmax)[0][0]
    pt.error("The layer {:d} in the atmospheric model has a higher "
             "temperature ({:.1f} K) than the highest allowed CS temperature "
             "({:.1f} K).".format(ihot, pyrat.atm.temp[ihot], pyrat.cs.tmax))

  # Allocate output extinction-coefficient array:
  pyrat.cs.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  for i in np.arange(pyrat.cs.nfiles):
    # Evaluate the spline:
    biv = sip.RectBivariateSpline(pyrat.cs.temp[i],
                                  pyrat.cs.wavenumber[i],
                                  pyrat.cs.absorption[i])

    # Interpolator requires sorted arrays:
    asort  = np.argsort(pyrat.atm.temp)
    # This one will de-sort to the original order:
    desort = np.argsort(asort)

    # Interpolate:
    sort_temp = pyrat.atm.temp[asort]
    cs_absorption = biv(sort_temp, pyrat.spec.wn)
    # Reverse sorting to the original order of the atmospheric layers:
    cs_absorption = cs_absorption[desort]

    # Apply density scale factor:
    dens = 1.0
    for j in np.arange(pyrat.cs.nmol):
      # Get index from the pyrat list of molecules:
      imol = np.where(pyrat.mol.name == pyrat.cs.molecules[i,j])[0][0]
      # Densities in amagat:
      dens *= pyrat.atm.d[:,imol] / pc.amagat

    # Compute CS absorption in cm-1 units (broadcasting):
    pyrat.cs.ec += (cs_absorption * np.expand_dims(dens, axis=1))

    pt.msg(pyrat.verb-40, "CS extinction: {}".format(pyrat.cs.ec[:,0]), 2)
  pt.msg(pyrat.verb, "Done.")
