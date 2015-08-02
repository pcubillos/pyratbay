import numpy as np
import scipy.constants   as sc
import scipy.interpolate as sip

import pconstants as pc
import ptools as pt

def read(pyrat):
  """
  Read a Collision-Induced Absorption (CIA) file.
  """

  pt.msg(pyrat.verb, "\nReading the Collision-Induced-Absorption files:", 0)
  # Number of CIA files:
  pyrat.cia.nfiles = len(pyrat.cia.files)
  # Allocate molecules array:
  pyrat.cia.molecules = np.zeros((pyrat.cia.nfiles, 2), pc.strfmt)
  # Allocate the number of temperature and wavenumber samples per file:
  pyrat.cia.ntemp     = np.zeros(pyrat.cia.nfiles, np.int)
  pyrat.cia.nwave     = np.zeros(pyrat.cia.nfiles, np.int)


  if pyrat.cia.nfiles == 0:
    pt.msg(pyrat.verb, "No CIA files to read.", 2)
  else:
    for i in np.arange(pyrat.cia.nfiles):
      pt.msg(pyrat.verb, "Reading CIA file: '{:s}'.".
                          format(pyrat.cia.files[i]), 2)
      f = open(pyrat.cia.files[i], "r")
      lines = f.readlines()
      f.seek(0)

      nheader = 0  # Number of lines in header
      # Read the file, extract the keywords info:
      while True:
        line = f.readline().strip()
        nheader += 1

        if line == "@DATA":
          break

        # Skip empty and comment lines:
        elif line == '' or line.startswith('#'):
          pass

        # Get the molecules involved:
        elif line == "@SPECIES":
          pyrat.cia.molecules[i] = f.readline().strip().split()
          # Check that CIA species are in the atmospheric file:
          absent = np.setdiff1d(pyrat.cia.molecules[i], pyrat.mol.name)
          if len(absent) > 0:
            pt.error("These species: {:s} are not listed in the atmospheric "
                     "file.\n".format(str(absent)))
          nheader += 1

        # Get the sampled temperatures:
        elif line == "@TEMPERATURES":
          pyrat.cia.temp.append(np.asarray(f.readline().strip().split(),
                                           np.double))
          pyrat.cia.ntemp[i] = len(pyrat.cia.temp[i])
          # Update temperature boundaries:
          pyrat.cia.tmin = np.amax((pyrat.cia.tmin, pyrat.cia.temp[i][ 0]))
          pyrat.cia.tmax = np.amin((pyrat.cia.tmax, pyrat.cia.temp[i][-1]))
          nheader += 1

        else:
          pt.error("CIA file has an unexpected line: \n'{:s}'".format(line))

      # Read the data:
      # Get number of wavenumber samples:
      pyrat.cia.nwave[i] = len(lines) - nheader
      # Allocate the wavenumber and absorption arrays:
      wavenumber = np.zeros(pyrat.cia.nwave[i], np.double)
      # Allocate the absorption (in cm-1 amagat-2):
      absorption = np.zeros((pyrat.cia.ntemp[i], pyrat.cia.nwave[i]), np.double)
      for j in np.arange(pyrat.cia.nwave[i]):
        data = f.readline().split()
        wavenumber[j]   = data[0]
        absorption[:,j] = data[1:]

      # Wavenumber range check:
      if (wavenumber[ 0] > pyrat.spec.wn[ 0] or
          wavenumber[-1] < pyrat.spec.wn[-1] ):
        pt.error("The wavenumber range [{:.3f}, {:.3f}] cm-1 of the CIA "
                 "file:\n  '{:s}',\ndoes not cover the Pyrat's wavenumber "
                 "range: [{:.3f}, {:.3f}] cm-1.".
                 format(wavenumber[0], wavenumber[-1], pyrat.cia.files[i],
                        pyrat.spec.wn[ 0], pyrat.spec.wn[-1]))

      # Add arrays to pyrat:
      pyrat.cia.wavenumber.append(wavenumber)
      pyrat.cia.absorption.append(absorption)

      # Screen output:
      pt.msg(pyrat.verb, "For {:s}-{:s} CIA, \nRead array of {:d} wavenumber "
           "and {:d} temperature samples.".format(pyrat.cia.molecules[i,0],
           pyrat.cia.molecules[i,1], pyrat.cia.nwave[i], pyrat.cia.ntemp[i]), 4)
      pt.msg(pyrat.verb, "Temperature sample limits: {:g}--{:g} K".format(
                         pyrat.cia.temp[i][0], pyrat.cia.temp[i][-1]), 4)
      pt.msg(pyrat.verb, "Wavenumber sample limits: {:.1f}--{:.1f} cm-1".format(
                     pyrat.cia.wavenumber[i][0], pyrat.cia.wavenumber[i][-1]),4)
      f.close()

  pt.msg(pyrat.verb, "Done.", 0)


def interpolate(pyrat):
  """
  Interpolate the CIA absorption to the planetary model temperature and
  wavenumber samples:

  Notes:
  ------

  Modification History:
  ---------------------
  2014-08-31  patricio  Initial Python implementation.
  """
  pt.msg(pyrat.verb, "\nBegin CIA interpolation.")

  # Check temperature boundaries:
  if np.any(pyrat.atm.temp < pyrat.cia.tmin):
    icold = np.where(pyrat.atm.temp < pyrat.cia.tmin)[0][0]
    pt.error("The layer {:d} in the atmospheric model has a lower "
             "temperature ({:.1f} K) than the lowest allowed CIA temperature "
             "({:.1f} K).".format(icold, pyrat.atm.temp[icold], pyrat.cia.tmin))
  if np.any(pyrat.atm.temp > pyrat.cia.tmax):
    ihot = np.where(pyrat.atm.temp > pyrat.cia.tmax)[0][0]
    pt.error("The layer {:d} in the atmospheric model has a higher "
             "temperature ({:.1f} K) than the highest allowed CIA temperature "
             "({:.1f} K).".format(ihot, pyrat.atm.temp[ihot], pyrat.cia.tmax))

  # Allocate output extinction-coefficient array:
  pyrat.cia.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  for i in np.arange(pyrat.cia.nfiles):
    # Get index from the pyrat list of molecules:
    imol1 = np.where(pyrat.mol.name == pyrat.cia.molecules[i,0])[0][0]
    imol2 = np.where(pyrat.mol.name == pyrat.cia.molecules[i,1])[0][0]

    # Evaluate the spline:
    biv = sip.RectBivariateSpline(pyrat.cia.temp[i],
                                  pyrat.cia.wavenumber[i],
                                  pyrat.cia.absorption[i])

    # Interpolator requires sorted arrays:
    asort  = np.argsort(pyrat.atm.temp)
    # This one will de-sort to the original order:
    desort = np.argsort(asort)

    # Interpolate:
    sort_temp = pyrat.atm.temp[asort]
    cia_absorption = biv(sort_temp, pyrat.spec.wn)
    # Reverse sorting to the original order of the atmospheric layers:
    cia_absorption = cia_absorption[desort]

    # Densities in amagat:
    dens1 = pyrat.atm.d[:,imol1]/(pyrat.mol.mass[imol1]*pc.amu) / pc.amagat
    dens2 = pyrat.atm.d[:,imol2]/(pyrat.mol.mass[imol2]*pc.amu) / pc.amagat

    # Compute CIA absorption in cm-1 units (broadcasting):
    pyrat.cia.ec += (cia_absorption * np.expand_dims(dens1*dens2, axis=1))

    pt.msg(pyrat.verb-40, "CIA extinction: {}".format(pyrat.cia.ec[:,0]), 2)
  pt.msg(pyrat.verb, "Done.")
