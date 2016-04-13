import sys, os
import numpy as np
import scipy.constants   as sc
import scipy.interpolate as sip

from .. import tools     as pt
from .. import constants as pc

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../lib')
import spline as sp

def read(pyrat):
  """
  Read a Cross-section (CS) file.
  """

  pt.msg(pyrat.verb-3, "\nReading cross-section files.", pyrat.log, 0)
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
    pt.msg(pyrat.verb-3, "No CS files to read.", pyrat.log, 2)
  else:
    for i in np.arange(pyrat.cs.nfiles):
      pt.msg(pyrat.verb-3, "Read CS file: '{:s}'.".format(pyrat.cs.files[i]),
             pyrat.log, 2)
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
            pt.error("These species: {:s} are not listed in the atmospheric "
                     "file.\n".format(str(absent)), pyrat.log)
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
          pt.error("CS file has an unexpected line: \n'{:s}'".format(line),
                   pyrat.log)

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
        pt.warning("The wavenumber range [{:.3f}, {:.3f}] cm-1 of the CS "
            "file:\n  '{:s}',\ndoes not cover the Pyrat's wavenumber "
            "range: [{:.3f}, {:.3f}] cm-1.".format(
            wavenumber[0], wavenumber[-1], pyrat.cs.files[i],
            pyrat.spec.wn[ 0], pyrat.spec.wn[-1]), pyrat.wlog, pyrat.log)

      # Add arrays to pyrat:
      pyrat.cs.wavenumber.append(wavenumber)
      pyrat.cs.absorption.append(absorption)

      # Screen output:
      pt.msg(pyrat.verb-4, "For {:s} CS,\nRead {:d} wavenumber and {:d} "
        "temperature samples.".format("-".join(pyrat.cs.molecules[i,0:nmol]),
                         pyrat.cs.nwave[i], pyrat.cs.ntemp[i]), pyrat.log, 4)
      pt.msg(pyrat.verb-4, "Temperature sample limits: {:g}--{:g} K".
                          format(pyrat.cs.temp[i][0], pyrat.cs.temp[i][-1]),
                                 pyrat.log, 4)
      pt.msg(pyrat.verb-4, "Wavenumber sample limits: {:.1f}--{:.1f} cm-1".
             format(pyrat.cs.wavenumber[i][0], pyrat.cs.wavenumber[i][-1]),
                    pyrat.log, 4)

      # Wavenumber-interpolated CS:
      iabsorp = np.zeros((pyrat.cs.ntemp[i], pyrat.spec.nwave), np.double)
      for j in np.arange(pyrat.cs.ntemp[i]):
        iabsorp[j] = sp.splinterp(absorption[j], wavenumber, pyrat.spec.wn, 0.0)
      pyrat.cs.iabsorp.append(iabsorp)
      # Array with second derivatives:
      # Note the shape is swapped w.r.t. iabsorp.
      iz      = np.zeros((pyrat.spec.nwave, pyrat.cs.ntemp[i]), np.double)
      wnlow = np.flatnonzero(pyrat.spec.wn >= wavenumber[ 0])[ 0]
      wnhi  = np.flatnonzero(pyrat.spec.wn <= wavenumber[-1])[-1]
      for j in np.arange(wnlow, wnhi):   # FINDME: wnhi+1 ??
        iz[j] = sp.spline_init(iabsorp[:,j], pyrat.cs.temp[i])
      pyrat.cs.iz.append(iz)

  pt.msg(pyrat.verb-3, "Done.", pyrat.log)


def tmp_interpolate(pyrat):
  """
  Interpolate the CS absorption to the planetary model temperature and
  wavenumber samples.
  """
  pt.msg(pyrat.verb-3, "\nBegin CS interpolation.", pyrat.log)

  # Check temperature boundaries:
  if np.any(pyrat.atm.temp < pyrat.cs.tmin):
    icold = np.where(pyrat.atm.temp < pyrat.cs.tmin)[0][0]
    pt.error("The layer {:d} in the atmospheric model has a lower temperature "
             "({:.1f} K) than the lowest allowed CS temperature ({:.1f} K).".
              format(icold, pyrat.atm.temp[icold], pyrat.cs.tmin), pyrat.log)
  if np.any(pyrat.atm.temp > pyrat.cs.tmax):
    ihot = np.where(pyrat.atm.temp > pyrat.cs.tmax)[0][0]
    pt.error("The layer {:d} in the atmospheric model has a higher temperature "
             "({:.1f} K) than the highest allowed CS temperature ({:.1f} K).".
              format(ihot, pyrat.atm.temp[ihot], pyrat.cs.tmax), pyrat.log)

  # Allocate output extinction-coefficient array:
  pyrat.cs.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  for i in np.arange(pyrat.cs.nfiles):
    cs_absorption = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))
    wnlow = np.flatnonzero(pyrat.spec.wn >= pyrat.cs.wavenumber[i][ 0])[ 0]
    wnhi  = np.flatnonzero(pyrat.spec.wn <= pyrat.cs.wavenumber[i][-1])[-1]
    for p in np.arange(pyrat.atm.nlayers):
      for wn in np.arange(wnlow, wnhi):
        cs_absorption[p,wn] = sp.splinterp_pt(pyrat.cs.iabsorp[i],
                                    pyrat.cs.temp[i], pyrat.cs.iz[i][wn],
                                    pyrat.cs.ntemp[i], pyrat.atm.temp[p])

    # Apply density scale factor:
    dens = 1.0
    for j in np.arange(pyrat.cs.nmol[i]):
      # Get index from the pyrat list of molecules:
      imol = np.where(pyrat.mol.name == pyrat.cs.molecules[i,j])[0][0]
      # Densities in amagat:
      dens *= pyrat.atm.d[:,imol] / pc.amagat

    # Compute CS absorption in cm-1 units (broadcasting):
    pyrat.cs.ec += (cs_absorption * np.expand_dims(dens, axis=1))

    pt.msg(pyrat.verb-6, "CS extinction: {}".format(pyrat.cs.ec[:,0]),
           pyrat.log, 2)
  pt.msg(pyrat.verb-3, "Done.", pyrat.log)


def interpolate(pyrat):
  """
  Interpolate the CS absorption to the planetary model temperature and
  wavenumber samples.
  """
  pt.msg(pyrat.verb, "\nBegin CS interpolation.", pyrat.log)

  # Check temperature boundaries:
  if np.any(pyrat.atm.temp < pyrat.cs.tmin):
    icold = np.where(pyrat.atm.temp < pyrat.cs.tmin)[0][0]
    pt.error("The layer {:d} in the atmospheric model has a lower temperature "
             "({:.1f} K) than the lowest allowed CS temperature ({:.1f} K).".
              format(icold, pyrat.atm.temp[icold], pyrat.cs.tmin), pyrat.log)
  if np.any(pyrat.atm.temp > pyrat.cs.tmax):
    ihot = np.where(pyrat.atm.temp > pyrat.cs.tmax)[0][0]
    pt.error("The layer {:d} in the atmospheric model has a higher temperature "
             "({:.1f} K) than the highest allowed CS temperature ({:.1f} K).".
              format(ihot, pyrat.atm.temp[ihot], pyrat.cs.tmax), pyrat.log)

  # Allocate output extinction-coefficient array:
  pyrat.cs.ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))

  # Interpolator will require sorted arrays:
  asort  = np.argsort(pyrat.atm.temp)
  # This one will de-sort to the original order:
  desort = np.argsort(asort)
  sort_temp = pyrat.atm.temp[asort]

  for i in np.arange(pyrat.cs.nfiles):
    # Evaluate the spline:
    biv = sip.RectBivariateSpline(pyrat.cs.temp[i],
                                  pyrat.cs.wavenumber[i],
                                  pyrat.cs.absorption[i])

    # Interpolate:
    cs_absorption = biv(sort_temp, pyrat.spec.wn)
    # Reverse sorting to the original order of the atmospheric layers:
    cs_absorption = cs_absorption[desort]

    # Apply density scale factor:
    dens = 1.0
    for j in np.arange(pyrat.cs.nmol[i]):
      # Get index from the pyrat list of molecules:
      imol = np.where(pyrat.mol.name == pyrat.cs.molecules[i,j])[0][0]
      # Densities in amagat:
      dens *= pyrat.atm.d[:,imol] / pc.amagat

    # Compute CS absorption in cm-1 units (broadcasting):
    pyrat.cs.ec += (cs_absorption * np.expand_dims(dens, axis=1))

    pt.msg(pyrat.verb-40, "CS extinction: {}".format(pyrat.cs.ec[:,0]),
           pyrat.log, 2)
  pt.msg(pyrat.verb, "Done.", pyrat.log)
