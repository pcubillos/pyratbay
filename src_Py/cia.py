import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.interpolate as sip

import pconstants as pc
import ptools as pt



def read(pyrat):
  """
  Read a CIA file as given by Borysow.

  Modification History:
  ---------------------
  2014-08-31  patricio  Initial python implementation.
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
      pt.msg(pyrat.verb, "Reading CIA file: '{:s}'.".format(
                                                      pyrat.cia.files[i]), 2)
      f = open(pyrat.cia.files[i], "r")
      lines = f.readlines()
      f.close()

      # Get the molecules involved:
      iline = 0
      while(lines[iline].startswith("#")):
        iline += 1
      pyrat.cia.molecules[i] = lines[iline].strip().split()
      iline += 1

      # Read the sampled temperatures:
      while(lines[iline].startswith("#")):
        iline += 1
      temps = lines[iline].split()[1:]
      pyrat.cia.ntemp[i] = len(temps)
      temp = np.zeros(pyrat.cia.ntemp[i])
      for j in np.arange(pyrat.cia.ntemp[i]):
        temp[j] = temps[j][:-1]  # Trim the 'k' from Kelvin.
      # Store the temperature array:
      pyrat.cia.temp.append(temp)
      iline += 1

      # Read the data:
      while(lines[iline].startswith("#")):
        iline += 1
      # Get number of wavenumber samples:
      pyrat.cia.nwave[i] = len(lines) - iline
      # Allocate the wavenumber and absorption arrays:
      wavenumber = np.zeros(pyrat.cia.nwave[i], np.double)
      # Allocate the absorption (in cm-1 amagat-2):
      absorption = np.zeros((pyrat.cia.nwave[i], pyrat.cia.ntemp[i]), np.double)
      for j in np.arange(pyrat.cia.nwave[i]):
        data = lines[iline].split()
        wavenumber[j] = data[0]
        absorption[j,:] = data[1:]
        iline += 1

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

  for i in np.arange(pyrat.cia.nfiles):
    # Get index from the pyrat list of molecules:
    imol1 = np.where(pyrat.mol.name == pyrat.cia.molecules[i,0])[0][0]
    imol2 = np.where(pyrat.mol.name == pyrat.cia.molecules[i,1])[0][0]
    #print(imol1, imol2)
    # FINDME: molecule not found exception

    # Evaluate the spline:
    #pt.msg(pyrat.verb, "{} {}".format(pyrat.wn[0], pyrat.wn[-1]))
    #pt.msg(pyrat.verb, "{} {}".format(pyrat.atm.temp[0], pyrat.atm.temp[-1]))

    interp = sip.interp2d(pyrat.cia.temp[i], pyrat.cia.wavenumber[i], 
                          pyrat.cia.absorption[i], kind='linear')

    # Whiny interp2d wants sorted arrays:
    asort  = np.argsort(pyrat.atm.temp)
    # This one will de-sort to the original order:
    desort = np.argsort(asort)

    # Interpolate:
    sort_temp = pyrat.atm.temp[asort]
    cia_absorption = interp(sort_temp, pyrat.wn)
    # Reverse sorting:
    cia_absorption = cia_absorption[:, desort]

    # Compute CIA absorption in cm-1 units:
    #print(np.shape(cia_absorption), np.shape(pyrat.atm.q[:,imol1]))
    pyrat.cia.extinc = (cia_absorption *           # Broadcasting here
                pyrat.atm.q[:,imol1] * pyrat.atm.q[:,imol2] / pc.amagat**2.0)

  pt.msg(pyrat.verb, "Done.")

