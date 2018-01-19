# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["readkurucz", "kunpack"]

import numpy as np
import scipy.constants   as sc


def readkurucz(kfile, temperature, logg):
  """
  Load a the Kurucz stellar spectrum with parameters closest to requested
  temperature and log(g).

  Parameters
  ----------
  kfile: String
     Path to the kurucz file.
  temperature: Scalar
     Surface temperature in K.
  logg: Scalar
     log10 of surface gravity (g in cgs units).

  Returns
  -------
  starfl: 1D ndarray
     Kurucz stellar flux in ergs s-1 cm-2 cm.
  starwn: 1D ndarray
     Array with wavenumber values in cm-1.
  tmodel: Scalar
     Surface temperature of the model in K.
  gmodel: Scalar
     log10 of the surface gravity for the model (g in cgs units).
  """

  inten, freq, grav, temp, nainten, head = kunpack(kfile, freq=True)

  # Wavenumber in cm^-1
  starwn = freq / sc.c * 1e-2

  # Find the model index with the nearest temp and log(g):
  # Nearest sampled temperature:
  tmodel = temp[np.argmin(np.abs(temp-temperature))]
  # Nearest sampled log(g):
  gmodel = grav[np.argmin(np.abs(grav-logg))]
  imodel = np.where((temp == tmodel) & (grav > gmodel))[0][0]

  # Get the stellar flux:
  starfl = inten[imodel]  # W m^-2 sr^-1 Hz^-1

  # Convert F_freq to F_wavenumber (Hz-1 --> m):
  #   multiply by c.
  # Convert units MKS to cgs:
  #   W m-2 = 1e3 ergs s-1 cm-2
  # Convert intensity (astrophysical flux) to flux:
  #   sr-1 = pi

  # Flux per wavenumber:  ergs s-1 cm-2 cm
  starfl = starfl * 1e3 * np.pi * (1e2 * sc.c)

  return starfl, starwn, tmodel, gmodel


def kunpack(filename, freq=False):
  """
  This function reads a file of stellar spectral intensity models
  from Bob Kurucz (Harvard) and returns its content.

  Parameters
  ----------
  filename: String
     Name of model file.  These come from http://kurucz.harvard.edu/grids.html
  freq: Boolean
     If True, reverse first dimension of model grid and return frequencies 
     instead of wavelengths in Wave.

  Returns
  -------
  inten: 2D ndarray
     Array of shape (nmod, nwavl) with the models brightnesses, with
     nwavl the number of wavelength samples, and nmod the number of models.
     These brightnesses in the file have been multiplied by 4, since they
     are Eddington fluxes (W m-2 sr-1 Hz-1).
  wave: 1D ndarray
     Array of size nwavl of the wavelengths (in meters) of inten. 
     If freq==True, wave contains the frequencies (in Hz).
  grav: 1D ndarray
     Array of size nmod with the log10 of the stellar surface
     gravities (g in cm s^-2) for the models.
  temp: 1D ndarray
     Array of size nmod with the temperature (in K) of the models.
  nainten: 2D ndarray
     Array of shape (nmod, nwavl) of the models brightnesses without
     line absorption.  Same units as inten.
  head: List of strings
     List of size nmod of the one-line header strings for the models.

  Example
  -------
  >>> import kurucz_inten as ki
  >>> import time
  >>> import numpy as np
  >>> import matplotlib.pyplot as plt

  >>> # Read the model file:
  >>> kfile = '/home/esp01/ancil/kurucz/fp00ak2odfnew.pck'
  >>> inten, wave, grav, temp, nainten, head = ki.read(kfile)

  >>> # Plot the intensities vs. frequency in Hz:
  >>> nmod  = len(head)
  >>> wait = 0.05
  >>> plt.figure(0)
  >>> for i in np.arange(nmod):
  >>>   plt.clf()
  >>>   a = plt.loglog(wave, inten[i],   "b")
  >>>   a = plt.loglog(wave, nainten[i], "r")
  >>>   a = plt.xlim(5e-9, 5e-4)
  >>>   a = plt.ylim(1e-15, 5e-5)
  >>>   a = plt.title( "model %d:  T=%d  log10(g)=%.1f"%(i, temp[i], grav[i]))
  >>>   plt.draw()
  >>>   plt.show()
  >>>   time.sleep(wait)

  >>> inten, wave, grav, temp, nainten, head = ki.read(kfile, freq=True)
  >>> # Estimate the luminosity of the sun (~4e26 W):
  >>> Tsun = 5770.0  # Sun's surface temperature
  >>> gsun = 4.44    # Sun's surface gravity
  >>> # Find the appropriate model according to Tsun and gsun:
  >>> isun = np.where((temp > Tsun) & (grav > gsun))[0][0]
  >>> fsun = inten[isun]
  >>> # The first pi converts from brightness to flux.
  >>> Lum = np.trapz(fsun, wave) * np.pi * 4 * np.pi * 6.96e8**2
  >>> print("Luminosity  L = %.4e  W/(m^2 sr Hz)"%(Lum))
  >>> Luminosity  L = 4.4718e+26  W/(m^2 sr Hz)
  >>> # The model is for a star with t=6000 and log g=4.5, so expect
  >>> # more than 4e26 W.
 
  Uncredited Developers
  ---------------------
  Joseph Harrington  (Cornell, UCF)
  Kevin Stevenson  (UCF)
  """
  # Read file into memory:
  f = open(filename, 'r')
  text = f.read()
  f.close()
  text = text.replace('\r','\n')
  filetxt = text.split('\n')

  # Get, parse, and count header lines:
  # Record effective temperatures and gravities:
  head      = []
  temp      = np.zeros(len(filetxt))
  grav      = np.zeros(len(filetxt))-1
  header    = np.zeros(len(filetxt), int)
  startwave = 0
  for i in np.arange(len(filetxt)):
    if filetxt[i].startswith("TEFF"):
      head.append(filetxt[i])
      temp[i]   = float(filetxt[i][ 5:12])
      grav[i]   = float(filetxt[i][22:29])
      header[i] = i
    elif filetxt[i].endswith("END"):
      startwave = i + 1

  temp   = temp  [np.where(temp   !=  0)]
  grav   = grav  [np.where(grav   != -1)]
  header = header[np.where(header !=  0)]
  nmod   = header.size
  nline  = (header[2] - header[1] - 1) / 2  # Omit the header line itself

  # Read and count wavelengths:
  wave = np.zeros(header[0]*len(filetxt[startwave])/10)
  k = 0
  string = ''.join(filetxt[startwave:header[0]])
  for j in np.arange(0, len(string), 10):
    wave[k] = float(string[j:j+10])
    k += 1

  wave = wave[np.where(wave != 0)] * 1e-9  # Convert nm to meters
  nwavl = wave.size

  # Allocate memory for models:
  inten   = np.zeros((nmod, nwavl))
  nainten = np.zeros((nmod, nwavl))

  #LOOP OVER MODELS
  for i in range(0, nmod):
    k = 0
    string1 = ''.join(filetxt[header[i]+1      :header[i]+nline+1  ])
    string2 = ''.join(filetxt[header[i]+nline+1:header[i]+2*nline+1])
    for j in range(0,len(string1),10):
      inten[i,k]   = float(string1[j:j+10])
      nainten[i,k] = float(string2[j:j+10])
      k += 1

  # Convert Eddington fluxes to brightnesses and erg cm-2 to J m-2:
  inten   *= 4.0 * 1e-3
  nainten *= 4.0 * 1e-3

  # Convert to frequency if requested:
  if freq:
    wave    = np.flipud(sc.c / wave)
    inten   = np.fliplr(inten)
    nainten = np.fliplr(nainten)
  
  return inten, wave, grav, temp, nainten, head
