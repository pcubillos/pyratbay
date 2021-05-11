# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'Vald',
    ]

import os
import re
import numpy as np

from ... import constants as pc
from .driver import Linelist


class Vald(Linelist):
  """
  Notes
  -----
    Download linelist from: http://vald.astro.uu.se/~vald/php/vald.php
       Selecting 'Extract Element' and 'Short format'.
    Download partition functions from:
  """
  def __init__(self, dbfile, pffile, log):
      """
      Initialize Basic data for the Database.

      Parameters
      ----------
      dbfile: String
          File with the Database line-transition info.
      pffile: String
          File with the partition function.
      log: File
          File object to store the log.
      """
      super(Vald, self).__init__(dbfile, pffile, log)

      # Molecule/atom properties:
      self.molecule, self.isotopes, self.mass, self.isoratio, \
          self.recsize, self.offset = self.getinfo()

      # Database name:
      self.name = f'VALD {self.molecule}'
      self.ion = ["I", "II", "III", "IV", "V", "VI"]


  def getinfo(self):
      """
      Doc me.
      """
      if not os.path.isfile(self.dbfile):
          self.log.error(f"VALD file '{self.dbfile}' does not exist.")

      # Extract name from database:
      with open(self.dbfile, "r") as data:
          offset  = len(data.readline())
          offset += len(data.readline())
          line = data.readline()

      recsize = len(line)
      name = line.split("'")[1]
      name = name.split()[0]      # Keep the species only

      # Read atomic info file from inputs folder:
      with open(f'{pc.ROOT}pyratbay/data/atoms.dat', 'r') as afile:
          atoms = afile.readlines()

      isotopes = []
      mass     = []
      isoratio = []

      # Get values for our molecule:
      for line in atoms:
          if line.startswith(name):
              line = line.split()
              mass.append(float(line[2]))
              isotopes.append("{:s}".format(line[0]))
              isoratio.append(1.0)

      return name, isotopes, mass, isoratio, recsize, offset


  def readwave(self, dbfile, irec):
      """
      Read irec-th wavenumber record from FILE dbfile.

      Parameters
      ----------
      dbfile: File object
          File where to extract the wavelength.
      irec: Integer
          Index of record.

      Returns
      -------
      wavenumber: Unsigned integer
          Wavenumber value in cm-1.
      """
      # Set pointer at required wavelength record location:
      dbfile.seek(self.offset + irec*self.recsize)
      # Read and extract the wavelength:
      wl = float(dbfile.read(self.recsize).split(',')[1])
      return 1.0/(wl*pc.A)


  def dbread(self, iwn, fwn, verb):
      """
      Read a VALD database.

      Parameters
      ----------
      iwn: Scalar
          Initial wavenumber limit (in cm-1).
      fwn: Scalar
          Final wavenumber limit (in cm-1).
      verb: Integer
          Verbosity threshold.

      Returns
      -------
      wnumber: 1D float ndarray
          Line-transition central wavenumber (cm-1).
      gf: 1D float ndarray
          gf value (unitless).
      elow: 1D float ndarray
          Lower-state energy (cm-1).
      isoID: 2D integer ndarray
        Isotope index.
      """
      # Open/read the file:
      data = open(self.dbfile, 'r')

      # Get the units from the header line:
      data.seek(0)
      data.readline()
      units = re.findall('\((.*?)\)', data.readline())
      # TBD: Do something with the units

      # Number of lines in the file:
      nlines = sum(line.startswith(f"'{self.molecule}")
                   for line in open(self.dbfile, 'r'))
      ifirst = 2
      ilast = ifirst + nlines

      # Check non-overlaping ranges:
      DBiwn = self.readwave(data, 0)
      DBfwn = self.readwave(data, nlines-1)
      if iwn > DBfwn or fwn < DBiwn:
          self.log.warning(f"The database ('{self.dbfile}') wavenumber "
             f"range ({DBiwn:.2f}--{DBfwn:.2f} cm-1) does not overlap with "
             f"the requested wavenumber range ({iwn:.2f}--{fwn:.2f} cm-1).")
          return None

      # Rewrite wavenumber limits as given in the VALD file (wl in Angstrom):
      fwl = 1.0 / (iwn * pc.A)
      iwl = 1.0 / (fwn * pc.A)

      # Find the positions of iwn and fwn:
      istart = self.binsearch(data, iwl, 0,      nlines-1, False)
      istop  = self.binsearch(data, fwl, istart, nlines-1, True)

      # Number of records to read
      nread = istop - istart + 1

      self.log.msg(f"Process {self.name} database between records "
                   f"{istart:,d} and {istop:,d}.", indent=2)

      interval = (istop - istart)//10  # Check-point interval
      if interval == 0:
          interval = 1

      # Line-transition data as given in database:
      wl    = np.zeros(nread, np.double)  # Wavelength (Angstrom)
      elo   = np.zeros(nread, np.double)  # Lower-state energy level (eV)
      loggf = np.zeros(nread, np.double)  # log10(gf)
      isoID = np.zeros(nread, int)

      i = 0
      while i < nread:
          # Read record:
          data.seek(self.offset + i*self.recsize)
          line = data.read(self.recsize)
          rec = line.split(',')

          name, ion = rec[0].strip("'").split()
          isoID[i]  = int(ion) - 1
          wl[i], elo[i], loggf[i] = rec[1:4]

          # Print a checkpoint statement every 10% interval:
          if i%interval == 0 and i != 0:
              self.log.msg(f"{10.*i/interval:5.1f}% completed.", indent=3)
              self.log.debug(
                  f"Wavenumber: {1.0/(wl[i]*pc.A):8.2f} cm-1   "
                  f"Wavelength: {wl[i]:6.3f} A\n"
                  f"Elow:     {elo[i]*pc.eV:.4e} cm-1   "
                  f"gf: {10**loggf[i]:.4e}   Iso ID: {isoID[i]:2d}", indent=6)
          i += 1
      data.close()

      wnumber = np.zeros(nread, np.double)
      elow    = np.zeros(nread, np.double)
      gf      = np.zeros(nread, np.double)

      # Wavenumber (in cm-1):
      wnumber[:] = 1.0 / (wl * pc.A)
      gf[:] = 10**(loggf)
      elow[:] = elo * pc.eV

      # FINDME: Hack for Na:
      gf[np.where(isoID>2)] = 0.0
      isoID[np.where(isoID>2)] = 0

      # Sort (increasingly) by wavenumber:
      return wnumber[::-1], gf[::-1], elow[::-1], isoID[::-1]
