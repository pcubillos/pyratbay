# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'Exomol',
    ]

import os
import numpy as np

from ... import constants as pc
from .driver import Linelist
from ... import tools as pt


class Exomol(Linelist):
  """Exomol database reader."""
  def __init__(self, dbfile, pffile, log):
      """
      Initialize Exomol database object.

      Parameters
      ----------
      dbfile: String
          File with the Database line-transition info.
      pffile: String
          File with the partition function.
      log: Log object
          An mc3.utils.Log instance to log screen outputs to file.
      """
      super(Exomol, self).__init__(dbfile, pffile, log)

      sfile = self.dbfile.replace('trans', 'states')
      if sfile.count('__') == 2:
          suffix = sfile[sfile.rindex('__'):sfile.index('.')]
          sfile = sfile.replace(suffix, '')
      else:
          suffix = ''
      # Check files exist:
      for dfile in [self.dbfile, sfile]:
          if not os.path.isfile(dfile):
              log.error(f"File '{dfile}' for Exomol database does not exist.")

      # Read states:
      if sfile.count('__') == 2:
          sfile = sfile.replace(sfile[sfile.rindex('__'):sfile.index('.')], '')
      with open(sfile, 'r') as f:
          lines = f.readlines()
      nstates = len(lines)
      self.E = np.zeros(nstates, np.double)  # State energy
      self.g = np.zeros(nstates, int)        # State degeneracy (incl. ns)
      for i in np.arange(nstates):
          self.E[i], self.g[i] = lines[i].split()[1:3]

      # Get info from file name:
      self.molecule, self.iso = pt.get_exomol_mol(dbfile)
      # Get isotopic info:
      ID, isotopes, mass, ratio = self.get_iso(self.molecule, dbtype='exomol')
      self.isotopes = isotopes
      self.mass     = mass
      self.isoratio = ratio
      # Database name:
      self.name = 'Exomol ' + self.molecule


  def readwave(self, dbfile, irec):
      """
      Read irec-th wavenumber record from FILE dbfile.

      Parameters
      ----------
      dbfile: File object
          File where to extract the wavenumber.
      irec: Integer
          Index of record.

      Returns
      -------
      wavenumber: Float
          Wavenumber value in cm-1.
      """
      dbfile.seek(irec*self.recsize)
      line = dbfile.readline()
      up  = int(line[ 0:12])
      low = int(line[13:25])
      wavenumber = self.E[up-1] - self.E[low-1]

      return wavenumber


  def dbread(self, iwn, fwn, verb):
      """
      Read line-transition info between wavenumbers iwn and fwn.

      Parameters
      ----------
      iwn: Float
          Lower wavenumber boundary in cm-1.
      fwn: Float
          Upper wavenumber boundary in cm-1.
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
      isoID: 1D integer ndarray
          Isotope index.
      """
      # Open file for reading:
      data = open(self.dbfile, 'r')

      # Read first line to get the record size:
      data.seek(0)
      line = data.readline()
      self.recsize = data.tell()

      # Get Total number of transitions in file:
      data.seek(0, 2)
      nlines = data.tell() // self.recsize

      # Check non-overlaping ranges:
      DBiwn = self.readwave(data, 0)
      DBfwn = self.readwave(data, nlines-1)
      if iwn > DBfwn or fwn < DBiwn:
          self.log.warning(
              f"Database ('{os.path.basename(self.dbfile)}') wavenumber "
              f"range ({DBiwn:.2f}--{DBfwn:.2f} cm-1) does not overlap with "
              f"the requested wavenumber range ({iwn:.2f}--{fwn:.2f} cm-1).")
          return None

      # Find the record index for iwn and fwn:
      istart = self.binsearch(data, iwn, 0,      nlines-1, False)
      istop  = self.binsearch(data, fwn, istart, nlines-1, True)
      # Number of records to read:
      nread = istop - istart + 1

      # Allocate arrays for values to extract:
      isoID = np.zeros(nread, int)
      A21   = np.zeros(nread, np.double)
      upID  = np.zeros(nread, int)
      loID  = np.zeros(nread, int)

      self.log.msg(f'Process {self.name} database between records '
                   f'{istart:,d} and {istop:,d}.', indent=2)

      interval = (istop - istart) // 10  # Check-point interval
      if interval == 0:
          interval = 1

      i = 0
      while (i < nread):
          # Read a record:
          data.seek((istart+i) * self.recsize)
          line = data.read(self.recsize)
          upID[i] = line[ 0:12]
          loID[i] = line[13:25]
          A21 [i] = line[26:36]
          # Print a checkpoint statement every 10% interval:
          if (i % interval) == 0.0  and  i != 0:
              wn = self.E[upID[i]-1] - self.E[loID[i]-1]
              gfval = self.g[upID[i]-1]*A21[i]*pc.C1 / (8.0*np.pi*pc.c) / wn**2
              self.log.msg(f'{10*i/interval:5.1f}% completed.', indent=3)
              self.log.debug(
                  f'Wavenumber: {wn:8.2f} cm-1   '
                  f'Wavelength: {1.0/(wn*pc.um):6.3f} um\n'
                  f'Elow:     {self.E[loID[i]-1]:.4e} cm-1   '
                  f'gf: {gfval:.4e}   '
                  f'Iso ID: {self.isotopes.index(self.iso):2d}', indent=6)
          i += 1
      data.close()

      wnumber = self.E[upID-1] - self.E[loID-1]
      # Equation (36) of Simeckova et al. (2006):
      gf = self.g[upID-1] * A21*pc.C1 / (8.0*np.pi*pc.c) / wnumber**2.0
      elow = self.E[loID-1]
      isoID[:] = self.isotopes.index(self.iso)

      return wnumber, gf, elow, isoID
