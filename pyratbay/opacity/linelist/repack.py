# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'Repack',
    ]

import os
import struct
import numpy as np

from ... import constants as pc
from .driver import Linelist


class Repack(Linelist):
  """Repack database reader."""
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
      super(Repack, self).__init__(dbfile, pffile, log)

      self.fmt = 'dddi'
      self.recsize = struct.calcsize(self.fmt)
      self.dsize   = struct.calcsize('d')

      # Get info from file name:
      self.molecule, self.dbtype = os.path.split(dbfile)[1].split('_')[0:2]
      # Get isotopic info:
      ID, isotopes, mass, ratio = self.get_iso(self.molecule, self.dbtype)
      self.isotopes = isotopes
      self.mass     = mass
      self.isoratio = ratio
      # Database name:
      self.name = f'repack {self.dbtype} {self.molecule}'


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
      wavenumber = struct.unpack('d', dbfile.read(self.dsize))[0]
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
      data = open(self.dbfile, 'rb')
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
      istart = self.binsearch(data, iwn, 0,      nlines-1, 0)
      istop  = self.binsearch(data, fwn, istart, nlines-1, 1)
      # Number of records to read:
      nread = istop - istart + 1

      # Allocate arrays for values to extract:
      wnumber = np.zeros(nread, np.double)
      elow    = np.zeros(nread, np.double)
      gf      = np.zeros(nread, np.double)
      iso     = np.zeros(nread,       int)

      self.log.msg(f'Process {self.name} database between records '
          f'{istart:,d} and {istop:,d}.', indent=2)

      interval = (istop - istart) // 10  # Check-point interval
      if interval == 0:
          interval = 1

      i = 0  # Stored record index
      while (i < nread):
          # Read a record:
          data.seek((istart+i) * self.recsize)
          wnumber[i], elow[i], gf[i], iso[i] = \
              struct.unpack(self.fmt, data.read(self.recsize))
          # Print a checkpoint statement every 10% interval:
          if (i % interval) == 0.0 and i != 0:
              self.log.msg(f'{10*i/interval:5.1f}% completed.', indent=3)
              self.log.debug(
                  f'Wavenumber: {wnumber[i]:8.2f} cm-1   '
                  f'Wavelength: {1.0/(wnumber[i]*pc.um):6.3f} um\n'
                  f'Elow:     {elow[i]:.4e} cm-1   gf: {gf[i]:.4e}   '
                  f'Iso ID: {iso[i]:2d}', indent=6)
          i += 1
      data.close()

      # Unique isotopes and inverse indices:
      uiso, inverse = np.unique(iso, return_inverse=True)
      isonamelen = len(str(np.amax(uiso)))  # Count how many digits
      idx = np.zeros(len(uiso), int)
      for i in np.arange(len(uiso)):
          idx[i] = self.isotopes.index(str(uiso[i]).zfill(isonamelen))
      isoID = idx[inverse]

      return wnumber, gf, elow, isoID
