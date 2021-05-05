# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'Hitran',
    ]

import os
import numpy as np

from ... import constants as pc
from .driver import Linelist
from .. import partitions as pf


class Hitran(Linelist):
  """HITRAN/HITEMP database reader."""
  def __init__(self, dbfile, pffile, log):
      """
      Initialize HITRAN database object.

      Parameters
      ----------
      dbfile: String
          File with the Database line-transition info.
      pffile: String
          File with the partition function.
      log: Log object
          An mc3.utils.Log instance to log screen outputs to file.
      """
      super(Hitran, self).__init__(dbfile, pffile, log)

      self.recsize      =   0 # Record length (will be set in self.dbread())
      self.rec_iso      =   2 # Isotope        position in record
      self.rec_wn       =   3 # Wavenumber     position in record
      self.rec_strength =  15 # Line intensity position in record
      self.rec_A21      =  25 # Einstein coef  position in record
      self.rec_air      =  35 # Air broadening position in record
      self.rec_elow     =  45 # Lower Energy   position in record
      self.rec_elow_end =  55 # Lower Energy   end position
      self.rec_g2       = 146 # Upper state weight position in record
      self.rec_g2_end   = 153 # Upper state weight end position
      self.rec_wn_len   =  12 # Wavenumber record length

      # Open DB file and read first two characters:
      if not os.path.isfile(self.dbfile):
          self.log.error(f"Input database file '{self.dbfile}' does not exist.")
      with open(self.dbfile, "r") as data:
          molID = int(data.read(2))  # Molecule ID is first two characters
      self.molecule = pf.get_tips_molname(molID)

      ID, isotopes, mass, ratio = self.get_iso(self.molecule, dbtype='hitran')

      self.molID = ID
      self.isotopes = isotopes
      self.mass = mass
      self.isoratio = ratio
      # Database name:
      self.name = 'HITRAN ' + self.molecule


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
      dbfile.seek(irec*self.recsize + self.rec_wn)
      wavenumber = float(dbfile.read(self.rec_wn_len))

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
      wnumber = np.zeros(nread, np.double)
      gf      = np.zeros(nread, np.double)
      elow    = np.zeros(nread, np.double)
      isoID   = np.zeros(nread,       int)
      A21     = np.zeros(nread, np.double)
      g2      = np.zeros(nread, np.double)

      self.log.msg(f'Process {self.name} database between records '
          f'{istart:,d} and {istop:,d}.', indent=2)

      interval = (istop - istart) // 10  # Check-point interval
      if interval == 0:
          interval = 1

      i = 0
      while i < nread:
          # Read a record:
          data.seek((istart+i) * self.recsize)
          line = data.read(self.recsize)
          isoID  [i] = float(line[self.rec_iso:self.rec_wn])
          wnumber[i] = float(line[self.rec_wn:self.rec_strength])
          elow   [i] = float(line[self.rec_elow:self.rec_elow_end])
          A21    [i] = float(line[self.rec_A21:self.rec_air])
          g2     [i] = float(line[self.rec_g2:self.rec_g2_end])
          # Print a checkpoint statement every 10% interval:
          if i%interval == 0.0 and i != 0:
              # Equation (36) of Simeckova et al. (2006)
              gfval = g2[i]*A21[i]*pc.C1 / (8.0*np.pi*pc.c) / wnumber[i]**2.0
              self.log.msg(f'{10*i/interval:5.1f}% completed.', indent=3)
              self.log.debug(
                  f'Wavenumber: {wnumber[i]:8.2f} cm-1   '
                  f'Wavelength: {1.0/(wnumber[i]*pc.um):6.3f} um\n'
                  f'Elow:     {elow[i]:.4e} cm-1   '
                  f'gf: {gfval:.4e}   Iso ID: {(isoID[i]-1)%10:2d}', indent=6)
          i += 1
      data.close()

      # Set isotopic index to start counting from 0:
      isoID -= 1
      isoID[np.where(isoID < 0)] = 9 # 10th isotope has index 0 --> 10-1=9

      # Equation (36) of Simeckova (2006):
      gf = g2*A21*pc.C1 / (8.0*np.pi*pc.c) / wnumber**2.0

      # Remove lines with unknown Elow, see Rothman et al. (1996):
      igood = np.where(elow > 0)
      return wnumber[igood], gf[igood], elow[igood], isoID[igood]
