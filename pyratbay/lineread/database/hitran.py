# Copyright (c) 2016-2020 Patricio Cubillos.
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE).

__all__ = ['hitran']

import os
import numpy as np

from ... import constants as pc
from .driver import dbdriver


class hitran(dbdriver):
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
      super(hitran, self).__init__(dbfile, pffile, log)

      self.recsize   =   0 # Record length (will be set in self.dbread())
      self.recisopos =   2 # Isotope        position in record
      self.recwnpos  =   3 # Wavenumber     position in record
      self.reclinpos =  15 # Line intensity position in record
      self.recApos   =  25 # Einstein coef  position in record
      self.recairpos =  35 # Air broadening position in record
      self.recelpos  =  45 # Low Energy     position in record
      self.recelend  =  55 # Low Energy     end position
      self.recg2pos  = 153 # Low stat weight position in record
      self.recg2end  = 160 # Low stat weight end position
      self.recmollen =   2 # Molecule   record length
      self.recwnlen  =  12 # Wavenumber record length

      ID, mol, isotopes, mass, ratio = self.getiso(fromfile=True)

      self.molID    = ID
      self.molecule = mol
      self.isotopes = isotopes
      self.mass     = mass
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
      dbfile.seek(irec*self.recsize + self.recwnpos)
      wavenumber = float(dbfile.read(self.recwnlen))

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
          self.log.warning("Database ('{:s}') wavenumber range ({:.2f}--{:.2f} "
              "cm-1) does not overlap with the requested wavenumber range "
              "({:.2f}--{:.2f} cm-1).".format(os.path.basename(self.dbfile),
                                              DBiwn, DBfwn, iwn, fwn))
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

      self.log.msg('Process {:s} database between records {:,d} and {:,d}.'.
          format(self.name, istart, istop), indent=2)

      interval = (istop - istart) // 10  # Check-point interval
      if interval == 0:
          interval = 1

      i = 0
      while i < nread:
          # Read a record:
          data.seek((istart+i) * self.recsize)
          line = data.read(self.recsize)
          isoID  [i] = float(line[self.recisopos:self.recwnpos ])
          wnumber[i] = float(line[self.recwnpos: self.reclinpos])
          elow   [i] = float(line[self.recelpos: self.recelend ])
          A21    [i] = float(line[self.recApos:  self.recairpos])
          g2     [i] = float(line[self.recg2pos: self.recg2end ])
          # Print a checkpoint statement every 10% interval:
          if i%interval == 0.0 and i != 0:
              gfval = A21[i]*g2[i]*pc.C1 / (8.0*np.pi*pc.c) / wnumber[i]**2.0
              self.log.msg('{:5.1f}% completed.'.format(10.*i/interval),
                  indent=3)
              self.log.debug(
                  'Wavenumber: {:8.2f} cm-1   Wavelength: {:6.3f} um\n'
                  'Elow:     {:.4e} cm-1   gf: {:.4e}   Iso ID: {:2d}'.
                  format(wnumber[i], 1.0/(wnumber[i]*pc.um), elow[i],
                         gfval, (isoID[i]-1)%10), indent=6)
          i += 1
      data.close()

      # Set isotopic index to start counting from 0:
      isoID -= 1
      isoID[np.where(isoID < 0)] = 9 # 10th isotope has index 0 --> 10-1=9

      # Equation (36) of Simekova (2006):
      gf = A21*g2*pc.C1 / (8.0*np.pi*pc.c) / wnumber**2.0

      # Remove lines with unknown Elow, see Rothman et al. (1996):
      igood = np.where(elow > 0)
      return wnumber[igood], gf[igood], elow[igood], isoID[igood]
