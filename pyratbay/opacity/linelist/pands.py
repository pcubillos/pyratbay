# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'Pands',
    ]

import os
import struct
import numpy as np

from ... import constants as pc
from .driver import Linelist


class Pands(Linelist):
  """Partridge & Schwenke (1997) H2O database reader."""
  def __init__(self, dbfile, pffile, log):
      """
      Initialize P&S database object.

      Parameters
      ----------
      dbfile: String
          File with the Database line-transition info.
      pffile: String
          File with the partition function.
      log: Log object
          An mc3.utils.Log instance to log screen outputs to file.
      """
      super(Pands, self).__init__(dbfile, pffile, log)

      # Isotopes names:
      self.isotopes = ['1H1H16O',   '1H1H17O',   '1H1H18O',   '1H2H16O']
      # Isotopes masses:
      self.mass     = [18.01056468, 19.01478156, 20.01481046, 19.01684143]
      # Isotopic abundance ratio:
      self.isoratio = [0.997000,    0.000508,    0.000508,    0.001984]

      # Molecule name:
      self.molecule = 'H2O'
      # Database name:
      self.name = 'Partridge & Schwenke (1997)'

      self.ratiolog = np.log(1 + 1/2e6)
      # Table of logarithms:
      self.tablog  = 10.0**(0.001*(np.arange(32769) - 16384))
      self.recsize = struct.calcsize('Ihh')
      self.wlsize  = struct.calcsize('I')


  def readwave(self, dbfile, irec):
      """
      Read irec-th wavelength record from FILE dbfile.

      Parameters
      ----------
      dbfile: File object
          File where to extract the wavelength.
      irec: Integer
          Index of record.

      Returns
      -------
      recwl: Unsigned integer
          Wavelength value as given in the P&S binary file.
      """
      dbfile.seek(irec*self.recsize)
      recwl = struct.unpack('I', dbfile.read(self.wlsize))[0]

      return recwl


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
      # Open the binary file:
      data = open(self.dbfile, 'rb')

      # Number of lines in the file:
      data.seek(0, 2)
      nlines = data.tell() // self.recsize

      # Wavelength limits as given in the P&S file:
      fwl = 1.0 / (iwn * pc.nm)
      iwl = 1.0 / (fwn * pc.nm)
      iwav = np.log(iwl) / self.ratiolog
      fwav = np.log(fwl) / self.ratiolog

      # Check non-overlaping ranges:
      DBfwn = 1.0/np.exp(self.ratiolog*self.readwave(data, 0))        / pc.nm
      DBiwn = 1.0/np.exp(self.ratiolog*self.readwave(data, nlines-1)) / pc.nm
      if iwn > DBfwn or fwn < DBiwn:
          self.log.warning("Database ('{:s}') wavenumber range ({:.2f}--{:.2f} "
              "cm-1) does not overlap with the requested wavenumber range "
              "({:.2f}--{:.2f} cm-1).".format(os.path.basename(self.dbfile),
                                              DBiwn, DBfwn, iwn, fwn))
          return None

      # Find the positions of iwav and fwav:
      istart = self.binsearch(data, iwav, 0,      nlines-1, False)
      istop  = self.binsearch(data, fwav, istart, nlines-1, True)
      # Number of records to read
      nread = istop - istart + 1

      # Allocate arrays:
      wnumber = np.zeros(nread, np.double)
      gf      = np.zeros(nread, np.double)
      elow    = np.zeros(nread, np.double)
      isoID   = np.zeros(nread, int)

      iw   = np.zeros(nread, int)
      ielo = np.zeros(nread, np.short)
      igf  = np.zeros(nread, np.short)

      self.log.msg('Process P&S H2O database between records '
          f'{istart:,d} and {istop:,d}.', indent=2)

      interval = (istop - istart)//10  # Check-point interval
      if interval == 0:
          interval = 1

      i = 0
      while i < nread:
          # Read record:
          data.seek((istart+i) * self.recsize)
          iw[i], ielo[i], igf[i] = struct.unpack('Ihh', data.read(self.recsize))
          # Print a checkpoint statement every 10% interval:
          if i%interval == 0 and i != 0:
              wl = np.exp(iw[i] * self.ratiolog) * pc.nm/pc.um
              self.log.msg('{10*i/interval:5.1f}% completed.', indent=3)
              self.log.debug(
                  f'Wavenumber: {1.0/(wl*pc.um):8.2f} cm-1   '
                  f'Wavelength: {wl:6.3f} um\n'
                  f'Elow:     {np.abs(ielo[i]):.4e} cm-1   '
                  f'gf: {4.0*self.tablog[np.abs(igf[i])]:.4e}   '
                  f'Iso ID: {2*(ielo[i] < 0) + 1*(igf[i] < 0):2d}', indent=6)
          i += 1
      data.close()

      # Calculate the wavenumber (in cm-1):
      wnumber[:] = 1.0 / (np.exp(iw * self.ratiolog) * pc.nm)
      gf[:]    = 4.0 * self.tablog[np.abs(igf)]
      elow[:]  = np.abs(ielo)
      # Assign indices for isotopes based on Kurucz's indices - 1:
      isoID[:] = 2*(ielo < 0) + 1*(igf < 0)

      # Sort by increasing wavenumber:
      return wnumber[::-1], gf[::-1], elow[::-1], isoID[::-1]
