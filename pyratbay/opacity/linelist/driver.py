# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import numpy as np

from ... import io as io
from ...constants import ROOT
from .. import partitions as pf


class Linelist(object):
  def __init__(self, dbfile, pffile, log):
      self.dbfile = dbfile
      self.pffile = pffile
      self.log    = log

  def getpf(self, verbose=0):
      """
      Compute partition function for specified source.

      Returns
      -------
      temp: 1D float ndarray
          Array with temperature sample.
      PF: 2D float ndarray
          The partition function data for each isotope at each temperature.
      isotopes: List of strings
          The names of the tabulated isotopes
      """
      # Calculate the partition-function from the CTIPS module:
      if self.pffile == 'tips':
          pf_data, isotopes, temp = pf.tips(self.molecule)
          return temp, pf_data, isotopes

      # Use polynomial expression:
      elif self.pffile == 'poly':
          temp = np.arange(1000.0, 7001.0, 50.0)
          ntemp = len(temp)
          niso  = len(self.isotopes)

          pf_data = np.zeros((niso, ntemp), np.double)
          for j in range(niso):
              for i in range(ntemp):
                  # Formula from Irwin 1981, ApJS 45, 621 (equation #2):
                  pf_data[j,i] = (self.PFcoeffs[j,0]                      +
                                  self.PFcoeffs[j,1]* np.log(temp[i])     +
                                  self.PFcoeffs[j,2]*(np.log(temp[i]))**2 +
                                  self.PFcoeffs[j,3]*(np.log(temp[i]))**3 +
                                  self.PFcoeffs[j,4]*(np.log(temp[i]))**4 +
                                  self.PFcoeffs[j,5]*(np.log(temp[i]))**5 )
          # Get the exponential of log(PF):
          pf_data = np.exp(pf_data)

          return temp, pf_data, self.isotopes

      # Extract the partition-function from the tabulated file:
      else:
          # TBD: Catch file not found error with self.log
          pf_data, iso, temp = io.read_pf(self.pffile)
          return temp, pf_data, iso


  def dbread(self, iwl, fwl, verb):
      """
      Read linelist values for specific database type.
      """
      pass


  def readwave(self, dbfile, irec):
      """
      Read the wavelength parameter as given in each database.
      """
      pass


  def binsearch(self, dbfile, wave, ilo, ihi, searchup=True):
      """
      Do a binary (and then linear) search for wavelength/wavenumber in
      file 'dbfile' between record positions ilo and ihi.

      Parameters
      ----------
      dbfile: File object
         File where to search.
      wave: Scalar
         Target wavelength/wavenumber (as given in each specific database).
      ilo: Integer
         Lowest index record to search.
      ihi: Integer
         highest index record to search.
      searchup: Boolean
         Search up (True) or down (False) the records for duplicate results
         after the binary search.

      Returns:
      --------
      irec:  Integer
         Record index for wave.
      """
      # Minimum and maximum boundaries where to search:
      imin, imax = ilo, ihi

      # Wavelength/wavenumber of record:
      recwave = 0

      # Binary search:
      while ihi - ilo > 1:
          # Middle record index:
          irec = (ihi + ilo)//2
          # Read wavelength/wavenumber, depending on linelist format:
          recwave = self.readwave(dbfile, irec)
          # Update search limits:
          if recwave > wave:
              ihi = irec
          else:
              ilo = irec

      # Linear search:
      if searchup:
          irec = ilo
      else:
          irec = ihi

      icheck  = irec  # Record index to check
      bounded = True  # Record's wavelength is still within boundaries

      while bounded:
          irec = icheck
          # Check index boundaries:
          if irec == imin or irec == imax:
              break
          if searchup:
              icheck += 1
              bounded = self.readwave(dbfile, icheck) < wave
          else:
              icheck -= 1
              bounded = self.readwave(dbfile, icheck) > wave

      return irec


  def get_iso(self, molname, dbtype):
      """
      Get isotopic info from isotopes.dat file.

      Parameters
      ----------
      mol: String
          If not None, extract data based on this molecule name.
      dbtype: String
          Database type (for isotope names).

      Returns
      -------
      molID: Integer
          HITRAN molecule ID.
      isotopes: List of strings
          Isotopes names.
      mass: List of floats
          Masses for each isotope.
      isoratio: List of integers
          Isotopic terrestrial abundance ratio.
      """
      ID, name, hit_iso, exo_iso, iso_ratio, iso_mass = \
          io.read_isotopes(ROOT + 'pyratbay/data/isotopes.dat')

      if dbtype == 'hitran':
          isotopes = hit_iso
      elif dbtype in ['exomol', 'kurucz']:
          isotopes = exo_iso
      else:
          self.log.error(f'Invalid database type: {dbtype}')

      isotopes = [iso   for mol,iso   in zip(name, isotopes) if mol==molname]
      mass     = [m     for mol,m     in zip(name, iso_mass) if mol==molname]
      isoratio = [ratio for mol,ratio in zip(name, iso_ratio) if mol==molname]

      molID = ID[name==molname][0]
      return molID, isotopes, mass, isoratio
