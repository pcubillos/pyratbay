# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys
import os
import numpy as np

from ... import tools     as pt

# Directory of db:
DBdir = os.path.dirname(os.path.realpath(__file__))
# Add path to ctips source code:
sys.path.append(DBdir + "/../../../modules/pytips")
import pytips as t


class dbdriver(object):
  def __init__(self, dbfile, pffile):
    self.dbfile = dbfile
    self.pffile = pffile


  def getpf(self, verbose=0):
    """
    Compute partition function for specified source.

    Returns
    -------
    temp: 1D float ndarray
       Array with temperature sample.
    PF:   2D float ndarray
       The partition function data for each isotope at each temperature.
    isotopes: List of strings
       The names of the tabulated isotopes
    """
    # Calculate the partition-function from the CTIPS module:
    if self.pffile == "ctips":
      # Number of isotopes:
      niso = len(self.isotopes)

      # Array of temperatures (TIPS range of temperature is 70K to 3000K):
      temp = np.arange(70.0, 3000.1, 10.0)
      ntemp = len(temp)
      # Output array for table of Temperature and PF values:
      PF = np.zeros((niso, ntemp), np.double)

      molID = np.repeat(int(t.molID(self.molecule)), ntemp)
      # Calculate the partition function for each isotope:
      for i in np.arange(niso):
        isoID = np.repeat(int(self.isotopes[i]), ntemp)
        PF[i] = t.tips(molID, isoID, temp)
      return temp, PF, self.isotopes

    # Use polynomial expression:
    elif self.pffile == "poly":
      # Temperature array:
      Temp = np.arange(1000.0, 7001.0, 50.0)
      Ntemp = len(Temp)

      # Number of isotopes:
      Niso = len(self.isotopes)

      # Intialize PF array:
      PF = np.zeros((Niso, Ntemp), np.double)

      # Calculate log(PF) at each temperature and isotope:
      for j in np.arange(Niso):
        for i in np.arange(Ntemp):
          # Formula from Irwin 1981, ApJS 45, 621 (equation #2):
          PF[j,i] = (self.PFcoeffs[j,0]                      +
                     self.PFcoeffs[j,1]* np.log(Temp[i])     +
                     self.PFcoeffs[j,2]*(np.log(Temp[i]))**2 +
                     self.PFcoeffs[j,3]*(np.log(Temp[i]))**3 +
                     self.PFcoeffs[j,4]*(np.log(Temp[i]))**4 +
                     self.PFcoeffs[j,5]*(np.log(Temp[i]))**5 )
      # Get the exponential of log(PF):
      PF = np.exp(PF)

      return Temp, PF, self.isotopes

    # Extract the partition-function from the tabulated file:
    else:
      return self.readpf()


  def dbread(self, iwl, fwl, verbose):
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

    Uncredited developers
    ---------------------
    Madison Stemm  (UCF)
    """
    # Minimum and maximum boundaries where to search:
    imin, imax = ilo, ihi

    # Wavelength/wavenumber of record:
    recwave = 0

    # Start binary search:
    while ihi - ilo > 1:
      # Middle record index:
      irec = (ihi + ilo)/2

      # Read wavelength/wavenumber, depending on linelist format:
      recwave = self.readwave(dbfile, irec)
      # Update search limits:
      if recwave > wave:
        ihi = irec
      else:
        ilo = irec

    # Start linear search:
    if searchup:
      irec = ilo
    else:
      irec = ihi

    # Check value of contiguous entries:
    icheck  = irec  # Record index to check
    bounded = True  # Record's wavelength is still within boundaries

    while bounded:
      # Update irec:
      irec = icheck
      # Check index boundaries:
      if irec == imin or irec == imax:
        break
      # Check wavelength boundaries:
      if searchup:
        icheck += 1
        bounded = self.readwave(dbfile, icheck) < wave
      else:
        icheck -= 1
        bounded = self.readwave(dbfile, icheck) > wave

    # Return record index:
    return irec


  def readpf(self):
    """
    Extract the partition-function and temperature from file.

    Returns
    -------
    temp: 1D float ndarray
       Array with temperature sample.
    PF:   2D float ndarray
       The partition function data for each isotope at each temperature.
    isotopes: List of strings
       The names of the tabulated isotopes
    """
    # Open-read file:
    if not os.path.isfile(self.pffile):
      pt.error("Partition-function file '{:s}' does not exist.".
               format(self.pffile), self.log)
    with open(self.pffile, "r") as f:
      lines = f.readlines()

    # Number of header lines (to skip when reading the tabulated data):
    nskip = 0
    while True:
      line = lines[nskip].strip()
      # Skip blank/empty lines:
      if line == "" or line.startswith('#'):
        pass
      # Read isotopes:
      elif line == "@ISOTOPES":
        isotopes = np.asarray(lines[nskip+1].strip().split())
      # Stop when the tabulated data begins:
      if line == "@DATA":
        nskip += 1
        break
      nskip += 1

    # Number of isotopes:
    niso = len(isotopes)
    # Number of temperature samples:
    ntemp = len(lines) - nskip

    # Allocate output arrays:
    temp = np.zeros(ntemp, np.double)
    PF   = np.zeros((niso, ntemp), np.double)

    # Read the data:
    for i in np.arange(ntemp):
      info = lines[nskip+i].strip().split()
      temp[i] = info[0]
      PF[:,i] = info[1:]

    return temp, PF, isotopes

  def getiso(self, fromfile=False, molname=None, dbtype="hitran"):
    """
    Get isotopic info from isotopes.dat file.

    Parameters
    ----------
    fromfile: String
       If True, extract data based on the database file info (for HITRAN).
    mol: String
       If not None, extract data based on this molecule name.
    dbtype: String
       Database type (for isotope names).

    Returns
    -------
    molID: Integer
       HITRAN molecule ID.
    molname: String
       Molecule's name.
    isotopes: List of strings
       Isotopes names.
    mass: List of floats
       Masses for each isotope.
    isoratio: List of integers
       Isotopic terrestrial abundance ratio.
    """
    if fromfile:
      # Open HITRAN DB file and read first two characters:
      if not os.path.isfile(self.dbfile):
        pt.error("HITRAN database file '{:s}' does not exist.".
                  format(self.dbfile), self.log)
      with open(self.dbfile, "r") as data:
        molID  = data.read(self.recmollen)
      molname = t.molname(int(molID))
    elif molname is None:
      pt.error("Neither fromfile nor mol were specified.", self.log)

    # Read isotopes info file:
    with open(DBdir + '/../../../inputs/isotopes.dat', 'r') as isofile:
      lines = isofile.readlines()

    isotopes = []
    mass     = []
    isoratio = []

    if dbtype   == "hitran":
      iiso = 2
    elif dbtype in ["exomol", "kurucz"]:
      iiso = 3
    # FINDME: Replace elif with else?

    # Get values for our molecule:
    for i in np.arange(len(lines)):
      if lines[i].startswith("#") or lines[i].strip() == "":
        continue
      info = lines[i].split()
      if info[1] == molname:
        molID = info[0]
        isotopes.append(info[iiso])
        isoratio.append(float(info[4]))
        mass.    append(float(info[5]))

    return molID, molname, isotopes, mass, isoratio
