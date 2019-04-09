# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["hitran"]

import os
import numpy as np

from ... import constants as pc
from .driver import dbdriver


class hitran(dbdriver):
  def __init__(self, dbfile, pffile, log):
    """
    Initialize the basic database info.

    Parameters
    ----------
    dbfile: String
       File with the Database info as given from HITRAN.
    pffile: String
       File with the partition function.
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

    # Get info from HITRAN configuration file:
    self.molID, self.molecule, self.isotopes, self.mass, self.isoratio = \
        self.getiso(fromfile=True)
    # Database name:
    self.name = "HITRAN " + self.molecule


  def readwave(self, dbfile, irec):
    """
    Read wavenumber from record irec in dbfile database.

    Parameters
    ----------
    dbfile: File object
       File where to extract the wavelength.
    irec: Integer
       Index of record.

    Returns
    -------
    wavenumber: Float
       Wavelength value in cm-1.
    """
    # Set pointer at required wavenumber record:
    dbfile.seek(irec*self.recsize + self.recwnpos)
    # Read:
    wavenumber = float(dbfile.read(self.recwnlen))

    return wavenumber


  def dbread(self, iwn, fwn, verb):
    """
    Read a HITRAN or HITEMP database (dbfile) between wavenumbers iwn and fwn.

    Parameters
    ----------
    dbfile: String
       A HITRAN or HITEMP database filename.
    iwn: Float
       Initial wavenumber limit (in cm-1).
    fwn: Float
       Final wavenumber limit (in cm-1).
    verb: Integer
       Verbosity threshold.

    Returns
    -------
    wnumber: 1D float ndarray
      Line-transition central wavenumber (centimeter-1).
    gf: 1D float ndarray
      gf value (unitless).
    elow: 1D float ndarray
      Lower-state energy (centimeter-1).
    isoID: 2D integer ndarray
      Isotope index (0, 1, 2, 3, ...).

    Notes
    -----
    The HITRAN data is provided in ASCII format.
    The line transitions are sorted in increasing wavenumber (cm-1) order.
    """
    # Open HITRAN file for reading:
    data = open(self.dbfile, "r")

    # Read first line to get the record size:
    data.seek(0)
    line = data.readline()
    self.recsize = data.tell()

    # Get Total number of transitions in file:
    data.seek(0, 2)
    nlines = data.tell() // self.recsize

    # Find the record index for iwn and fwn:
    istart = self.binsearch(data, iwn, 0,      nlines-1, False)
    istop  = self.binsearch(data, fwn, istart, nlines-1, True)

    # Non-overlaping wavenumber ranges:
    data.seek(0)
    line = data.read(self.recsize)
    DBiwn = float(line[self.recwnpos: self.reclinpos])
    data.seek((nlines-1) * self.recsize)
    line = data.read(self.recsize)
    DBfwn = float(line[self.recwnpos: self.reclinpos])
    if iwn > DBfwn or fwn < DBiwn:
      self.log.warning("Database ('{:s}') wavenumber range ({:.2f}--{:.2f} "
        "cm-1) does not overlap with the requested wavenumber range "
        "({:.2f}--{:.2f} cm-1).".format(os.path.basename(self.dbfile),
                                        DBiwn, DBfwn, iwn, fwn))
      return None

    # Number of records to read:
    nread = istop - istart + 1

    # Allocate arrays for values to extract:
    wnumber = np.zeros(nread, np.double)
    gf      = np.zeros(nread, np.double)
    elow    = np.zeros(nread, np.double)
    isoID   = np.zeros(nread,       int)
    A21     = np.zeros(nread, np.double)  # Einstein A coefficient
    g2      = np.zeros(nread, np.double)  # Lower statistical weight

    self.log.msg("Process {:s} database between records {:,d} and {:,d}.".
                 format(self.name, istart, istop), verb=2, indent=2)
    interval = (istop - istart)//10  # Check-point interval

    i = 0  # Stored record index
    while (i < nread):
      # Read a record:
      data.seek((istart+i) * self.recsize)
      line = data.read(self.recsize)
      # Extract values:
      isoID  [i] = float(line[self.recisopos:self.recwnpos ])
      wnumber[i] = float(line[self.recwnpos: self.reclinpos])
      elow   [i] = float(line[self.recelpos: self.recelend ])
      A21    [i] = float(line[self.recApos:  self.recairpos])
      g2     [i] = float(line[self.recg2pos: self.recg2end ])
      # Print a checkpoint statement every 10% interval:
      if (i % interval) == 0.0 and i != 0:
        gfval = A21[i]*g2[i]*pc.C1/(8.0*np.pi*pc.c)/wnumber[i]**2.0
        self.log.msg("{:5.1f}% completed.".format(10.*i/interval),
                     verb=2, indent=3)
        self.log.msg("Wavenumber: {:8.2f} cm-1   Wavelength: {:6.3f} um\n"
                     "Elow:     {:.4e} cm-1   gf: {:.4e}   Iso ID: {:2d}".
                     format(wnumber[i], 1.0/(wnumber[i]*pc.um),
                            elow[i], gfval, (isoID[i]-1)%10), verb=3, indent=6)
      i += 1

    # Set isotopic index to start counting from 0:
    isoID -= 1
    isoID[np.where(isoID < 0)] = 9 # 10th isotope had index 0 --> 10-1=9

    # Calculate gf using Equation (36) of Simekova (2006):
    gf = A21 * g2 * pc.C1 / (8.0 * np.pi * pc.c) / wnumber**2.0

    data.close()

    # Remove lines with unknown Elow, see Rothman et al. (1996):
    igood = np.where(elow > 0)
    return wnumber[igood], gf[igood], elow[igood], isoID[igood]
