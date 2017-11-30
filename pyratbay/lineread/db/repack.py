# Copyright (c) 2016-2017 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["repack"]

import os
import re
import struct
import numpy as np

from ... import tools     as pt
from ... import constants as pc
from .driver import dbdriver


class repack(dbdriver):
  def __init__(self, dbfile, pffile, log):
    """
    Initialize the basic database info.

    Parameters
    ----------
    dbfile: String
       File with the Database info as given from Exomol.
    pffile: String
       File with the partition function.
    log: FILE
       A log file.
    """
    super(repack, self).__init__(dbfile, pffile)

    # Log file:
    self.log = log

    self.fmt = "dddi"
    self.recsize = struct.calcsize(self.fmt)
    self.dsize   = struct.calcsize("d")

    # Get info from file name:
    self.molecule, self.dbtype = os.path.split(dbfile)[1].split("_")[0:2]
    # Database name:
    self.name = "repack {:s} {:s}".format(self.dbtype, self.molecule)

    # Get isotopic info:
    ID, mol, isotopes, mass, ratio = self.getiso(molname=self.molecule,
                                                 dbtype=self.dbtype)
    self.isotopes = isotopes
    self.mass     = mass
    self.isoratio = ratio


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
    dbfile.seek(irec*self.recsize)
    # Read:
    wavenumber = struct.unpack("d", dbfile.read(self.dsize))[0]

    return wavenumber


  def dbread(self, iwn, fwn, verb, *args):
    """
    Read a repack database (dbfile) between wavenumbers iwn and fwn.

    Parameters
    ----------
    dbfile: String
       An Exomol line-list database filename.
    iwn: Float
       Initial wavenumber limit (in cm-1).
    fwn: Float
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

    Notes
    -----
      The line transitions are sorted in increasing wavenumber (cm-1) order.
    """

    # Open file for reading:
    data = open(self.dbfile, "rb")
    # Get Total number of transitions in file:
    data.seek(0, 2)
    nlines   = data.tell() / self.recsize

    # Find the record index for iwn and fwn:
    istart = self.binsearch(data, iwn, 0,      nlines-1, 0)
    istop  = self.binsearch(data, fwn, istart, nlines-1, 1)

    # Data-base wavenumber ranges:
    DBiwn = self.readwave(data, 0)
    DBfwn = self.readwave(data, nlines-1)

    if iwn > DBfwn or fwn < DBiwn:
      pt.warning(verb-2, "Database ('{:s}') wavenumber range ({:.2f}--{:.2f} "
        "cm-1) does not overlap with the requested wavenumber range "
        "({:.2f}--{:.2f} cm-1).".format(os.path.basename(self.dbfile),
                                        DBiwn, DBfwn, iwn, fwn), self.log, [])
      return None

    # Number of records to read:
    nread = istop - istart + 1
    # Allocate arrays for values to extract:
    wnumber = np.zeros(nread, np.double)
    elow    = np.zeros(nread, np.double)
    gf      = np.zeros(nread, np.double)
    iso     = np.zeros(nread,       int)

    pt.msg(verb-4, "Process repack database between records {:,d} and {:,d}.".
                   format(istart, istop), self.log, 2)
    interval = (istop - istart)/10  # Check-point interval

    i = 0  # Stored record index
    while (i < nread):
      # Read a record:
      data.seek((istart+i) * self.recsize)
      # Extract values:
      wnumber[i], elow[i], gf[i], iso[i] = \
                             struct.unpack(self.fmt, data.read(self.recsize))
      # Print a checkpoint statement every 10% interval:
      if (i % interval) == 0.0  and  i != 0:
        pt.msg(verb-4, "{:5.1f}% completed.".format(10.*i/interval),
               self.log, 3)
        pt.msg(verb-5, "Wavenumber: {:8.2f} cm-1   Wavelength: {:6.3f} um\n"
                       "Elow:     {:.4e} cm-1   gf: {:.4e}   Iso ID: {:2d}".
                         format(wnumber[i], 1.0/(wnumber[i]*pc.um), elow[i],
                                gf[i], iso[i], self.log, 6))
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
