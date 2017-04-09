# Copyright (c) 2016-2017 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["hitran"]

import os
import re
import numpy as np

from ... import tools     as pt
from ... import constants as pc
from .driver import dbdriver

# Directory of db:
DBdir = os.path.dirname(os.path.realpath(__file__))


class exomol(dbdriver):
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
    super(exomol, self).__init__(dbfile, pffile)

    self.recwnpos  =  37 # Wavenumber position in record
    self.recwnlen  =  12 # Wavenumber record length

    # Log file:
    self.log = log

    # Read states:
    with open(dbfile.replace(".trans", ".states"), "r") as f:
      lines = f.readlines()
    nstates = len(lines)
    self.stateID = np.zeros(nstates, int)
    self.E       = np.zeros(nstates, np.double)
    self.J       = np.zeros(nstates, int)
    for i in np.arange(nstates):
      self.stateID[i], self.E[i], self.J[i] = lines[i].split()[:3]

    # Get info from file name:
    s = os.path.split(dbfile)[1].split("_")[0].split("-")
    self.molecule = ""
    isotopes = ""
    mass = 0.0
    for i in np.arange(len(s)):
      match = re.match(r"([0-9]+(?:.[0-9]+)?)([a-z]+)", s[i], re.I)
      self.molecule += match.group(2)
      isotopes += match.group(1)[-1:]
      mass += float(match.group(1))
    self.iso = isotopes
    self.isotopes = [isotopes]
    self.mass     = [mass]  # get somehow
    self.isoratio = [1.0]  # ditto
    # Database name:
    self.name = "Exomol " + self.molecule


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


  def dbread(self, iwn, fwn, verb, *args):
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
    pffile: String
       Partition function filename.

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
    - The HITRAN data is provided in ASCII format.
    - The line transitions are sorted in increasing wavenumber (cm-1) order.
    """

    # Open HITRAN file for reading:
    data = open(self.dbfile, "r")

    # Read first line to get the record size:
    data.seek(0)
    line = data.readline()
    self.recsize = len(line)

    # Get Total number of transitions in file:
    data.seek(0, 2)
    nlines   = data.tell() / self.recsize

    # Find the record index for iwn and fwn:
    istart = self.binsearch(data, iwn, 0,      nlines-1, 0)
    istop  = self.binsearch(data, fwn, istart, nlines-1, 1)
    #print(istart, istop)

    # Non-overlaping wavenumber ranges:
    data.seek(self.recwnpos)
    DBiwn = float(data.read(self.recwnlen))

    data.seek((nlines-1) * self.recsize + self.recwnpos)
    DBfwn = float(data.read(self.recwnlen))
    #print(DBiwn, DBfwn)
    if iwn > DBfwn or fwn < DBiwn:
      pt.warning(verb-2, "Database ('{:s}') wavenumber range ({:.2f}--{:.2f} "
        "cm-1) does not overlap with the requested wavenumber range "
        "({:.2f}--{:.2f} cm-1).".format(os.path.basename(self.dbfile),
                                        DBiwn, DBfwn, iwn, fwn), self.log, [])
      return None

    # Number of records to read:
    nread = istop - istart + 1
    print(nread)
    # Allocate arrays for values to extract:
    wnumber = np.zeros(nread, np.double)
    gf      = np.zeros(nread, np.double)
    elow    = np.zeros(nread, np.double)
    isoID   = np.zeros(nread,       int)
    A21     = np.zeros(nread, np.double)  # Einstein A coefficient
    upID    = np.zeros(nread,       int)
    loID    = np.zeros(nread,       int)

    pt.msg(verb-4, "Process Exomol database between records {:,d} and {:,d}.".
                   format(istart, istop), self.log, 2)
    interval = (istop - istart)/10  # Check-point interval

    i = 0  # Stored record index
    while (i < nread):
      # Read a record:
      data.seek((istart+i) * self.recsize)
      line = data.read(self.recsize)
      # Extract values:
      upID[i], loID[i], A21[i], wnumber[i] = line.split()
      # Print a checkpoint statement every 10% interval:
      if (i % interval) == 0.0  and  i != 0:
        #gfval = A21[i]*g2[i]*pc.C1/(8.0*np.pi*pc.c)/wnumber[i]**2.0
        pt.msg(verb-4, "{:5.1f}% completed.".format(10.*i/interval),
               self.log, 3)
        #pt.msg(verb-5,"Wavenumber: {:8.2f} cm-1   Wavelength: {:6.3f} um\n"
        #                "Elow:     {:.4e} cm-1   gf: {:.4e}   Iso ID: {:2d}".
        #                 format(wnumber[i], 1.0/(wnumber[i]*pc.um),
        #                        elow[i], gfval, (isoID[i]-1)%10), self.log, 6)
      i += 1


    # Calculate gf using Equation (4) of Harris et al. (2006):
    gf = (2*self.J[upID]+1) * A21 * pc.C1 / (8.0 * np.pi * pc.c) / wnumber**2.0
    # Elow:
    elow = self.E[loID]
    # Isotopic index:
    isoID[:] = self.isotopes.index(self.iso)
    data.close()

    # Remove lines with unknown Elow, see Rothman et al. (1996):
    return wnumber, gf, elow, isoID
