# ****************************** START LICENSE ******************************
# ******************************* END LICENSE ******************************

__all__ = ["pands"]

import os
import struct
import numpy as np

from ... import tools     as pt
from ... import constants as pc
from .driver import dbdriver

# Directory of db:
DBdir = os.path.dirname(os.path.realpath(__file__))


class vald(dbdriver):
  """
  Notes
  -----
    Download linelist from: http://vald.astro.uu.se/~vald/php/vald.php
       Selecting 'Extract Element' and 'Short format'.
    Download partition functions from:
  """
  def __init__(self, dbfile, pffile, log):
    """
    Initialize Basic data for the Database.

    Parameters
    ----------
    dbfile: String
       File with the Database line-transition info.
    pffile: String
       File with the partition function.
    log: File
       File object to store the log.
    """
    super(vald, self).__init__(dbfile, pffile)

    # log file:
    self.log = log

    # Molecule/atom properties:
    self.molecule, self.isotopes, self.mass, self.isoratio, \
                   self.recsize, self.offset = self.getinfo()

    # Database name:
    self.name = "VALD {:s}".format(self.molecule)
    self.ion = ["I", "II", "III", "IV", "V", "VI"]


  def getinfo(self):
    """
    Doc me.
    """
    if not os.path.isfile(self.dbfile):
      pt.error("VALD database file '{:s}' does not exist.".
                format(self.dbfile), self.log)
    data = open(self.dbfile, "r")

    # Extract name from database:
    offset  = len(data.readline())
    offset += len(data.readline())
    line = data.readline()

    recsize = len(line)
    name    = line.split("'")[1]
    name    = name.split()[0]      # Keep the species only

    # Read atomic info file from inputs folder:
    afile = open(DBdir + '/../../../inputs/atoms.dat', 'r')
    lines = afile.readlines()
    afile.close()

    isotopes = []
    mass     = []
    isoratio = []

    # Get values for our molecule:
    for i in np.arange(len(lines)):
      if lines[i].startswith(name):
        line     = lines[i].split()
        mass    .append(float(line[2]))
        isotopes.append("{:s}".format(line[0]))
        isoratio.append(1.0)

    return name, isotopes, mass, isoratio, recsize, offset


  def readwave(self, dbfile, irec):
    """
    Read wavelength parameter from irec record in dbfile database.
 
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

    Notes
    -----
    To convert to wavelength in nanometers do: exp(recwl * ratiolog)
    """
    # Set pointer at required wavelength record location:
    dbfile.seek(self.offset + irec*self.recsize)
    # Read and extract the wavelength:
    recwl = float(dbfile.read(self.recsize).split(",")[1])
 
    return recwl


  def dbread(self, iwn, fwn, verb, *args):
    """
    Read a VALD database.
 
    Parameters
    ----------
    iwn: Scalar
       Initial wavenumber limit (in cm-1).
    fwn: Scalar
       Final wavenumber limit (in cm-1).
    verb: Integer
       Verbosity threshold.
    args:
       Additional arguments, not needed?.
 
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
    The line transitions are sorted in increasing wavenlength order.
    """
 
    # Open/read the file:
    if not os.path.isfile(self.dbfile):
      pt.error("VALD database file '{:s}' does not exist.".
               format(self.dbfile), log)
    data = open(self.dbfile, "r")
    lines = data.readlines()

    # Get the number of lines in the file:
    ifirst =  2  # Line index of first transition
    ilast  = -1  # Line index of last transition
    while not lines[ilast].startswith("'{:s}".format(self.molecule)):
      ilast -= 1
    nlines = len(lines) - ifirst + ilast + 1

    # Rewrite wavenumber limits as given in the VALD file (wl in Angstrom):
    fwl = 1.0 / (iwn * pc.A)  # cm to Angstrom
    iwl = 1.0 / (fwn * pc.A)
 
    # Find the positions of iwav and fwav:
    istart = self.binsearch(data, iwl, 0,      nlines-1, 0)
    istop  = self.binsearch(data, fwl, istart, nlines-1, 1)

    # Number of records to read
    nread = istop - istart + 1
 
 
    pt.msg(verb-4, "Starting to read VALD database between records {:d} and "
                   "{:d}.".format(istart, istop), self.log, 2)

    interval = (istop - istart)/10  # Check-point interval

    # Line-transition data as given in database:
    wl    = np.zeros(nread, np.double)  # Wavelength (Angstrom)
    elo   = np.zeros(nread, np.double)  # Lower-state energy level (eV)
    loggf = np.zeros(nread, np.double)  # log10(gf)
    isoID = np.zeros(nread, int)

    i   = 0  # Stored record index
    data.seek(istart*self.recsize)
    while (i < nread):
      # Read a record:
      data.seek(self.offset + i*self.recsize)
      line = data.read(self.recsize)
      name, ion = line.split("'")[1].split()
      isoID[i]  = int(ion) - 1
      rec = line.split(",")
      wl[i], elo[i], loggf[i] = rec[1:4]
 
      # Print a checkpoint statement every 10% interval:
      if (i % interval) == 0 and i != 0:
        pt.msg(verb-4, "{:5.1f}% completed.".format(10.*i/interval),
               self.log, 3)
        pt.msg(verb-5,"Wavenumber: {:8.2f} cm-1   Wavelength: {:6.3f} A\n"
                        "Elow:     {:.4e} cm-1   gf: {:.4e}   Iso ID: {:2d}".
                        format(1.0/(wl[i] * pc.A), wl[i], elo[i]*pc.eV,
                               10**loggf[i], isoID[i]), self.log, 6)
      i += 1

    # Store data in two arrays for doubles and integers:
    wnumber = np.zeros(nread, np.double)
    elow    = np.zeros(nread, np.double)
    gf      = np.zeros(nread, np.double)

    # Wavenumber (in cm-1):
    wnumber[:] = 1.0 / (wl * pc.A)
    # Get gf fom log:
    gf[:]      = 10**(loggf)
    # Energy of lowest transition level:
    elow[:]    = elo * pc.eV

    # FINDME: Hack for Na:
    gf[np.where(isoID>2)] = 0.0
    isoID[np.where(isoID>2)] = 0

    data.close()
    # Sort (increasingly) by wavenumber:
    return wnumber[::-1], gf[::-1], elow[::-1], isoID[::-1]
