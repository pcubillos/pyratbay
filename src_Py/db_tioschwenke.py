# ****************************** START LICENSE ******************************
# ******************************* END LICENSE ******************************

import struct, time
import numpy as np

import ptools     as pt
import pconstants as pc
from db_driver import dbdriver

class tioschwenke(dbdriver):
  """
  Notes:
  ------
  Linelist and partition function downloaded from:
    http://kurucz.harvard.edu/molecules/tio/tioschwenke.bin
    http://kurucz.harvard.edu/molecules/tio/tiopart.dat

  There might be a problem with the linebreak character of the partition
  function.  One way to fix is, on vim do: :%s/\r/\r/g
  """
  def __init__(self, dbfile, pffile):
    super(tioschwenke, self).__init__(dbfile, pffile)

    # Database name:
    self.name ="Schwenke TiO (1998)"
    # Isotopes names:
    self.isotopes = [" 46TiO", "47TiO", "48TiO", "49TiO", "50TiO"]
    # Isotopes mass:
    self.mass = [61.94754403, 62.94667863, 63.94286193, 64.94278573,
                 65.93970673]
    # Isotopic abundance ratio:
    self.isoratio = [0.080, 0.073, 0.738, 0.055, 0.054]
    # Molecule name:
    self.molecule = "TiO"

    self.recsize = 16 # Record size (bytes)
    self.recdata = 10 # Useful record data size (bytes)
    self.ratiolog  = np.log(1.0 + 1.0/2000000)
    # Table of logarithms:
    self.tablog    = 10.0**(0.001*(np.arange(32769) - 16384))
    self.pf_isonames = 0 # PF line with isotopes names 

  def readwave(self, dbfile, irec):
    """
    Read wavelength parameter from irec record in dbfile database.

    Parameters:
    -----------
    dbfile: File object
       File where to extract the wavelength.
    irec: Integer
       Index of record.

    Returns:
    --------
    rec_wl: integer
       Wavelength value at record irec, as given in dbfile database.
    """
    # Set pointer at required wavelength record location:
    dbfile.seek(irec*self.recsize)
    # Read and extract the wavelength:
    recwl = struct.unpack('ihhh', dbfile.read(self.recdata))[0]

    return recwl


  def getinfo(self):
    """
    Get isotopes names from partition function file
    """
    # Open and read the partition function file:
    partDB = open(self.pffile)
    PFlines = partDB.readlines()
    partDB.close()

    # Get isotopes names from first line:
    return PFlines[self.pf_isonames].split()[1:]  # Skip first word


  def dbread(self, iwn, fwn, verbose, *args):
    """
    Read the Schwenke TiO database.
 
    Parameters:
    -----------
    iwn: Scalar
       Initial wavenumber limit (in cm-1).
    fwn: Scalar
       Final wavenumber limit (in cm-1).
    verbose: Integer
       Verbosity threshold.
    args:
       Additional arguments, not needed.
 
    Returns:
    --------
    wnumber: 1D float ndarray
      Line-transition central wavenumber (centimeter-1).
    gf: 1D float ndarray
      gf value (unitless).
    elow: 1D float ndarray
      Lower-state energy (centimeter-1).
    isoID: 2D integer ndarray
      Isotope index (0, 1, 2, 3, ...).
    """
 
    # Open the file:
    data = open(self.dbfile, "rb")
 
    # Get the number of lines in the file:
    data.seek(0, 2)                     # Set pointer at the file's end
    nlines = data.tell() / self.recsize # Number of lines (bytes/record_size)
 
    # Rewrite wavelength limits as given in the Database file:
    iwl = 1.0/(fwn * pc.NTC)         # cm to nanometer
    fwl = 1.0/(iwn * pc.NTC)
    iwav = np.log(iwl) / self.ratiolog
    fwav = np.log(fwl) / self.ratiolog
 
    # Find the positions of iwl and fwl, then jump to wl_i position:
    istart = self.binsearch(data, iwav, 0,      nlines-1, 0)
    istop  = self.binsearch(data, fwav, istart, nlines-1, 1)

    # Number of records to read:
    nread = istop - istart + 1

    # Store data in two arrays for doubles and integers:
    # For Wavenumber, Elow, gf, and isotope ID:
    wnumber = np.zeros(nread, np.double)
    gf      = np.zeros(nread, np.double)
    elow    = np.zeros(nread, np.double)
    isoID   = np.zeros(nread,     int)
 
    pt.msg(verbose, "Starting to read Schwenke database between "
                    "records {:d} and {:d}.".format(istart, istop))

    interval = (istop - istart)/10  # Check-point interval


    iw   = np.zeros(nread, int)
    ieli = np.zeros(nread, np.short)
    ielo = np.zeros(nread, np.short)
    igf  = np.zeros(nread, np.short)

    i   = 0  # Stored record index
    while (i < nread):
      # Read a record:
      data.seek((istart+i)*self.recsize)
      iw[i], ieli[i], ielo[i], igf[i] = struct.unpack('ihhh',
                                                      data.read(self.recdata))
 
      # Print a checkpoint statement every 1/20th interval:
      if verbose > 1:
        if (i % interval) == 0 and i != 0:
          wl = np.exp(iw[i] * self.ratiolog) * pc.NTC/pc.MTC
          pt.msg(verbose-1,"Checkpoint {:5.1f}%".format(10.*i/interval), 2)
          pt.msg(verbose-2,"Wavenumber: {:8.2f} cm-1   Wavelength: {:6.3f} um\n"
                          "Elow:     {:.4e} cm-1   gf: {:.4e}   Iso ID: {:2d}".
                           format(1.0/ (wl * pc.MTC), wl,
                                  self.tablog[ielo[i]], self.tablog[igf[i]],
                                  np.abs(ieli[i]) - 8950))
      i += 1

    # Calculate the wavenumber (in cm-1):
    wnumber[:] = 1.0/(np.exp(iw * self.ratiolog) * pc.NTC)
    # Get gf from log table:
    gf[:]      = self.tablog[igf]
    # Get lowest state energy from log table:
    elow[:]    = self.tablog[ielo]
    # Get isotopic index:
    isoID[:]   = np.abs(ieli) - 8950

    pt.msg(verbose, "Done.\n")
    data.close()

    # Sort by wavenumber (in ascending order):
    return wnumber[::-1], gf[::-1], elow[::-1], isoID[::-1]
