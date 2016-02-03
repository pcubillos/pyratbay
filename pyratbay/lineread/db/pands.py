# ****************************** START LICENSE ******************************
# ******************************* END LICENSE ******************************

__all__ = ["pands"]

import struct
import numpy as np

from ... import tools     as pt
from ... import constants as pc
from .driver import dbdriver


class pands(dbdriver):
  """
  Notes:
  ------
  Linelist and partition function files downloaded from:
    http://kurucz.harvard.edu/molecules/h2o/h2ofastfix.bin
    http://kurucz.harvard.edu/molecules/h2o/h2opartfn.dat
  """
  def __init__(self, dbfile, pffile):
    """
    Initialize Basic data for the Database.

    Parameters:
    -----------
    dbfile: String
       File with the Database line-transition info.
    pffile: String
       File with the partition function.
    """
    super(pands, self).__init__(dbfile, pffile)

    # Database name:
    self.name = "Partridge & Schwenke (1997)"
    # Isotopes names:
    self.isotopes  = ['1H1H16O',   '1H1H17O',   '1H1H18O',   '1H2H16O'] 
    # Isotopes masses:
    self.mass      = [18.01056468, 19.01478156, 20.01481046, 19.01684143]
    # Isotopic abundance ratio:
    self.isoratio  = [0.997000,    0.000508,    0.000508,    0.001984]

    # Molecule name:
    self.molecule  = "H2O"

    self.ratiolog  = np.log(1 + 1/2e6)
    # Table of logarithms: 
    self.tablog = 10.0**(0.001*(np.arange(32769) - 16384))
    self.recsize     = 8 # Record size


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
    recwl: Unsigned integer
       Wavelength value as given in the P&S binary file.

    Notes:
    ------
    To convert to wavelength in nanometers do: exp(recwl * ratiolog)
    """
    # Set pointer at required wavelength record location:
    dbfile.seek(irec*self.recsize)
    # Read and extract the wavelength:
    recwl = struct.unpack('Ihh', dbfile.read(self.recsize))[0]
 
    return recwl


  def dbread(self, iwn, fwn, verbose, *args):
    """
    Read the Partridge and Schwenke H2O database.
 
    Parameters:
    -----------
    iwn: Scalar
       Initial wavenumber limit (in cm-1).
    fwn: Scalar
       Final wavenumber limit (in cm-1).
    verbose: Integer
       Verbosity threshold.
    args:
       Additional arguments, not needed for pands.
 
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

    Developers:
    -----------
    Madison Stemm      madison.stemm@ucf edu
    Patricio Cubillos  pcubillos@fulbrightmail.org

    Notes:
    ------
    The Partridge & Schwenke database is a binary format.
    The line transitions are sorted in increasing wavenlength order.
    """
 
    # Open the binary file:
    data = open(self.dbfile, "rb")
 
    # Get the number of lines in the file:
    data.seek(0, 2)                      # Set pointer at the file's end
    nlines   = data.tell()/ self.recsize # Number of lines (8 bytes per line)
 
    # Rewrite wavelength limits as given in the P&S file:
    fwl = 1.0 / (iwn * pc.nm)          # cm to nanometer
    iwl = 1.0 / (fwn * pc.nm)
    iwav = np.log(iwl) / self.ratiolog  # As given in file
    fwav = np.log(fwl) / self.ratiolog
 
    # Find the positions of iwav and fwav:
    istart = self.binsearch(data, iwav, 0,      nlines-1, 0)
    istop  = self.binsearch(data, fwav, istart, nlines-1, 1)

    # Number of records to read
    nread = istop - istart + 1
 
    # Store data in two arrays for doubles and integers:
    wnumber = np.zeros(nread, np.double)
    gf      = np.zeros(nread, np.double)
    elow    = np.zeros(nread, np.double)
    isoID   = np.zeros(nread, int)
 
    pt.msg(verbose, "Starting to read P&S database between "
                    "records {:d} and {:d}.".format(istart, istop))

    interval = (istop - istart)/10  # Check-point interval

    iw   = np.zeros(nread, int)
    ielo = np.zeros(nread, np.short)
    igf  = np.zeros(nread, np.short)

    i   = 0  # Stored record index
    data.seek(istart*self.recsize)
    while (i < nread):
      # Read a record:
      iw[i], ielo[i], igf[i] = struct.unpack('Ihh', data.read(self.recsize))
 
      # Print a checkpoint statement every 10% interval:
      if verbose > 1:
        if (i % interval) == 0 and i != 0:
          wl = np.exp(iw[i] * self.ratiolog) * pc.nm/pc.um
          pt.msg(verbose-1,"Checkpoint {:5.1f}%".format(10.*i/interval), 2)
          pt.msg(verbose-2,"Wavenumber: {:8.2f} cm-1   Wavelength: {:6.3f} um\n"
                          "Elow:     {:.4e} cm-1   gf: {:.4e}   Iso ID: {:2d}".
                             format(1.0/ (wl * pc.um), wl,
                                  np.abs(ielo[i]), self.tablog[np.abs(igf[i])],
                                  2*(ielo[i] < 0) + 1*(igf[i] < 0)), 4)
      i += 1

    # Calculate the wavenumber (in cm-1):
    wnumber[:] = 1.0 / (np.exp(iw * self.ratiolog) * pc.nm)
    # Get gf fom log table:
    gf[:]      = self.tablog[np.abs(igf)]
    # Energy of lowest transition level:
    elow[:]    = np.abs(ielo)
    # Assign indices for isotopes based on Kurucz's indices-1:
    isoID[:]   = 2*(ielo < 0) + 1*(igf < 0)

    data.close()
    pt.msg(verbose, "Done.\n")
    # Sort (increasingly) by wavenumber:
    return wnumber[::-1], gf[::-1], elow[::-1], isoID[::-1]
