# ****************************** START LICENSE ******************************

# ******************************* END LICENSE ******************************

import struct
import numpy as np

import ptools     as pt
import pconstants as pc
from driver import dbdriver

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


  def readwl(self, dbfile, irec):
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
    rec_wl: Float
       Wavelength value at record irec, as given in dbfile database.
    """
    # Set pointer at required wavelength record location:
    dbfile.seek(irec*self.recsize)
    # Read and extract the wavelength:
    rec_wl = struct.unpack('Ihh', dbfile.read(self.recsize))[0]
 
    return rec_wl


  def dbread(self, iwn, fwn, verbose, *args):
    """
    Read the Partridge and Schwenke H2O database (dbfile) between the
    wavelengths iwl and fwl.
 
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
    wnumber: 1D ndarray (double)
      Line-transition central wavenumber (cm-1).
    gf: 1D ndarray (double)
      gf value (unitless).
    elow: 1D ndarray (double)
      Lower-state energe (centimeter^-1).
    isoID: 2D ndarray (integer)
      Isotope index (1, 2, 3, ...).

    Developers:
    -----------
    Madison Stemm      madison.stemm@ucf edu
    Patricio Cubillos  pcubillos@fulbrightmail.org
    """
 
    # Open the binary file:
    data = open(self.dbfile, "rb")
 
    # Get the number of lines in the file:
    data.seek(0, 2)                      # Set pointer at the file's end
    nlines   = data.tell()/ self.recsize # Number of lines (8 bytes per line)
 
    # Rewrite wavelength limits as given in the P&S file:
    #print(iwn, fwn)
    fwav = 1.0 / (iwn * pc.NTC)         # cm to nanometer
    iwav = 1.0 / (fwn * pc.NTC)
    #print(iwav, fwav)
    iwav = np.log(iwav) / self.ratiolog
    fwav = np.log(fwav) / self.ratiolog
 
    #print(iwav, fwav)
    # Find the positions of iwl and fwl, then jump to wl_i position:
    #print("FLAG 01")
    irec = self.binsearch(data, iwav, 0,    nlines, 0)
    #print("FLAG 02")
    frec = self.binsearch(data, fwav, irec, nlines, 1)
    #print("FLAG 03")
    nread = frec - irec + 1  # Number of records to read
 
    # Store data in two arrays for doubles and integers:
    # For Wavelength, Elow, and log(gf):
    wlength = np.zeros(nread, np.double)
    gf      = np.zeros(nread, np.double)
    elow    = np.zeros(nread, np.double)
    # For Isotope index:
    isoID   = np.zeros(nread,     int)
 
    pt.msg(verbose, "Beginning to read P&S database, between "
                    "records {:d} and {:d}.".format(irec, frec))

    interval = (frec - irec)/20  # Check-point interval

    iw   = np.zeros(nread, int)
    ielo = np.zeros(nread, np.short)
    igf  = np.zeros(nread, np.short)

    i   = 0  # Stored record index
    data.seek(irec*self.recsize) 
    while (i < nread):
      # Read a record:
      iw[i], ielo[i], igf[i] = struct.unpack('Ihh', data.read(self.recsize))
 
      # Print a checkpoint statement every 1/20th interval
      if verbose > 1:
        pos = float(data.tell()/self.recsize)
        if (i % interval) == 0 and i != 0:
          pt.msg(verbose-1, "Checkpoint {:2d}/20...".format(i/interval))
          pt.msg(verbose-20, "iwl: {:d}, ielow: {:5d}, igf: {:6d}".
                             format(iw[i], ielo[i], igf[i]))
          pt.msg(verbose-3, "Wavelength:  {:.3f} um,  IsoID: {:d},  Elow: "
                            "{:.2e} cm-1,  gf: {:.2e}".
                         format(np.exp(iw[i] * self.ratiolog) * pc.NTC/pc.MTC,
                                2*(ielo[i] < 0) + 1*(igf[i] < 0),
                                np.abs(ielo[i]), self.tablog[np.abs(igf[i])]))
      i += 1

    # Compute wavelength in microns:
    wlength[:] = np.exp(iw * self.ratiolog) * pc.NTC/pc.MTC
    # Calculate the wavenumber for TLI (cm-1):
    wnumber = 1.0/ (wlength * pc.MTC)
    # Get gf fom log table:
    gf[:]      = self.tablog[np.abs(igf)]
    # Set energy of lowest transition level:
    elow[:]    = np.abs(ielo)
    # Assign indices for isotopes based on Kurucz's indices-1:
    isoID[:]   = 2*(ielo < 0) + 1*(igf < 0)

    pt.msg(verbose, "Done.\n")
    data.close()
    return wnumber, gf, elow, isoID
