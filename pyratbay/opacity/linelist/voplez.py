# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    "Voplez",
    ]

import os
import numpy as np

from ... import constants as pc
from .driver import Linelist


class Voplez(Linelist):
  """
  Download the linelist from:
  """
  def __init__(self, dbfile, pffile, log):
    """
    Initializer.
    """
    super(Voplez, self).__init__(dbfile, pffile, log)

    # Database name:
    self.name = "Bertrand Plez VO"
    # Isotopic names:
    self.isotopes = ["16"]  # I'm Using AFGL naming convention
    # Isotopic masses (amu):
    self.mass     = [66.941]
    # Isotopic abundance ratio:
    self.isoratio = [1.0]

    # Molecule name:
    self.molecule = "VO"

    # Partition-function polynomial coefficients:
    # (from communication with B. Pelz):
    self.PFcoeffs = np.array([[ 6.62090157e+02, -4.03350494e+02,
                                9.82836218e+01, -1.18526504e+01,
                                7.08429905e-01, -1.67235124e-02]])

    # Other utilities:
    self.recsize  = 53  # Record length
    self.recwnpos = 33  # Wavenumber position in record
    self.recelpos = 44  # Elow       position in record
    self.recgfpos = 21  # gf         position in record
    self.recwnlen = 10  # Record lengths
    self.recwnend = 43  # Record lengths
    self.recelend = 50
    self.recgfend = 32


  def readwave(self, dbfile, irec):
    """
    Extract the wavelength from record irec.

    Parameters:
    -----------
    dbfile: File pointer
       Pointer to file being read.
    irec: Integer
       Index of record to read.

    Returns:
    --------
    wl: Float
       The wavelength in microns for record irec.
    """
    # Set pointer at required record:
    dbfile.seek(irec*self.recsize + self.recwnpos)
    # Read record (wavenumber in cm-1):
    wave = float(dbfile.read(self.recwnlen))
    # Convert to wavelength (micron) and return:
    return 1.0/(wave*pc.um)


  def dbread(self, iwn, fwn, verb):
    """
    Read the B. Plez VO database between the wavelengths iwl and fwl.

    Parameters:
    -----------
    iwn: Scalar
       Initial wavenumber limit (in cm-1).
    fwn: Scalar
       Final wavenumber limit (in cm-1).
    verb: Integer
       Verbosity threshold.

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
    Patricio Cubillos (UCF).
    Sarah Blumenthal (UCF).

    Notes:
    ------
    The Plez VO database is an ASCII format.
    The line transitions are sorted in increasing wavelength (micron) order.
    """
    # Open the file:
    if not os.path.isfile(self.dbfile):
        self.log.error(f"Plez VO database file '{self.dbfile}' does not exist.")
    data = open(self.dbfile, "r")
    # Get the total number of transitions:
    data.seek(0, 2)
    nlines   = data.tell() / self.recsize

    # Conver input wavenumber to database units (wavelength):
    iwl = 1.0 / (fwn * pc.um)
    fwl = 1.0 / (iwn * pc.um)

    # Find the record index for iwl and fwl:
    istart = self.binsearch(data, iwl, 0,      nlines-1, 0)
    istop  = self.binsearch(data, fwl, istart, nlines-1, 1)

    # Number of records to read:
    nread = istop - istart + 1

    # Store data in two arrays for doubles and integers:
    wnumber = np.zeros(nread, np.double)
    gf      = np.zeros(nread, np.double)
    elow    = np.zeros(nread, np.double)
    isoID   = np.zeros(nread, int)

    self.log.msg(f"Starting to read Plez VO database between records "
        f"{istart:,d} and {istop:,d}.", indent=2)

    interval = (istop - istart)/10  # Check-point interval
    if interval == 0:
        interval = 1

    i = 0  # Record index counter
    while (i < nread):
        # Read a record:
        data.seek((istop-i) * self.recsize)
        line = data.read(self.recsize)
        # Store values:
        wnumber[i] = float(line[self.recwnpos:self.recwnend])
        gf     [i] = float(line[self.recgfpos:self.recgfend])
        elow   [i] = float(line[self.recelpos:self.recelend])

        # Print a checkpoint statement every 10% interval:
        if (i % interval) == 0.0  and  i != 0:
            self.log.msg(f"{10*i/interval:5.1f}% completed.", indent=3)
            self.log.debug(
                f"Wavenumber: {wnumber[i]:8.2f} cm-1   "
                f"Wavelength: {1.0/(wnumber[i]*pc.um):6.3f} um\n"
                f"Elow:     {elow[i]*pc.eV:.4e} cm-1   "
                f"gf: {gf[i]:.4e}   Iso ID: {isoID[i]:2d}", indent=6)
        i += 1

    # Convert Elow from eV to cm-1:
    elow[:] = elow * pc.eV
    data.close()

    return wnumber, gf, elow, isoID
