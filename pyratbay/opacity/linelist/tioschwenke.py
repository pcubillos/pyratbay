# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    "Tioschwenke",
    ]

import struct
import numpy as np

from ... import constants as pc
from .driver import Linelist


class Tioschwenke(Linelist):
  """
  Notes:
  ------
  Linelist and partition function downloaded from:
    http://kurucz.harvard.edu/molecules/tio/tioschwenke.bin
    http://kurucz.harvard.edu/molecules/tio/tiopart.dat

  There might be a problem with the linebreak character of the partition
  function.  One way to fix is, on vim do: :%s/\r/\r/g
  """
  def __init__(self, dbfile, pffile, log):
    super(Tioschwenke, self).__init__(dbfile, pffile, log)

    # Database name:
    self.name ="Schwenke TiO (1998)"
    # Molecule name:
    self.molecule = "TiO"
    ID, isotopes, mass, ratio = self.get_iso(self.molecule, dbtype="exomol")
    self.isotopes = isotopes
    self.mass     = mass
    self.isoratio = ratio

    self.recsize = 16  # Record size (bytes)
    self.recdata = 10  # Useful record data size (bytes)
    self.ratiolog = np.log(1.0 + 1.0/2000000)
    # Table of logarithms:
    self.tablog   = 10.0**(0.001*(np.arange(32769) - 16384))
    self.pf_isonames = 0  # PF line with isotopes names


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
    rec_wl: integer
       Wavelength value at record irec, as given in dbfile database.
    """
    # Set pointer at required wavelength record location:
    dbfile.seek(irec*self.recsize)
    # Read and extract the wavelength:
    recwl = struct.unpack('ihhh', dbfile.read(self.recdata))[0]

    return recwl


  def dbread(self, iwn, fwn, verb):
    """
    Read the Schwenke TiO database.

    Parameters
    ----------
    iwn: Scalar
        Initial wavenumber limit (in cm-1).
    fwn: Scalar
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
    """
    # Open the file:
    data = open(self.dbfile, "rb")

    # Get the number of lines in the file:
    data.seek(0, 2)                     # Set pointer at the file's end
    nlines = data.tell() / self.recsize # Number of lines (bytes/record_size)

    # Rewrite wavelength limits as given in the Database file:
    iwl = 1.0/(fwn * pc.nm)         # cm to nanometer
    fwl = 1.0/(iwn * pc.nm)
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

    self.log.msg("Starting to read Schwenke database between records "
        f"{istart:,d} and {istop:,d}.", indent=2)

    interval = (istop - istart)/10  # Check-point interval
    if interval == 0:
        interval = 1


    iw   = np.zeros(nread, int)
    ieli = np.zeros(nread, np.short)
    ielo = np.zeros(nread, np.short)
    igf  = np.zeros(nread, np.short)

    i   = 0  # Stored record index
    while (i < nread):
        # Read a record:
        data.seek((istart+i)*self.recsize)
        iw[i], ieli[i], ielo[i], igf[i] = struct.unpack(
            'ihhh', data.read(self.recdata))

        # Print a checkpoint statement every 10% interval:
        if (i % interval) == 0 and i != 0:
            wl = np.exp(iw[i] * self.ratiolog) * pc.nm/pc.um
            self.log.msg("{10*i/interval:5.1f}% completed.", indent=3)
            self.log.debug(
                f"Wavenumber: {1.0/(wl*pc.um):8.2f} cm-1   "
                f"Wavelength: {wl:6.3f} um\n"
                f"Elow:     {self.tablog[ielo[i]]:.4e} cm-1   "
                f"gf: {self.tablog[igf[i]]:.4e}   "
                f"Iso ID: {np.abs(ieli[i]) - 8950:2d}", indent=6)
        i += 1

    # Calculate the wavenumber (in cm-1):
    wnumber[:] = 1.0/(np.exp(iw * self.ratiolog) * pc.nm)
    # Get gf from log table:
    gf[:] = self.tablog[igf]
    # Get lowest state energy from log table:
    elow[:] = self.tablog[ielo]
    # Get isotopic index:
    isoID[:] = np.abs(ieli) - 8950

    data.close()

    # Sort by wavenumber (in ascending order):
    return wnumber[::-1], gf[::-1], elow[::-1], isoID[::-1]
