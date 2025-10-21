# Copyright (c) 2021-2025 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Vald',
]

import os
import numpy as np

from ... import constants as pc
from ... import io as io
from .driver import Linelist


class Vald(Linelist):
    """
    Notes
    -----
        Download linelist from: http://vald.astro.uu.se/~vald/php/vald.php
           Selecting 'Extract Element' and 'Short format'.
        Download partition functions from:

    import os
    import re
    import numpy as np
    import mc3.utils as mu
    import pyratbay.constants as pc
    import pyratbay.tools as pt
    from pyratbay.opacity.linelist.driver import Linelist

    dbfile = 'VALD_Fe.dat'
    pffile = 'PF_barklem_Fe.dat'
    log = mu.Log('delete_me_Fe.log')
    self = vald = Vald(dbfile, pffile, log)
    wn_init = 1e4/33.0
    wn_end = 1e4/0.15
    """
    def __init__(self, dbfile, ion, pffile, log):
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
        super(Vald, self).__init__(dbfile, pffile, log)

        # Open/read the file:
        if not os.path.isfile(self.dbfile):
            self.log.error(f"VALD file '{self.dbfile}' does not exist.")
        with open(self.dbfile, 'r') as f:
            self._data = f.readlines()

        # Get the units from the header line:
        #units_line = self._data[1]
        #units = re.findall('\((.*?)\)', data.readline())
        # TBD: check units are the right ones?

        # Molecule/atom properties:
        self.molecule = ion
        self.isotopes = [ion]
        self.isoratio = [1.0]
        self.mass = self.getinfo()

        atom = ion.replace('+', '')
        ion_count = 1 + ion.count('+')
        ion_label = f"'{atom} {ion_count}'"

        # Number of lines in the file:
        self._data = [
            line
            for line in self._data
            if line.startswith(ion_label)
        ]

        # Database name:
        self.name = f'VALD {self.molecule}'


    def getinfo(self):
        """
        Doc me.
        """
        # Read atomic info file from inputs folder:
        with open(f'{pc.ROOT}pyratbay/data/atoms.dat', 'r') as afile:
            atoms = afile.readlines()

        mass = []
        name = self.molecule.replace('+', '')
        # Get values for our molecule:
        for line in atoms:
            if line.startswith(name):
                line = line.split()
                mass.append(float(line[2]))

        return mass


    def readwave(self, dbfile, irec):
        """
        Read irec-th wavenumber record from FILE dbfile.

        Parameters
        ----------
        dbfile: File object
            File where to extract the wavelength.
        irec: Integer
            Index of record.

        Returns
        -------
        wavenumber: Unsigned integer
            Wavenumber value in cm-1.
        """
        line = self._data[irec]
        wn = float(line.split(',')[1])
        return wn


    def dbread(self, wn_init, wn_end, verb):
        """
        Read a VALD database.

        Parameters
        ----------
        wn_init: Scalar
            Initial wavenumber limit (in cm-1).
        wn_end: Scalar
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
        iso_id: 2D integer ndarray
          Isotope index.
        """
        nlines = len(self._data)
        # Check non-overlaping ranges:
        db_wn_init = self.readwave(self.dbfile, 0)
        db_wn_end = self.readwave(self.dbfile, nlines-1)
        if wn_init > db_wn_end or wn_end < db_wn_init:
            self.log.warning(
               f"The database ('{self.dbfile}') wavenumber range "
               f"({db_wn_init:.2f}--{db_wn_end:.2f} cm-1) does not overlap "
               "with the requested wavenumber range "
               f"({wn_init:.2f}--{wn_end:.2f} cm-1)."
            )
            return None

        # Find the positions of wn_init and wn_end:
        istart = self.binsearch(self._data, wn_init, 0, nlines-1, False)
        istop = self.binsearch(self._data, wn_end, istart, nlines-1, True)

        # Number of records to read
        nread = istop - istart + 1

        self.log.msg(
            f"Process {self.name} database between records "
            f"{istart:,d} and {istop:,d}.",
            indent=2,
        )

        interval = (istop - istart)//10  # Check-point interval
        if interval == 0:
            interval = 1

        # Line-transition data as given in database:
        wn = np.zeros(nread, np.double)  # Wavenumber (cm-1)
        elow = np.zeros(nread, np.double)  # Lower-state energy level (cm-1)
        loggf = np.zeros(nread, np.double)  # log10(gf)
        # Use the iso_id label for ionic states instead:
        iso_id = np.zeros(nread, int)

        for i in range(nread):
            line = self._data[i+istart]
            record = line.split(',')

            name, ion = record[0].strip("'").split()
            iso_id[i]  = int(ion) - 1
            iso_id[i]  = 0
            wn[i], elow[i], loggf[i] = record[1:4]

            # Print a checkpoint statement every 10% interval:
            if i%interval == 0 and i != 0:
                self.log.msg(f"{10.*i/interval:5.1f}% completed.", indent=3)
                self.log.debug(
                    f"Wavenumber: {wn[i]:6.3f} A\n"
                    f"Wavelength: {1.0/(wn[i]*pc.um):8.2f} cm-1   "
                    f"Elow:     {elow[i]:.4e} cm-1   "
                    f"gf: {10**loggf[i]:.4e}   Iso ID: {iso_id[i]:2d}",
                    indent=6,
                )

        gf = 10**loggf
        # Mask out ions not found in the pf_file?

        return wn, gf, elow, iso_id
