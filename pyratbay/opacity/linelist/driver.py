# Copyright (c) 2021-2025 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import numpy as np

from ... import io as io
from ...constants import ROOT
from .. import partitions as pf


class Linelist():
    def __init__(self, dbfile, pffile, log):
        self.dbfile = dbfile
        self.pffile = pffile
        self.log = log


    def getpf(self, verbose=0):
        """
        Compute partition function for specified source.

        Returns
        -------
        temp: 1D float ndarray
            Array with temperature sample.
        PF: 2D float ndarray
            The partition function data for each isotope at each temperature.
        isotopes: List of strings
            The names of the tabulated isotopes
        """
        # Calculate the partition-function from the TIPS module:
        if self.pffile == 'tips':
            pf_data, isotopes, temp = pf.tips(self.molecule)
            return temp, pf_data, isotopes

        # Use polynomial expression:
        elif self.pffile == 'poly':
            temp = np.arange(1000.0, 7001.0, 50.0)
            ntemp = len(temp)
            niso  = len(self.isotopes)

            pf_data = np.zeros((niso, ntemp), np.double)
            for j in range(niso):
                for i in range(ntemp):
                    # Formula from Irwin 1981, ApJS 45, 621 (equation #2):
                    pf_data[j,i] = (
                        self.PFcoeffs[j,0] +
                        self.PFcoeffs[j,1]* np.log(temp[i]) +
                        self.PFcoeffs[j,2]*(np.log(temp[i]))**2 +
                        self.PFcoeffs[j,3]*(np.log(temp[i]))**3 +
                        self.PFcoeffs[j,4]*(np.log(temp[i]))**4 +
                        self.PFcoeffs[j,5]*(np.log(temp[i]))**5
                    )
            # Get the exponential of log(PF):
            pf_data = np.exp(pf_data)
            return temp, pf_data, self.isotopes

        # Extract the partition-function from a tabulated file
        else:
            # TBD: Catch file not found error with self.log
            pf_data, iso, temp = io.read_pf(self.pffile)
            return temp, pf_data, iso.tolist()


    def dbread(self, iwl, fwl, verb):
        """
        Read linelist values for specific database type.
        """
        pass


    def readwave(self, dbfile, irec):
        """
        Read the wavelength parameter as given in each database.
        """
        pass


    def binsearch(self, dbfile, wave, ilo, ihi, searchup=True):
        """
        Do a binary (and then linear) search for wavelength/wavenumber in
        file 'dbfile' between record positions ilo and ihi.

        Parameters
        ----------
        dbfile: File object
           File where to search.
        wave: Scalar
           Target wavelength/wavenumber (as given in each specific database).
        ilo: Integer
           Lowest index record to search.
        ihi: Integer
           highest index record to search.
        searchup: Boolean
           Search up (True) or down (False) the records for duplicate results
           after the binary search.

        Returns:
        --------
        irec:  Integer
           Record index for wave.
        """
        # Minimum and maximum boundaries where to search:
        imin, imax = ilo, ihi

        # Binary search:
        recwave = 0
        while ihi - ilo > 1:
            irec = (ihi + ilo)//2
            recwave = self.readwave(dbfile, irec)
            if recwave > wave:
                ihi = irec
            else:
                ilo = irec

        # Then, linear search (values can be repeated):
        irec = ilo if searchup else ihi
        icheck  = irec
        bounded = True

        while bounded:
            irec = icheck
            if irec == imin or irec == imax:
                break
            if searchup:
                icheck += 1
                bounded = self.readwave(dbfile, icheck) < wave
            else:
                icheck -= 1
                bounded = self.readwave(dbfile, icheck) > wave

        return irec


    def get_iso(self, molname):
        """
        Get isotopic info from isotopes.dat file.

        Parameters
        ----------
        mol: String
            If not None, extract data based on this molecule name.
        dbtype: String
            Database type (for isotope names).

        Returns
        -------
        isotopes: List of strings
            Isotopes names.
        mass: List of floats
            Masses for each isotope.
        isoratio: List of integers
            Isotopic terrestrial abundance ratio.
        """
        iso_data = io.read_isotopes(ROOT + 'pyratbay/data/isotopes.dat')
        name, hit_iso, exo_iso, iso_ratio, iso_mass = iso_data
        isotopes = exo_iso

        isotopes = [str(iso) for mol,iso in zip(name, isotopes) if mol==molname]
        mass = [m for mol,m in zip(name, iso_mass) if mol==molname]
        isoratio = [ratio for mol,ratio in zip(name, iso_ratio) if mol==molname]

        return isotopes, mass, isoratio
