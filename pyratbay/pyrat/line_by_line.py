# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import ctypes
import sys
import multiprocessing as mp

import numpy as np
import scipy.interpolate as sip

from .. import constants as pc
from .. import tools as pt
from . import extinction as ex


class Database():
    """A Line-by-line database"""
    def __init__(self, tli, log):
        # Database name:
        len_name = pt.unpack(tli, 1, "h")
        self.name = pt.unpack(tli, len_name, "s")
        log.msg(f"Data base name: '{self.name}'", indent=2)
        # Molecule name:
        len_name = pt.unpack(tli, 1, "h")
        self.molname = pt.unpack(tli, len_name, "s")
        log.msg(f"Molecule name: '{self.molname}'", indent=2)
        # temperature array:
        self.ntemp = pt.unpack(tli, 1, "h")
        self.niso = pt.unpack(tli, 1, "h")
        self.temp = np.asarray(pt.unpack(tli, self.ntemp, "d"))
        log.msg(
            f"Temperature range: {self.temp[0]:4.1f}--{self.temp[-1]:4.1f} K.",
            indent=2,
        )
        self.iso_pf = np.zeros((self.niso, self.ntemp))

        # Isotopic info:
        self.iso_name = np.zeros(self.niso, 'U20')
        self.iso_mass = np.zeros(self.niso)
        self.iso_ratio = np.zeros(self.niso)

        # Store per-isotope info:
        for j in range(self.niso):
            len_name = pt.unpack(tli, 1, "h")
            self.iso_name[j] = pt.unpack(tli, len_name, "s")
            self.iso_mass[j] = pt.unpack(tli, 1, "d")
            self.iso_ratio[j] = pt.unpack(tli, 1, "d")
            self.iso_pf[j] = np.asarray(pt.unpack(tli, self.ntemp, "d"))

            log.msg(
                f'Isotope: {self.iso_name[j]},  '
                f'mass: {self.iso_mass[j]:.4f} u,  '
                f'isotopic ratio: {self.iso_ratio[j]:.5g}',
                indent=4,
            )

    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write('Database name (name): {:s}', self.name)
        fw.write('Species name (molname):  {:s}', self.molname)
        fw.write('Number of isotopes (niso): {:d}', self.niso)
        fw.write('Number of temperature samples (ntemp): {:d}', self.ntemp)
        fw.write('Temperature (temp, K):\n    {}', self.temp, prec=3, edge=3)
        fw.write('Partition function for each isotope (z):')
        for z in self.iso_pf:
            fw.write('    {}', z, fmt={'float':'{: .3e}'.format}, edge=3)
        return fw.text



class Line_By_Line():
    """All Line-by-line data"""
    def __init__(self, inputs, species, wn_low, wn_high, log, pyrat):
        self.name = 'line by line'
        self.tmin = -np.inf
        self.tmax =  np.inf
        self.db = []
        self.species = []
        self.nspec = 0
        self.ntransitions = 0
        self.ndb = 0
        # LBL data
        self.wn = np.array([], np.double)
        self.elow = np.array([], np.double)
        self.gf = np.array([], np.double)
        self.isoid = np.array([], int)

        self.pyrat = pyrat
        self.nwave = pyrat.spec.nwave
        self.nlayers = pyrat.atm.nlayers
        # Line-transition data file
        self.tlifile = None
        self.ethresh = inputs.ethresh

        # Count number of TLI files:
        if inputs.tlifile is None:
            log.head("\nNo line transition file to read.")
            return

        self.ec = np.zeros((self.nlayers, self.nwave))
        sm_ext = mp.Array(
            ctypes.c_double,
            np.zeros(self.nlayers*self.nwave, np.double),
        )
        self.ec = np.ctypeslib.as_array(
            sm_ext.get_obj()).reshape((self.nlayers, self.nwave))

        log.head("\nReading line-by-line info.")
        with pt.log_error(log):
            pt.file_exists('tlifile', 'TLI', inputs.tlifile)
        self.tlifile = inputs.tlifile

        # Collect all LBL databases from TLI files:
        for tli_file in self.tlifile:
            log.head(f"Read TLI file: '{tli_file}'.", indent=2)
            lbl_data = read_tli_file(tli_file, wn_low, wn_high, log)
            databases, wn, gf, elow, iso_id = lbl_data

            # Add the pre-existing number of isotopes:
            niso = np.sum([db.niso for db in self.db])
            self.isoid = np.concatenate((self.isoid, iso_id+niso))
            self.db += databases
            self.wn = np.concatenate((self.wn, wn))
            self.gf = np.concatenate((self.gf, gf))
            self.elow = np.concatenate((self.elow, elow))
        self.isoid = np.asarray(self.isoid, int)

        # Temperature boundaries:
        self.tmin = np.amax([np.amin(db.temp) for db in self.db])
        self.tmax = np.amin([np.amax(db.temp) for db in self.db])

        # Sizes:
        self.ntransitions = len(self.wn)
        self.ndb = len(self.db)

        # Sort out isotopic info:
        niso = self.niso = np.sum([db.niso for db in self.db])
        total_niso = 0

        self.iso_name = []
        self.iso_mass = np.zeros(niso)
        self.iso_ratio = np.zeros(niso)
        self.iso_atm_index = np.zeros(niso, int)
        self.iso_pf_interp = []
        for j,db in enumerate(self.db):
            iso_mask = np.arange(total_niso, total_niso+db.niso)
            self.iso_name += list(db.iso_name)
            self.iso_mass[iso_mask] = db.iso_mass
            self.iso_ratio[iso_mask] = db.iso_ratio
            if db.molname not in species:
                log.error(
                    f"The species '{db.molname}' for isotopes "
                    f"{db.iso_name} is not present in the atmosphere"
                )
            self.species.append(db.molname)
            self.iso_atm_index[iso_mask] = species.index(db.molname)
            for j in range(db.niso):
                pf_interp = sip.interp1d(db.temp, db.iso_pf[j], kind='slinear')
                self.iso_pf_interp.append(pf_interp)
            total_niso += db.niso


        self.species = np.unique(self.species)
        self.nspec = len(self.species)
        self.mol_index = [
            species.index(mol)
            for mol in self.species
        ]
        # Get species indices in opacity table for each isotope:
        self.iso_mol_index = np.array([
            list(self.species).index(species[i])
            for i in self.iso_atm_index
        ])

        self.iso_name = np.array(self.iso_name)
        log.msg(
            f'Number of isotopes: {niso}\n'
            f"Read a total of {self.ntransitions:,d} line transitions.",
            indent=2,
        )
        log.debug(f"Isotope's molecule indices:\n  {self.iso_atm_index}", indent=2)
        log.head("Read LBL transitions done.\n")


    def calc_extinction_coefficient(
        self, temperature, density, layer=None, skip_mol=[],
    ):
        """
        Calculate the extinction coefficient on the spot over
        temperature and number density profiles.

        Parameters
        ----------
        temperature: 1D float array
            Atmospheric temperature (K)
        density: 2D float array
            Atmospheric number density (gr cm-3).
        layer: Integer
            If not None, compute the extinction coefficient at a single
            layer set by the given index.
        skip_mol: 1D iterable of strings
            Species listed here will be flagged to neglect their opacity.
        """
        # Update partition functions:
        self.iso_pf = np.zeros((self.niso, self.nlayers))
        for i in range(self.niso):
            self.iso_pf[i] = self.iso_pf_interp[i](temperature)

        # Single layer:
        if layer is not None:
            ec = ex.extinction(self.pyrat, [layer], grid=False, add=False)
            for i in range(self.nspec):
                ec[i] *= density[layer]
            return ec

        self.ec[:] = 0.0
        processes = []
        indices = np.arange(self.nlayers) % self.pyrat.ncpu
        for i in range(self.pyrat.ncpu):
            subproc_indices = np.where(indices==i)[0]
            grid = False
            add = True
            args = (self.pyrat, subproc_indices, grid, add, skip_mol)
            proc = mp.get_context('fork').Process(
                target=ex.extinction,
                args=args,
            )
            processes.append(proc)
            proc.start()
        for proc in processes:
            proc.join()

        return self.ec


    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write('Line-transition information:')
        if self.tlifile is None:
            fw.write('No input TLI files.')
            return fw.text
        fw.write('Input TLI files (tlifile): {}', self.tlifile)
        fw.write('Number of databases (ndb): {:d}', self.ndb)
        for db in self.db:
            fw.write('\n'+str(db))
        fw.write(
            '\nTotal number of line transitions (ntransitions): {:,d}\n'
            'Minimum and maximum temperatures (tmin, tmax): [{:.1f}, {:.1f}] K',
            self.ntransitions, self.tmin, self.tmax)
        fw.write('Line-transition isotope IDs (isoid):\n    {}', self.isoid,
            edge=7)
        fw.write('Line-transition wavenumbers (wn, cm-1):\n    {}',
            self.wn,   fmt={'float':'{:.3f}'.format}, edge=3)
        fw.write('Line-transition lower-state energy (elow, cm-1):\n    {}',
            self.elow, fmt={'float':'{: .3e}'.format}, edge=3)
        fw.write('Line-transition gf (gf, cm-1):\n    {}',
            self.gf,   fmt={'float':'{: .3e}'.format}, edge=3)
        fw.write('Line-transition strength threshold (ethresh): {:.2e}',
            self.ethresh)

        fw.write('Isotopes information:')
        fw.write('Number of isotopes (niso): {:d}', self.niso)
        fw.write(
            '\nIsotope  Molecule      Mass    Isotopic   Database'
            '\n            index     g/mol       ratio'
            '\n (name)    (imol)    (mass)     (ratio)')

        iso_index = 0
        db_index = 0
        for i in range(self.niso):
            fw.write(
                '{:>7s}  {:8d}  {:8.4f}   {:.3e}   {}',
                self.iso_name[i], self.iso_atm_index[i], self.iso_mass[i],
                self.iso_ratio[i], self.db[db_index].name,
            )
            iso_index += 1
            if iso_index == self.db[db_index].niso:
                db_index += 1
                iso_index = 0
        return fw.text


def read_tli_file(tli_file, wn_low, wn_high, log):
    """
    Read a TLI file, extract LBL info between given wavenumber ranges.

    Parameters
    ----------
    tli_file: String
        Path to a TLI file.
    wn_low: Float
        Lower wavenumber boundary to extract (cm-1)
    wn_high: Float
        Highest wavenumber boundary to extract (cm-1)
    log: mc3.utils.Log
        Log

    Returns
    -------
    databases: List of Database objects
    wn: 1D float array
        Line transition central wavenumber.
    gf: 1D float array
        Line transition oscillator strength.
    elow: 1D float array
        Line transition lowest-state energy.
    iso_id: 1D short array
        Isotope ID of line transition.
    """
    tli = open(tli_file, "rb")

    tli.seek(0)
    # Read header:
    endian = tli.read(1).decode()
    if endian == 'b':
        log.msg("TLI data storage: Big endian", indent=2)
    elif endian == 'l':
        log.msg("TLI data storage: Little endian", indent=2)

    # Compare endianness:
    endianness = sys.byteorder
    if endianness[0:1] != endian:
        log.error(
            f"Incompatible endianness between TLI file ({endian}) and "
            f"Pyrat ({endianness[0:1]})"
        )

    # Read TLI version:
    TLI_ver, TLI_min, TLI_rev = pt.unpack(tli, 3, "h")
    log.msg(f"TLI version: {TLI_ver}.{TLI_min}.{TLI_rev}.", indent=2)
    if TLI_ver != 6 or TLI_min not in [1,2,3,4,5]:
        log.error(
            "Incompatible TLI version.  The TLI file must be created "
            "with Lineread version 6.1-6.5."
        )

    # Read initial and final wavenumber from TLI:
    lbl_wn_low, lbl_wn_high = pt.unpack(tli, 2, "d")
    log.msg(
        "TLI wavenumber range (cm-1): "
        f"[{lbl_wn_low:.1f}, {lbl_wn_high:.1f}]",
        indent=2,
    )
    # Check TLI and pyrat wavelength ranges:
    if lbl_wn_low > wn_high or lbl_wn_high < wn_low:
        log.warning(
            "TLI wavenumber range "
            f"({lbl_wn_low:.1f}--{lbl_wn_high:.1f} cm-1) does "
            f"not overlap with Pyrat wavenumber range "
            f"({wn_low:.1f}--{wn_high:.1f} cm-1)."
        )
    # TLI range is smaller than the pyrat required range:
    elif lbl_wn_low > wn_low or lbl_wn_high < wn_high:
        log.warning(
            f"TLI wavenumber range "
            f"({lbl_wn_low:.1f}--{lbl_wn_high:.2f} cm-1) does "
            f"not cover the full Pyrat wavenumber range "
            f"({wn_low:.1f}--{wn_high:.1f} cm-1)."
        )

    # Read number of data bases:
    n_db = pt.unpack(tli, 1, "h")
    log.msg(f"Number of data bases in TLI file: {n_db}", indent=2)

    databases = []
    for i in range(n_db):
        db = Database(tli, log)
        databases.append(db)

    # number of line transitions:
    n_transitions = pt.unpack(tli, 1, "i")
    # number of isotopes in line-transition array:
    n_iso = pt.unpack(tli, 1, "i")
    # number of transitions per isotope:
    niso_tran = np.atleast_1d(pt.unpack(tli, n_iso, "i"))

    # Position where the line-transition data begins:
    init_wl = tli.tell()
    init_iso = init_wl + n_transitions*pc.dreclen
    init_el = init_iso + n_transitions*pc.sreclen
    init_gf = init_el + n_transitions*pc.dreclen

    # Count the number of transitions:
    tli.seek(0, 2)
    endrec = tli.tell()
    nrec = (endrec - init_wl)*1.0 / pc.tlireclen
    if nrec != n_transitions:
        log.error(
            f'The remaining data file size ({nrec:.1f}) does not '
            f'correspond to the number of transitions ({n_transitions})'
        )
    log.msg(
        f'There are {n_transitions:,d} line transitions in TLI file.',
        indent=2,
    )

    # Allocate arrays:
    wn = np.zeros(n_transitions)
    isoid = np.zeros(n_transitions, np.short)
    elow = np.zeros(n_transitions)
    gf = np.zeros(n_transitions)

    # Number of line-transitions offset for a given isotope:
    offset = 0
    start = init_wl
    nlt = 0  # Total number of line-transitions read
    for i in range(n_iso):
        # Search lower and higher line-transition indices to read:
        ifirst = pt.binsearch(tli, wn_low,  start, niso_tran[i], upper=False)
        ilast  = pt.binsearch(tli, wn_high, start, niso_tran[i], upper=True)
        # Add offset for this isotope:
        ifirst += offset
        ilast  += offset

        if ifirst < 0 or ilast < 0:
            start  += niso_tran[i]*pc.dreclen
            offset += niso_tran[i]
            continue

        # Print boundaries:
        tli.seek(ifirst*pc.dreclen + init_wl, 0)
        first = pt.unpack(tli, 1, 'd')
        tli.seek(ilast*pc.dreclen  + init_wl, 0)
        last = pt.unpack(tli, 1, 'd')
        log.debug(
            f"Found initial transition ({ifirst:8d}):  {first:13.4f} cm-1\n"
            f"Found final transition ({ilast:8d}):  {last:13.4f} cm-1",
            indent=2,
        )

        # Number of transitions to read:
        nread = ilast - ifirst + 1
        # Get isotope ID:
        i0 = pt.binsearch(tli, 0, start, niso_tran[i], upper=False)
        tli.seek((i0+offset)*pc.sreclen + init_iso, 0)
        isoID = pt.unpack(tli, 1, 'h')

        # Read data into arrays:
        tli.seek(ifirst*pc.dreclen + init_wl,  0)
        wn[nlt:nlt+nread] = pt.unpack(tli, nread, "d")

        tli.seek(ifirst*pc.sreclen + init_iso, 0)
        isoid[nlt:nlt+nread] = pt.unpack(tli, nread, 'h')

        tli.seek(ifirst*pc.dreclen + init_el,  0)
        elow[nlt:nlt+nread] = pt.unpack(tli, nread, 'd')

        tli.seek(ifirst*pc.dreclen + init_gf,  0)
        gf[nlt:nlt+nread] = pt.unpack(tli, nread, 'd')
        log.msg(
            f"Read {nread:11,d} transitions for isotope {isoID:2d}.",
            indent=4,
        )

        start += niso_tran[i]*pc.dreclen
        offset += niso_tran[i]
        nlt += nread

    tli.close()

    return (
        databases,
        wn[0:nlt],
        gf[0:nlt],
        elow[0:nlt],
        isoid[0:nlt],
    )

