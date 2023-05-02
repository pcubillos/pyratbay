# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import ctypes
import importlib
import multiprocessing as mp

import numpy as np

from .. import tools as pt
from .. import constants as pc
from .. import io as io
from ..lib import _extcoeff as ec


def read_opacity(pyrat, wn_mask=None):
    """
    Read opacity table(s), which are actually cross-section tables
    with units of (cm2 molecule-1).
    """
    ex = pyrat.ex
    log = pyrat.log
    ex.extfile = pyrat.inputs.extfile

    # No need to read anything:
    if pyrat.runmode == 'opacity' or ex.extfile is None:
        return

    if pt.isfile(ex.extfile) == 0:
         missing_files = [
             extfile
             for extfile in ex.extfile
             if pt.isfile(extfile) == 0
         ]
         log.error(f'Missing cross-section files: {missing_files}')

    log.head("\nReading cross-section table file(s):")
    for extfile in ex.extfile:
        log.head(f"  '{extfile}'.")

    # Load opacities into shared memory if and only if possible and needed:
    use_shared_memory = False
    mpi_exists = importlib.util.find_spec('mpi4py') is not None
    if mpi_exists:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        use_shared_memory = comm.size > 1

    # Get dimensions first:
    species, ex.temp, ex.press, ex.wn = io.read_opacity(
        ex.extfile[0], extract='arrays',
    )
    ex.ntemp = len(ex.temp)
    ex.nlayers = len(ex.press)
    ex.nwave = len(ex.wn)

    ex.species = []
    species_per_file = []
    for extfile in ex.extfile:
        efile = os.path.basename(extfile)
        species, temp, press, wn = io.read_opacity(extfile, extract='arrays')
        species_per_file.append(list(species))
        ntemp = len(temp)
        nlayers = len(press)
        nwave = len(wn)

        if ntemp != ex.ntemp or nlayers != ex.nlayers or nwave != ex.nwave:
            log.error(
                f"Shape of the cross-section file '{efile}' "
                "does not match with previous file shapes."
            )
        # Species, temperature (K), pressure (barye), and wavenumber (cm-1):
        if ex.temp is not None:
            value_mismatch = [
                np.any(np.abs(1.0-temp/ex.temp) > 0.01),
                np.any(np.abs(1.0-press/ex.press) > 0.01),
                np.any(np.abs(1.0-wn/ex.wn) > 0.01),
            ]
            if np.any(value_mismatch):
                vals = np.array(['temperature', 'pressure', 'wavenumber'])
                mismatch = ', '.join(vals[value_mismatch])
                log.error(
                    f"Tabulated {mismatch} values in file '{efile}' "
                    "do not match with previous arrays"
                )
        # Add new species:
        ex.species += [
            spec
            for spec in species
            if spec not in ex.species
        ]

    spec_indices = []
    for species in species_per_file:
        spec_indices.append([
            ex.species.index(spec) for spec in species
        ])
    ex.species = np.array(ex.species)
    ex.nspec = len(ex.species)

    if wn_mask is None:
        wn_mask = np.ones(len(ex.wn), bool)
    ex.wn = ex.wn[wn_mask]
    ex.nwave = len(ex.wn)

    # Cross-sections table (cm2 molecule-1):
    cs_shape = (ex.nspec, ex.ntemp, ex.nlayers, ex.nwave)
    if not use_shared_memory:
        ex.etable = np.zeros(cs_shape)
        for idx,extfile in zip(spec_indices, ex.extfile):
            ex.etable[idx] += \
                io.read_opacity(extfile, extract='opacity')[:,:,:,wn_mask]
    else:
        itemsize = MPI.DOUBLE.Get_size()
        if comm.rank == 0:
            nbytes = np.prod(cs_shape) * itemsize
        else:
            nbytes = 0
        # on rank 0, create the shared block
        # else get a handle to it
        win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)

        buf, itemsize = win.Shared_query(0)
        assert itemsize == MPI.DOUBLE.Get_size()
        ex.etable = np.ndarray(buffer=buf, dtype='d', shape=cs_shape)
        if comm.rank == 0:
            for idx,extfile in zip(spec_indices, ex.extfile):
                ex.etable[idx] += \
                    io.read_opacity(extfile, extract='opacity')[:,:,:,wn_mask]
        comm.Barrier()

    log.msg(
        f"File(s) have {ex.nspec} species, {ex.ntemp} temperature "
        f"samples, {ex.nlayers} layers, and {ex.nwave} wavenumber samples.",
        indent=2,
    )

    # Set tabulated temperature extrema:
    ex.tmin = np.amin(ex.temp)
    ex.tmax = np.amax(ex.temp)

    with np.printoptions(precision=1):
        str_temp = str(ex.temp)
    with np.printoptions(formatter={'float':'{: .2f}'.format}, threshold=100):
        str_wn = str(ex.wn)
    with np.printoptions(formatter={'float':'{:.3e}'.format}):
        str_press = str(ex.press/pc.bar)
    log.msg(
        f"Species names: {ex.species}\n"
        f"Temperatures (K):\n   {str_temp}\n"
        f"Pressure layers (bar):\n{str_press}\n"
        f"Wavenumber array (cm-1):\n   {str_wn}",
        indent=2,
    )


def compute_opacity(pyrat):
    """
    Compute the cross-sections spectum (cm2 molecule-1) over a tabulated
    grid of temperature, pressure, and wavenumber.
    """
    ex = pyrat.ex
    spec = pyrat.spec
    log = pyrat.log

    if ex.extfile is None:
        log.error('Undefined extinction-coefficient file (extfile)')
    if len(ex.extfile) > 1:
        log.error(
            'Computing cross-section table, but there is more than one'
            'cross-section file set (extfile)'
        )

    extfile = ex.extfile[0]
    log.head(f"\nGenerating new cross-section table file:\n  '{extfile}'")
    # Temperature boundaries check:
    if ex.tmin < pyrat.lt.tmin:
        log.error(
            'Requested cross-section table temperature '
            f'(tmin={ex.tmin:.1f} K) below the lowest available TLI '
            f'temperature ({pyrat.lt.tmin:.1f} K)'
        )
    if ex.tmax > pyrat.lt.tmax:
        log.error(
            'Requested cross-section table temperature '
            f'(tmax={ex.tmax:.1f} K) above the highest available TLI '
            f'temperature ({pyrat.lt.tmax:.1f} K)'
        )

    # Create the temperature array:
    ex.ntemp = int((ex.tmax-ex.tmin)/ex.tstep) + 1
    ex.temp = np.linspace(ex.tmin, ex.tmin + (ex.ntemp-1)*ex.tstep, ex.ntemp)

    with np.printoptions(formatter={'float':'{:.1f}'.format}):
        log.msg(f"Temperature sample (K):\n {ex.temp}", indent=2)

    # Evaluate the partition function at the given temperatures:
    log.msg("Interpolate partition function.", indent=2)
    ex.z = np.zeros((pyrat.iso.niso, ex.ntemp), np.double)
    for i in range(pyrat.iso.niso):
        ex.z[i] = pyrat.iso.zinterp[i](ex.temp)

    # Allocate wavenumber, pressure, and isotope arrays:
    ex.wn = spec.wn
    ex.nwave = spec.nwave
    ex.press = pyrat.atm.press
    ex.nlayers = pyrat.atm.nlayers

    # Allocate extinction-coefficient array:
    log.msg("Calculate cross-sections.", indent=2)
    size = ex.nspec * ex.ntemp * ex.nlayers * ex.nwave
    sm_ect = mp.Array(ctypes.c_double, np.zeros(size, np.double))
    ex.etable = np.ctypeslib.as_array(
        sm_ect.get_obj()).reshape((ex.nspec, ex.ntemp, ex.nlayers, ex.nwave))

    # Multi-processing extinction calculation (in C):
    processes = []
    indices = np.arange(ex.ntemp*ex.nlayers) % pyrat.ncpu  # CPU indices
    grid = True
    add = False
    for i in range(pyrat.ncpu):
        args = (pyrat, np.where(indices==i)[0], grid, add)
        proc = mp.get_context('fork').Process(target=extinction, args=args)
        processes.append(proc)
        proc.start()
    for proc in processes:
        proc.join()

    # Store values in file:
    io.write_opacity(extfile, ex.species, ex.temp, ex.press, ex.wn, ex.etable)
    log.head(
        f"Cross-section table written to file: '{extfile}'.",
        indent=2,
    )


def extinction(pyrat, indices, grid=False, add=False):
    """
    Python multiprocessing wrapper for the extinction-coefficient (EC)
    calculation function for the atmospheric layers or EC grid.

    Parameters
    ----------
    pyrat: Pyrat Object
    indices: 1D integer list
       The indices of the atmospheric layer or EC grid to calculate.
    grid: Bool
       If True, compute EC per species for EC grid.
       If False, compute EC for atmospheric layer.
    add: Bool
       If True, co-add EC contribution (cm-1) from all species.
       If False, keep EC contribution (cm2 molec-1) from each species separated.
    """
    atm = pyrat.atm
    if add:  # Total extinction coefficient spectrum (cm-1)
        extinct_coeff = np.zeros((1, pyrat.spec.nwave))
    else:  # Cross-section spectra for each species (cm2 molecule-1)
        extinct_coeff = np.zeros((pyrat.ex.nspec, pyrat.spec.nwave))

    if pyrat.iso.iext is None:
        # Get species indices in opacity table for the isotopes:
        pyrat.iso.iext = np.zeros(pyrat.iso.niso, int)
        for i in range(pyrat.iso.niso):
            if pyrat.iso.imol[i] != -1:
                pyrat.iso.iext[i] = np.where(
                    pyrat.ex.species == pyrat.atm.species[pyrat.iso.imol[i]])[0][0]
            else:
                pyrat.iso.iext[i] = -1  # Isotope not in atmosphere
                # FINDME: find patch for this case in ec.extinction()

    # Turn off verb of all processes except the first:
    verb = pyrat.log.verb
    pyrat.log.verb = (0 in indices) * verb

    for i,index in enumerate(indices):
        ilayer = index % atm.nlayers  # Layer index
        pressure = atm.press[ilayer]  # Layer pressure
        molq = atm.vmr[ilayer]  # Molecular abundance
        density  = atm.d[ilayer]  # Molecular density
        if grid:  # Take from grid
            itemp = int(index / atm.nlayers)  # Temp. index in EC table
            temp = pyrat.ex.temp[itemp]
            ziso = pyrat.ex.z[:,itemp]      # Isotopic ratio
            pyrat.log.msg(
                "Extinction-coefficient table: "
                f"layer {ilayer+1:3d}/{atm.nlayers}, "
                f"iteration {i+1:3d}/{len(indices)}.",
                indent=2,
            )
        else: # Take from atmosphere
            temp = atm.temp [ilayer]
            ziso = pyrat.iso.z  [:,ilayer]  # Isotopic ratio
            pyrat.log.msg(
                f"Calculating extinction at layer {ilayer+1:3d}/{atm.nlayers} "
                f"(T={temp:6.1f} K, p={pressure/pc.bar:.1e} bar).",
                indent=2,
            )

        # Calculate extinction-coefficient in C:
        extinct_coeff[:] = 0.0
        interpolate_flag = int(
            pyrat.spec.resolution is not None
            or pyrat.spec.wnstep is not None
        )
        ec.extinction(
            extinct_coeff,
            pyrat.voigt.profile, pyrat.voigt.size, pyrat.voigt.index,
            pyrat.voigt.lorentz, pyrat.voigt.doppler,
            pyrat.spec.wn, pyrat.spec.own, pyrat.spec.odivisors,
            density, molq, pyrat.atm.mol_radius, pyrat.atm.mol_mass,
            pyrat.iso.imol, pyrat.iso.mass, pyrat.iso.ratio,
            ziso, pyrat.iso.iext,
            pyrat.lt.wn, pyrat.lt.elow, pyrat.lt.gf, pyrat.lt.isoid,
            pyrat.voigt.cutoff, pyrat.ex.ethresh, pressure, temp,
            verb-10, int(add),
            interpolate_flag,
        )
        # Store output:
        if grid:   # Into grid
            pyrat.ex.etable[:, itemp, ilayer] = extinct_coeff
        elif add:  # Into ex.ec array for atmosphere
            pyrat.ex.ec[ilayer:ilayer+1] = extinct_coeff
        else:      # return single-layer EC of given layer
            return extinct_coeff


def get_ec(pyrat, layer):
    """
    Compute per-species extinction coefficient at requested layer.
    """
    # Interpolate:
    if pyrat.ex.extfile is not None:
        exc = np.zeros((pyrat.ex.nspec, pyrat.spec.nwave))
        label = []
        temp = pyrat.atm.temp[layer]
        itemp = np.where(pyrat.ex.temp <= temp)[0][-1]
        if itemp == len(pyrat.ex.temp):
            itemp -= 1
        for i in range(pyrat.ex.nspec):
            imol = np.where(pyrat.atm.species == pyrat.ex.species[i])[0][0]
            label.append(pyrat.atm.species[imol])
            etable = pyrat.ex.etable[i,:,layer,:]
            exc[i] = ((etable[itemp  ] * (pyrat.ex.temp[itemp+1] - temp) +
                       etable[itemp+1] * (temp - pyrat.ex.temp[itemp]  ) ) /
                      (pyrat.ex.temp[itemp+1] - pyrat.ex.temp[itemp])      )
            exc[i] *= pyrat.atm.d[layer, imol]
    # Line-by-line:
    else:
        exc = extinction(pyrat, [layer], grid=False, add=False)
        label = []
        for i in range(pyrat.ex.nspec):
            imol = np.where(pyrat.atm.species == pyrat.ex.species[i])[0][0]
            exc[i] *= pyrat.atm.d[layer,imol]
            label.append(pyrat.atm.species[imol])
    return exc, label
