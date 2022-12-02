# Copyright (c) 2021-2022 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import ctypes
import multiprocessing as mp

import numpy as np
import scipy.interpolate as sip

from .. import tools as pt
from .. import constants as pc
from .. import io as io
from ..lib import _extcoeff as ec


def exttable(pyrat):
    """
    Handle extinction-coefficient table (read/calculate/write file).
    """
    # ID list of species with isotopes:
    if np.size(np.where(pyrat.iso.imol>=0)[0]) == 0:
        pyrat.ex.nspec = 0
    else:
        imols = pyrat.iso.imol[np.where(pyrat.iso.imol>=0)]
        pyrat.ex.species = pyrat.mol.name[np.unique(imols)]
        pyrat.ex.nspec = len(pyrat.ex.species)

    # If the extinction file was not defined, skip this step:
    if pyrat.ex.extfile is None:
        pyrat.log.head("\nNo extinction-coefficient table requested.")
        return

    # If the extinction file exists, read it:
    if pt.isfile(pyrat.ex.extfile) == 1:
        pyrat.log.head("\nReading extinction-coefficient table file(s):")
        for extfile in pyrat.ex.extfile:
            pyrat.log.head(f"  '{extfile}'.")
        read_extinction(pyrat)
    # If the extinction file doesn't exist, calculate it:
    else:
        pyrat.log.head(
            "\nGenerating new extinction-coefficient table file:"
            f"\n  '{pyrat.ex.extfile[0]}'."
        )
        calc_extinction(pyrat)


def read_extinction(pyrat):
    """
    Read extinction-coefficient tables from files.
    """
    ex = pyrat.ex
    log = pyrat.log
    spec = pyrat.spec
    ex.species, ex.etable = [], []

    ex.ntemp = 0
    ex.temp = ex.press = ex.wn = None
    for extfile in ex.extfile:
        edata = io.read_opacity(extfile)
        efile = os.path.basename(extfile)
        # Array sizes:
        nspec, ntemp, nlayers, nwave = edata[0]
        size_mismatch = (
            ex.ntemp > 0 and
            (ntemp != ex.ntemp or nlayers != ex.nlayers or nwave != ex.nwave)
        )
        if size_mismatch:
            log.error(
                f"Shape of the extinction-coefficient file '{efile}' "
                "does not match with previous file shapes."
            )
        ex.ntemp, ex.nlayers, ex.nwave = ntemp, nlayers, nwave

        # Species, temperature (K), pressure (barye), and wavenumber (cm-1):
        species, temp, press, wn = edata[1]
        if ex.temp is not None and np.any(np.abs(1-temp/ex.temp) > 0.01):
            log.error(
                f"Tabulated temperature array of file '{efile}' "
                "does not match with the previous temperature arrays."
            )
        if ex.press is not None and np.any(np.abs(1-press/ex.press) > 0.01):
            log.error(
                f"Tabulated pressure array of file '{efile}' "
                "does not match with the previous pressure arrays."
            )
        if ex.wn is not None and np.any(np.abs(1-wn/ex.wn) > 0.01):
            log.error(
                f"Tabulated wavelength array of file '{efile}' "
                "does not match with the previous wavelength arrays."
            )
        ex.temp, ex.press, ex.wn = temp, press, wn
        ex.species += list(species)
        # Extinction-coefficient table (cm2 molecule-1):
        ex.etable.append(edata[2])

    ex.species = np.array(ex.species)
    ex.nspec = len(ex.species)
    ex.etable = np.vstack(ex.etable)
    if np.ndim(ex.etable) == 3:
        ex.etable = np.expand_dims(ex.etable, axis=0)

    log.msg(
        f"File(s) have {ex.nspec} molecules, {ex.ntemp} "
        f"temperature samples, {ex.nlayers} layers, and {ex.nwave} "
        "wavenumber samples.",
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

    # New wavelength sampling:
    if ex.nwave != spec.nwave or np.sum(np.abs(ex.wn-spec.wn)) > 0:
        wn_mask = (ex.wn >= spec.wnlow) & (ex.wn <= spec.wnhigh)

        if spec.wnlow <= ex.wn[0]:
            spec.wnlow = ex.wn[0]
        if spec.wnhigh >= ex.wn[-1]:
            spec.wnhigh = ex.wn[-1]

        if len(set(np.ediff1d(ex.wn))) == 1:
            spec.wnstep = ex.wn[1] - ex.wn[0]
            spec.resolution = None
            #spec.wnlow = ex.wn[0]
            sampling_text = f'sampling rate = {spec.wnstep:.2f} cm-1'
        else:
            g = ex.wn[-2]/ex.wn[-1]
            # Assuming no one would care to set a R with more than 5 decimals:
            spec.resolution = np.round(0.5*(1+g)/(1-g), decimals=5)
            #spec.wnhigh = 2 * ex.wn[-1]/(1+g)
            sampling_text = f'R = {spec.resolution:.1f}'

        # Update wavenumber sampling:
        spec.wn = ex.wn = ex.wn[wn_mask]
        ex.etable = ex.etable[:,:,:,wn_mask]
        spec.nwave = ex.nwave = len(ex.wn)
        # Keep wavenumber oversampling factor:
        spec.ownstep = spec.wnstep / spec.wnosamp
        spec.onwave = (spec.nwave - 1) * spec.wnosamp + 1
        spec.own = np.linspace(spec.wn[0], spec.wn[-1], spec.onwave)

        # Update interpolated stellar spectrum:
        if pyrat.phy.starflux is not None:
            fill_value = (pyrat.phy.starflux[0], pyrat.phy.starflux[-1])
            sinterp = sip.interp1d(
                pyrat.phy.starwn, pyrat.phy.starflux,
                bounds_error=False,
                fill_value=fill_value,
            )
            spec.starflux = sinterp(spec.wn)
        # Update observational variables:
        pyrat.set_filters()
        log.warning(
            "Wavenumber sampling from extinction-coefficient "
            "table does not match the input wavenumber sampling.  Adopting "
            f"tabulated array with {ex.nwave} samples, {sampling_text}, and "
            f"ranges [{spec.wnlow:.2f}, {spec.wnhigh:.2f}] cm-1."
        )


def calc_extinction(pyrat):
    """
    Compute the extinction-coefficient (per species) for a tabulated grid of
    temperatures and pressures, over a wavenumber array.
    """
    ex = pyrat.ex
    spec = pyrat.spec
    log = pyrat.log

    # Temperature boundaries check:
    if ex.tmin < pyrat.lt.tmin:
        log.error(
            'Requested extinction-coefficient table temperature '
            f'(tmin={ex.tmin:.1f} K) below the lowest available TLI '
            f'temperature ({pyrat.lt.tmin:.1f} K).'
        )
    if ex.tmax > pyrat.lt.tmax:
        log.error(
            'Requested extinction-coefficient table temperature '
            f'(tmax={ex.tmax:.1f} K) above the highest available TLI '
            f'temperature ({pyrat.lt.tmax:.1f} K).'
        )

    # Create the temperature array:
    ex.ntemp = int((ex.tmax-ex.tmin)/ex.tstep) + 1
    ex.temp = np.linspace(ex.tmin, ex.tmin + (ex.ntemp-1)*ex.tstep, ex.ntemp)

    with np.printoptions(formatter={'float':'{:.1f}'.format}):
        log.msg(f"Temperature sample (K):\n {ex.temp}", indent=2)

    # Evaluate the partition function at the given temperatures:
    log.msg("Interpolate partition function.", indent=2)
    ex.z = np.zeros((pyrat.iso.niso, ex.ntemp), np.double)
    for i in range(pyrat.lt.ndb):           # For each Database
        for j in range(pyrat.lt.db[i].niso):  # For each isotope in DB
            zinterp = sip.interp1d(
                pyrat.lt.db[i].temp, pyrat.lt.db[i].z[j], kind='slinear',
            )
            ex.z[pyrat.lt.db[i].iiso+j] = zinterp(ex.temp)

    # Allocate wavenumber, pressure, and isotope arrays:
    ex.wn = spec.wn
    ex.nwave = spec.nwave

    ex.press = pyrat.atm.press
    ex.nlayers = pyrat.atm.nlayers

    # Allocate extinction-coefficient array:
    log.msg("Calculate extinction coefficient.", indent=2)
    size = ex.nspec * ex.ntemp * ex.nlayers * ex.nwave
    sm_ect = mp.Array(ctypes.c_double, np.zeros(size, np.double))
    ex.etable = np.ctypeslib.as_array(
        sm_ect.get_obj()).reshape((ex.nspec, ex.ntemp, ex.nlayers, ex.nwave))

    # Multi-processing extinction calculation (in C):
    processes = []
    indices = np.arange(ex.ntemp*ex.nlayers) % pyrat.ncpu  # CPU indices
    for i in range(pyrat.ncpu):
        grid = True
        add = False
        args = (pyrat, np.where(indices==i)[0], grid, add)
        proc = mp.get_context('fork').Process(target=extinction, args=args)
        processes.append(proc)
        proc.start()
    for proc in processes:
        proc.join()

    # Store values in file:
    io.write_opacity(
        ex.extfile[0], ex.species, ex.temp, ex.press, ex.wn, ex.etable)
    log.head(
        f"Extinction-coefficient table written to file: '{ex.extfile[0]}'.",
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
    else:  # Opacity spectrum per molecule (cm2 molecule-1)
        extinct_coeff = np.zeros((pyrat.ex.nspec, pyrat.spec.nwave))

    if pyrat.iso.iext is None:
        # Get species indices in opacity table for the isotopes:
        pyrat.iso.iext = np.zeros(pyrat.iso.niso, int)
        for i in range(pyrat.iso.niso):
            if pyrat.iso.imol[i] != -1:
                pyrat.iso.iext[i] = np.where(
                    pyrat.ex.species == pyrat.mol.name[pyrat.iso.imol[i]])[0][0]
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
            density, molq, pyrat.mol.radius, pyrat.mol.mass,
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
            imol = np.where(pyrat.mol.name == pyrat.ex.species[i])[0][0]
            label.append(pyrat.mol.name[imol])
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
            imol = np.where(pyrat.mol.name == pyrat.ex.species[i])[0][0]
            exc[i] *= pyrat.atm.d[layer,imol]
            label.append(pyrat.mol.name[imol])
    return exc, label
