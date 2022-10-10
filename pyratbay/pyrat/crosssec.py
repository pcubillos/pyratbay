# Copyright (c) 2021-2022 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import numpy as np

from .. import constants as pc
from .. import io as io
from ..lib import _spline as sp


def read(pyrat):
    """
    Read a Cross-section (CS) file.
    """
    log = pyrat.log
    log.head('\nReading cross-section files.')
    cs = pyrat.cs = pyrat.cs.clone_new(pyrat)

    if cs.files is None:
        log.head('No CS files to read.', indent=2)
        return

    cs.nfiles = len(cs.files)

    for csfile in cs.files:
        log.head(f"Read CS file: '{csfile}'.", indent=2)
        absorption, molecules, temp, wn = io.read_cs(csfile)

        cs.absorption.append(absorption)
        cs.molecules.append(molecules)
        cs.temp.append(temp)
        cs.wavenumber.append(wn)

        ntemp = len(temp)
        nwave = len(wn)

        # Check that CS species are in the atmospheric file:
        absent = np.setdiff1d(molecules, pyrat.mol.name)
        if len(absent) > 0:
            log.error(
                f'These cross-section species {absent} are not listed in '
                'the atmospheric file\n'
            )

        # Update temperature boundaries:
        cs.tmin = np.amax((cs.tmin, temp[0]))
        cs.tmax = np.amin((cs.tmax, temp[-1]))

        # Wavenumber range check:
        if wn[0] > pyrat.spec.wn[0] or wn[-1] < pyrat.spec.wn[-1]:
            log.warning(
                f"The wavenumber range [{wn[0]:.3f}, {wn[-1]:.3f}] cm-1 "
                f"of the CS file:\n  '{csfile}',"
                "\ndoes not cover the Pyrat's wavenumber range: "
                f"[{pyrat.spec.wn[0]:.3f}, {pyrat.spec.wn[-1]:.3f}] cm-1."
            )

        # Screen output:
        molecs = '-'.join(molecules)
        log.msg(
            f'Cross-section opacity for {molecs}:\n'
            f'Read {nwave} wavenumber and {ntemp} temperature samples.',
            indent=4,
        )
        log.msg(
            f'Temperature sample limits: {temp[0]:g}--{temp[-1]:g} K',
            indent=4,
        )
        log.msg(
            f'Wavenumber sample limits: {wn[0]:.1f}--{wn[-1]:.1f} cm-1',
            indent=4,
        )

        # Wavenumber-interpolated CS:
        iabsorp = np.zeros((ntemp, pyrat.spec.nwave), np.double)
        for j in range(ntemp):
            z = sp.second_deriv(absorption[j], wn)
            iabsorp[j] = sp.splinterp_1D(
                absorption[j], wn, z, pyrat.spec.wn, 0.0,
            )
        cs.iabsorp.append(iabsorp)
        # Array with second derivatives:
        iz   = np.zeros((pyrat.spec.nwave, ntemp), np.double)
        wnlo = np.flatnonzero(pyrat.spec.wn >= wn[ 0])[ 0]
        wnhi = np.flatnonzero(pyrat.spec.wn <= wn[-1])[-1] + 1
        for j in range(wnlo, wnhi):
            iz[j] = sp.second_deriv(iabsorp[:,j], temp)
        cs.iz.append(iz.T)
        cs.iwnlo.append(wnlo)
        cs.iwnhi.append(wnhi)

    log.head('Cross-section read done.')


def interpolate(pyrat, layer=None):
    """
    Interpolate the CS absorption into the planetary model temperature.
    """
    pyrat.log.head('\nBegin CS interpolation.')

    # Allocate output extinction-coefficient array:
    if layer is None:   # Take a single layer
        ec = np.zeros((pyrat.atm.nlayers, pyrat.spec.nwave))
        li, lf = 0, pyrat.atm.nlayers
    else: # Take whole atmosphere
        ec = np.zeros((pyrat.cs.nfiles, pyrat.spec.nwave))
        li, lf = layer, layer+1
        label = []

    for i in range(pyrat.cs.nfiles):
        cs_absorption = np.zeros((lf-li, pyrat.spec.nwave))
        sp.splinterp_2D(
            pyrat.cs.iabsorp[i], pyrat.cs.temp[i], pyrat.cs.iz[i],
            pyrat.atm.temp[li:lf], cs_absorption,
            pyrat.cs.iwnlo[i], pyrat.cs.iwnhi[i],
        )

        # Get density scale factor in amagat:
        dens = 1.0
        for mol in pyrat.cs.molecules[i]:
            imol = np.where(pyrat.mol.name == mol)[0][0]
            dens *= pyrat.atm.d[li:lf,imol] / pc.amagat

        # Compute CS absorption in cm-1 units:
        if layer is None:
            ec += cs_absorption * np.expand_dims(dens, axis=1)
        else:
            ec[i] = cs_absorption * dens
            if len(pyrat.cs.molecules[i]) == 2:
                label.append('CIA ' + '-'.join(pyrat.cs.molecules[i]))
            else:
                label.append(pyrat.cs.molecules[i][0])

    # Return per-database EC if single-layer run:
    if layer is not None:
        return ec, label
    # Else, store cumulative result into pyrat object:
    pyrat.cs.ec = ec
    pyrat.log.head('Cross-section interpolate done.')
