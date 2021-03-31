# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import numpy as np
import scipy.interpolate as sip

from .. import tools as pt
from .. import constants as pc


def make_wavenumber(pyrat):
    """
    Make the wavenumber sample from user inputs.
    """
    # Alias for Pyrat's Spectrum object:
    spec = pyrat.spec
    log  = pyrat.log

    log.head('\nGenerating wavenumber array.')
    # Low wavenumber boundary:
    if pyrat.inputs.wnlow is None:
        if pyrat.inputs.wlhigh is None:
            log.error('Low wavenumber boundary is undefined.  Either set '
                      'wnlow or wlhigh.')
        else:
            spec.wnlow = 1.0 / spec.wlhigh
    elif pyrat.inputs.wlhigh is not None:
        log.warning(
            f'Both wnlow ({spec.wnlow:.2e} cm-1) and wlhigh ({spec.wlhigh:.2e} '
             'cm) were defined.  Pyrat will take wnlow and ignore wlhigh.')

    # High wavenumber boundary:
    if pyrat.inputs.wnhigh is None:
        if pyrat.inputs.wllow is None:
            log.error('High wavenumber boundary is undefined.  Either set '
                      'wnhigh or wllow.')
        else:
            spec.wnhigh = 1.0 / spec.wllow
    elif pyrat.inputs.wllow is not None:
        log.warning(
            f'Both wnhigh ({spec.wnhigh:.2e} cm-1) and wllow ({spec.wllow:.2e}'
             ' cm) were defined.  Pyrat will take wnhigh and ignore wllow.')

    # Consistency check (wnlow < wnhigh):
    if spec.wnlow > spec.wnhigh:
      log.error(f'Wavenumber low boundary ({spec.wnlow:.1f} cm-1) must be '
                f'larger than the high boundary ({spec.wnhigh:.1f} cm-1).')

    if spec.wnstep is None:
        log.error('Undefined wavenumber sampling step size (wnstep).')

    if spec.wnosamp is None:
        log.error('Undefined wavenumber oversampling factor (wnosamp).')

    # Set wavelength limits based on the wavenumber limits:
    spec.wlhigh = 1.0 / spec.wnlow
    spec.wllow  = 1.0 / spec.wnhigh

    if spec.resolution is not None:
        # Constant-resolving power wavenumber sampling:
        f = 0.5 / spec.resolution
        g = (1.0+f) / (1.0-f)
        spec.nwave = int(np.ceil(-np.log(spec.wnlow/spec.wnhigh) / np.log(g)))
        spec.wn = spec.wnlow * g**np.arange(spec.nwave)
    else:
        # Constant-sampling rate wavenumber sampling:
        spec.nwave = int((spec.wnhigh-spec.wnlow)/spec.wnstep) + 1
        spec.wn = spec.wnlow + np.arange(spec.nwave) * spec.wnstep

    # Fine-sampled wavenumber array (constant sampling rate):
    spec.ownstep = spec.wnstep / spec.wnosamp
    spec.onwave = int(np.ceil((spec.wn[-1]-spec.wnlow)/spec.ownstep)) + 1
    spec.own = spec.wnlow + np.arange(spec.onwave) * spec.ownstep

    # Get list of divisors:
    spec.odivisors = pt.divisors(spec.wnosamp)

    # Re-set final boundary (stay inside given boundaries):
    if spec.wn[-1] != spec.wnhigh:
        log.warning(
            f'Final wavenumber modified from {spec.wnhigh:.4f} cm-1 (input)\n'
            f'                            to {spec.wn[-1]:.4f} cm-1 (Pyrat).')

    # Screen output:
    log.msg(f'Initial wavenumber boundary:  {spec.wnlow:.5e} cm-1  '
            f'({spec.wlhigh/pt.u(spec.wlunits):.3e} {spec.wlunits})', indent=2)
    log.msg(f'Final   wavenumber boundary:  {spec.wnhigh:.5e} cm-1  '
            f'({spec.wllow/pt.u(spec.wlunits):.3e} {spec.wlunits})', indent=2)
    if spec.resolution is None:
        log.msg(f'Wavenumber sampling interval: {spec.wnstep:.2g} cm-1',
            indent=2)
    else:
        log.msg(f'Spectral resolving power: {spec.resolution:.1f}', indent=2)
    log.msg(f'Wavenumber sample size:      {spec.nwave:8d}\n'
            f'Wavenumber fine-sample size: {spec.onwave:8d}\n', indent=2)
    log.head('Wavenumber sampling done.')


def make_atmprofiles(pyrat):
    """
    Define atmospheric-profile layers sampling.

    Notes
    -----
    There are multiple things going on in this functions, here's a summary:
    - Compute hydrostatic-equilibrium radius profile if necessary.
    - Reset (pressure/radius) boundaries if requested.
    - Resample (pressure/radius/abundance/temperature) profiles if requested.
    - Check whether top of atmosphere crosses the Hill radius.
    - Compute partition-function at layers temperatures.
    """
    log = pyrat.log
    log.head('\nGenerating atmospheric profile sample.')

    # Pyrat and user-input atmospheric-data objects:
    atm    = pyrat.atm
    atm_in = pyrat.inputs.atm

    # Check that the layers are sorted from the top to the bottom of
    #  the atmosphere:
    sort    = np.all(np.ediff1d(atm_in.press) > 0)  # Top to bottom
    reverse = np.all(np.ediff1d(atm_in.press) < 0)  # Bottom to top
    if atm_in.radius is not None:
        sort    *= np.all(np.ediff1d(atm_in.radius) < 0)
        reverse *= np.all(np.ediff1d(atm_in.radius) > 0)

    if sort:       # Layers are in the correct order
        pass
    elif reverse:  # Layers in reverse order
        log.warning('The atmospheric layers are in reversed order '
                    '(bottom-top).  Resorting to be from the top down.')
        if atm_in.radius is not None:
            atm_in.radius = np.flipud(atm_in.radius)
        atm_in.press = np.flipud(atm_in.press)
        atm_in.temp  = np.flipud(atm_in.temp)
        atm_in.mm    = np.flipud(atm_in.mm)
        atm_in.q     = np.flipud(atm_in.q)
        atm_in.d     = np.flipud(atm_in.d)
    else:
        log.error('The atmospheric layers are neither sorted from the '
                  'bottom up, nor from the top down.')

    missing = [
        pyrat.phy.mplanet is None,
        pyrat.phy.gplanet is None,
        pyrat.phy.rplanet is None,
        atm.refpressure is None,
        ]

    if pyrat.atm.rmodelname is None and atm_in.radius is None \
       and pyrat.runmode != "opacity":
        log.error('Cannot compute radius profile.  Need to set a radius '
                  'model (radmodel) or provide an input radius array in '
                  'the atmospheric file.')

    if pyrat.atm.rmodelname is not None and pyrat.runmode != "opacity":
        if np.any(missing[0:3]):
            log.error('Cannot compute hydrostatic-equilibrium radius profile.  '
                'Must define at least two of mplanet, rplanet, or gplanet.')
        if missing[3]:
            log.error('Cannot compute hydrostatic-equilibrium radius profile.  '
                'Undefined reference pressure level (refpressure).')

    punits = pt.u(atm.punits)
    runits = pt.u(atm.runits)
    if pyrat.atm.rmodelname is not None and not np.any(missing):
        # Calculate hydostatic-equilibrium radius profile:
        log.msg(
            f'Reference pressure: {atm.refpressure/punits:.3e} {atm.punits}.',
            indent=2)
        log.msg(
            f'Reference radius: {pyrat.phy.rplanet/runits:8.1f} {atm.runits}.',
            indent=2)
        if not np.isinf(pyrat.phy.rhill):
            log.msg(f'Hill radius:      {pyrat.phy.rhill/runits:8.1f} '
                    f'{atm.runits}.', indent=2)
        atm_in.radius = pyrat.hydro(
            atm_in.press, atm_in.temp, atm_in.mm,
            pyrat.phy.gplanet, pyrat.phy.mplanet,
            atm.refpressure, pyrat.phy.rplanet)

    # Check if Hydrostatic Eq. breaks down:
    ibreak = 0  # Index and flag at the same time
    if atm_in.radius is not None and np.any(np.isinf(atm_in.radius)):
        ibreak = np.where(np.isfinite(atm_in.radius))[0][0]

    # Set the interpolating function (for use later):
    try:
        radinterp = sip.interp1d(
            atm_in.press[ibreak:], atm_in.radius[ibreak:], kind='slinear')
        pressinterp = sip.interp1d(
            np.flipud(atm_in.radius[ibreak:]),
            np.flipud(atm_in.press [ibreak:]), kind='slinear')
    except:
        pass

    # Set radius/pressure boundaries if exist:
    if atm.ptop is not None:
        if atm.ptop >= atm_in.press[ibreak]:
            ibreak = 0   # Turn-off break flag
    elif atm.radhigh is not None:
        if atm.radhigh <= atm_in.radius[ibreak]:
            ibreak = 0   # Turn-off break flag
            atm.ptop = pressinterp(atm.radhigh)[0]
        #else:
        #  out-of-bounds error
    else:
        atm.ptop = np.amin(atm_in.press)

    if atm.radlow is not None:
        atm.pbottom = pressinterp(atm.radlow)[0]
    elif atm.pbottom is None:
        atm.pbottom = np.amax(atm_in.press)

    if ibreak != 0 and np.isinf(pyrat.phy.rhill):
        log.error(
            'Unbounded atmosphere.  Hydrostatic-equilibrium radius solution '
           f'diverges at pressure {atm_in.press[ibreak]/pc.bar:.3e} bar.  '
            'Set mstar and smaxis to define a Hill radius (top boundary) '
            'and avoid error.')

    # Out of bounds errors:
    pmin = np.amin(atm_in.press)
    pmax = np.amax(atm_in.press)
    if atm.ptop < pmin or atm.ptop > pmax:
        log.error(
           f'Top-pressure boundary (ptop={atm.ptop/punits:.2e} {atm.punits}) '
            'lies outside of the atmospheric-file range '
           f'{pmin/punits:.2e}--{pmax/punits:.2e} {atm.punits}.')
    if atm.pbottom < pmin or atm.pbottom > pmax:
        log.error(
           f'Bottom-pressure boundary (pbottom={atm.pbottom/punits:.2e} '
           f'{atm.punits}) lies outside of the atmospheric-file range '
           f'{pmin/punits:.2e}--{pmax/punits:.2e} {atm.punits}.')
    if atm.pbottom <= atm.ptop:
        log.error(
           f'Bottom-layer pressure ({atm.pbottom/punits:.2e} {atm.punits}) '
            'must be higher than the top-layer pressure '
           f'({atm.ptop/punits:.2e} {atm.punits}).')

    log.msg('User pressure boundaries: '
           f'{atm.ptop/pc.bar:.2e}--{atm.pbottom/pc.bar:.2e} bar.', indent=2)

    # Resample to equispaced log-pressure array if requested:
    if pyrat.inputs.nlayers is not None:
        atm.press = np.logspace(
            np.log10(atm.ptop), np.log10(atm.pbottom), atm.nlayers)
        atm.radius = radinterp(atm.press)
        resample = True

    # Resample to equispaced radius array if requested:
    elif atm.radstep is not None:
        # Make equispaced radius array:
        if atm.radlow is None:
            atm.radlow  = radinterp(atm.pbottom)[0]
        if atm.radhigh is None:
            atm.radhigh = radinterp(atm.ptop)[0]

        atm.radius = np.arange(atm.radlow, atm.radhigh, atm.radstep)
        atm.nlayers = len(atm.radius)
        # Interpolate to pressure array:
        atm.press = pressinterp(atm.radius)
        resample = True

    else:  # Take the atmospheric-file sampling:
        # Get top-bottom indices:
        ilow  = np.where(atm_in.press >= atm.ptop)   [0][ 0]
        ihigh = np.where(atm_in.press <= atm.pbottom)[0][-1]
        # Take values within the boundaries:
        atm.press = atm_in.press[ilow:ihigh+1]
        atm.temp  = atm_in.temp [ilow:ihigh+1]
        atm.mm    = atm_in.mm   [ilow:ihigh+1]
        atm.q     = atm_in.q    [ilow:ihigh+1]
        atm.d     = atm_in.d    [ilow:ihigh+1]
        if atm_in.radius is not None:
            atm.radius = atm_in.radius[ilow:ihigh+1]
        atm.nlayers = len(atm.press)
        resample = False

    # Check the radii that lie within Hill radius:
    if atm.radius is not None:
        atm.rtop = pt.ifirst(atm.radius < pyrat.phy.rhill, default_ret=0)
    if atm.rtop > 0:
        rhill = pyrat.phy.rhill/runits
        log.warning(
            'The atmospheric pressure array extends beyond the Hill radius '
           f'({rhill:.5f} {atm.runits}) at pressure '
           f'{atm_in.press[atm.rtop]/pc.bar:.2e} bar (layer {atm.rtop}).  '
            'Extinction beyond this layer will be neglected.')

    # Print radius array:
    if atm.radius is not None:
        radstr = ', '.join(f'{k:.3f}' for k in atm.radius/runits)
        log.msg(f'Radius array ({atm.runits}) =   [{radstr}]', indent=2, si=4)
        log.msg(
            'Valid upper/lower radius boundaries:    '
           f'{atm.radius[atm.rtop]/runits:.5f} - {atm.radius[-1]/runits:.5f} '
           f'{atm.runits}.', indent=2)

    log.msg(
        'Valid lower/higher pressure boundaries: '
       f'{atm.press[atm.rtop]/punits:.2e} - {atm.pbottom/punits:.2e} '
       f'{atm.punits}.', indent=2)
    log.msg(f'Number of valid model layers: {atm.nlayers-atm.rtop}.', indent=2)

    # Interpolate to new atm-layer sampling if necessary:
    if resample:
        tempinterp = sip.interp1d(atm_in.press, atm_in.temp, kind='slinear')
        mminterp   = sip.interp1d(atm_in.press, atm_in.mm,   kind='slinear')
        atm.temp = tempinterp(atm.press)
        atm.m    =   mminterp(atm.press)
        # Interpolate abundance profiles:
        atm.q = np.zeros((atm.nlayers, pyrat.mol.nmol))
        atm.d = np.zeros((atm.nlayers, pyrat.mol.nmol))
        for i in np.arange(pyrat.mol.nmol):
            qinterp = sip.interp1d(atm_in.press, atm_in.q[:, i], kind='slinear')
            dinterp = sip.interp1d(atm_in.press, atm_in.d[:, i], kind='slinear')
            atm.q[:,i] = qinterp(atm.press)
            atm.d[:,i] = dinterp(atm.press)

    # Interpolate isotopes partition function:
    log.msg(f'Number of isotopes: {pyrat.iso.niso}', indent=2)
    # Initialize the partition-function array for pyrat.iso:
    pyrat.iso.z = np.zeros((pyrat.iso.niso, atm.nlayers))
    for db in pyrat.lt.db:            # For each Database
        for j in np.arange(db.niso):  # For each isotope in DB
            log.debug(
                f'Interpolating (isotope ID {db.iiso+j}) partition function.',
                indent=4)
            zinterp = sip.interp1d(db.temp, db.z[j], kind='slinear')
            pyrat.iso.z[db.iiso+j] = zinterp(atm.temp)

    # Base abundance profiles:
    atm.qbase = np.copy(atm.q)

    log.head('Make atmosphere done.')
