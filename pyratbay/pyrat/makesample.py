# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import numpy as np

from .. import constants as pc
from .. import spectrum as ps
from .. import tools as pt


def make_wavenumber(pyrat):
    """
    Make the wavenumber sample from user inputs.
    """
    spec = pyrat.spec
    ex = pyrat.ex
    log = pyrat.log

    log.head('\nGenerating wavenumber array.')
    wl_units = pt.u(spec.wlunits)
    # Low wavenumber boundary:
    if pyrat.inputs.wnlow is None:
        if pyrat.inputs.wlhigh is None:
            log.error(
                'Low wavenumber boundary is undefined.  Either set '
                'wnlow or wlhigh'
            )
        else:
            spec.wnlow = 1.0 / spec.wlhigh
    elif pyrat.inputs.wlhigh is not None:
        log.warning(
            f'Both wnlow ({spec.wnlow:.2e} cm-1) and wlhigh ({spec.wlhigh:.2e} '
             'cm) were defined.  Pyrat will take wnlow and ignore wlhigh'
        )

    # High wavenumber boundary:
    if pyrat.inputs.wnhigh is None:
        if pyrat.inputs.wllow is None:
            log.error(
                'High wavenumber boundary is undefined.  Either set '
                'wnhigh or wllow'
            )
        else:
            spec.wnhigh = 1.0 / spec.wllow
    elif pyrat.inputs.wllow is not None:
        log.warning(
            f'Both wnhigh ({spec.wnhigh:.2e} cm-1) and wllow ({spec.wllow:.2e}'
             ' cm) were defined.  Pyrat will take wnhigh and ignore wllow'
        )

    # Consistency check (wnlow < wnhigh):
    if spec.wnlow > spec.wnhigh:
        log.error(
            f'Wavenumber low boundary ({spec.wnlow:.1f} cm-1) must be '
            f'larger than the high boundary ({spec.wnhigh:.1f} cm-1)'
        )

    if spec.wnstep is None:
        log.error('Undefined wavenumber sampling step size (wnstep)')

    if spec.wnosamp is None:
        log.error('Undefined wavenumber oversampling factor (wnosamp)')

    # Set wavelength limits based on the wavenumber limits:
    spec.wlhigh = 1.0 / spec.wnlow
    spec.wllow  = 1.0 / spec.wnhigh

    if spec.resolution is not None:
        # Constant-resolving power wavenumber sampling:
        print(spec.wnlow, spec.wnhigh, spec.resolution)
        spec.wn = ps.constant_resolution_spectrum(
            spec.wnlow, spec.wnhigh, spec.resolution,
        )
        spec.nwave = len(spec.wn)
        spec.wlstep = None
    elif spec.wlstep is not None:
        # Constant-sampling rate wavelength sampling:
        wl = np.arange(spec.wllow, spec.wlhigh, spec.wlstep)
        spec.nwave = len(wl)
        spec.wn = 1.0/np.flip(wl)
        spec.wnlow = spec.wn[0]
        spec.resolution = None
    else:
        # Constant-sampling rate wavenumber sampling:
        spec.nwave = int((spec.wnhigh-spec.wnlow)/spec.wnstep) + 1
        spec.wn = spec.wnlow + np.arange(spec.nwave) * spec.wnstep

    # Fine-sampled wavenumber array (constant sampling rate):
    spec.ownstep = spec.wnstep / spec.wnosamp
    spec.onwave = int(np.ceil((spec.wn[-1]-spec.wnlow)/spec.ownstep)) + 1
    spec.own = spec.wnlow + np.arange(spec.onwave) * spec.ownstep
    if spec.spectrum is None:
        spec.spectrum = np.zeros(spec.nwave, np.double)

    # Get list of divisors:
    spec.odivisors = pt.divisors(spec.wnosamp)

    # Re-set final boundary (stay inside given boundaries):
    if spec.wn[-1] != spec.wnhigh:
        log.warning(
            f'Final wavenumber modified from {spec.wnhigh:.4f} cm-1 (input)\n'
            f'                            to {spec.wn[-1]:.4f} cm-1 (Pyrat)'
        )

    # Screen output:
    log.msg(
        f'Initial wavenumber boundary:  {spec.wnlow:.5e} cm-1  '
        f'({spec.wlhigh/wl_units:.3e} {spec.wlunits})\n'
        f'Final   wavenumber boundary:  {spec.wnhigh:.5e} cm-1  '
        f'({spec.wllow/wl_units:.3e} {spec.wlunits})',
        indent=2,
    )

    if spec.resolution is not None:
        msg = f'Spectral resolving power: {spec.resolution:.1f}'
    elif spec.wlstep is not None:
        wl_step = spec.wlstep / wl_units
        msg = f'Wavelength sampling interval: {wl_step:.2g} {spec.wlunits}'
    else:
        msg = f'Wavenumber sampling interval: {spec.wnstep:.2g} cm-1'
    log.msg(
        f'{msg}\n'
        f'Wavenumber sample size:      {spec.nwave:8d}\n'
        f'Wavenumber fine-sample size: {spec.onwave:8d}\n',
        indent=2,
    )
    log.head('Wavenumber sampling done.')


    # If there are opacity tables:
    modify_from_opacity = (
        ex.wn is not None and
        (ex.nwave != spec.nwave or np.sum(np.abs(ex.wn-spec.wn)) > 0)
    )
    if modify_from_opacity:
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

        log.warning(
            "Wavenumber sampling from extinction-coefficient "
            "table does not match the input wavenumber sampling.  Adopting "
            f"tabulated array with {ex.nwave} samples, {sampling_text}, and "
            f"ranges [{spec.wnlow:.2f}, {spec.wnhigh:.2f}] cm-1."
        )


def make_atmprofiles(pyrat):
    """
    ** NOT IN USE **
    Define atmospheric-profile layers sampling.

    Notes
    -----
    There are multiple things going on in this functions, here's a summary:
    - Reset (pressure/radius) boundaries if requested.
    - Resample (radius) profiles if requested.
    - Check whether top of atmosphere crosses the Hill radius.
    - Compute partition-function at layers temperatures.
    """
    log = pyrat.log

    # Pyrat and user-input atmospheric-data objects:
    atm = pyrat.atm
    atm_in = pyrat.inputs.atm

    punits = pt.u(atm.punits)
    runits = pt.u(atm.runits)
    if pyrat.atm.rmodelname is not None and not np.any(missing):
        if not np.isinf(pyrat.phy.rhill):
            log.msg(
                f'Hill radius: {pyrat.phy.rhill/runits:8.1f} {atm.runits}.',
                indent=2,
            )
        atm.radius = pyrat.hydro(
            atm_in.press, atm_in.temp, atm_in.mm,
            pyrat.phy.gplanet, pyrat.phy.mplanet,
            atm.refpressure, pyrat.phy.rplanet,
        )

    # Check if Hydrostatic Eq. breaks down:
    if atm_in.radius is not None and np.any(np.isinf(atm_in.radius)):
        ibreak = np.where(np.isfinite(atm_in.radius))[0][0]

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
           f'Top-pressure boundary (ptop={atm.ptop/punits:.2e} '
           f'{atm.punits}) lies outside of the atmospheric file range '
           f'{pmin/punits:.2e}--{pmax/punits:.2e} {atm.punits}'
        )
    if atm.pbottom < pmin or atm.pbottom > pmax:
        log.error(
           f'Bottom-pressure boundary (pbottom={atm.pbottom/punits:.2e} '
           f'{atm.punits}) lies outside of the atmospheric file range '
           f'{pmin/punits:.2e}--{pmax/punits:.2e} {atm.punits}'
        )
    if atm.pbottom <= atm.ptop:
        log.error(
           f'Bottom-layer pressure ({atm.pbottom/punits:.2e} {atm.punits}) '
            'must be higher than the top-layer pressure '
           f'({atm.ptop/punits:.2e} {atm.punits})'
        )

    log.msg(
        'User pressure boundaries: '
        f'{atm.ptop/pc.bar:.2e}--{atm.pbottom/pc.bar:.2e} bar.',
        indent=2,
    )

    # Resample to equispaced radius array if requested:
    if atm.radstep is not None:
        # Make equispaced radius array:
        if atm.radlow is None:
            atm.radlow  = radinterp(atm.pbottom)[0]
        if atm.radhigh is None:
            atm.radhigh = radinterp(atm.ptop)[0]

        atm.radius = np.arange(atm.radlow, atm.radhigh, atm.radstep)
        atm.nlayers = len(atm.radius)
        # Interpolate to pressure array:
        atm.press = pressinterp(atm.radius)

    # Check that the radii lie within Hill radius:
    if atm.radius is not None:
        atm.rtop = pt.ifirst(atm.radius < pyrat.phy.rhill, default_ret=0)
    if atm.rtop > 0:
        rhill = pyrat.phy.rhill/runits
        log.warning(
            'The atmospheric pressure array extends beyond the Hill radius '
           f'({rhill:.5f} {atm.runits}) at pressure '
           f'{atm_in.press[atm.rtop]/pc.bar:.2e} bar (layer {atm.rtop}).  '
            'Extinction beyond this layer will be neglected.')


    log.head('Make atmosphere done.')
