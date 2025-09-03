# Copyright (c) 2021-2025 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import multiprocessing as mp

import numpy as np
import scipy.interpolate as si

from .. import constants as pc


def check_spectrum(pyrat):
    """
    Check that user input arguments make sense.
    """
    # Shortcuts:
    log = pyrat.log
    spec = pyrat.spec
    atm = pyrat.atm

    if pyrat.runmode == 'spectrum' and spec.specfile is None:
        log.error('Undefined output spectrum file (specfile).')

    if atm.vmr is None and atm.chemistry is None:
        log.error(
            'Missing atmospheric volume mixing ratios. Need to either read '
            'an input profile or compute one via a chemistry model (chemistry)'
        )

    # Check that the radius profile exists or can be computed:
    if atm.radius is None and pyrat.runmode != 'opacity':
        log.error(
            'Missing atmospheric radius profile.  Need to either read an '
            'input profile or compute one via the radmodel argument'
        )

    # Not needed for f_lambda
    missing_radius_ratio = (
        pyrat.od.rt_path in pc.eclipse_rt and
        (atm.rplanet is None or atm.rstar is None)
    )
    if missing_radius_ratio:
        msg = "Undefined radius ratio, need to define both rplanet and rstar"
        log.error(msg)

    # Accept ray-path argument:
    if pyrat.runmode in ['spectrum', 'retrieval'] and pyrat.od.rt_path is None:
        log.error(
            "Undefined radiative-transfer observing geometry (rt_path)."
            f"  Select from {pc.rt_paths}"
        )

    # Check system arguments:
    if pyrat.od.rt_path in pc.transmission_rt and atm.rstar is None:
        log.error(
            'Undefined stellar radius (rstar), required for transmission '
            'calculation'
        )

    if pyrat.ncpu >= mp.cpu_count():
        n_cpu = mp.cpu_count()
        log.warning(
            f'Number of requested CPUs ({pyrat.ncpu}) is >= than the number '
            f'of available CPUs ({n_cpu}).  Enforced ncpu to {n_cpu-1}'
        )
        pyrat.ncpu = n_cpu - 1
    log.head('Check spectrum done.')


def setup(pyrat):
    """
    Process stellar spectrum.
    Process the oberving filter bands.
    """
    atm = pyrat.atm
    log = pyrat.log

    # Not needed for f_lambda
    if pyrat.od.rt_path == 'f_lambda':
        return

    # Store interpolated stellar spectrum:
    if atm.starflux is not None:
        # 1D spectra
        if np.ndim(atm.starflux) == 1:
            sinterp = si.interp1d(atm.starwn, atm.starflux)
            pyrat.spec.starflux = sinterp(pyrat.spec.wn)
            # Band-integrate the stellar flux
            pyrat.obs.bandflux_star = np.array([
                band(pyrat.spec.starflux)
                for band in pyrat.obs.filters
            ])
        # 2D spectra
        else:
            sinterp = si.interp1d(atm.starwn, atm.starflux, axis=1)
            starflux = sinterp(pyrat.spec.wn)
            pyrat.spec.flux_interp = si.interp1d(atm.sed_temps, starflux, axis=0)
            pyrat.spec.starflux = pyrat.spec.flux_interp(atm.tstar)

    is_eclipse = pyrat.od.rt_path in pc.eclipse_rt
    if is_eclipse and pyrat.spec.starflux is None:
        log.error(
            'Undefined stellar flux model, required for eclipse calculations.  '
            'Set starspec, kurucz, or tstar (for a blackbody spectrum)'
        )

