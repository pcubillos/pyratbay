# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import multiprocessing as mp

import numpy as np
import scipy.interpolate as si

from .. import tools as pt
from .. import constants as pc
from .. import spectrum as ps
from .. import io as io


def check_spectrum(pyrat):
    """
    Check that user input arguments make sense.
    """
    # Shortcuts:
    log = pyrat.log
    phy = pyrat.phy
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

    missing_radius_ratio = (
        pyrat.od.rt_path in pc.emission_rt and
        (atm.rplanet is None or phy.rstar is None)
    )
    if missing_radius_ratio:
        log.error("Undefined radius ratio, need to define rplanet and rstar")

    # Accept ray-path argument:
    if pyrat.runmode in ['spectrum', 'retrieval'] and pyrat.od.rt_path is None:
        log.error(
            "Undefined radiative-transfer observing geometry (rt_path)."
            f"  Select from {pc.rt_paths}"
        )

    # Check system arguments:
    if pyrat.od.rt_path in pc.transmission_rt and phy.rstar is None:
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
    phy = pyrat.phy
    log = pyrat.log

    # Read stellar spectrum model: starspec, kurucz, or blackbody
    if phy.starspec is not None:
        starflux, starwn, star_temps = io.read_spectra(phy.starspec)
    elif phy.kurucz is not None:
        if phy.tstar is None:
            log.error(
                'Undefined stellar temperature (tstar), required for '
                'Kurucz model'
            )
        if phy.log_gstar is None:
            log.error(
                'Undefined stellar gravity (log_gstar) needed for Kurucz model'
            )
        starflux, starwn, kuruczt, kuruczg = ps.read_kurucz(
            phy.kurucz, phy.tstar, phy.log_gstar,
        )
        log.msg(
            f'Input stellar params: T={phy.tstar:7.1f} K, log(g)={phy.log_gstar:4.2f}\n'
            f'Best Kurucz match:    T={kuruczt:7.1f} K, log(g)={kuruczg:4.2f}'
        )
    elif phy.marcs is not None:
        pass
    elif phy.phoenix is not None:
        pass
    elif phy.tstar is not None:
        starwn = pyrat.spec.wn
        starflux = ps.bbflux(starwn, phy.tstar)
    else:
        starflux, starwn = None, None

    phy.starflux = starflux
    phy.starwn = starwn

    # Store interpolated stellar spectrum:
    if phy.starflux is not None:
        # 1D spectra
        if np.ndim(phy.starflux) == 1:
            sinterp = si.interp1d(phy.starwn, phy.starflux)
            pyrat.spec.starflux = sinterp(pyrat.spec.wn)
            # Band-integrate the stellar flux
            pyrat.obs.bandflux_star = np.array([
                band(pyrat.spec.starflux)
                for band in pyrat.obs.filters
            ])
        # 2D spectra
        else:
            sinterp = si.interp1d(phy.starwn, phy.starflux, axis=1)
            starflux = sinterp(pyrat.spec.wn)
            pyrat.spec.flux_interp = si.interp1d(star_temps, starflux, axis=0)
            pyrat.spec.starflux = pyrat.spec.flux_interp(phy.tstar)

    is_emission = pyrat.od.rt_path in pc.emission_rt
    if is_emission and starflux is None:
        log.error(
            'Undefined stellar flux model.  Set starspec, kurucz, or '
            'tstar (for a blackbody spectrum)'
        )

