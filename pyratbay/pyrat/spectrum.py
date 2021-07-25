# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import numpy as np
from scipy.interpolate import interp1d
import scipy.special as ss

from .. import constants as pc
from .. import io as io
from .. import spectrum as ps
from ..lib import _trapz as t


def spectrum(pyrat):
    """
    Spectrum calculation driver.
    """
    pyrat.log.head('\nCalculate the planetary spectrum.')

    # Initialize the spectrum array:
    pyrat.spec.spectrum = np.empty(pyrat.spec.nwave, np.double)
    if pyrat.cloud.fpatchy is not None:
        pyrat.spec.clear  = np.empty(pyrat.spec.nwave, np.double)
        pyrat.spec.cloudy = np.empty(pyrat.spec.nwave, np.double)

    # Call respective function depending on the RT/geometry:
    if pyrat.od.rt_path in pc.transmission_rt:
        modulation(pyrat)

    elif pyrat.od.rt_path == 'emission':
        intensity(pyrat)
        flux(pyrat)

    elif pyrat.od.rt_path == 'emission_two_stream':
        two_stream(pyrat)

    # Print spectrum to file:
    if pyrat.od.rt_path in pc.transmission_rt:
        spec_type = 'transit'
    elif pyrat.od.rt_path in pc.emission_rt:
        spec_type = 'emission'

    io.write_spectrum(
        1.0/pyrat.spec.wn, pyrat.spec.spectrum, pyrat.spec.specfile, spec_type)
    pyrat.log.head('Done.')


def modulation(pyrat):
    """Calculate modulation spectrum for transit geometry."""
    rtop = pyrat.atm.rtop
    radius = pyrat.atm.radius
    depth = pyrat.od.depth

    # Get Delta radius (and simps' integration variables):
    h = np.ediff1d(radius[rtop:])
    # The integrand:
    integ = (np.exp(-depth[rtop:,:]) * np.expand_dims(radius[rtop:],1))

    if pyrat.cloud.fpatchy is not None:
        h_clear = np.copy(h)
        integ_clear = (
            np.exp(-pyrat.od.depth_clear[rtop:,:]) *
                  np.expand_dims(radius[rtop:],1))

    if 'deck' in (m.name for m in pyrat.cloud.models):
        # Replace (interpolating) last layer with cloud top:
        deck = pyrat.cloud.models[pyrat.cloud.model_names.index('deck')]
        if deck.itop > rtop:
            h[deck.itop-rtop-1] = deck.rsurf - radius[deck.itop-1]
            integ[deck.itop-rtop] = interp1d(
                radius[rtop:], integ, axis=0)(deck.rsurf)

    # Number of layers for integration at each wavelength:
    nlayers = pyrat.od.ideep - rtop + 1
    spectrum = t.trapz2D(integ, h, nlayers-1)
    pyrat.spec.spectrum = (radius[rtop]**2 + 2*spectrum) / pyrat.phy.rstar**2

    if pyrat.cloud.fpatchy is not None:
        nlayers = pyrat.od.ideep_clear - rtop + 1
        pyrat.spec.clear = t.trapz2D(integ_clear, h_clear, nlayers-1)

        pyrat.spec.clear = ((radius[rtop]**2 + 2*pyrat.spec.clear)
                             / pyrat.phy.rstar**2)
        pyrat.spec.cloudy = pyrat.spec.spectrum
        pyrat.spec.spectrum = (   pyrat.cloud.fpatchy  * pyrat.spec.cloudy +
                               (1-pyrat.cloud.fpatchy) * pyrat.spec.clear  )

    if pyrat.spec.specfile is not None:
        specfile = f": '{pyrat.spec.specfile}'"
    else:
        specfile = ""
    pyrat.log.head(f"Computed transmission spectrum{specfile}.", indent=2)


def intensity(pyrat):
    """
    Calculate the intensity spectrum [units] for eclipse geometry.
    """
    spec = pyrat.spec
    pyrat.log.msg('Computing intensity spectrum.', indent=2)
    if spec.quadrature is not None:
        spec.raygrid = np.arccos(np.sqrt(spec.qnodes))

    # Allocate intensity array:
    spec.nangles = len(spec.raygrid)
    spec.intensity = np.empty((spec.nangles, spec.nwave), np.double)

    # Calculate the Planck Emission:
    pyrat.od.B = np.zeros((pyrat.atm.nlayers, spec.nwave), np.double)
    ps.blackbody_wn_2D(spec.wn, pyrat.atm.temp, pyrat.od.B, pyrat.od.ideep)

    if 'deck' in (m.name for m in pyrat.cloud.models):
        deck = pyrat.cloud.models[pyrat.cloud.model_names.index('deck')]
        pyrat.od.B[deck.itop] = ps.blackbody_wn(pyrat.spec.wn, deck.tsurf)

    # Plane-parallel radiative-transfer intensity integration:
    spec.intensity = t.intensity(
        pyrat.od.depth, pyrat.od.ideep, pyrat.od.B, np.cos(spec.raygrid),
        pyrat.atm.rtop)


def flux(pyrat):
    """
    Calculate the hemisphere-integrated flux spectrum [units] for eclipse
    geometry.
    """
    spec = pyrat.spec
    # Calculate the projected area:
    boundaries = np.linspace(0, 0.5*np.pi, spec.nangles+1)
    boundaries[1:spec.nangles] = 0.5 * (spec.raygrid[:-1] + spec.raygrid[1:])
    area = np.pi * (np.sin(boundaries[1:])**2 - np.sin(boundaries[:-1])**2)

    if spec.quadrature is not None:
        area = spec.qweights * np.pi
    # Weight-sum the intensities to get the flux:
    spec.spectrum[:] = np.sum(spec.intensity * np.expand_dims(area,1), axis=0)
    pyrat.log.head(f"Computed emission spectrum: '{spec.specfile}'.", indent=2)


def two_stream(pyrat):
    """
    Two-stream approximation radiative transfer.

    This function defines downward (flux_down) and uppward fluxes
    (flux_up) into pyrat.spec, and sets the emission spectrum as the
    uppward flux at the top of the atmosphere (flux_up[0]):

    flux_up: 2D float ndarray
        Upward flux spectrum through each layer under the two-stream
        approximation (erg s-1 cm-2 cm).
    flux_down: 2D float ndarray
        Downward flux spectrum through each layer under the two-stream
        approximation (erg s-1 cm-2 cm).
    """
    pyrat.log.msg('Computing two-stream flux spectrum.', indent=2)
    spec = pyrat.spec
    phy = pyrat.phy
    nlayers = pyrat.atm.nlayers

    # Set internal net bolometric flux to sigma*Tint**4:
    spec.f_int = ps.blackbody_wn(spec.wn, phy.tint)
    total_f_int = np.trapz(spec.f_int, spec.wn)
    if total_f_int > 0:
        spec.f_int *= pc.sigma * phy.tint**4 / total_f_int

    # Diffusivity factor (Eq. B5 of Heng et al. 2014):
    dtau0 = np.diff(pyrat.od.depth, n=1, axis=0)
    trans = (1-dtau0)*np.exp(-dtau0) + dtau0**2 * ss.exp1(dtau0)

    B = pyrat.od.B = ps.blackbody_wn_2D(spec.wn, pyrat.atm.temp)
    Bp = np.diff(pyrat.od.B, n=1, axis=0) / dtau0

    # Diffuse approximation to compute downward and upward fluxes:
    spec.flux_down = np.zeros((nlayers, spec.nwave))
    spec.flux_up = np.zeros((nlayers, spec.nwave))

    # Top boundary condition:
    spec.flux_down[0] = \
        phy.beta_irr * (0.5*phy.rstar/phy.smaxis)**2 * spec.starflux
    # (Corrected) Eqs. 10 of Malik et al. (2017),
    # setting the diffusivity factor to D=2 (epsilon=0.5):
    for i in range(nlayers-1):
        spec.flux_down[i+1] = \
            trans[i] * spec.flux_down[i] \
            + np.pi*(B[i+1] - trans[i]*B[i] - 0.5*Bp[i]*(1-trans[i]))

    spec.flux_up[nlayers-1] = spec.flux_down[nlayers-1] + spec.f_int
    for i in reversed(range(nlayers-1)):
        spec.flux_up[i] = \
            trans[i] * spec.flux_up[i+1] \
            + np.pi*(B[i] - trans[i]*B[i+1] + 0.5*Bp[i]*(1-trans[i]))

    spec.spectrum = spec.flux_up[0]
