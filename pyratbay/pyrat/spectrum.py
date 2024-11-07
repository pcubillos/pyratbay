# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Spectrum',
    'spectrum',
    'modulation',
    'intensity',
    'flux',
    'two_stream',
]

import numpy as np
import scipy.constants as sc
import scipy.special as ss
from scipy.interpolate import interp1d

from .. import constants as pc
from .. import io
from .. import spectrum as ps
from .. import tools as pt
from ..lib import _trapz as t


class GetWavelength:
    def __get__(self, obj, objtype=None):
        return 1.0 / (obj.wn * pc.um)


class Spectrum():
    wl = GetWavelength()
    def __init__(self, inputs, log):
        """
        Make the wavenumber sample from user inputs.
        """
        log.head('\nGenerating wavenumber array.')

        self.specfile = inputs.specfile # Transmission/Emission spectrum file
        self.intensity = None  # Intensity spectrum array
        self.clear     = None  # Clear modulation spectrum for patchy model
        self.cloudy    = None  # Cloudy modulation spectrum for patchy model
        self.starflux  = None  # Stellar flux spectrum
        self.wnosamp = None

        # Gaussian-quadrature flux integration over hemisphere (for emission)
        if inputs.quadrature is not None:
            self.quadrature = inputs.quadrature
            qnodes, qweights = ss.p_roots(self.quadrature)
            qnodes = 0.5*(qnodes + 1.0)
            self.raygrid = np.arccos(np.sqrt(qnodes))
            self.quadrature_mu = np.sqrt(qnodes)
            weights = 0.5 * np.pi * qweights
        else:
            # Custom defined mu-angles:
            raygrid = inputs.raygrid * sc.degree
            if raygrid[0] != 0:
                log.error('First angle in raygrid must be 0.0 (normal)')
            if np.any(raygrid < 0) or np.any(raygrid > 0.5*np.pi):
                log.error('raygrid angles must lie between 0 and 90 deg')
            if np.any(np.ediff1d(raygrid) <= 0):
                log.error('raygrid angles must be monotonically increasing')
            self.quadrature_mu = np.cos(raygrid)
            # Weights are the projected area for each mu range:
            bounds = np.linspace(0, 0.5*np.pi, len(raygrid)+1)
            bounds[1:-1] = 0.5 * (raygrid[:-1] + raygrid[1:])
            weights = np.pi * (np.sin(bounds[1:])**2 - np.sin(bounds[:-1])**2)

        self.quadrature_weights = np.expand_dims(weights, axis=1)
        self.nangles = len(self.quadrature_mu)
        self.f_dilution = inputs.f_dilution

        # Radiative-transfer path:
        self._rt_path = inputs.rt_path

        if inputs.wlunits is None:
            self.wlunits = 'um'
        else:
            self.wlunits = inputs.wlunits
        wl_units = pt.u(inputs.wlunits)

        # Low wavenumber boundary:
        if inputs.wnlow is None and inputs.wlhigh is None:
            log.error(
                'Undefined low wavenumber boundary.  Either set wnlow or wlhigh'
            )
        if inputs.wnlow is not None and inputs.wlhigh is not None:
            log.warning(
                f'Both wnlow ({self.wnlow:.2e} cm-1) and wlhigh '
                 '({self.wlhigh:.2e} cm) were defined.  wlhigh will be ignored'
            )

        if inputs.wnlow is not None:
            self.wnlow = inputs.wnlow
            self.wlhigh = 1.0 / self.wnlow
        else:
            self.wlhigh = inputs.wlhigh
            self.wnlow = 1.0 / self.wlhigh


        # High wavenumber boundary:
        if inputs.wnhigh is None and inputs.wllow is None:
            log.error(
                'Undefined high wavenumber boundary. Either set wnhigh or wllow'
            )
        if inputs.wnhigh is not None and inputs.wllow is not None:
            log.warning(
                f'Both wnhigh ({self.wnhigh:.2e} cm-1) and wllow '
                 '({self.wllow:.2e} cm) were defined.  wllow will be ignored'
            )

        if inputs.wnhigh is not None:
            self.wnhigh = inputs.wnhigh
            self.wllow = 1.0 / self.wnhigh
        else:
            self.wllow = inputs.wllow
            self.wnhigh = 1.0 / self.wllow

        # Consistency check (wnlow < wnhigh):
        if self.wnlow > self.wnhigh:
            log.error(
                f'Wavenumber low boundary ({self.wnlow:.1f} cm-1) must be '
                f'larger than the high boundary ({self.wnhigh:.1f} cm-1)'
            )

        self.resolution = None
        self.wlstep = None

        # If there are cross-section tables, take sampling from there:
        if pt.isfile(inputs.extfile) == 1 and inputs.runmode != 'opacity':
            wn = io.read_opacity(inputs.extfile[0], extract='arrays')[3]

            # Update wavenumber sampling:
            wn_mask = (wn >= self.wnlow) & (wn <= self.wnhigh)
            self.wn = wn[wn_mask][::inputs.wn_thinning]
            self.nwave = len(self.wn)
            self.spectrum = np.zeros(self.nwave, np.double)

            if self.wnlow <= self.wn[0]:
                self.wnlow = self.wn[0]
            if self.wnhigh >= self.wn[-1]:
                self.wnhigh = self.wn[-1]

            # Guess sampling by looking at std of sampling rates:
            dwn = np.ediff1d(np.abs(self.wn))
            dwl = np.ediff1d(np.abs(1.0/self.wn))
            res = np.abs(self.wn[1:]/dwn)

            std_dwn = np.std(dwn/np.mean(dwn))
            std_dwl = np.std(dwl/np.mean(dwl))
            std_res = np.std(res/np.mean(res))
            if std_dwn < std_dwl and std_dwn < std_res:
                self.wnstep = self.wn[1] - self.wn[0]
                sampling_text = f'sampling rate = {self.wnstep:.2f} cm-1'
            elif std_dwl < std_dwn and std_dwl < std_res:
                self.wlstep = np.abs(1/self.wn[0] - 1/self.wn[1]) / pc.um
                sampling_text = f'sampling rate = {self.wlstep:.6f} um'
            else:
                g = self.wn[-2]/self.wn[-1]
                # Assuming no one would care for a R with more than 5 decimals:
                self.resolution = np.round(0.5*(1+g)/(1-g), decimals=5)
                #self.wnhigh = 2 * ex.wn[-1]/(1+g)
                sampling_text = f'R = {self.resolution:.1f}'

            log.msg(
                "Reading spectral sampling from extinction-coefficient "
                f"table.  Adopting array with {sampling_text}, "
                f"and {self.nwave} samples between "
                f"[{self.wnlow:.2f}, {self.wnhigh:.2f}] cm-1."
            )
            return


        # At least one sampling mode must be defined:
        undefined_sampling_rate = (
            inputs.wnstep is None and
            inputs.wlstep is None and
            inputs.resolution is None
        )
        if undefined_sampling_rate:
            log.error(
                'Undefined spectral sampling rate, either set resolution, '
                'wnstep, or wlstep'
            )

        # highly composite numbers
        hcn = np.array([
            1, 2, 4, 6, 12, 24, 36, 48, 60, 120, 180, 240, 360, 720, 840,
            1260, 1680, 2160, 2520, 5040, 7560, 10080, 15120, 20160, 25200,
            27720, 45360, 50400, 55440, 83160, 110880, 221760, 277200,
        ])

        if inputs.wnstep is not None:
            self.wnstep = inputs.wnstep
        # Default to a wavenumber supersampling ~0.0004 (R ~2e7 at 1.0 um)
        if inputs.wnosamp is None:
            if inputs.wnstep is None:
                self.wnstep = 1.0
            self.wnosamp = hcn[self.wnstep/hcn <= 0.0004][0]
        else:
            self.wnosamp = inputs.wnosamp

        self.resolution = inputs.resolution
        self.wlstep = inputs.wlstep

        if inputs.resolution is not None:
            # Constant-resolving power wavenumber sampling:
            self.wn = ps.constant_resolution_spectrum(
                self.wnlow, self.wnhigh, self.resolution,
            )
            self.wlstep = None
        elif self.wlstep is not None:
            # Constant-sampling rate wavelength sampling:
            wl = np.arange(self.wllow, self.wlhigh, self.wlstep)
            self.wn = 1.0/np.flip(wl)
            self.wnlow = self.wn[0]
            self.resolution = None
        else:
            # Constant-sampling rate wavenumber sampling:
            nwave = int((self.wnhigh-self.wnlow)/self.wnstep) + 1
            self.wn = self.wnlow + np.arange(nwave) * self.wnstep
        self.nwave = len(self.wn)

        # Fine-sampled wavenumber array:
        self.ownstep = self.wnstep / self.wnosamp
        self.onwave = int(np.ceil((self.wn[-1]-self.wnlow)/self.ownstep)) + 1
        self.own = self.wnlow + np.arange(self.onwave) * self.ownstep
        self.spectrum = np.zeros(self.nwave, np.double)

        # Get list of divisors:
        self.odivisors = pt.divisors(self.wnosamp)

        # Re-set final boundary (stay inside given boundaries):
        if self.wn[-1] != self.wnhigh:
            log.warning(
                f'Final wavenumber modified from {self.wnhigh:.4f} cm-1 (input)'
                f'\n                            to {self.wn[-1]:.4f} cm-1'
            )

        # Screen output:
        log.msg(
            f'Initial wavenumber boundary:  {self.wnlow:.5e} cm-1  '
            f'({self.wlhigh/wl_units:.3e} {self.wlunits})\n'
            f'Final   wavenumber boundary:  {self.wnhigh:.5e} cm-1  '
            f'({self.wllow/wl_units:.3e} {self.wlunits})',
            indent=2,
        )

        if self.resolution is not None:
            msg = f'Spectral resolving power: {self.resolution:.1f}'
        elif self.wlstep is not None:
            wl_step = self.wlstep / wl_units
            msg = f'Wavelength sampling interval: {wl_step:.2g} {self.wlunits}'
        else:
            msg = f'Wavenumber sampling interval: {self.wnstep:.2g} cm-1'
        log.msg(
            f'{msg}\n'
            f'Wavenumber sample size:      {self.nwave:8d}\n'
            f'Wavenumber fine-sample size: {self.onwave:8d}\n',
            indent=2,
        )
        log.head('Wavenumber sampling done.')


    def __str__(self):
        fmt = {'float': '{: .3e}'.format}
        fw = pt.Formatted_Write()
        fw.write('Spectral information:')
        fw.write('Wavenumber internal units: cm-1')
        fw.write('Wavelength internal units: cm')
        fw.write('Wavelength display units (wlunits): {:s}', self.wlunits)
        fw.write(
            'Low wavenumber boundary (wnlow):   {:10.3f} cm-1  '
            '(wlhigh = {:6.2f} {})',
            self.wnlow, self.wlhigh/pt.u(self.wlunits), self.wlunits,
        )
        fw.write(
            'High wavenumber boundary (wnhigh): {:10.3f} cm-1  '
            '(wllow  = {:6.2f} {})',
            self.wnhigh, self.wllow/pt.u(self.wlunits), self.wlunits,
        )
        fw.write('Number of samples (nwave): {:d}', self.nwave)
        if self.resolution is None:
            fw.write('Sampling interval (wnstep): {:.3f} cm-1', self.wnstep)
        else:
            fw.write(
                'Spectral resolving power (resolution): {:.1f}',
                self.resolution,
            )
        fw.write(
            'Wavenumber array (wn, cm-1):\n    {}',
            self.wn,
            fmt={'float': '{: .3f}'.format},
        )
        fw.write('Oversampling factor (wnosamp): {:d}', self.wnosamp)


        fw.write(
            '\nGaussian quadrature cos(theta) angles (quadrature_mu):\n    {}',
            self.quadrature_mu,
            prec=3,
        )
        fw.write(
            'Gaussian quadrature weights (quadrature_weights):\n    {}',
            self.quadrature_weights.flatten(),
            prec=3,
        )
        if self.intensity is not None:
            fw.write('Intensity spectra (intensity, erg s-1 cm-2 sr-1 cm):')
            for intensity in self.intensity:
                fw.write('    {}', intensity, fmt=fmt, edge=3)
        if self._rt_path in pc.emission_rt:
            fw.write(
                'Emission spectrum (spectrum, erg s-1 cm-2 cm):\n    {}',
                self.spectrum, fmt=fmt, edge=3,
            )
        elif self._rt_path in pc.transmission_rt:
            fw.write(
                '\nTransmission spectrum, (Rp/Rs)**2 (spectrum):\n    {}',
                self.spectrum, fmt=fmt, edge=3,
            )
        return fw.text


def spectrum(pyrat):
    """
    Spectrum calculation driver.
    """
    pyrat.log.head('\nCalculate the planetary spectrum.')

    # Initialize the spectrum array:
    pyrat.spec.spectrum = np.empty(pyrat.spec.nwave, np.double)
    if pyrat.opacity.is_patchy:
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

    if pyrat.spec.f_dilution is not None and pyrat.od.rt_path in pc.emission_rt:
        pyrat.spec.spectrum *= pyrat.spec.f_dilution

    # Print spectrum to file:
    if pyrat.od.rt_path in pc.transmission_rt:
        spec_type = 'transit'
    elif pyrat.od.rt_path in pc.emission_rt:
        spec_type = 'emission'

    io.write_spectrum(
        pyrat.spec.wl,
        pyrat.spec.spectrum,
        pyrat.spec.specfile,
        spec_type,
    )
    if pyrat.spec.specfile is not None:
        specfile = f": '{pyrat.spec.specfile}'"
    else:
        specfile = ""
    pyrat.log.head(f"Computed {spec_type} spectrum{specfile}.", indent=2)
    pyrat.log.head('Done.')


def modulation(pyrat):
    """Calculate transmission spectrum for transit geometry"""
    rtop = pyrat.atm.rtop
    radius = pyrat.atm.radius
    depth = pyrat.od.depth

    # Get Delta radius (and simps' integration variables):
    h = np.ediff1d(radius[rtop:])
    # The integrand:
    integ = (np.exp(-depth[rtop:,:]) * np.expand_dims(radius[rtop:],1))

    if pyrat.opacity.is_patchy:
        depth_clear = pyrat.od.depth_clear
        h_clear = np.copy(h)
        integ_clear = (
            np.exp(-depth_clear[rtop:,:]) * np.expand_dims(radius[rtop:],1)
        )

    for model in pyrat.opacity.models:
        if model.name == 'deck':
            # Replace (by interpolating) last layer with cloud top:
            if model.itop > rtop:
                h[model.itop-rtop-1] = model.rsurf - radius[model.itop-1]
                integ[model.itop-rtop] = interp1d(
                    radius[rtop:], integ, axis=0)(model.rsurf)
            break

    # Number of layers for integration at each wavelength:
    nlayers = pyrat.od.ideep - rtop + 1
    spectrum = t.trapz2D(integ, h, nlayers-1)
    pyrat.spec.spectrum = (radius[rtop]**2 + 2*spectrum) / pyrat.phy.rstar**2

    if pyrat.opacity.is_patchy:
        nlayers = pyrat.od.ideep_clear - rtop + 1
        pyrat.spec.clear = t.trapz2D(integ_clear, h_clear, nlayers-1)

        pyrat.spec.clear = (
            (radius[rtop]**2 + 2*pyrat.spec.clear) / pyrat.phy.rstar**2
        )
        pyrat.spec.cloudy = pyrat.spec.spectrum
        pyrat.spec.spectrum = (
            pyrat.spec.cloudy * pyrat.opacity.fpatchy +
            pyrat.spec.clear * (1-pyrat.opacity.fpatchy)
        )


def intensity(pyrat):
    """
    Calculate the intensity spectrum (erg s-1 cm-2 sr-1 cm) for
    eclipse geometry.
    """
    spec = pyrat.spec
    pyrat.log.msg('Computing intensity spectrum.', indent=2)

    # Allocate intensity array:
    spec.intensity = np.empty((spec.nangles, spec.nwave), np.double)

    # Calculate the Planck Emission:
    pyrat.od.B = np.zeros((pyrat.atm.nlayers, spec.nwave), np.double)
    ps.blackbody_wn_2D(spec.wn, pyrat.atm.temp, pyrat.od.B, pyrat.od.ideep)

    for model in pyrat.opacity.models:
        if model.name == 'deck':
            pyrat.od.B[model.itop] = ps.blackbody_wn(pyrat.spec.wn, model.tsurf)

    # Plane-parallel radiative-transfer intensity integration:
    spec.intensity = t.intensity(
        pyrat.od.depth, pyrat.od.ideep, pyrat.od.B, spec.quadrature_mu,
        pyrat.atm.rtop,
    )


def flux(pyrat):
    """
    Calculate the hemisphere-integrated flux spectrum (erg s-1 cm-2 cm)
    for eclipse geometry.
    """
    # Weight-sum the intensities to get the flux:
    pyrat.spec.spectrum[:] = np.sum(
        pyrat.spec.intensity * pyrat.spec.quadrature_weights,
        axis=0,
    )


def two_stream(pyrat):
    """
    Two-stream approximation radiative transfer
    following Heng et al. (2014)

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
    pyrat.log.msg('Compute two-stream flux spectrum.', indent=2)
    spec = pyrat.spec
    phy = pyrat.phy
    nlayers = pyrat.atm.nlayers

    # Set internal net bolometric flux to sigma*Tint**4:
    spec.f_int = ps.blackbody_wn(spec.wn, pyrat.atm.tint)
    total_f_int = np.trapezoid(spec.f_int, spec.wn)
    if total_f_int > 0:
        spec.f_int *= pc.sigma * pyrat.atm.tint**4 / total_f_int

    # Diffusivity factor (Eq. B5 of Heng et al. 2014):
    dtau0 = np.diff(pyrat.od.depth, n=1, axis=0)
    trans = (1-dtau0)*np.exp(-dtau0) + dtau0**2 * ss.exp1(dtau0)

    B = pyrat.od.B = ps.blackbody_wn_2D(spec.wn, pyrat.atm.temp)
    Bp = np.diff(pyrat.od.B, n=1, axis=0) / dtau0

    # Diffuse approximation to compute downward and upward fluxes:
    spec.flux_down = np.zeros((nlayers, spec.nwave))
    spec.flux_up = np.zeros((nlayers, spec.nwave))

    # TBD: Make this an error if any data is not there
    is_irradiation = (
        spec.starflux is not None
        and pyrat.atm.smaxis is not None
        and phy.rstar is not None
    )
    # Top boundary condition:
    if is_irradiation:
        spec.flux_down[pyrat.atm.rtop] = \
            pyrat.atm.beta_irr * (phy.rstar/pyrat.atm.smaxis)**2 * spec.starflux
    # Eqs. (B6) of Heng et al. (2014):
    for i in range(nlayers-1):
    # TBD: Can I make this work of rtop is below the top layer?
    #for i in range(pyrat.atm.rtop, nlayers-1):
        spec.flux_down[i+1] = (
            trans[i] * spec.flux_down[i]
            + np.pi * B[i] * (1-trans[i])
            + np.pi * Bp[i] * (
                  -2/3 * (1-np.exp(-dtau0[i])) + dtau0[i]*(1-trans[i]/3))
        )

    spec.flux_up[nlayers-1] = spec.flux_down[nlayers-1] + spec.f_int
    for i in reversed(range(nlayers-1)):
    #for i in reversed(range(pyrat.atm.rtop, nlayers-1)):
        spec.flux_up[i] = (
            trans[i] * spec.flux_up[i+1]
            + np.pi * B[i+1] * (1-trans[i])
            + np.pi * Bp[i] * (
                  2/3 * (1-np.exp(-dtau0[i])) - dtau0[i]*(1-trans[i]/3))
        )

    spec.spectrum = spec.flux_up[0]

