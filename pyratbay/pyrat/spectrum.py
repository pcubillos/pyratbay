# Copyright (c) 2021-2025 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Spectrum',
    'spectrum',
    'two_stream',
]

import os
import numpy as np
import scipy.constants as sc
import scipy.special as ss

from .. import constants as pc
from .. import io
from .. import spectrum as ps
from .. import tools as pt


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
        self.clear = None  # Clear modulation spectrum for patchy model
        self.cloudy = None  # Cloudy modulation spectrum for patchy model
        self.starflux = None  # Stellar flux spectrum
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
        if pt.isfile(inputs.sampled_cs) == 1 and inputs.runmode != 'opacity':
            wn = io.read_opacity(inputs.sampled_cs[0], extract='arrays')[3]

            # Update wavenumber sampling:
            wn_mask = ps.wn_mask(wn, self.wnlow, self.wnhigh)
            self.wn = wn[wn_mask][::inputs.wn_thinning]
            self.nwave = len(self.wn)
            self.spectrum = np.zeros(self.nwave, np.double)
            if self._rt_path not in pc.transmission_rt:
                self.fplanet = np.zeros(self.nwave, np.double)

            wn_min = self.wn[0]
            wn_max = self.wn[-1]

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
                f"[{wn_min:.2f}, {wn_max:.2f}] cm-1."
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
        if self._rt_path not in pc.transmission_rt:
            self.fplanet = np.zeros(self.nwave, np.double)

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
        wn_min = np.amin(self.wn)
        wn_max = np.amax(self.wn)
        wl_min = 1/(wn_max*pt.u(self.wlunits))
        wl_max = 1/(wn_min*pt.u(self.wlunits))

        fw.write(
            f'Low wavenumber boundary (wnlow):   {wn_min:10.3f} cm-1  '
            f'(wlhigh = {wl_max:6.2f} {self.wlunits})',
        )
        fw.write(
            f'High wavenumber boundary (wnhigh): {wn_max:10.3f} cm-1  '
            f'(wllow  = {wl_min:6.2f} {self.wlunits})',
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


def _get_cloud_deck(pyrat):
    """Wrapper to get cloud deck parameters if they exist"""
    for model in pyrat.opacity.models:
        if model.name == 'deck':
            return model.rsurf, model.tsurf, model.itop
    return None, None, None


def spectrum(pyrat):
    """
    Spectrum calculation driver.
    """
    pyrat.log.head('\nCalculate the planetary spectrum.')
    spec = pyrat.spec

    # Initialize the spectrum array:
    spec.spectrum = np.empty(spec.nwave, np.double)
    if pyrat.opacity.is_patchy:
        spec.clear  = np.empty(spec.nwave, np.double)
        spec.cloudy = np.empty(spec.nwave, np.double)

    deck_rsurf, deck_tsurf, deck_itop = _get_cloud_deck(pyrat)
    f_patchy = pyrat.opacity.fpatchy


    # Transmission spectroscopy:
    if pyrat.od.rt_path in pc.transmission_rt:
        spec.spectrum = ps.transmission(
            pyrat.od.depth, pyrat.atm.radius, pyrat.atm.rstar,
            pyrat.od.ideep, pyrat.atm.rtop,
            deck_rsurf, deck_itop,
        )
        if pyrat.opacity.is_patchy:
            spec.cloudy = spec.spectrum
            spec.clear = ps.transmission(
                pyrat.od.depth_clear, pyrat.atm.radius, pyrat.atm.rstar,
                pyrat.od.ideep_clear, pyrat.atm.rtop,
            )
            spec.spectrum = f_patchy*spec.cloudy + (1.0-f_patchy)*spec.clear


    # Plane-parallel emission
    elif pyrat.od.rt_path in ['emission', 'eclipse', 'f_lambda']:
        pyrat.od.B = np.zeros((pyrat.atm.nlayers, spec.nwave), np.double)
        ps.blackbody_wn_2D(spec.wn, pyrat.atm.temp, pyrat.od.B)
        weights = spec.quadrature_weights

        spec.intensity, spec.spectrum = ps.plane_parallel_rt(
            pyrat.od.depth, pyrat.od.B, spec.wn, spec.quadrature_mu,
            spec.quadrature_weights, pyrat.od.ideep, pyrat.atm.rtop,
            deck_tsurf, deck_itop,
        )
        spec.spectrum = np.sum(spec.intensity * weights, axis=0)
        if pyrat.opacity.is_patchy:
            spec.cloudy = spec.spectrum
            intensity_clear, spec.clear = ps.plane_parallel_rt(
                pyrat.od.depth_clear, pyrat.od.B, spec.wn, spec.quadrature_mu,
                spec.quadrature_weights, pyrat.od.ideep_clear, pyrat.atm.rtop,
            )
            spec.spectrum = f_patchy*spec.cloudy + (1.0-f_patchy)*spec.clear
        spec.fplanet = spec.spectrum


    # Two stream emission
    elif pyrat.od.rt_path in ['emission_two_stream', 'eclipse_two_stream']:
        two_stream(pyrat)


    # Scaling factors
    is_emission = pyrat.od.rt_path in pc.eclipse_rt + pc.emission_rt
    if is_emission and spec.f_dilution is not None:
        spec.fplanet *= spec.f_dilution
        if pyrat.opacity.is_patchy:
            spec.clear *= spec.f_dilution
            spec.cloudy *= spec.f_dilution

    if pyrat.od.rt_path in pc.eclipse_rt:
        fstar_rprs = 1/spec.starflux * (pyrat.atm.rplanet/pyrat.atm.rstar)**2
        spec.fplanet = np.copy(spec.spectrum)
        spec.spectrum = spec.eclipse = spec.fplanet * fstar_rprs
        if pyrat.opacity.is_patchy:
            spec.fplanet_clear = np.copy(spec.clear)
            spec.fplanet_cloudy = np.copy(spec.cloudy)
            spec.clear = spec.fplanet_clear * fstar_rprs
            spec.cloudy = spec.fplanet_cloudy * fstar_rprs


    # Print spectra to file:
    if pyrat.od.rt_path in pc.transmission_rt:
        spec_type = 'transit'
    elif pyrat.od.rt_path in pc.emission_rt:
        spec_type = 'emission'
    elif pyrat.od.rt_path in pc.eclipse_rt:
        spec_type = 'eclipse'
    io.write_spectrum(
        spec.wl,
        spec.spectrum,
        spec.specfile,
        spec_type,
    )

    # Also save fstar and fplanet when possible
    if spec.specfile is not None:
        file, extension = os.path.splitext(spec.specfile)
        if is_emission and spec.starflux is not None:
            fstar_file = f'{file}_fstar{extension}'
            io.write_spectrum(
                spec.wl,
                spec.starflux,
                fstar_file,
                'emission',
            )
        if pyrat.od.rt_path in pc.eclipse_rt:
            fplanet_file = f'{file}_fplanet{extension}'
            io.write_spectrum(
                spec.wl,
                spec.fplanet,
                fplanet_file,
                'emission',
            )

    if spec.specfile is not None:
        specfile = f": '{spec.specfile}'"
    else:
        specfile = ""
    pyrat.log.head(f"Computed {spec_type} spectrum{specfile}.", indent=2)
    pyrat.log.head('Done.')


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
        and pyrat.atm.rstar is not None
    )
    # Top boundary condition:
    if is_irradiation:
        spec.flux_down[pyrat.atm.rtop] = pyrat.atm.beta_irr * \
            (pyrat.atm.rstar/pyrat.atm.smaxis)**2 * spec.starflux
    # Eqs. (B6) of Heng et al. (2014):
    # TBD: Can I make this work if rtop is below the top layer?
    #for i in range(pyrat.atm.rtop, nlayers-1):
    for i in range(nlayers-1):
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

    pyrat.spec.fplanet = pyrat.spec.spectrum = spec.flux_up[0]

