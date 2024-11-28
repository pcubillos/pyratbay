# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import multiprocessing as mp
from collections import OrderedDict
import os
import subprocess

import numpy as np
import mc3

from .. import atmosphere as pa
from .. import constants as pc
from .. import io as io
from .. import plots as pp
from .. import spectrum as ps
from .. import tools as pt

from .atmosphere import Atmosphere
from .observation import Observation
from .opacity import Opacity
from .retrieval import Retrieval
from .voigt import Voigt
from . import spectrum as sp
from .  import extinction as ex
from .  import optical_depth as od
from .  import objects as ob
from .  import argum as ar


class Pyrat():
    """
    Main Pyrat object.
    """
    def __init__(self, cfg_file, log=True, mute=False):
        """
        Parse the command-line arguments into the pyrat object.

        Parameters
        ----------
        cfg_file: String
            A Pyrat Bay configuration file.
        log: Bool
            Flag to save screen outputs to file (True) or not (False)
            (e.g., to prevent overwritting log of a previous run).
        mute: Bool
            If True, enforce verb to take a value of -1.

        Examples
        --------
        >>> import pyratbay as pb
        >>> pyrat = pb.run('spectrum_transmission.cfg')
        """
        # Setup time tracker:
        timer = pt.Timer()
        self.timestamps = OrderedDict()
        self.timestamps['init'] = timer.clock()

        # Parse config file inputs:
        if isinstance(cfg_file, str):
            self.inputs, self.log = pt.parse(cfg_file, log, mute)
        else:
            # If cfg file was already parsed (pb.run), sneakily use log
            # argument as the logging object instead of a boolean flag
            self.inputs = cfg_file
            self.log = log

        self.ncpu = self.inputs.ncpu
        self.runmode = self.inputs.runmode

        self.phy = ob.Physics(self.inputs)
        # TBD: Remove self.ex entirely?
        self.ex = ob.Extinction(self.inputs, self.log)
        self.od = ob.Optdepth(self.inputs, self.log)

        # Initialize Atmosphere:
        self.atm = Atmosphere(self.inputs, self.log)
        self.timestamps['atmosphere'] = timer.clock()

        # Initialize wavenumber sampling:
        self.spec = sp.Spectrum(self.inputs, self.log)
        self.timestamps['spectrum'] = timer.clock()

        self.obs = Observation(self.inputs, self.spec.wn, self.log)
        ar.check_spectrum(self)

        # Setup opacity models:
        self.opacity = Opacity(
            self.inputs,
            self.spec.wn,
            self.atm.species,
            self.atm.press,
            self.log,
            self,
        )
        self.timestamps['read opacities'] = timer.clock()

        if 'lbl' in self.opacity.models_type:
            i_lbl = self.opacity.models_type.index('lbl')
            lbl = self.opacity.models[i_lbl]
            self.voigt = Voigt(
                self.inputs,
                lbl,
                self.ex,
                self.spec,
                self.atm,
                self.log,
            )
            self.timestamps['voigt'] = timer.clock()

        # Setup more retrieval parameters:
        ar.setup(self)
        self.ret = Retrieval(
            self.inputs,
            self.atm,
            self.phy,
            self.obs,
            self.opacity,
            self.log,
        )


    def compute_opacity(self):
        """
        Calculate opacity (cm2 molecule-1) tabulated over
        temperature, pressure, and wavenumber arrays
        """
        ex.compute_opacity(self)


    def run(self, temp=None, vmr=None, radius=None, skip=[]):
        """
        Evaluate a Pyrat spectroscopic model

        Parameters
        ----------
        temp: 1D float ndarray
            Updated atmospheric temperature profile in Kelvin, of size nlayers.
        abund: 2D float ndarray
            Updated atmospheric abundances profile by number density, of
            shape [nlayers, nmol].
        radius: 1D float ndarray
            Updated atmospheric altitude profile in cm, of size nlayers.
        skip: List of strings
            If set, the opacity from the model names or line-sample species
            listed here will be neglected.
        """
        timer = pt.Timer()

        # Re-calculate atmospheric properties if required:
        self.atm.calc_profiles(temp, vmr, radius, self.phy.mstar, self.log)

        out_of_bounds = self.opacity.check_temp_bounds(self.atm.temp)
        good_status = len(out_of_bounds) == 0
        if not good_status:
            self.log.warning(
                "Temperature values lie out of the cross-section "
                f"boundaries for: {out_of_bounds}"
            )
            self.spec.spectrum[:] = 0.0
            return

        # Calculate extinction coefficient:
        self.opacity.calc_extinction_coefficient(
            self.atm.temp, self.atm.radius, self.atm.d, skip=skip,
        )
        self.timestamps['extinction'] = timer.clock()

        # Calculate the optical depth:
        od.optical_depth(self)
        self.timestamps['odepth'] = timer.clock()

        # Calculate the spectrum:
        sp.spectrum(self)
        self.timestamps['spectrum'] = timer.clock()

        self.log.msg(
            "\nTimestamps (s):\n" +
            "\n".join(
                f"{key:10s}: {val:10.6f}"
                for key,val in self.timestamps.items()
            )
        )


    def eval(self, params, retmodel=True, skip=[]):
        """
        Fitting routine for atmospheric retrieval

        Parameters
        ----------
        params: 1D float iterable
            Array of fitting parameters that define the atmosphere.
        retmodel: Bool
            Flag to include the model spectra in the return.
        skip: List of strings
            If set, the opacity from the model names or line-sample species
            listed here will be neglected.

        Returns
        -------
        spectrum: 1D float ndarray
            The output model spectra.  Returned only if retmodel=True.
        bandflux: 1D float ndarray
            The waveband-integrated spectrum values.
        """
        atm = self.atm
        ret = self.ret
        obs = self.obs
        params = np.asarray(params)

        if len(params) != ret.nparams:
            self.log.warning(
                f'The number of input fitting parameters ({len(params)}) does '
                f'not match\nthe number of required parameters ({ret.nparams})'
            )
            return None, None if retmodel else None

        # Update models parameters:
        if ret.itemp is not None:
            ifree = ret.map_pars['temp']
            atm.tpars[ifree] = params[ret.itemp]

        if ret.imol is not None:
            for j,imol in enumerate(ret.imol):
                imodel,idx = ret.map_pars['mol'][j]
                atm.vmr_pars[imodel][idx] = params[imol]

        if ret.irad is not None:
            self.atm.rplanet = params[ret.irad][0] * pt.u(atm.runits)
        elif ret.ipress is not None:
            self.atm.refpressure = 10.0**params[ret.ipress][0]

        if ret.imass is not None:
            self.atm.mplanet = params[ret.imass][0] * pt.u(self.atm.mass_units)

        ifree = ret.map_pars['opacity']
        for j,model in enumerate(self.opacity.models):
            if ifree[j] == []:
                continue
            idx = ifree[j]
            ipar = ret.iopacity[j]
            model.pars[idx] = params[ipar]

        if ret.ipatchy is not None:
            self.opacity.fpatchy = params[ret.ipatchy][0]

        if ret.itstar is not None:
            self.phy.tstar = params[ret.itstar][0]
            self.spec.starflux = self.spec.flux_interp(self.phy.tstar)
            self.obs.bandflux_star = np.array([
                band(self.spec.starflux)
                for band in self.obs.filters
            ])

        if ret.idilut is not None:
            self.spec.f_dilution = params[ret.idilut][0]

        # Calculate atmosphere and spectrum:
        self.run(skip=skip)

        reject_flag = False
        # Turn-on reject flag if temperature is out-of-bounds:
        temp = atm.temp
        if np.any(temp < ret.tlow) or np.any(temp > ret.thigh):
            temp[:] = 0.0
            reject_flag = True
            self.log.warning(
                "Input temperature profile runs out of "
                f"boundaries ({ret.tlow:.1f}--{ret.thigh:.1f} K)"
            )
        # Check abundaces stay within bounds:
        if pa.qcapcheck(atm.vmr, ret.qcap, atm.ibulk):
            reject_flag = True
            self.log.warning(
                "The sum of trace abundances' VMRs exceeds "
                f"the cap of {ret.qcap:.3f}"
            )

        # Band-integrate spectrum:
        obs.bandflux = self.band_integrate()

        # Instrumental offset:
        if ret.ioffset is not None:
            ifree = ret.map_pars['offset']
            obs.offset_pars[ifree] = params[ret.ioffset]
            obs.data = obs.depth.offset_data(obs.offset_pars, obs.units)

        # Uncertainty scaling:
        if ret.ierror is not None:
            ifree = ret.map_pars['error']
            obs.uncert_pars[ifree] = params[ret.ierror]
            obs.uncert = obs.depth.scale_errors(obs.uncert_pars, obs.units)

        # Invalid model:
        if not np.any(obs.bandflux):
            reject_flag = True

        # Reject this iteration if there are invalid temperatures or radii:
        if obs.bandflux is not None and reject_flag:
            obs.bandflux[:] = np.inf

        ret.params = np.copy(params)
        if retmodel:
            return self.spec.spectrum, obs.bandflux

        return obs.bandflux


    def retrieval(self):
        """
        Run an MCMC or nested-sampling atmospheric retrieval.
        """
        ret = self.ret
        obs = self.obs
        log = self.log

        if ret.mcmcfile is None:
            log.error('Undefined retrieval file (mcmcfile)')
        if ret.sampler is None:
            log.error(
                'Undefined retrieval algorithm (sampler).  '
                f'Select from {pc.samplers}'
            )
        if ret.params is None:
            log.error(
                'Undefined retrieval fitting parameters (retrieval_params)'
            )
        if ret.pstep is None:
            log.error('Missing pstep argument, required for retrieval runs')

        if obs.data is None:
            log.error("Undefined transit/eclipse data (data)")
        if obs.uncert is None:
            log.error("Undefined data uncertainties (uncert)")
        if obs.nfilters == 0:
            log.error("Undefined transmission filters (filters)")

        # Create output folder if needed:
        pt.mkdir(ret.mcmcfile)
        # Basename of the output files (no extension):
        output, extension = os.path.splitext(ret.mcmcfile)

        # MultiNest wrapper call:
        if ret.sampler == 'multinest':
            sampler_output = pt.multinest_run(self, output)
            if pt.get_mpi_rank() != 0:
                return
            posterior = sampler_output['posterior']

        # mc3 MCMC wrapper call:
        if ret.sampler == 'snooker':
            if ret.nsamples is None:
                log.error('Undefined number of retrieval samples (nsamples)')
            if ret.burnin is None:
                log.error('Undefined number of retrieval burn-in samples (burnin)')
            if ret.nchains is None:
                log.error('Undefined number of retrieval parallel chains (nchains)')

            # TBD: Fix resuming
            ret.resume = False
            # Mute logging in pyrat object, but not in mc3:
            self.log = mc3.utils.Log(verb=-1, width=80)
            self.spec.specfile = None  # Avoid writing spectrum file during MCMC
            retmodel = False  # Return only the band-integrated spectrum
            # Run MCMC:
            sampler_output = mc3.sample(
                data=self.obs.data, uncert=self.obs.uncert,
                func=self.eval, indparams=[retmodel], params=ret.params,
                pmin=ret.pmin, pmax=ret.pmax, pstep=ret.pstep,
                prior=ret.prior, priorlow=ret.priorlow, priorup=ret.priorup,
                sampler=ret.sampler, nsamples=ret.nsamples,
                nchains=ret.nchains, burnin=ret.burnin, thinning=ret.thinning,
                grtest=True, grbreak=ret.grbreak, grnmin=ret.grnmin,
                log=log, ncpu=self.ncpu,
                plots=True, showbp=True, theme=ret.theme,
                pnames=ret.pnames, texnames=ret.texnames,
                resume=ret.resume, savefile=ret.mcmcfile,
            )
            if sampler_output is None:
                log.error("Error in mc3")
            posterior, zchain, zmask = mc3.utils.burn(sampler_output)


        # Post processing (can be done directly from posterior outputs)
        ret.bestp = bestp = sampler_output['bestp']
        ret.posterior = posterior

        # Best-fitting model:
        self.spec.specfile = f"{output}_bestfit_spectrum.dat"
        ret.spec_best, ret.bestbandflux = self.eval(bestp)

        atm = self.atm
        header = "# Retrieval best-fitting atmospheric model.\n\n"
        bestatm = f"{output}_bestfit_atmosphere.atm"
        io.write_atm(
            bestatm, atm.press, atm.temp, atm.species,
            atm.vmr, radius=atm.radius,
            punits=atm.punits, runits=atm.runits, header=header,
        )
        filename = f'{output}_bestfit_spectrum.png'
        self.plot_spectrum(spec='best', filename=filename)

        # Temperature profiles:
        if atm.temp_model is not None:
            tparams = atm.tpars
            tparams[ret.map_pars['temp']] = bestp[ret.itemp]
            ret.temp_best = atm.temp_model(tparams)

            nsamples, nfree = np.shape(posterior)
            t_posterior = np.tile(tparams, (nsamples,1))
            # Map temperature free parameters from posterior to tparams:
            ifree = np.where(self.ret.pstep>0)[0]
            for j, imap in zip(ret.itemp, ret.map_pars['temp']):
                if j in ifree:
                    ipost = list(ifree).index(j)
                    t_posterior[:,imap] = posterior[:,ipost]
            tpost = pa.temperature_posterior(t_posterior, atm.temp_model)
            ret.temp_median = tpost[0]
            ret.temp_post_boundaries = tpost[1:]
            self.plot_temperature(
                filename=f'{output}_posterior_temperature_profile.png',
            )

        is_emission = self.od.rt_path in pc.emission_rt
        is_transmission = self.od.rt_path in pc.transmission_rt

        if is_emission:
            contrib = ps.contribution_function(
                self.od.depth, atm.press, self.od.B,
            )
        elif is_transmission:
            contrib = ps.transmittance(self.od.depth, self.od.ideep)
        bands_idx = [band.idx for band in self.obs.filters]
        bands_response = [band.response for band in self.obs.filters]
        band_cf = ps.band_cf(contrib, bands_response, self.spec.wn, bands_idx)

        path = 'transit' if is_transmission else 'emission'
        pp.contribution(
            band_cf, 1.0/(self.obs.bandwn*pc.um),
            path, atm.press,
            filename=f'{output}_bestfit_cf.png',
        )

        self.log = log  # Un-mute
        root_output = os.path.split(output)[0]
        log.msg(f"\nOutput retrieval files located at {root_output}")

        if self.inputs.post_processing:
            os.environ['PBAY_NO_MPI'] = "1"
            subprocess.call(
                f'pbay --post {self.inputs.config_file} &',
                shell=True,
            )


    def radiative_equilibrium(
            self, nsamples=None, continue_run=False, convection=False,
        ):
        """
        Compute radiative-thermochemical equilibrium atmosphere.
        Currently there is no convergence criteria implemented,
        some 100--300 iterations are typically sufficient to converge
        to a stable temperature-profile solution.

        Parameters
        ----------
        nsamples: Integer
            Number of radiative-equilibrium iterations to run.
        continue_run: Bool
            If True, continue from a previous radiative-equilibrimu run.
        convection: Bool
            If True, skip convective flux calculation in the radiative
            equilibrium calculation.

        Returns
        -------
        There are no returned values, but this method updates the
        temperature profile (self.atm.temp) and abundances (self.atm.vmr)
        with the values from the last radiative-equilibrium iteration.

        This method also defines self.atm.radeq_temps, a 2D array
        containing all temperature-profile iterations.
        """
        atm = self.atm

        if nsamples is None:
            nsamples = self.inputs.nsamples

        # No outputs while iterating
        tmp_verb = self.log.verb
        self.log.verb = 0
        basename, extension = os.path.splitext(self.spec.specfile)
        self.spec.specfile = None

        # Enforce two-stream RT:
        rt_path = self.od.rt_path
        self.od.rt_path = 'emission_two_stream'
        tmin = np.amax(list(self.opacity.tmin.values()))
        tmax = np.amin(list(self.opacity.tmax.values()))

        # Initial temperature scale factor
        if not hasattr(atm, '_dt_scale') or not continue_run:
            atm._dt_scale = np.tile(1.0e5, atm.nlayers)

        if hasattr(atm, 'radeq_temps') and continue_run:
            radeq_temps = atm.radeq_temps
        else:
            radeq_temps = np.atleast_2d(atm.temp)

        print("\nRadiative-thermochemical equilibrium calculation:")
        radeq_temps = ps.radiative_equilibrium(
            atm.press,
            radeq_temps,
            nsamples,
            atm.chem_model,
            self.run,
            self.spec.wn,
            self.spec,
            atm,
            convection,
            tmin, tmax,
        )

        # Update last tempertature iteration and save to file:
        self.atm.radeq_temps = radeq_temps
        atm.temp = radeq_temps[-1]
        io.write_atm(
            f'{basename}.atm',
            atm.press, atm.temp, atm.species, atm.vmr,
            punits="bar",
            header="# Radiative-thermochemical equilibrium profile.\n\n",
        )
        self.od.rt_path = rt_path
        np.savez(f'{basename}.npz', pressure=atm.press, temps=radeq_temps)
        self.spec.specfile = f'{basename}.dat'
        spec_type = 'emission'
        io.write_spectrum(
            self.spec.wl,
            self.spec.spectrum,
            self.spec.specfile,
            spec_type,
        )
        self.log.verb = tmp_verb


    def band_integrate(self, spectrum=None):
        """
        Band-integrate transmission spectrum (transit) or planet-to-star
        flux ratio (eclipse) over transmission band passes.
        """
        if self.obs.filters is None:
            return None
        if spectrum is None:
            spectrum = self.spec.spectrum

        bandflux = np.array([band(spectrum) for band in self.obs.filters])
        if self.od.rt_path in pc.emission_rt:
            rprs_square = (self.atm.rplanet/self.phy.rstar)**2.0
            bandflux = bandflux / self.obs.bandflux_star * rprs_square

        self.obs.bandflux = bandflux
        return self.obs.bandflux



    def get_ec(self, layer):
        """
        Extract extinction-coefficient contribution (in cm-1) from each
        component of the atmosphere at the requested layer.

        Parameters
        ----------
        layer: Integer
           The index of the atmospheric layer where to extract the EC.

        Returns
        -------
        ec: 2D float ndarray
           An array of shape [ncomponents, nwave] with the EC spectra
           (in cm-1) from each component of the atmosphere.
        label: List of strings
           The names of each atmospheric component that contributed to EC.
        """
        if len(self.opacity.models) > 0:
            return self.opacity.get_ec(self.atm.temp, self.atm.d, layer)
        return None, []


    def percentile_spectrum(self, nmax=None):
        """Compute spectrum posterior percentiles."""
        if self.ret.posterior is None:
            print('pyrat objec does not have a posterior distribution.')
            return

        nsamples = np.shape(self.ret.posterior)[0]
        draws = np.arange(nsamples)
        if nmax is not None:
            nmax = np.clip(nmax, 0, nsamples)
            draws = np.random.choice(draws, nmax, replace=False)

        # Unique MCMC samples:
        u, uind, uinv = np.unique(self.ret.posterior[draws,0],
            return_index=True, return_inverse=True)
        print('Computing {:d} models.'.format(len(u)))

        # Array of all model parameters (with unique samples)
        posterior = np.repeat([self.ret.params], len(u), axis=0)
        ifree = np.where(self.ret.pstep >0)[0]
        posterior[:,ifree] = self.ret.posterior[uind]
        # Need to keep FILE objects out of pool:
        logfile, self.log.file = self.log.file, None
        verb, self.log.verb = self.log.verb, -1

        with mp.get_context('fork').Pool(self.ncpu) as pool:
            models = pool.map(self.eval, posterior)
        models = np.array([model for model, bandm in models])

        self.log.file = logfile
        self.log.verb = verb

        nwave = len(self.spec.wn)
        low1   = np.zeros(nwave)
        low2   = np.zeros(nwave)
        median = np.zeros(nwave)
        high1  = np.zeros(nwave)
        high2  = np.zeros(nwave)
        for i in range(nwave):
            msample = models[uinv,i]
            low2[i]   = np.percentile(msample,  2.275)
            low1[i]   = np.percentile(msample, 15.865)
            median[i] = np.percentile(msample, 50.000)
            high1[i]  = np.percentile(msample, 84.135)
            high2[i]  = np.percentile(msample, 97.725)

        self.ret.spec_median = median
        self.ret.spec_low1 = low1
        self.ret.spec_low2 = low2
        self.ret.spec_high1 = high1
        self.ret.spec_high2 = high2


    def plot_spectrum(self, spec='model', **kwargs):
        """
        Plot spectrum.

        Parameters
        ----------
        spec: String
            Flag indicating which model to plot.  By default plot the
            latest evaulated model (spec='model').  Other options are
            'best' or 'median' to plot the posterior best-fit or median
            model, in which case, the code will plot the 1- and 2-sigma
            boundaries if they have been computed (see
            self.percentile_spectrum).
        kwargs: dict
            Dictionary of arguments to pass into plots.spectrum().
            See help(pyratbay.plots.spectrum).

        Returns
        -------
        ax: AxesSubplot instance
            The matplotlib Axes of the figure.
        """
        args = {
            'wavelength': 1.0/(self.spec.wn*pc.um),
            'data': self.obs.data,
            'uncert': self.obs.uncert,
            'logxticks': self.inputs.logxticks,
            'yran': self.inputs.yran,
            'theme': self.ret._default_theme,
        }

        obs = self.obs
        if obs.nfilters > 0:
            args['bands_wl0'] = [band.wl0 for band in obs.filters]
            args['bands_wl'] = [band.wl for band in obs.filters]
            args['bands_response'] = [band.response for band in obs.filters]
            args['bands_flux'] = obs.bandflux

        if self.ret.spec_low2 is not None:
            args['bounds'] = [
                self.ret.spec_low2,  self.ret.spec_low1,
                self.ret.spec_high1, self.ret.spec_high2,
            ]

        if spec == 'model':
            args['label'] = 'model'
            spectrum = np.copy(self.spec.spectrum)
        elif spec == 'best':
            args['label'] = 'best-fit model'
            spectrum = np.copy(self.ret.spec_best)
            args['bands_flux'] = self.ret.bestbandflux
        elif spec == 'median':
            args['label'] = 'median model'
            spectrum = np.copy(self.ret.spec_median)
            args['bands_flux'] = self.band_integrate(spectrum)
        else:
            print(
                "Invalid 'spec'.  Select from 'model' (default), 'best', "
                "or 'median'."
            )
            return

        # kwargs can overwite any of the previous value:
        args.update(kwargs)

        is_eclipse = (
            self.od.rt_path in pc.emission_rt and
            self.spec.starflux is not None and
            self.atm.rplanet is not None and
            self.phy.rstar is not None
        )

        if self.od.rt_path in pc.transmission_rt:
            args['rt_path'] = 'transit'
        elif is_eclipse:
            args['rt_path'] = 'eclipse'
            rprs = self.atm.rplanet/self.phy.rstar
            spectrum = spectrum/self.spec.starflux * rprs**2.0
            if 'bounds' in args:
                args['bounds'] = [
                    bound/self.spec.starflux * rprs**2.0
                    for bound in args['bounds']
                ]
        else:
            args['rt_path'] = 'emission'
        ax = pp.spectrum(spectrum, **args)
        return ax


    def plot_temperature(self, **kwargs):
        """
        Plot temperature profile.
        If self.ret.posterior exitst, plot the best fit, median, and
        the '1sigma/2sigma' boundaries of the temperature posterior
        distribution.

        Parameters
        ----------
        kwargs: dict
            Dictionary of arguments to pass into plots.temperature().
            See help(pyratbay.plots.temperature).

        Returns
        -------
        ax: AxesSubplot instance
            The matplotlib Axes of the figure.
        """
        kwargs['pressure'] = self.atm.press
        kwargs['theme'] = self.ret.theme
        if self.ret.posterior is None:
            kwargs['profiles'] = [self.atm.temp]
        else:
            kwargs['profiles'] = [self.ret.temp_median, self.ret.temp_best]
            kwargs['labels'] = ['median', 'best-fit']
            kwargs['bounds'] = self.ret.temp_post_boundaries

        ax = pp.temperature(**kwargs)
        return ax


    def __str__(self):
        if self.spec.resolution is not None:
            wave = f"R={self.spec.resolution:.1f}"
        elif self.spec.wlstep is not None:
            wave = f'delta-wl={self.spec.wlstep:.2f}'
        else:
            wave = f"delta-wn={self.spec.wnstep:.3f} cm-1"

        opacities = []
        if len(self.opacity.models) > 0:
            for i,model in enumerate(self.opacity.models):
                if self.opacity.models_type[i] in ['line_sample', 'lbl']:
                    opacities += model.species.tolist()
                elif self.opacity.models_type[i] in ['alkali']:
                    opacities.append(model.species)
                else:
                    opacities.append(model.name)

        pmin = self.atm.press[ 0]
        pmax = self.atm.press[-1]
        wlmin = 1.0/(self.spec.wn[-1]*pc.um)
        wlmax = 1.0/(self.spec.wn[ 0]*pc.um)
        return (
            "Pyrat atmospheric model\n"
            f"configuration file:  '{self.inputs.config_file}'\n"
            f"Pressure profile:  {pmin:.2e} -- {pmax:.2e} bar "
            f"({self.atm.nlayers:d} layers)\n"
            f"Wavelength range:  {wlmin:.2f} -- {wlmax:.2f} um "
            f"({self.spec.nwave:d} samples, {wave})\n"
            f"Composition:\n  {self.atm.species}\n"
            f"Opacity sources:\n  {opacities}"
        )

