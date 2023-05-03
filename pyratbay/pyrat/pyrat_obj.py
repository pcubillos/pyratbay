# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import multiprocessing as mp
from collections import OrderedDict

import numpy as np

from .. import atmosphere as pa
from .. import constants as pc
from .. import io as io
from .. import plots as pp
from .. import spectrum as ps
from .. import tools as pt

from .alkali import Alkali
from .atmosphere import Atmosphere
from .crosssec import CIA
from .h_ion import H_Ion
from .line_by_line import Line_By_Line
from .observation import Observation
from . import spectrum as sp
from .  import extinction as ex
from .  import clouds as cl
from .  import optical_depth as od
from .  import objects as ob
from .  import argum as ar
from .  import voigt as v


class Pyrat(object):
  """
  Main Pyrat object.
  """
  def __init__(self, cfile, no_logfile=False, mute=False):
      """
      Parse the command-line arguments into the pyrat object.

      Parameters
      ----------
      cfile: String
          A Pyrat Bay configuration file.
      no_logfile: Bool
          If True, enforce not to write outputs to a log file
          (e.g., to prevent overwritting log of a previous run).
      mute: Bool
          If True, enforce verb to take a value of -1.

      Examples
      --------
      >>> import pyratbay as pb
      >>> # Initialize and execute task:
      >>> pyrat = pb.run('spectrum_transmission.cfg')

      >>> # Initialize only:
      >>> pyrat = pb.Pyrat('spectrum_transmission.cfg')
      >>> # Then, setup internal varible for spectra evaluation:
      >>> pyrat.set_atmosphere()
      >>> pyrat.set_spectrum()
      """
      # Sub-classes:
      self.voigt = ob.Voigt()         # Voigt profile
      self.ex = ob.Extinction()       # Extinction-coefficient
      self.od = ob.Optdepth()         # Optical depth

      # Parse config file inputs:
      pt.parse(self, cfile, no_logfile, mute)

      self.phy = ob.Physics(self.inputs)
      self.ret = ob.Retrieval(self.inputs, self.log)

      # Setup time tracker:
      timer = pt.Timer()
      self.timestamps = OrderedDict()
      self.timestamps['init'] = timer.clock()


  def set_atmosphere(self):
      timer = pt.Timer()
      self.atm = Atmosphere(self.inputs, self.log)
      self.timestamps['read atm'] = timer.clock()


  def set_spectrum(self):
      timer = pt.Timer()
      # Initialize wavenumber sampling:
      self.spec = sp.Spectrum(self.inputs, self.log)
      self.timestamps['wn sample'] = timer.clock()

      # Read opacity tables (if needed):
      ex.read_opacity(self, self.spec._wn_mask)
      self.timestamps['read opacity'] = timer.clock()

      self.obs = Observation(self.inputs, self.spec.wn, self.log)
      ar.check_spectrum(self)

      # Read line-by-line data:
      self.lt = Line_By_Line(
          self.inputs,
          list(self.atm.species),
          self.spec.wnlow, self.spec.wnhigh,
          self.log,
      )

      self.timestamps['read tli'] = timer.clock()

      # Setup more observational/retrieval parameters:
      ar.setup(self)

      # Extinction Voigt grid:
      v.voigt(self)
      self.timestamps['voigt'] = timer.clock()

      self.alkali = Alkali(
          self.inputs.model_names,
          self.atm.press,
          self.spec.wn,
          self.inputs.alkali_cutoff,
          self.atm.species,
          self.log,
      )

      self.h_ion = H_Ion(
          self.inputs,
          self.spec.wn,
          self.atm.species,
          self.log,
      )

      self.cs = CIA(
          self.inputs.cia_files,
          self.spec.wn,
          self.atm.species,
          self.log,
      )
      self.timestamps['read cs'] = timer.clock()


  def compute_opacity(self):
      """
      Calculate opacity (cm2 molecule-1) tabulated over
      temperature, pressure, and wavenumber arrays
      """
      ex.compute_opacity(self)


  def run(self, temp=None, vmr=None, radius=None):
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
      """
      timer = pt.Timer()

      # Re-calculate atmospheric properties if required:
      self.atm.calc_profiles(temp, vmr, radius, self.phy.mstar, self.log)

      # Interpolate CIA absorption:
      good_status = self.cs.calc_extinction_coefficient(
          self.atm.temp, self.atm.d, self.log,
      )
      self.timestamps['interp cs'] = timer.clock()

      # Calculate cloud, Rayleigh, and H- absorption:
      cl.absorption(self)
      self.rayleigh.calc_extinction_coefficient(self.atm.d)
      self.h_ion.calc_extinction_coefficient(self.atm.temp, self.atm.d)
      self.timestamps['cloud+ray'] = timer.clock()

      # Calculate the alkali absorption:
      self.alkali.calc_extinction_coefficient(self.atm.temp, self.atm.d)
      self.timestamps['alkali'] = timer.clock()

      # Calculate the optical depth:
      good_status &= od.optical_depth(self)
      self.timestamps['odepth'] = timer.clock()

      if not good_status:
          self.spec.spectrum[:] = 0.0
          return

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


  def eval(self, params, retmodel=True):
      """
      Fitting routine for atmospheric retrieval

      Parameters
      ----------
      params: 1D float iterable
          Array of fitting parameters that define the atmosphere.
      retmodel: Bool
          Flag to include the model spectra in the return.

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


      # Update temperature parameters (if necessary):
      if ret.itemp is not None:
          ifree = ret.map_pars['temp']
          atm.tpars[ifree] = params[ret.itemp]

      # Update abundance parameters:
      if ret.imol is not None:
          ifree = ret.map_pars['mol']
          atm.molpars[ifree] = params[ret.imol]

      # Update reference radius/pressure if requested:
      if ret.irad is not None:
          self.atm.rplanet = params[ret.irad][0] * pt.u(atm.runits)
      elif ret.ipress is not None:
          p_ref = 10.0**params[ret.ipress][0] * pc.bar
          self.atm.refpressure = p_ref

      # Update planetary mass if requested:
      if ret.imass is not None:
          self.atm.mplanet = params[ret.imass][0] * pt.u(self.atm.mass_units)

      # Update Rayleigh parameters if requested:
      if ret.iray is not None:
          ifree = ret.map_pars['ray']
          self.rayleigh.pars[ifree] = params[ret.iray]

      # Update cloud parameters if requested:
      if ret.icloud is not None:
          ifree = ret.map_pars['cloud']
          self.cloud.pars[ifree] = params[ret.icloud]
          j = 0
          for model in self.cloud.models:
              model.pars = self.cloud.pars[j:j+model.npars]
              j += model.npars

      # Update patchy-cloud fraction if requested:
      if ret.ipatchy is not None:
          self.cloud.fpatchy = params[ret.ipatchy][0]

      # Update stellar effective temperature:
      if ret.itstar is not None:
          self.phy.tstar = params[ret.itstar][0]
          self.spec.starflux = self.spec.flux_interp(self.phy.tstar)

      # Calculate atmosphere and spectrum:
      self.run()

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
              "The sum of trace abundances' fraction exceeds "
              f"the cap of {ret.qcap:.3f}"
          )

      # Band-integrate spectrum:
      obs.bandflux = self.band_integrate()

      # Instrumental offset:
      if obs.offset_instruments is not None:
          if ret.ioffset is not None:
              ifree = ret.map_pars['offset']
              obs.offset_pars[ifree] = params[ret.ioffset]
          for j, offset in enumerate(obs.offset_pars):
              obs.bandflux[obs.offset_indices[j]] -= offset*obs._dunits

      # update_atm() in self.run() broke:
      if not np.any(obs.bandflux):
          reject_flag = True

      # Reject this iteration if there are invalid temperatures or radii:
      if obs.bandflux is not None and reject_flag:
          obs.bandflux[:] = np.inf

      ret.params = np.copy(params)
      if retmodel:
          return self.spec.spectrum, obs.bandflux

      return obs.bandflux


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

      This method also defines pyrat.atm.radeq_temps, a 2D array
      containing all temperature-profile iterations.
      """
      atm = self.atm

      if nsamples is None:
          nsamples = self.inputs.nsamples

      self.log.verb = 0  # Mute it
      # Enforce two-stream RT:
      rt_path = self.od.rt_path
      self.od.rt_path = 'emission_two_stream'
      tmin = np.amax((self.cs.tmin, self.ex.tmin))
      tmax = np.amin((self.cs.tmax, self.ex.tmax))

      # Initial temperature scale factor
      f_scale = atm._fscale if hasattr(atm,'_fscale') else None

      if hasattr(atm, 'radeq_temps') and continue_run:
          radeq_temps = atm.radeq_temps
          f_scale = None
      else:
          radeq_temps = None

      print("\nRadiative-thermochemical equilibrium calculation:")
      radeq_temps, f_scale = ps.radiative_equilibrium(
          atm.press, atm.temp, nsamples,
          atm.chem_model,
          self.run,
          self.spec.wn,
          self.spec,
          self.atm,
          radeq_temps,
          convection,
          tmin, tmax,
          f_scale,
          self.atm.mplanet, self.atm.mol_mass,
      )
      print("\nDone.")

      # Update last tempertature iteration and save to file:
      atm.temp = radeq_temps[-1]
      io.write_atm(
          self.spec.specfile.replace('.dat','.atm'),
          atm.press, atm.temp, self.atm.species, atm.vmr,
          punits="bar",
          header="# Radiative-thermochemical equilibrium profile.\n\n",
      )
      self.atm._fscale = f_scale
      self.od.rt_path = rt_path
      self.log.verb = self.verb


  def band_integrate(self, spectrum=None):
      """
      Band-integrate transmission spectrum (transit) or planet-to-star
      flux ratio (eclipse) over transmission band passes.
      """
      if self.obs.filters is None:
          return None

      if spectrum is None:
          spectrum = self.spec.spectrum

      rprs_square = (self.atm.rplanet/self.phy.rstar)**2.0
      if self.od.rt_path in pc.emission_rt:
          spectrum = spectrum / self.spec.starflux * rprs_square

      self.obs.bandflux = np.array([
          np.trapz(spectrum[band.idx]*band.response, band.wn)
          for band in self.obs.filters
      ])

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
      # Allocate outputs:
      ec = np.empty((0, self.spec.nwave))
      label = []
      # Line-by-line extinction coefficient:
      if self.ex.nspec != 0:
          e, lab = ex.get_ec(self, layer)
          ec = np.vstack((ec, e))
          label += lab
      # CIA extinction coefficient:
      if self.cs.nfiles != 0:
          e, lab = self.cs.get_ec(self.atm.temp, self.atm.d, layer)
          ec = np.vstack((ec, e))
          label += lab
      # Rayleigh scattering extinction coefficient:
      if self.rayleigh.models != []:
          e, lab = self.rayleigh.get_ec(self.atm.d, layer)
          ec = np.vstack((ec, e))
          label += lab
      # Haze/clouds extinction coefficient:
      if self.cloud.models != []:
          e, lab = cl.get_ec(self, layer)
          ec = np.vstack((ec, e))
          label += lab
      # H- extinction coefficient:
      if self.h_ion.model is not None:
          e, lab = self.h_ion.get_ec(self.atm.temp, self.atm.d, layer)
          ec = np.vstack((ec, e))
          label += [lab]
      # Alkali resonant lines extinction coefficient:
      if self.alkali.models != []:
          e, lab = self.alkali.get_ec(self.atm.temp, self.atm.d, layer)
          ec = np.vstack((ec, e))
          label += lab
      return ec, label


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
      }

      if self.obs.nfilters > 0:
          args['bands_wl0'] = [band.wl0 for band in self.obs.filters]
          args['bands_wl'] = [band.wl for band in self.obs.filters]
          args['bands_response'] = [band.response for band in self.obs.filters]
          args['bands_flux'] = self.obs.bandflux

      if self.ret.spec_low2 is not None and spec != 'model':
          args['bounds'] = [
              self.ret.spec_low2,  self.ret.spec_low1,
              self.ret.spec_high1, self.ret.spec_high2]

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
      if self.ret.posterior is None:
          ax = pp.temperature(
              self.atm.press, profiles=[self.atm.temp], **kwargs)
          return ax

      ax = pp.temperature(
          self.atm.press,
          profiles=[self.ret.temp_median, self.ret.temp_best],
          labels=['median', 'best'],
          bounds=self.ret.temp_post_boundaries, **kwargs)
      return ax


  def __str__(self):
      if self.spec.resolution is not None:
         wave = "R={:.1f}".format(self.spec.resolution)
      elif self.spec.wlstep is not None:
         wave = f'delta-wl={self.spec.wlstep:.2f}'
      else:
         wave = f"delta-wn={self.spec.wnstep:.3f} cm-1"

      opacities = []
      if self.ex.nspec != 0:
          for mol in self.ex.species:
              imol = np.where(self.atm.species == mol)[0][0]
              opacities.append(self.atm.species[imol])
      if self.cs.nfiles > 0:
          for cia in self.cs.models:
              opacities.append(cia.name)
      for rmodel in self.rayleigh.models:
          opacities.append(rmodel.name)
      for cloud in self.cloud.models:
          opacities.append(cloud.name)
      for alkali in self.alkali.models:
          opacities.append(alkali.mol)
      if self.h_ion.model is not None:
          opacities.append(self.h_ion.model.name)

      pmin = self.atm.press[ 0]/pc.bar
      pmax = self.atm.press[-1]/pc.bar
      wlmin = 1.0/(self.spec.wn[-1]*pc.um)
      wlmax = 1.0/(self.spec.wn[ 0]*pc.um)
      return (
          "Pyrat atmospheric model\n"
          f"configuration file:  '{self.inputs.configfile}'\n"
          f"Pressure profile:  {pmin:.2e} -- {pmax:.2e} bar "
          f"({self.atm.nlayers:d} layers)\n"
          f"Wavelength range:  {wlmin:.2f} -- {wlmax:.2f} um "
          f"({self.spec.nwave:d} samples, {wave})\n"
          f"Composition:\n  {self.atm.species}\n"
          f"Opacity sources:\n  {opacities}"
      )

