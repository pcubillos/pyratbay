# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import multiprocessing as mp
from collections import OrderedDict

import numpy as np

from .. import atmosphere as pa
from .. import constants as pc
from .. import io as io
from .. import opacity as op
from .. import plots as pp
from .. import spectrum as ps
from .. import tools as pt

from .  import extinction as ex
from .  import crosssec as cs
from .  import clouds as cl
from .  import alkali as al
from .  import read_atm as ra
from .  import optical_depth as od
from .  import spectrum as sp
from .  import objects as ob
from .  import argum as ar
from .  import makesample as ms
from .  import voigt as v
from .  import read_tli as rtli


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
      self.spec = ob.Spectrum()       # Spectrum data
      self.atm = ob.Atm()             # Modeling atmospheric model
      self.lt = ob.Linetransition()   # Line-transition data
      self.mol = ob.Molecules()       # Molecules data
      self.iso = ob.Isotopes()        # Isotopes data
      self.voigt = ob.Voigt()         # Voigt profile
      self.ex = ob.Extinction()       # Extinction-coefficient
      self.od = ob.Optdepth()         # Optical depth
      self.obs = ob.Observation()     # Observational data
      self.phy = ob.Physics()         # System physical parameters
      self.ret = ob.Retrieval()       # Retrieval variables
      self.timestamps = OrderedDict()

      # Parse config file inputs:
      pt.parse(self, cfile, no_logfile, mute)
      self.inputs.atm = ob.Atm()


  def set_atmosphere(self):
      # Setup time tracker:
      timer = pt.Timer()
      self.timestamps['init'] = timer.clock()

      # Check that user input arguments make sense:
      ar.check_atmosphere(self)  # TBD: Need this here?

      # Read the atmospheric file:
      ra.make_atmosphere(self)
      self.timestamps['read atm'] = timer.clock()

      # TBD: Revise and repurpose this function:
      #ms.make_atmprofiles(self)
      #self.timestamps['atm sample'] = timer.clock()


  def set_spectrum(self):
      timer = pt.Timer()
      # Initialize wavenumber sampling:
      ms.make_wavenumber(self)
      self.timestamps['wn sample'] = timer.clock()

      # Read opacity tables (if needed):
      ex.read_opacity(self, self.spec._wn_mask)
      self.timestamps['read opacity'] = timer.clock()

      ar.check_spectrum(self)

      # Read line database:
      rtli.read_tli(self)
      self.timestamps['read tli'] = timer.clock()

      # Setup more observational/retrieval parameters:
      ar.setup(self)

      # Extinction Voigt grid:
      v.voigt(self)

      # Alkali opacity models:
      self.alkali = al.Alkali(
          self.inputs.model_names,
          self.atm.press,
          self.spec.wn,
          self.inputs.alkali_cutoff,
          self.mol.name,
          self.log,
      )
      self.timestamps['voigt'] = timer.clock()

      # Hydrogen ion opacity:
      self.h_ion = op.Hydrogen_Ion_Opacity(
          1.0/self.spec.wn/pc.um, self.mol.name,
      )
      # At the moment work as an on/off flag, as there's only one model
      self.h_ion.has_opacity &= self.od.h_ion_models is not None

      # CIA opacity:
      self.cs = cs.CIA(
          self.inputs.cia_files,
          self.spec.wn,
          self.mol.name,
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
      status = ra.update_atm(self, temp, vmr, radius)
      if status == 0:
          self.spec.spectrum[:] = 0.0
          return
      # Interpolate CIA absorption:
      self.cs.calc_extinction_coefficient(self.atm.temp, self.atm.d)
      self.timestamps['interp cs'] = timer.clock()

      # Calculate cloud, Rayleigh, and H- absorption:
      cl.absorption(self)
      self.rayleigh.absorption(self.atm.d)
      self.h_ion.absorption(self.atm.temp, self.atm.d)
      self.timestamps['cloud+ray'] = timer.clock()

      # Calculate the alkali absorption:
      self.alkali.calc_extinction_coefficient(self.atm.temp, self.atm.d)
      self.timestamps['alkali'] = timer.clock()

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
          self.phy.rplanet = params[ret.irad][0] * pt.u(atm.runits)
      elif ret.ipress is not None:
          p_ref = 10.0**params[ret.ipress][0] * pc.bar
          self.atm.refpressure = p_ref

      # Update planetary mass if requested:
      if ret.imass is not None:
          self.phy.mplanet = params[ret.imass][0] * pt.u(self.phy.mpunits)

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
          obs.starflux = [
              self.spec.starflux[band.idx]
              for band in obs.filters
          ]

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
          self.phy.mplanet, self.mol.mass,
      )
      print("\nDone.")

      # Update last tempertature iteration and save to file:
      atm.temp = radeq_temps[-1]
      io.write_atm(
          self.spec.specfile.replace('.dat','.atm'),
          atm.press, atm.temp, self.mol.name, atm.vmr,
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
      if self.obs.bandtrans is None:
          return None

      if spectrum is None:
          spectrum = self.spec.spectrum
      specwn = self.spec.wn
      bandidx = self.obs.bandidx
      rprs_square = (self.phy.rplanet/self.phy.rstar)**2.0

      if self.od.rt_path in pc.transmission_rt:
          bandtrans = self.obs.bandtrans
      elif self.od.rt_path in pc.emission_rt:
          bandtrans = [
              btrans/sflux * rprs_square
              for btrans, sflux in zip(self.obs.bandtrans, self.obs.starflux)
          ]

      self.obs.bandflux = np.array([
          np.trapz(spectrum[idx]*btrans, specwn[idx])
          for btrans, idx in zip(bandtrans, bandidx)
      ])

      return self.obs.bandflux


  def hydro(self, pressure, temperature, mu, g, mass, p0, r0):
      """
      Hydrostatic-equilibrium driver.
      Depending on self.atm.rmodelname, select between the g=GM/r**2
      (hydro_m) or constant-g (hydro_g) formula to compute
      the hydrostatic-equilibrium radii of the planet layers.

      Parameters
      ----------
      pressure: 1D float ndarray
         Atmospheric pressure for each layer (in barye).
      temperature: 1D float ndarray
         Atmospheric temperature for each layer (in K).
      mu: 1D float ndarray
         Mean molecular mass for each layer (in g mol-1).
      g: Float
         Atmospheric gravity (in cm s-2).
      mass: Float
         Planetary mass (in g).
      p0: Float
         Reference pressure level (in barye) where radius(p0) = r0.
      r0: Float
         Reference radius level (in cm) corresponding to p0.
      """
      if self.atm.rmodelname is None:
          print('No hydrostatic-equilibrium model defined.')
          return None
      # H.E. with  g=GM/r**2:
      elif self.atm.rmodelname == 'hydro_m':
          return pa.hydro_m(pressure, temperature, mu, mass, p0, r0)
      # H.E. with constant g:
      elif self.atm.rmodelname == 'hydro_g':
          return pa.hydro_g(pressure, temperature, mu, g, p0, r0)


  def set_filters(self):
      """
      Set observational variables (pyrat.obs) based on given parameters.
      """
      if self.obs.filters is None:
          return

      bandidx = []  # Filter wavenumber indices
      starflux = []  # Interpolated stellar flux
      bandtrans = []  # Normalized interpolated filter transmission
      bandwn = []  # Band's mean wavenumber
      for passband in self.obs.filters:
          # Resample the filters into the planet wavenumber array:
          passband(wn=self.spec.wn)
          bandidx.append(passband.idx)
          bandtrans.append(passband.response)
          bandwn.append(passband.wn0)
          if self.phy.starflux is not None:
              starflux.append(self.spec.starflux[passband.idx])

      # Per-band variables:
      self.obs.starflux = starflux
      self.obs.bandidx = bandidx
      self.obs.bandtrans = bandtrans
      self.obs.bandwn = np.asarray(bandwn)
      self.obs.bandflux = np.zeros(self.obs.nfilters, np.double)


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
      # Cross-section extinction coefficient:
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
      # H- opacity:
      if self.h_ion.has_opacity:
          e, lab = self.h_ion.get_ec(self.atm.temp[layer], self.atm.d[layer])
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
      pyrat_args = {
          'data': self.obs.data,
          'uncert': self.obs.uncert,
          'bandtrans': self.obs.bandtrans,
          'bandidx': self.obs.bandidx,
          'starflux': self.spec.starflux,
          'logxticks': self.inputs.logxticks,
          'yran': self.inputs.yran,
      }

      wavelength = 1.0/(self.spec.wn*pc.um)
      if self.obs.bandwn is not None:
          pyrat_args['bandwl'] = 1.0/(self.obs.bandwn*pc.um)

      if self.obs.bandflux is not None:
          pyrat_args['bandflux'] = self.obs.bandflux

      if self.ret.spec_low2 is not None and spec != 'model':
          pyrat_args['bounds'] = [
              self.ret.spec_low2,  self.ret.spec_low1,
              self.ret.spec_high1, self.ret.spec_high2]

      if spec == 'model':
          pyrat_args['label'] = 'model'
          spectrum = self.spec.spectrum
      elif spec == 'best':
          pyrat_args['label'] = 'best-fit model'
          pyrat_args['bandflux'] = self.ret.bestbandflux
          spectrum = self.ret.spec_best
      elif spec == 'median':
          pyrat_args['label'] = 'median model'
          pyrat_args['bandflux'] = self.band_integrate(spectrum)
          spectrum = self.ret.spec_median
      else:
          print(
              "Invalid 'spec'.  Select from 'model' (default), 'best', "
              "or 'median'."
          )
          return

      if self.phy.rplanet is not None and self.phy.rstar is not None:
          pyrat_args['rprs'] = self.phy.rplanet/self.phy.rstar

      # kwargs can overwite any of the previous value:
      pyrat_args.update(kwargs)

      if self.od.rt_path in pc.transmission_rt:
          rt_path = 'transit'
      elif self.od.rt_path in pc.emission_rt:
          rt_path = 'eclipse'
      ax = pp.spectrum(spectrum, wavelength, rt_path, **pyrat_args)
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
              imol = np.where(self.mol.name == mol)[0][0]
              opacities.append(self.mol.name[imol])
      if self.cs.nfiles > 0:
          for cia in self.cs.models:
              opacities.append(cia.name)
      for rmodel in self.rayleigh.models:
          opacities.append(rmodel.name)
      for cloud in self.cloud.models:
          opacities.append(cloud.name)
      for alkali in self.alkali.models:
          opacities.append(alkali.mol)
      if self.h_ion.has_opacity:
          opacities.append(self.h_ion.name)

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
          f"Composition:\n  {self.mol.name}\n"
          f"Opacity sources:\n  {opacities}"
      )

