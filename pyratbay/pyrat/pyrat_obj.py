# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import multiprocessing as mp
from collections import OrderedDict

import numpy  as np

from .. import constants  as pc
from .. import atmosphere as pa
from .. import tools      as pt
from .. import plots      as pp
from .. import io         as io

from .  import extinction as ex
from .  import crosssec   as cs
from .  import rayleigh   as ray
from .  import clouds     as cl
from .  import alkali     as al
from .  import readatm    as ra
from .  import optdepth   as od
from .  import spectrum   as sp
from .  import objects    as ob
from .  import argum      as ar
from .  import makesample as ms
from .  import voigt      as v
from .  import readlinedb as rl


class Pyrat(object):
  """
  Main Pyrat object.
  """
  def __init__(self, cfile):
      """
      Parse the command-line arguments into the pyrat object.

      Parameters
      ----------
      args: Namespace
          Object storing user-input attributes to initialize Pyrat.
      log: Log object
          An mc3.utils.Log instance to log screen outputs to file.

      Examples
      --------
      >>> import pyratbay as pb
      >>> # There are two ways to initialize a Pyrat object:
      >>> # Initialize only:
      >>> pyrat = pb.init('spectrum_transmission.cfg')
      >>> # This is equivalent to:
      >>> args, log = pb.tools.parse('spectrum_transmission.cfg')
      >>> pyrat = pb.Pyrat(args, log)

      >>> # Initialize and compute a spectrum:
      >>> pyrat = pb.run('spectrum_transmission.cfg')
      """
      # Sub-classes:
      self.spec     = ob.Spectrum()        # Spectrum data
      self.atm      = ob.Atm()             # Modeling atmospheric model
      self.lt       = ob.Linetransition()  # Line-transition data
      self.mol      = ob.Molecules()       # Molecules data
      self.iso      = ob.Isotopes()        # Isotopes data
      self.voigt    = ob.Voigt()           # Voigt profile
      self.ex       = ob.Extinction()      # Extinction-coefficient
      self.cs       = ob.Cross()           # Cross-section extinction
      self.od       = ob.Optdepth()        # Optical depth
      self.cloud    = ob.Cloud()           # Cloud models
      self.rayleigh = ob.Rayleigh()        # Rayleigh models
      self.alkali   = ob.Alkali()          # Alkali opacity models
      self.obs      = ob.Observation()     # Observational data
      self.phy      = ob.Physics()         # System physical parameters
      self.ret      = ob.Retrieval()       # Retrieval variables
      self.timestamps = OrderedDict()

      # Parse config file inputs:
      pt.parse(self, cfile)
      self.inputs.atm = ob.Atm()


  def setup_spectrum(self):
      # Setup time tracker:
      timer = pt.clock()

      # Check that user input arguments make sense:
      ar.check_spectrum(self)
      self.timestamps['init'] = next(timer)

      # Initialize wavenumber sampling:
      ms.make_wavenumber(self)
      self.timestamps['wn sample'] = next(timer)

      # Read the atmospheric file:
      ra.read_atm(self)
      self.timestamps['read atm'] = next(timer)

      # Read line database:
      rl.read_tli(self)
      self.timestamps['read tli'] = next(timer)

      # Make atmospheric profiles (pressure, radius, temperature, abundances):
      ms.make_atmprofiles(self)
      self.timestamps['atm sample'] = next(timer)

      # Setup more observational/retrieval parameters:
      ar.setup(self)

      # Extinction Voigt grid:
      v.voigt(self)
      # Alkali Voigt grid:
      al.init(self)
      self.timestamps['voigt'] = next(timer)

      # Calculate extinction-coefficient table:
      ex.exttable(self)
      self.timestamps['ext table'] = next(timer)

      # Read CIA files:
      cs.read(self)
      self.timestamps['read cs'] = next(timer)


  def run(self, temp=None, abund=None, radius=None):
      """
      Evaluate a Pyrat spectroscopic model

      Parameters
      ----------
      pyrat: A Pyrat instance
      temp: 1D float ndarray
          Updated atmospheric temperature profile in Kelvin, of size nlayers.
      abund: 2D float ndarray
          Updated atmospheric abundances profile by number density, of
          shape [nlayers, nmol].
      radius: 1D float ndarray
          Updated atmospheric altitude profile in cm, of size nlayers.
      """
      timer = pt.clock()

      # Re-calculate atmospheric properties if required:
      status = ra.reloadatm(self, temp, abund, radius)
      if status == 0:
          return
      # Interpolate CIA absorption:
      cs.interpolate(self)
      self.timestamps['interp cs'] = next(timer)

      # Calculate cloud and Rayleigh absorption:
      cl.absorption(self)
      ray.absorption(self)
      self.timestamps['cloud+ray'] = next(timer)

      # Calculate the alkali absorption:
      al.absorption(self)
      self.timestamps['alkali'] = next(timer)

      # Calculate the optical depth:
      od.opticaldepth(self)
      self.timestamps['odepth'] = next(timer)

      # Calculate the spectrum:
      sp.spectrum(self)
      self.timestamps['spectrum'] = next(timer)

      self.log.msg("\nTimestamps (s):\n" +
                   "\n".join("{:10s}: {:10.6f}".format(key,val)
                             for key,val in self.timestamps.items()))

      if len(self.log.warnings) > 0:
          # Write all warnings to file:
          wpath, wfile = os.path.split(self.log.logname)
          wfile = "{:s}/warnings_{:s}".format(wpath, wfile)
          with open(wfile, "w") as f:
              f.write("Warnings log:\n\n{:s}\n".format(self.log.sep))
              f.write("\n\n{:s}\n".format(self.log.sep).join(self.log.warnings))
          # Report it:
          self.log.head("\n{:s}\n  There were {:d} warnings raised.  "
              "See '{:s}'.\n{:s}".format(self.log.sep, len(self.log.warnings),
                                         wfile, self.log.sep))


  def eval(self, params, retmodel=True, verbose=False):
      """
      Fitting routine for MCMC.

      Parameters
      ----------
      params: 1D float iterable
         Array of fitting parameters that define the atmosphere.
      retmodel: Bool
         Flag to include the model spectra in the return.
      verbose: Bool
         Flag to print out if a run failed.

      Returns
      -------
      spectrum: 1D float ndarray
         The output model spectra.  Returned only if retmodel=True.
      bandflux: 1D float ndarray
         The waveband-integrated spectrum values.
      """
      params = np.asarray(params)
      q0 = np.copy(self.atm.qbase)

      if self.ret.imass is not None:
          self.phy.mplanet = params[self.ret.imass][0] * pt.u(self.phy.mpunits)

      rejectflag = False
      # Update temperature profile if requested:
      if self.ret.itemp is not None:
          temp = self.atm.tmodel(params[self.ret.itemp], *self.atm.targs)
      else:
          temp = self.atm.temp
      # Turn-on reject flag if temperature is out-of-bounds:
      if np.any(temp < self.ret.tlow) or np.any(temp > self.ret.thigh):
          temp[:] = 0.5*(self.ret.tlow + self.ret.thigh)
          rejectflag = True
          if verbose:
              self.log.warning("Input temperature profile runs out of "
                               "boundaries ({:.1f--{:.1f}} K)".
                               format(self.ret.tlow, self.ret.thigh))

      # Update abundance profiles if requested:
      if self.ret.imol is not None:
          q2 = pa.qscale(q0, self.mol.name, self.atm.molmodel,
                         self.atm.molfree, params[self.ret.imol],
                         self.atm.bulk,
                         iscale=self.atm.ifree, ibulk=self.atm.ibulk,
                         bratio=self.atm.bulkratio, invsrat=self.atm.invsrat)
      else:
          q2 = self.atm.q

      # Check abundaces stay within bounds:
      if pa.qcapcheck(q2, self.ret.qcap, self.atm.ibulk):
          rejectflag = True
          if verbose:
              self.log.warning("The sum of trace abundances' fraction exceeds "
                               "the cap of {:.3f}.".format(self.ret.qcap))

      # Update reference radius if requested:
      if self.ret.irad is not None:
          self.phy.rplanet = params[self.ret.irad][0] * pt.u(self.atm.runits)

      # Update Rayleigh parameters if requested:
      if self.ret.iray is not None:
          j = 0
          rpars = params[self.ret.iray]
          for rmodel in self.rayleigh.models:
              rmodel.pars = rpars[j:j+rmodel.npars]
              j += rmodel.npars

      # Update cloud parameters if requested:
      if self.ret.icloud is not None:
          j = 0
          pars = params[self.ret.icloud]
          for model in self.cloud.models:
              model.pars = pars[j:j+model.npars]
              j += model.npars

      # Update patchy-cloud fraction if requested:
      if self.ret.ipatchy is not None:
          self.cloud.fpatchy = params[self.ret.ipatchy]

      # Calculate spectrum:
      self.run(temp=temp, abund=q2)

      # Band-integrate spectrum:
      self.obs.bandflux = self.band_integrate()

      # Reject this iteration if there are invalid temperatures or radii:
      if self.obs.bandflux is not None and rejectflag:
          self.obs.bandflux[:] = np.inf

      if retmodel:
          return self.spec.spectrum, self.obs.bandflux

      return self.obs.bandflux


  def band_integrate(self, spectrum=None):
      """
      Band-integrate transmission spectrum (transit) or planet-to-star
      flux ratio (eclipse) over transmission band passes.
      """
      if self.obs.bandtrans is None:
          return None

      if spectrum is None:
          spectrum = self.spec.spectrum
      specwn   = self.spec.wn
      bandidx  = self.obs.bandidx

      if self.od.path == 'transit':
          bandtrans = self.obs.bandtrans
      elif self.od.path == 'eclipse':
          bandtrans = [btrans/sflux * self.phy.rprs**2
               for btrans, sflux in zip(self.obs.bandtrans, self.obs.starflux)]

      self.obs.bandflux = np.array([np.trapz(spectrum[idx]*btrans, specwn[idx])
                                    for btrans, idx in zip(bandtrans, bandidx)])

      return self.obs.bandflux


  def hydro(self, pressure, temperature, mu, g, mass, p0, r0):
      """
      Hydrostatic-equilibrium driver.
      Depending on the self.atm.hydrom flag, select between the g=GM/r**2
      (hydrom=True) or constant-g (hydrom=False) formula to compute
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
      # H.E. with  g=GM/r**2:
      if self.atm.hydrom:
          return pa.hydro_m(pressure, temperature, mu, mass, p0, r0)
      # H.E. with constant g:
      return pa.hydro_g(pressure, temperature, mu, g, p0, r0)


  def set_filters(self):
      """
      Set observational variables (pyrat.obs) based on given parameters.
      """
      if self.obs.filters is None:
          return

      bandidx   = []  # Filter wavenumber indices
      starflux  = []  # Interpolated stellar flux
      bandtrans = []  # Normalized interpolated filter transmission
      bandwn    = []  # Band's mean wavenumber
      for filter in self.obs.filters:
          # Read filter wavenumber and transmission curves:
          filterwn, filtertr = io.read_spectrum(filter)
          # Resample the filters into the planet wavenumber array:
          btrans, bidx = pt.resample(filtertr, filterwn, self.spec.wn,
              normalize=True)
          bandidx.append(bidx)
          bandtrans.append(btrans)
          bandwn.append(np.sum(filterwn*filtertr)/np.sum(filtertr))
          if self.phy.starflux is not None:
              starflux.append(self.spec.starflux[bidx])

      # Per-band variables:
      self.obs.bandidx   = bandidx
      self.obs.bandtrans = bandtrans
      self.obs.starflux  = starflux
      self.obs.bandwn    = np.asarray(bandwn)
      self.obs.bandflux  = np.zeros(self.obs.nfilters, np.double)


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
      if self.ex.nmol != 0:
          e, lab = ex.get_ec(self, layer)
          ec = np.vstack((ec, e))
          label += lab
      # Cross-section extinction coefficient:
      if self.cs.nfiles != 0:
          e, lab = cs.interpolate(self, layer)
          ec = np.vstack((ec, e))
          label += lab
      # Rayleigh scattering extinction coefficient:
      if self.rayleigh.models != []:
          e, lab = ray.get_ec(self, layer)
          ec = np.vstack((ec, e))
          label += lab
      # Haze/clouds extinction coefficient:
      if self.cloud.models != []:
          e, lab = cl.get_ec(self, layer)
          ec = np.vstack((ec, e))
          label += lab
      # Alkali resonant lines extinction coefficient:
      if self.alkali.models != []:
          e, lab = al.get_ec(self, layer)
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

      with mp.Pool(self.ncpu) as pool:
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


  def plot_spectrum(self, spec='model', logxticks=None, gaussbin=2.0,
                    yran=None, filename=None):
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
      """
      wavelength = 1.0/(self.spec.wn*pc.um)
      if self.obs.bandwn is not None:
          bandwl = 1.0/(self.obs.bandwn*pc.um)
      else:
          bandwl = None
      if self.obs.bandtrans is not None and np.all(self.obs.bandflux==0):
          bandflux = self.band_integrate()
      else:
          bandflux = None

      if logxticks is None:
          logxticks = self.inputs.logxticks
      if yran is None:
          yran = self.inputs.yran

      bounds = None
      if self.ret.spec_low2 is not None and spec != 'model':
          bounds = [self.ret.spec_low2,  self.ret.spec_low1,
                    self.ret.spec_high1, self.ret.spec_high2]

      if spec == 'model':
          label = 'model'
          spectrum = self.spec.spectrum
      elif spec == 'best':
          label = 'best-fit model'
          spectrum = self.ret.spec_best
          bandflux = self.ret.bestbandflux
      elif spec == 'median':
          label = 'median model'
          spectrum = self.ret.spec_median
          bandflux = self.band_integrate(spectrum)
      else:
          print("Invalid 'spec'.  Select from 'model' (default), 'best', "
                "or 'median'.")
          return
      pp.spectrum(spectrum, wavelength, self.od.path,
          self.obs.data, self.obs.uncert, bandwl, bandflux,
          self.obs.bandtrans, self.obs.bandidx,
          self.spec.starflux, self.phy.rprs, label, bounds,
          logxticks, gaussbin, yran, filename)


  def plot_posterior_pt(self, filename=None):
      """Plot posterior distribution of PT profile."""
      if self.ret.posterior is None:
          print('pyrat objec does not have a posterior distribution.')
          return

      posterior = self.ret.posterior
      ifree = self.ret.pstep[self.ret.itemp] > 0
      itemp = np.arange(np.sum(ifree))
      if filename is None:
          outfile = os.path.splitext(os.path.basename(self.ret.mcmcfile))[0]
          filename = '{:s}_posterior_PT_profile.png'.format(outfile)
      pp.posterior_pt(posterior[:,itemp], self.atm.tmodel, self.atm.targs,
          self.ret.params[self.ret.itemp], ifree, self.atm.press,
          self.ret.bestp[self.ret.itemp], filename)


  def __str__(self):
      if self.spec.resolution is not None:
         wave = "R={:.1f}".format(self.spec.resolution)
      else:
         wave = "dwn={:.3f} cm-1".format(self.spec.wnstep)

      opacities = []
      if self.ex.nmol != 0:
          for molID in self.ex.molID:
              imol = np.where(self.mol.ID == molID)[0][0]
              opacities.append(self.mol.name[imol])
      if self.cs.nfiles != 0:
          for molecs in self.cs.molecules:
              if len(molecs) == 2:
                  opacities.append('CIA ' + '-'.join(molecs))
              else:
                  opacities.append(molecs[0])
      for rmodel in self.rayleigh.models:
          opacities.append(rmodel.name)
      for cloud in self.cloud.models:
          opacities.append(cloud.name)
      for alkali in self.alkali.models:
          opacities.append(alkali.mol)

      return ("Pyrat atmospheric model\n"
          "configuration file:  '{:s}'\n"
          "Pressure profile (bar):  {:.2e} -- {:.2e} ({:d} layers)\n"
          "Wavelength range (um):  {:.2f} -- {:.2f} ({:d} samples, {:s})\n"
          "Composition:  {}\n"
          "Opacity sources:  {}".format(
          self.inputs.configfile,
          self.atm.press[ 0]/pc.bar,
          self.atm.press[-1]/pc.bar,
          self.atm.nlayers,
          1.0/(self.spec.wn[-1]*pc.um),
          1.0/(self.spec.wn[ 0]*pc.um),
          self.spec.nwave,
          wave,
          self.mol.name,
          opacities))

