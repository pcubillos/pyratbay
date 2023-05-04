# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'Extinction',
    'Optdepth',
    'Cloud',
    'Physics',
    'Retrieval',
]

import numpy as np

from .. import tools as pt
from .. import constants as pc
from .. import atmosphere as pa


class Extinction(object):
  def __init__(self):
      self.ec      = None # Molecular line-transition extinction coefficient
                          #  in cm-1 [nlayers, nwave]
      self.ethresh = None # Extinction-coefficient threshold
      self.extfile = None # Extinction-coefficient table filename
      self.etable  = None # Tabulated extinction coefficient (cm-2 molecule-1)
                          # with shape [nmol, nlayer, ntemp, nwave]
      self.tmin    = None # Minimum temperature to sample
      self.tmax    = None # Maximum temperature to sample
      self.tstep   = None # Temperature-sample step interval
      self.z       = None # Partition function at tabulated temperatures
                          #   [niso, ntemp]
      self.nspec   = 0    # Number of species
      self.ntemp   = 0    # Number of temperature samples
      self.nlayers = 0    # Number of pressure layers
      self.nwave   = 0    # Number of wavenumber spectral samples

      self.temp    = None # Tabulated temperatures
      self.press   = None # Tabulated pressures
      self.wn      = None # Tabulated wavenumber

  def __str__(self):
      fmt = {'float': '{:.2e}'.format}
      fw = pt.Formatted_Write()
      fw.write('Extinction-coefficient information:')
      fw.write('Line-transition strength threshold (ethresh): {:.2e}',
          self.ethresh)
      if self.ec is not None:
          fw.write(
              '\nLBL extinction coefficient for the atmospheric model '
              '(ec, cm-1) [layer, wave]:\n{}', self.ec, fmt=fmt)
      extfile = ['None'] if self.extfile is None else self.extfile
      fw.write("Extinction-coefficient table filename(s) (extfile): {}",
          '\n    '.join(extfile))
      if self.extfile is None:
          return fw.text
      fw.write('Minimum temperature (tmin, K): {:6.1f}', self.tmin)
      fw.write('Maximum temperature (tmax, K): {:6.1f}', self.tmax)
      fw.write('Temperature sampling interval (tstep, K): {:6.1f}', self.tstep)
      fw.write('\nNumber of species (nspec):          {:5d}', self.nspec)
      fw.write('Number of temperatures (ntemp):     {:5d}', self.ntemp)
      fw.write('Number of layers (nlayers):         {:5d}', self.nlayers)
      fw.write('Number of spectral samples (nwave): {:5d}', self.nwave)
      fw.write('\nSpecies array (species): {}', self.species)
      fw.write('Temperature array (temp, K):\n   {}', self.temp)
      fw.write('Partition function (z): {}', self.z)
      fw.write('Pressure array (press, bar):\n   {}', self.press/pc.bar,
               fmt=fmt, lw=8000)
      fw.write('Wavenumber array (wn, cm-1):\n    {}', self.wn, prec=4)
      fw.write('Tabulated extinction coefficient (etable, cm2 molecule-1) '
               'of shape\n    [nmol, ntemp, nlayers, nwave]:\n{}', self.etable,
               fmt=fmt, edge=2)
      return fw.text


class Cloud(object):
    """Interface to collect all Cloud opacity models"""
    def __init__(self, model_names, pars, fpatchy=None, log=None):
        self.model_names = model_names
        self.models = []    # List of cloud models
        self.pnames = []
        self.texnames = []
        self.npars = 0
        self.pars = None
        self.ec = None  # Cloud extinction coefficient
        self.fpatchy = fpatchy

        if model_names is None:
            return
        for name in model_names:
            model = pa.clouds.get_model(name)
            self.models.append(model)
            self.npars += model.npars
            self.pnames += model.pnames
            self.texnames += model.texnames

        # Parse parameters:
        if pars is None:
            return
        self.pars = pars
        input_npars = len(self.pars)
        if self.npars != input_npars:
            log.error(
                f'Number of input cloud parameters ({input_npars}) '
                'does not match the number of required '
                f'model parameters ({self.npars})'
            )
        j = 0
        for model in self.models:
            model.pars = self.pars[j:j+model.npars]
            j += model.npars

    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write('Cloud-opacity models (models):')
        for model in self.models:
            fw.write('\n' + str(model))
        fw.write('\nPatchiness fraction (fpatchy): {:.3f}', self.fpatchy)
        fw.write('Total atmospheric cloud extinction-coefficient '
                 '(ec, cm-1):\n{}', self.ec, fmt={'float':' {:.3e}'.format})
        return fw.text


class Optdepth():
  def __init__(self):
      self.maxdepth = None  # Maximum optical depth to calculate
      self.rt_path  = None  # Radiative-transfer observing geometry
      self.ec       = None  # Total extinction coefficient [nlayers, nwave]
      self.epatchy  = None  # Cloudy extinction coefficient for patchy model
      self.raypath  = []    # Distance along ray path  [nlayers]
      self.depth    = None  # Optical depth at raypath [nlayers, nwave]
      self.pdepth   = None  # Cloudy optical depth for patchy model
      self.B        = None  # Blackbody Planck emission [nlayers, nwave]
      self.ideep    = None  # Layer index where depth reached maxdepth [nwave]


  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write('Optical depth information:')
      fw.write('Observing geometry (rt_path): {}', self.rt_path)
      if self.ec is not None:
          fw.write(
              'Total atmospheric extinction coefficient (ec, cm-1) [layer'
              ', wave]:\n{}',
              self.ec,
              fmt={'float':'{: .3e}'.format},
          )
      if self.depth is None:
          fw.write(
              '\nMaximum optical depth to calculate (maxdepth): {:.2f}',
              self.maxdepth,
          )
          return fw.text

      ideepest = np.amax(self.ideep)
      if self.rt_path in pc.transmission_rt:
          fw.write(
              '\nDistance along the ray path across each layer '
              '(outside-in) at each impact parameter (raypath, km):',
          )
          with np.printoptions(formatter={'float':'{:.1f}'.format},threshold=6):
              fw.write('    IP[{:3d}]: {}', 1, self.raypath[1]/pc.km)
              fw.write('    IP[{:3d}]: {}', 2, self.raypath[2]/pc.km)
              fw.write('    IP[{:3d}]: {}', 3, self.raypath[3]/pc.km)
              fw.write('    ...')
              fw.write('    IP[{:3d}]: {}', len(self.raypath)-1,
                       self.raypath[len(self.raypath)-1]/pc.km)
          od_text = (
              '\nOptical depth at each impact parameter, down to '
              'max(ideep) (depth):'
          )
      elif self.rt_path in pc.emission_rt:
          fw.write(
              '\nDistance across each layer along a normal ray path '
              '(raypath, km):\n    {}',
              self.raypath/pc.km,
              fmt={'float':'{:.1f}'.format},
              edge=4,
          )
          od_text = (
              '\nOptical depth at each layer along a normal ray '
              'path into the planet, down to max(ideep) (depth):'
          )

      fw.write(
          '\nMaximum optical depth to calculate (maxdepth): {:.2f}',
          self.maxdepth,
      )
      fw.write(
          'Layer index where the optical depth reaches maxdepth (ideep):'
          '\n    {}',
          self.ideep,
          fmt={'int': '{:3d}'.format},
          edge=7,
      )
      fw.write('Maximum ideep (deepest layer reaching maxdepth): {}', ideepest)

      if self.rt_path in pc.emission_rt:
          fw.write(
              '\nPlanck emission down to max(ideep) (B, erg s-1 cm-2 '
              'sr-1 cm):\n{}',
              self.B[0:ideepest+1],
              fmt={'float':'{: .3e}'.format},
          )

      fw.write(
          '{}\n{}', od_text, self.depth[0:ideepest+1],
          fmt={'float':'{: .3e}'.format},
      )
      return fw.text


class Retrieval():
  def __init__(self, inputs, log):
      self.retflag = None
      self.nparams = 0     # Number of free parameters
      self.itemp   = None  # Temperature-model parameter indices
      self.irad    = None  # Reference-radius model parameter index
      self.ipress  = None  # Reference-pressuew model parameter index
      self.imol    = None  # Abundance-model parameter indices
      self.iray    = None  # Rayleigh-model parameter indices
      self.icloud  = None  # Cloud-model parameter indices
      self.ipatchy = None  # Patchy-model parameter index
      self.imass   = None
      self.itstar = None
      self.ioffset = None
      self.posterior = None
      self.bestp     = None
      self.spec_best = None
      self.spec_low1 = None
      self.spec_low2 = None
      self.spec_high1 = None
      self.spec_high2 = None
      self.pnames   = []   # Model parameter names (screen)
      self.texnames = []   # Model parameter names (figures)

      self.mcmcfile = inputs.mcmcfile
      self.retflag = inputs.retflag
      self.qcap = inputs.qcap
      # Lower and upper temperature retrieval boundaries
      self.tlow = inputs.tlow
      self.thigh = inputs.thigh

      self.sampler = inputs.sampler
      # MCMC options
      self.nsamples = inputs.nsamples
      self.burnin = inputs.burnin
      self.thinning = inputs.thinning
      self.nchains = inputs.nchains
      self.grbreak = inputs.grbreak
      self.grnmin = inputs.grnmin

      # Overrides retflag. At some point this will be the only way.
      if inputs.retrieval_params is not None:
          pars = [
              par for par in inputs.retrieval_params.splitlines()
              if par != ''
          ]
          nparams = len(pars)
          pnames = []
          params = np.zeros(nparams)
          pmin = np.tile(-np.inf, nparams)
          pmax = np.tile(np.inf, nparams)
          pstep = np.zeros(nparams)
          prior = np.zeros(nparams)
          priorlow = np.zeros(nparams)
          priorup = np.zeros(nparams)
          for i,par in enumerate(pars):
              fields = par.split()
              nfields = len(fields)
              if nfields not in [2, 5, 7, 8]:
                  log.error(
                      "Invalid number of fields for retrieval_params entry\n"
                      f"'{par}'"
                  )
              pnames.append(fields[0])
              params[i] = fields[1]
              if nfields == 2:
                  continue
              pmin[i] = fields[2]
              pmax[i] = fields[3]
              pstep[i] = fields[4]
              if nfields == 5:
                  continue
              prior[i] = fields[5]
              priorlow[i] = fields[6]
              if nfields == 7:
                  priorup[i] = fields[6]
              else:
                  priorup[i] = fields[7]
          self.pnames = pnames
          self.params = params
          self.pstep = pstep
          self.pmin = pmin
          self.pmax = pmax
          self.prior = prior
          self.priorlow = priorlow
          self.priorup = priorup
      else:
          self.params = inputs.params
          self.pstep = inputs.pstep
          self.pmin = inputs.pmin
          self.pmax = inputs.pmax
          self.prior = inputs.prior
          self.priorlow = inputs.priorlow
          self.priorup = inputs.priorup


  def __str__(self):
      fw = pt.Formatted_Write()
      fw.write('Retrieval information:')
      if self.params is None:
          fw.write('No retrieval parameters set.')
          return fw.text

      pmin = [None for _ in self.params] if self.pmin is None else self.pmin
      pmax = [None for _ in self.params] if self.pmax is None else self.pmax
      psteps = [None for _ in self.params] if self.pstep is None else self.pstep

      fw.write(
          '  Parameter name        value        pmin        pmax       pstep',
      )
      fw.write(
          '  {:15}  {:>10}  {:>10}  {:>10}  {:>10}',
          '(pnames)', '(params)', '(pmin)', '(pmax)', '(pstep)',
      )
      n_obs = len(self.pnames)

      for i in range(n_obs):
          fw.write(
              '  {:15s}  {:10.3e}  {:10.3e}  {:10.3e}  {:10.3e}',
              self.pnames[i], self.params[i], pmin[i], pmax[i], psteps[i],
          )

      if self.prior is not None:
          fw.write('\nParameter name     Prior')
          for i in range(n_obs):
              pname = self.pnames[i]
              val = self.params[i]
              if psteps[i] == 0.0:
                  fw.write('  {:15s}  Fixed at  {:10.3e}', pname, val)
              elif psteps[i] < 0.0:
                  j = -int(psteps[i]) - 1
                  fw.write('  {:15s}  Shared with  {}', pname, self.pnames[j])
              elif self.priorlow[i]==0 and self.priorup[i]==0:
                  fw.write(
                      '  {:15s}  Uniform between     [{:10.3e}, {:10.3e}]',
                      pname, pmin[i], pmax[i],
                  )
              else:
                  fw.write(
                      '  {:15s}  Gaussian  {:10.3e} -{:.3e}  {:+.3e}',
                      pname, self.prior[i], self.priorlow[i], self.priorup[i],
                  )

      fw.write('\nRetrieval algorithm (sampler): {}', self.sampler)
      if self.sampler is None:
          return fw.text
      fw.write('Number of retrieval samples (nsamples): {:,}', self.nsamples)
      # if self.sampler == 'snooker':
      fw.write('Number of parallel chains (nchains):   {}', self.nchains)
      fw.write('Number of burned-in samples (burnin):  {:,}', self.burnin)
      fw.write('Thinning factor (thinning): {}', self.thinning)
      if self.grbreak > 0.0:
          fw.write(
              'Gelman-Rubin convergence criterion to stop (grbreak): '
              f'{self.grbreak}',
          )
          if self.grnmin > 1:
              fw.write(
                  'Minimum number of samples before GR stop (grnmin): '
                  f'{int(self.grnmin)}'
              )
          else:
              fw.write(
                  'Minimum fraction of samples before GR stop (grnmin): '
                  f'{self.grnmin}'
              )

      fw.write(
          f'\nUpper boundary for sum of metal abundances (qcap): {self.qcap}',
      )
      fw.write(f'Temperature upper boundary (tlow, K):  {self.tlow:6.1f}')
      fw.write(f'Temperature lower boundary (thigh, K): {self.thigh:6.1f}')

      fw.write('\nRetrieval posterior file (mcmcfile): {}', self.mcmcfile)
      if self.posterior is not None:
          nsamples, nparams = np.shape(self.posterior)
          fw.write(
              '\nParameter name     Best-fit   Posterior distribution '
              f'of shape [{nsamples},{nparams}]\n'
              '                   (bestp)    (posterior)',
          )
          post = iter(self.posterior.T)
          for pname, bestp, pstep in zip(self.pnames, self.bestp, psteps):
              if pstep > 0:
                  fw.write(
                      '  {:15} {:10.3e}  {}', pname, bestp, next(post),
                      fmt={'float':'{: .3e}'.format}, edge=2,
                  )
              else:
                  fw.write('  {:15} {:10.3e}', pname, bestp)
          fw.write(
              '\nBest-fit spectrum (spec_best):\n    {}', self.spec_best,
              fmt={'float':'{: .3e}'.format},
          )
      return fw.text


class Physics(object):
    """Physical properties about the planet and star"""
    def __init__(self, inputs):
        # Stellar properties
        self.tstar = inputs.tstar
        self.rstar = inputs.rstar
        self.mstar = inputs.mstar
        self.log_gstar = inputs.log_gstar

        # Stellar spectrum filename
        self.starspec = inputs.starspec
        self.kurucz = inputs.kurucz
        self.marcs = inputs.marcs
        self.phoenix = inputs.phoenix

        self.starwn = None  # Input stellar wavenumber array
        self.starflux = None  # Input stellar flux spectrum in  FINDME units

    def __str__(self):
        fw = pt.Formatted_Write()
        fw.write('Physical properties information:')

        rstar = pt.none_div(self.rstar, pc.rsun)
        mstar = pt.none_div(self.mstar, pc.msun)
        fw.write(
            '\nStellar effective temperature (tstar, K): {:.1f}',
            self.tstar,
        )
        fw.write('Stellar radius (rstar, Rsun): {:.3f}', rstar)
        fw.write('Stellar mass (mstar, Msun):   {:.3f}', mstar)
        fw.write(
            'Stellar surface gravity (log_gstar, cm s-2): {:.2f}',
            self.log_gstar,
        )
        #fw.write('Planet-to-star radius ratio (rprs):   {:.5f}', rprs)
        if self.starspec is not None:
            fw.write(f"Input stellar spectrum (starspec): '{self.starspec}'")
        elif self.kurucz is not None:
            fw.write(f"Input Kurucz stellar spectrum (kurucz): '{self.kurucz}'")
        elif self.marcs is not None:
            fw.write(f"Input MARCS stellar spectrum (marcs): '{self.marcs}'")
        elif self.phoenix is not None:
            fw.write(
                f"Input PHOENIX stellar spectrum (phoenix): '{self.phoenix}'",
            )
        elif self.starflux is not None:
            fw.write(
                "Input stellar spectrum is a blackbody at Teff = {:.1f} K.",
                self.tstar,
            )
        fw.write('Stellar spectrum wavenumber (starwn, cm-1):\n    {}',
            self.starwn, fmt={'float': '{:10.3f}'.format})
        fw.write('Stellar flux spectrum (starflux, erg s-1 cm-2 cm):\n    {}',
            self.starflux, fmt={'float': '{: .3e}'.format})
        return fw.text
