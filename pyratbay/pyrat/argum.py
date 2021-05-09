# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import multiprocessing as mp

import numpy as np
import scipy.constants   as sc
import scipy.interpolate as si
import scipy.special     as ss

from .. import tools      as pt
from .. import constants  as pc
from .. import spectrum   as ps
from .. import atmosphere as pa
from .. import io         as io


def check_spectrum(pyrat):
  """
  Check that user input arguments make sense.
  """
  # Shortcuts:
  log  = pyrat.log
  phy  = pyrat.phy
  spec = pyrat.spec
  atm  = pyrat.atm
  obs  = pyrat.obs

  # Check that input files exist:
  if pyrat.mol.molfile is None:
      pyrat.mol.molfile = pc.ROOT + 'pyratbay/data/molecules.dat'

  with pt.log_error(log):
      pt.file_exists('atmfile', 'Atmospheric',    pyrat.atm.atmfile)
      pt.file_exists('tlifile', 'TLI',            pyrat.lt.tlifile)
      pt.file_exists('molfile', 'Molecular-data', pyrat.mol.molfile)

  if pyrat.runmode == 'spectrum' and spec.specfile is None:
      log.error('Undefined output spectrum file (specfile).')

  # Compute the Hill radius for the planet:
  if (phy.mstar is not None and phy.mplanet is not None
      and phy.smaxis is not None):
      phy.rhill = phy.smaxis * (phy.mplanet/(3*phy.mstar))**(1.0/3.0)

  # Check Voigt-profile arguments:
  if (pyrat.voigt.dmin is not None and pyrat.voigt.dmax is not None
      and pyrat.voigt.dmax <= pyrat.voigt.dmin):
      log.error('dmax ({:g} cm-1) must be > dmin ({:g} cm-1).'.
                format(pyrat.voigt.dmax, pyrat.voigt.dmin))

  if (pyrat.voigt.lmin is not None and pyrat.voigt.lmax is not None
      and pyrat.voigt.lmax <= pyrat.voigt.lmin):
      log.error('lmax ({:g} cm-1) must be > lmin ({:g} cm-1).'.
                format(pyrat.voigt.lmax, pyrat.voigt.lmin))

  if pyrat.runmode == 'opacity' or pt.isfile(pyrat.ex.extfile) == 0:
      if pyrat.ex.tmin is None:
          log.error('Undefined lower temperature boundary (tmin) for '
                    'extinction-coefficient grid.')
      if pyrat.ex.tmax is None:
          log.error('Undefined upper temperature boundary (tmax) for '
                    'extinction-coefficient grid.')
      if pyrat.ex.tstep is None:
          log.error('Undefined temperature sampling step (tstep) for '
                    'extinction-coefficient grid.')
      if pyrat.lt.tlifile is None:
          log.error('Requested extinction-coefficient table, but there '
                    'are no input TLI files.')

  if pyrat.runmode == 'mcmc':
      if pyrat.od.rt_path in pc.emission_rt:
          if pyrat.phy.rplanet is None or pyrat.phy.rstar is None:
              log.error("Undefined radius ratio (need rplanet and rstar).")
      if pyrat.obs.data is None:
          log.error("Undefined transit/eclipse data (data).")
      if pyrat.obs.uncert is None:
          log.error("Undefined data uncertainties (uncert).")
      if pyrat.obs.filters is None:
          log.error("Undefined transmission filters (filters).")
      if pyrat.ret.retflag == []:
          log.error('Undefined retrieval model flags.  Select from {}.'.
                    format(pc.retflags))
      if pyrat.ret.sampler is None:
          log.error('Undefined retrieval algorithm (sampler).  Select from '
                    '[snooker].')
      if pyrat.ret.nsamples is None:
          log.error('Undefined number of retrieval samples (nsamples).')
      if pyrat.ret.burnin is None:
          log.error('Undefined number of retrieval burn-in samples (burnin).')
      if pyrat.ret.nchains is None:
          log.error('Undefined number of retrieval parallel chains (nchains).')
      if pyrat.ret.params is None:
          log.error('Undefined retrieval fitting parameters (params).')

  # Check cloud models:
  if pyrat.cloud.model_names is not None:
      pyrat.cloud.models = []
      npars = 0
      for name in pyrat.cloud.model_names:
          model  = pa.clouds.get_model(name)
          npars += model.npars
          pyrat.cloud.models.append(model)
      # Parse the cloud parameters:
      if pyrat.cloud.pars is not None:
          if npars != len(pyrat.cloud.pars):
              log.error('Number of input cloud parameters ({:d}) does not '
                        'match the number of required model parameters ({:d}).'.
                        format(len(pyrat.cloud.pars), npars))
          j = 0
          for model in pyrat.cloud.models:
              npars = model.npars
              model.pars = pyrat.cloud.pars[j:j+npars]
              j += npars

  # Check Rayleigh models:
  if pyrat.rayleigh.model_names is not None:
      pyrat.rayleigh.models = []
      npars = 0
      for name in pyrat.rayleigh.model_names:
          model = pa.rayleigh.get_model(name)
          npars += model.npars
          pyrat.rayleigh.models.append(model)
      # Process the Rayleigh parameters:
      if npars == 0 and pyrat.rayleigh.pars is None:
          pyrat.rayleigh.pars = []
      if pyrat.rayleigh.pars is not None:
          if npars != len(pyrat.rayleigh.pars):
              log.error('Number of input Rayleigh parameters ({:d}) does not '
                        'match the number of required model parameters ({:d}).'.
                        format(len(pyrat.rayleigh.pars), npars))
          j = 0
          for model in pyrat.rayleigh.models:
              npars = model.npars
              model.pars = pyrat.rayleigh.pars[j:j+npars]
              j += npars

  # Check alkali arguments:
  if pyrat.alkali.model_names is not None:
      pyrat.alkali.models = [
          pa.alkali.get_model(name, pyrat.alkali.cutoff)
          for name in pyrat.alkali.model_names]

  # Accept ray-path argument:
  print(pyrat.od)
  if pyrat.runmode in ['spectrum', 'mcmc'] and pyrat.od.rt_path is None:
      log.error(
          "Undefined radiative-transfer observing geometry (rt_path)."
          f"  Select from {pc.rt_paths}.")

  if 'temp' in pyrat.ret.retflag and atm.tmodelname is None:
      log.error('Requested temp in retflag, but there is no tmodel.')
  if 'mol' in pyrat.ret.retflag:
      if atm.molmodel is None:
          log.error("Requested mol in retflag, but there is no 'molmodel'.")
      if atm.bulk is None:
          log.error('Requested mol in retflag, but there are no bulk species.')
  if 'ray' in pyrat.ret.retflag and pyrat.rayleigh.models == []:
      log.error('Requested ray in retflag, but there are no rayleigh models.')
  if 'cloud' in pyrat.ret.retflag and pyrat.cloud.models == []:
      log.error('Requested cloud in retflag, but there are no cloud models.')

  # Check system arguments:
  if pyrat.od.rt_path in pc.transmission_rt and phy.rstar is None:
      log.error('Undefined stellar radius (rstar), required for transmission '
                'calculation.')

  # Check raygrid:
  if spec.raygrid[0] != 0:
      log.error('First angle in raygrid must be 0.0 (normal to surface).')
  if np.any(spec.raygrid < 0) or np.any(spec.raygrid > 90):
      log.error('raygrid angles must lie between 0 and 90 deg.')
  if np.any(np.ediff1d(spec.raygrid) <= 0):
      log.error('raygrid angles must be monotonically increasing.')
  # Store raygrid values in radians:
  spec.raygrid *= sc.degree

  # Gauss quadrature integration variables:
  if spec.quadrature is not None:
      qnodes, qweights = ss.p_roots(spec.quadrature)
      spec.qnodes   = 0.5*(qnodes + 1.0)
      spec.qweights = 0.5 * qweights

  # Number of datapoints and filters:
  if obs.data is not None:
      obs.ndata = len(obs.data)
  if obs.filters is not None:
      obs.nfilters = len(obs.filters)
  # Number checks:
  if pyrat.obs.uncert is not None and pyrat.obs.ndata != len(pyrat.obs.uncert):
      log.error('Number of data uncertainty values ({:d}) does not match '
                'the number of data points ({:d}).'.
                format(len(pyrat.obs.uncert), pyrat.obs.ndata))
  if (obs.filters is not None and obs.ndata > 0 and obs.ndata != obs.nfilters):
      log.error('Number of filter bands ({:d}) does not match the '
                'number of data points ({:d}).'.
                format(obs.nfilters, obs.ndata))

  if pyrat.ncpu >= mp.cpu_count():
      log.warning('Number of requested CPUs ({:d}) is >= than the number '
                  'of available CPUs ({:d}).  Enforced ncpu to {:d}.'.
                  format(pyrat.ncpu, mp.cpu_count(), mp.cpu_count()-1))
      pyrat.ncpu = mp.cpu_count() - 1
  log.head('Check spectrum done.')


def setup(pyrat):
  """
  Process retrieval variables: bulk, molmodel.
  Process stellar spectrum.
  Process the oberving filter bands.
  """
  # Shortcuts:
  phy = pyrat.phy
  ret = pyrat.ret
  atm = pyrat.atm
  obs = pyrat.obs
  log = pyrat.log

  # Setup bulk and variable-abundance species:
  species = pyrat.mol.name
  # Non-overlapping species:
  if atm.bulk is not None and len(np.setdiff1d(atm.bulk, species)) > 0:
      log.error('These bulk species are not present in the atmosphere: {}.'.
                format(np.setdiff1d(atm.bulk, species)))

  if atm.molmodel is not None:
      if atm.molfree is None:
          log.error('molmodel is set, but there are no molfree.')
      if len(atm.molmodel) != len(atm.molfree):
          log.error('There should be one molfree for each molmodel:\n'
              'molmodel: {}\nmolfree: {}'.format(atm.molmodel, atm.molfree))
      if len(np.setdiff1d(atm.molfree, species)) > 0:
          log.error('These species are not present in the atmosphere: {}.'.
                    format(np.setdiff1d(atm.molfree, species)))

  # Overlapping species:
  if (atm.bulk is not None and atm.molfree is not None
      and len(np.intersect1d(atm.bulk, atm.molfree)) > 0):
      log.error('These species were marked as both bulk and variable-'
          'abundance: {}.'.format(np.intersect1d(atm.bulk, atm.molfree)))

  # Obtain abundance ratios between the bulk species:
  spec = list(species)
  if atm.bulk is not None:
      atm.ibulk = [spec.index(mol) for mol in atm.bulk]
      atm.bulkratio, atm.invsrat = pa.ratio(pyrat.atm.q, atm.ibulk)
  if atm.molmodel is not None:
      atm.ifree = [spec.index(mol) for mol in atm.molfree]
      nabund = len(atm.ifree)
      # Abundance free-parameter names:
      mpnames = [f'log({mol})' for mol in atm.molfree]
      mtexnames = [fr'$\log_{{10}}(X_{{\rm {mol}}})$' for mol in atm.molfree]
  else:
      nabund = 0
      mpnames, mtexnames = [], []
      atm.ibulk = None

  # Read stellar spectrum model (starspec, kurucz, or blackbody(tstar)):
  if phy.starspec is not None:
      starwn, starflux = io.read_spectrum(phy.starspec)
  elif phy.kurucz is not None:
      if phy.tstar is None:
          log.error('Undefined stellar temperature (tstar), required for '
                    'Kurucz model.')
      if phy.gstar is None:
          log.error('Undefined stellar gravity (gstar), required for '
                    'Kurucz model.')
      starflux, starwn, kuruczt, kuruczg = ps.read_kurucz(
          phy.kurucz, phy.tstar, np.log10(phy.gstar))
      log.msg('Input stellar params: T={:7.1f} K, log(g)={:4.2f}\n'
              'Best Kurucz match:    T={:7.1f} K, log(g)={:4.2f}'.
              format(phy.tstar, np.log10(phy.gstar), kuruczt, kuruczg))
  elif phy.marcs is not None:
      pass
  elif phy.phoenix is not None:
      pass
  elif phy.tstar is not None:
      starwn   = pyrat.spec.wn
      starflux = ps.bbflux(starwn, phy.tstar)
  else:
      starflux, starwn = None, None

  phy.starflux = starflux
  phy.starwn   = starwn

  # Store interpolated stellar spectrum:
  if phy.starflux is not None:
      sinterp = si.interp1d(phy.starwn, phy.starflux)
      pyrat.spec.starflux = sinterp(pyrat.spec.wn)

  is_emission = pyrat.od.rt_path in pc.emission_rt
  if pyrat.runmode=='mcmc' and is_emission and starflux is None:
      log.error('Undefined stellar flux model.  Set starspec, kurucz, or '
                'tstar (for a blackbody spectrum).')

  # Set observational variables (for given filters and other parameters):
  if obs.filters is not None:
      pyrat.set_filters()

  # Temperature models and arguments:
  if atm.tmodelname not in pc.tmodels:
      ntemp = 0
      tpnames, ttexnames = [], []
  else:
      atm.tmodel = pa.tmodels.get_model(
          atm.tmodelname, pressure=pyrat.atm.press, nlayers=pyrat.atm.nlayers)
      ntemp = atm.tmodel.npars
      tpnames = atm.tmodel.pnames
      ttexnames = atm.tmodel.texnames

  # Rayleigh models:
  nray = 0
  rpnames, rtexnames = [], []
  for rmodel in pyrat.rayleigh.models:
      rpnames   += rmodel.pnames
      rtexnames += rmodel.texnames
      nray      += rmodel.npars

  # Cloud models:
  ncloud = 0
  cpnames, ctexnames = [], []
  for model in pyrat.cloud.models:
      cpnames   += model.pnames
      ctexnames += model.texnames
      ncloud    += model.npars

  # TeX unit conversions:
  utex = {
      'mjup':r'$M_{\rm Jup}$',
      'mearth':r'$M_{\oplus}$',
      'kg':'kg',
      'gram':'g',
      'rjup':r'$R_{\rm Jup}$',
      'rearth':r'$R_{\oplus}$',
      'km':'km',
      'm':'m',
      'cm':'cm',
  }

  # Indices to parse the array of fitting parameters:
  ret.nparams = 0
  ret.pnames, ret.texnames = [], []
  if 'temp' in ret.retflag:
      ret.itemp  = np.arange(ret.nparams, ret.nparams + ntemp)
      ret.pnames   += tpnames
      ret.texnames += ttexnames
      ret.nparams += ntemp
  if 'rad' in ret.retflag:
      ret.irad   = np.arange(ret.nparams, ret.nparams + 1)
      ret.pnames   += [f'Rp ({atm.runits})']
      ret.texnames += [fr'$R_{{\rm planet}}$ ({utex[atm.runits]})']
      ret.nparams += 1
  if 'mol' in ret.retflag:
      ret.imol = np.arange(ret.nparams, ret.nparams + nabund)
      ret.pnames   += mpnames
      ret.texnames += mtexnames
      ret.nparams += nabund
  if 'ray' in ret.retflag:
      ret.iray   = np.arange(ret.nparams, ret.nparams + nray)
      ret.pnames   += rpnames
      ret.texnames += rtexnames
      ret.nparams += nray
  if 'cloud' in ret.retflag:
      ret.icloud  = np.arange(ret.nparams, ret.nparams + ncloud)
      ret.pnames   += cpnames
      ret.texnames += ctexnames
      ret.nparams += ncloud
  if 'patchy' in ret.retflag:
      ret.ipatchy = np.arange(ret.nparams, ret.nparams + 1)
      ret.pnames   += ['f_patchy']
      ret.texnames += [r'$\phi_{\rm patchy}$']
      ret.nparams += 1
  if 'mass' in ret.retflag:
      ret.imass = np.arange(ret.nparams, ret.nparams + 1)
      ret.pnames   += [f'Mp ({phy.mpunits})']
      ret.texnames += [fr'$M_{{\rm planet}}$ ({utex[phy.mpunits]})']
      ret.nparams += 1

  # Retrieval variables:
  if ret.params is not None and len(ret.params) != ret.nparams:
      nparams = len(ret.params)
      log.error(
          f'The number of input fitting parameters (params, {nparams}) does '
          f'not match\nthe number of required parameters ({ret.nparams}).')
  if pyrat.runmode == 'mcmc':
      if ret.pstep is None:
          log.error('Missing pstep argument, required for MCMC runs.')

  # Check for non-retrieval model/parameters:
  if (pyrat.rayleigh.models != []
      and (pyrat.runmode != 'mcmc' or 'ray' not in ret.retflag)
      and pyrat.rayleigh.pars is None):
      log.error('Rayleigh parameters (rpars) have not been specified.')
  if (pyrat.cloud.models != []
      and (pyrat.runmode != 'mcmc' or 'cloud' not in ret.retflag)
      and pyrat.cloud.pars is None):
      log.error('Cloud parameters (cpars) have not been specified.')
