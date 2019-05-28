# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import multiprocessing as mp

import numpy as np
import scipy.constants   as sc
import scipy.interpolate as si
import scipy.special     as ss

from .. import tools      as pt
from .. import constants  as pc
from .. import starspec   as ps
from .. import atmosphere as pa
from .. import io         as io

from .  import alkali as al


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
      pyrat.mol.molfile = pc.ROOT + 'inputs/molecules.dat'

  with pt.log_error(log):
      pt.file_exists('atmfile', 'Atmospheric',    pyrat.atm.atmfile)
      pt.file_exists('tlifile', 'TLI',            pyrat.lt.tlifile)
      pt.file_exists('molfile', 'Molecular-data', pyrat.mol.molfile)

  if pyrat.runmode == 'spectrum' and spec.outspec is None:
      log.error('Undefined output spectrum file (outspec).')

  # Hydrostatic by constant g or g(M,R):
  if pyrat.inputs.mplanet is not None:
      atm.hydrom = True

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
      if pyrat.od.path == 'eclipse':
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
      if pyrat.ret.walk is None:
          log.error('Undefined retrieval algorithm (walk).  Select from '
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
      pyrat.alkali.models = []
      for aname in pyrat.alkali.model_names:
          ialkali = np.where(al.mnames == aname)[0][0]
          pyrat.alkali.models.append(al.models[ialkali])

  # Accept ray-path argument:
  if pyrat.runmode in ['spectrum', 'mcmc'] and pyrat.od.path is None:
      log.error("Undefined observing geometry (path).  Select between "
                "'transit' or 'eclipse'.")

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
  if pyrat.od.path == 'transit' and phy.rstar is None:
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

  # Retrieval variables:
  if pyrat.ret.params is not None:
      pyrat.ret.nparams = len(pyrat.ret.params)

  if atm.tmodelname == 'tcea':
      if phy.rstar is None:
          log.error('Undefined stellar radius (rstar), required for '
                    'temperature model.')
      if phy.tstar is None:
          log.error('Undefined stellar temperature (tstar), required for '
                    'temperature model.')
      if phy.smaxis is None:
          log.error('Undefined orbital semi-major axis (smaxis), required for '
                    'temperature model.')
      if phy.gplanet is None:
          log.error('Undefined planetary surface gravity (gplanet), required '
                    'for temperature model.')

  if pyrat.ncpu >= mp.cpu_count():
      log.warning('Number of requested CPUs ({:d}) is >= than the number '
                  'of available CPUs ({:d}).  Enforced ncpu to {:d}.'.
                  format(pyrat.ncpu, mp.cpu_count(), mp.cpu_count()-1))
      pyrat.ncpu = mp.cpu_count() - 1
  log.msg('Check spectrum done.')


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
      mpnames   = ['log({:s})'.format(mol) for mol in atm.molfree]
      mtexnames = [r'$\log_{{10}}(f_{{\rm {:s}}})$'.format(mol)
                   for mol in atm.molfree]
  else:
      nabund = 0
      mpnames, mtexnames = [], []

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
      starflux, starwn, kuruczt, kuruczg = ps.read_kurucz(phy.kurucz,
          phy.tstar, np.log10(phy.gstar))
      log.msg('Input stellar params: T={:7.1f} K, log(g)={:4.2f}\n'
              'Best Kurucz match:    T={:7.1f} K, log(g)={:4.2f}'.
              format(phy.tstar, np.log10(phy.gstar), kuruczt, kuruczg), verb=2)
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

  if pyrat.runmode=='mcmc' and pyrat.od.path=='eclipse' and starflux is None:
      log.error('Undefined stellar flux model.  Set starspec, kurucz, or '
                'tstar (for a blackbody spectrum).')

  # Set observational variables (for given filters and other parameters):
  if obs.filters is not None:
      bandidx   = []  # Filter wavenumber indices
      starflux  = []  # Interpolated stellar flux
      bandtrans = []  # Normalized interpolated filter transmission
      bandwn    = []  # Band's mean wavenumber
      for filter in obs.filters:
          # Read filter wavenumber and transmission curves:
          filterwn, filtertr = io.read_spectrum(filter)
          # Resample the filters into the planet wavenumber array:
          btrans, bidx = pt.resample(filtertr, filterwn, pyrat.spec.wn,
              normalize=True)
          bandidx.append(bidx)
          bandtrans.append(btrans)
          bandwn.append(np.sum(filterwn*filtertr)/np.sum(filtertr))
          if phy.starflux is not None:
              starflux.append(pyrat.spec.starflux[bidx])

      # Per-band variables:
      obs.bandidx   = bandidx
      obs.bandtrans = bandtrans
      obs.starflux  = starflux
      obs.bandwn    = np.asarray(bandwn)
      obs.bandflux  = np.zeros(obs.nfilters, np.double)

  # Planet-to-star radius ratio:
  if phy.rplanet is not None and phy.rstar is not None:
      phy.rprs = phy.rplanet/phy.rstar

  # Temperature models and arguments:
  if atm.tmodelname == 'tcea':
      ntemp = 5
      atm.tmodel = pa.tmodels.tcea
      atm.targs  = [pyrat.atm.press, phy.rstar, phy.tstar, phy.tint,
                    phy.gplanet, phy.smaxis]
      tpnames   = ['log(kappa)', 'log(gamma1)', 'log(gamma2)', 'alpha', 'beta']
      ttexnames = [r'$\log_{10}(\kappa)$', r'$\log_{10}(\gamma_1)$',
                   r'$\log_{10}(\gamma2)$', r'$\alpha$', r'$\beta$']
  elif atm.tmodelname == 'isothermal':
      ntemp = 1
      atm.tmodel = pa.tmodels.isothermal
      atm.targs = [pyrat.atm.nlayers]
      tpnames   = ['T (K)']
      ttexnames = [r'$T\ ({\rm K})$']
  elif atm.tmodelname == 'madhu_noinv':
      ntemp = 5
      atm.tmodel = pa.tmodels.madhu_noinv
      atm.targs  = [pyrat.atm.press*1e-6]
      tpnames    = ['a1', 'a2', 'p1', 'p3', 'T3']
      ttexnames  = [r'$a_1$', r'$a_2$', r'$p_1$', r'$p_3$', r'$T_3$']
  elif atm.tmodelname == 'madhu_inv':
      ntemp = 6
      atm.tmodel = pa.tmodels.madhu_inv
      atm.targs  = [pyrat.atm.press*1e-6]
      tpnames    = ['a1', 'a2', 'p1', 'p2', 'p3', 'T3']
      ttexnames  = [r'$a_1$', r'$a_2$', r'$p_1$', r'$p_2$', r'$p_3$', r'$T_3$']
  else:
      ntemp = 0
      tpnames, ttexnames = [], []

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

  # Indices to parse the array of fitting parameters:
  nparams = 0
  ret.pnames, ret.texnames = [], []
  if 'temp' in ret.retflag:
      ret.itemp  = np.arange(nparams, nparams + ntemp)
      ret.pnames   += tpnames
      ret.texnames += ttexnames
      nparams += ntemp
  if 'rad' in ret.retflag:
      ret.irad   = np.arange(nparams, nparams + 1)  # nrad is always 1
      ret.pnames   += ['Radius (km)']
      ret.texnames += [r'${\rm Radius\ (km)}$']
      nparams += 1
  if 'mol' in ret.retflag:
      ret.imol = np.arange(nparams, nparams + nabund)
      ret.pnames   += mpnames
      ret.texnames += mtexnames
      nparams += nabund
  if 'ray' in ret.retflag:
      ret.iray   = np.arange(nparams, nparams + nray)
      ret.pnames   += rpnames
      ret.texnames += rtexnames
      nparams += nray
  if 'cloud' in ret.retflag:
      ret.icloud  = np.arange(nparams, nparams + ncloud)
      ret.pnames   += cpnames
      ret.texnames += ctexnames
      nparams += ncloud
  if 'patchy' in ret.retflag:
      ret.ipatchy = np.arange(nparams, nparams + 1)  # npatchy is always 1
      ret.pnames   += ['f_patchy']
      ret.texnames += [r'$f_{\rm patchy}$']
      nparams += 1

  if pyrat.runmode == 'mcmc':
      if ret.nparams != nparams:
          log.error('The input number of fitting parameters ({:d}) does not '
                    'match the number of required parameters ({:d}).'.
                    format(ret.nparams, nparams))

  # Check for non-retrieval model/parameters:
  if (pyrat.rayleigh.models != []
      and (pyrat.runmode != 'mcmc' or 'ray' not in ret.retflag)
      and pyrat.rayleigh.pars is None):
      log.error('Rayleigh parameters (rpars) have not been specified.')
  if (pyrat.cloud.models != []
      and (pyrat.runmode != 'mcmc' or 'cloud' not in ret.retflag)
      and pyrat.cloud.pars is None):
      log.error('Cloud parameters (cpars) have not been specified.')

