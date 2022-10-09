# Copyright (c) 2021-2022 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import multiprocessing as mp

import numpy as np
import scipy.constants as sc
import scipy.interpolate as si
import scipy.special as ss

from .. import tools as pt
from .. import constants as pc
from .. import spectrum as ps
from .. import atmosphere as pa
from .. import io as io


def check_spectrum(pyrat):
    """
    Check that user input arguments make sense.
    """
    # Shortcuts:
    log = pyrat.log
    phy = pyrat.phy
    spec = pyrat.spec
    atm = pyrat.atm
    obs = pyrat.obs

    # Check that input files exist:
    if pyrat.mol.molfile is None:
        pyrat.mol.molfile = pc.ROOT + 'pyratbay/data/molecules.dat'

    with pt.log_error(log):
        pt.file_exists('atmfile', 'Atmospheric', atm.atmfile)
        pt.file_exists('tlifile', 'TLI', pyrat.lt.tlifile)
        pt.file_exists('molfile', 'Molecular-data', pyrat.mol.molfile)

    if pyrat.runmode == 'spectrum' and spec.specfile is None:
        log.error('Undefined output spectrum file (specfile).')

    # Compute the Hill radius for the planet:
    if (phy.mstar is not None and phy.mplanet is not None
        and phy.smaxis is not None):
        phy.rhill = phy.smaxis * (phy.mplanet/(3*phy.mstar))**(1.0/3.0)

    # Check Voigt-profile arguments:
    dmin = pyrat.voigt.dmin
    dmax = pyrat.voigt.dmax
    if dmin is not None and dmax is not None and dmax <= dmin:
        log.error(f'dmax ({dmax:g} cm-1) must be > dmin ({dmin:g} cm-1)')

    lmin = pyrat.voigt.lmin
    lmax = pyrat.voigt.lmax
    if lmin is not None and lmax is not None and lmax <= lmin:
        log.error(f'lmax ({lmax:g} cm-1) must be > lmin ({lmin:g} cm-1)')

    if pyrat.runmode == 'opacity' or pt.isfile(pyrat.ex.extfile) == 0:
        if pyrat.ex.tmin is None:
            log.error(
                'Undefined lower temperature boundary (tmin) for '
                'extinction-coefficient grid',
            )
        if pyrat.ex.tmax is None:
            log.error(
                'Undefined upper temperature boundary (tmax) for '
                'extinction-coefficient grid',
            )
        if pyrat.ex.tstep is None:
            log.error(
                'Undefined temperature sampling step (tstep) for '
                'extinction-coefficient grid',
            )
        if pyrat.lt.tlifile is None:
            log.error(
                'Requested extinction-coefficient table, but there '
                'are no input TLI files',
            )

    if pyrat.runmode == 'mcmc':
        if pyrat.od.rt_path in pc.emission_rt:
            if phy.rplanet is None or phy.rstar is None:
                log.error("Undefined radius ratio (need rplanet and rstar)")
        if obs.data is None:
            log.error("Undefined transit/eclipse data (data)")
        if obs.uncert is None:
            log.error("Undefined data uncertainties (uncert)")
        if obs.filters is None:
            log.error("Undefined transmission filters (filters)")
        if pyrat.ret.retflag == []:
            flags = pc.retflags
            log.error(f'Undefined retrieval model flags.  Select from {flags}')
        if pyrat.ret.sampler is None:
            log.error(
                'Undefined retrieval algorithm (sampler).  '
                'Select from [snooker]'
            )
        if pyrat.ret.nsamples is None:
            log.error('Undefined number of retrieval samples (nsamples)')
        if pyrat.ret.burnin is None:
            log.error('Undefined number of retrieval burn-in samples (burnin)')
        if pyrat.ret.nchains is None:
            log.error('Undefined number of retrieval parallel chains (nchains)')
        if pyrat.ret.params is None:
            log.error('Undefined retrieval fitting parameters (params)')

    # Check cloud models:
    if pyrat.cloud.model_names is not None:
        pyrat.cloud.models = []
        npars = 0
        for name in pyrat.cloud.model_names:
            model  = pa.clouds.get_model(name)
            npars += model.npars
            pyrat.cloud.models.append(model)
        # Parse the cloud parameters:
        input_npars = len(pyrat.cloud.pars)
        if pyrat.cloud.pars is not None:
            if npars != input_npars:
                log.error(
                    f'Number of input cloud parameters ({input_npars}) '
                    'does not match the number of required '
                    f'model parameters ({npars})'
                )
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
            input_npars = len(pyrat.rayleigh.pars)
            if npars != input_npars:
                log.error(
                    f'Number of input Rayleigh parameters ({input_npars}) '
                    'does not match the number of required '
                    f'model parameters ({npars})'
                )
            j = 0
            for model in pyrat.rayleigh.models:
                npars = model.npars
                model.pars = pyrat.rayleigh.pars[j:j+npars]
                j += npars

    # Check alkali arguments:
    if pyrat.alkali.model_names is not None:
        pyrat.alkali.models = [
            pa.alkali.get_model(name, pyrat.alkali.cutoff)
            for name in pyrat.alkali.model_names
        ]

    # Accept ray-path argument:
    if pyrat.runmode in ['spectrum', 'mcmc'] and pyrat.od.rt_path is None:
        log.error(
            "Undefined radiative-transfer observing geometry (rt_path)."
            f"  Select from {pc.rt_paths}"
        )

    if 'temp' in pyrat.ret.retflag and atm.tmodelname is None:
        log.error('Requested temp in retflag, but there is no tmodel')
    if 'mol' in pyrat.ret.retflag:
        if atm.molmodel == []:
            log.error("Requested mol in retflag, but there is no 'molmodel'")
        # TBD: This will break for pure eq-chem runs
        if atm.bulk is None:
            log.error('Requested mol in retflag, but there are no bulk species')
    if 'ray' in pyrat.ret.retflag and pyrat.rayleigh.models == []:
        log.error('Requested ray in retflag, but there are no rayleigh models')
    if 'cloud' in pyrat.ret.retflag and pyrat.cloud.models == []:
        log.error('Requested cloud in retflag, but there are no cloud models')

    # Check system arguments:
    if pyrat.od.rt_path in pc.transmission_rt and phy.rstar is None:
        log.error(
            'Undefined stellar radius (rstar), required for transmission '
            'calculation'
        )

    # Check raygrid:
    if spec.raygrid[0] != 0:
        log.error('First angle in raygrid must be 0.0 (normal to surface)')
    if np.any(spec.raygrid < 0) or np.any(spec.raygrid > 90):
        log.error('raygrid angles must lie between 0 and 90 deg')
    if np.any(np.ediff1d(spec.raygrid) <= 0):
        log.error('raygrid angles must be monotonically increasing')
    # Store raygrid values in radians:
    spec.raygrid *= sc.degree

    # Gauss quadrature integration variables:
    if spec.quadrature is not None:
        qnodes, qweights = ss.p_roots(spec.quadrature)
        spec.qnodes = 0.5*(qnodes + 1.0)
        spec.qweights = 0.5 * qweights

    # Number of datapoints and filters:
    if obs.data is not None:
        obs.ndata = len(obs.data)
    if obs.filters is not None:
        obs.nfilters = len(obs.filters)
    # Number checks:
    if obs.uncert is not None:
        n_uncert = len(obs.uncert)
        if obs.ndata != n_uncert:
            log.error(
                f'Number of data uncertainty values ({n_uncert}) does not '
                f'match the number of data points ({obs.ndata})'
            )

    if obs.filters is not None and obs.ndata > 0 and obs.ndata != obs.nfilters:
        log.error(
            f'Number of filter bands ({obs.nfilters}) does not match the '
            f'number of data points ({obs.ndata})'
        )

    if pyrat.ncpu >= mp.cpu_count():
        n_cpu = mp.cpu_count()
        log.warning(
            f'Number of requested CPUs ({pyrat.ncpu}) is >= than the number '
            f'of available CPUs ({n_cpu}).  Enforced ncpu to {n_cpu-1}'
        )
        pyrat.ncpu = n_cpu - 1
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
      missing = np.setdiff1d(atm.bulk, species)
      log.error(
          f'These bulk species are not present in the atmosphere: {missing}.')

  if atm.molmodel != []:
      if atm.molfree is None:
          log.error('molmodel is set, but there are no molfree.')
      if len(atm.molmodel) != len(atm.molfree):
          log.error(
              'There should be one molfree for each molmodel:\n'
              f'molmodel: {atm.molmodel}\nmolfree: {atm.molfree}')
      free_molecules = [
          mol for mol,model in zip(atm.molfree,atm.molmodel)
          if model != 'equil']
      missing_molecules = np.setdiff1d(free_molecules, species)
      if len(missing_molecules) > 0:
          log.error(
              'These molfree species are not present in the '
              f'atmosphere: {missing_molecules}.')

  # Overlapping species:
  if atm.bulk is not None and atm.molfree is not None:
      bulk_free_species = np.intersect1d(atm.bulk, atm.molfree)
      if len(bulk_free_species) > 0:
          log.error(
              'These species were marked as both bulk and '
              f'variable-abundance: {bulk_free_species}.'
          )

  # Obtain abundance ratios between the bulk species:
  spec = list(species)
  if atm.bulk is not None:
      atm.ibulk = [spec.index(mol) for mol in atm.bulk]
      atm.bulkratio, atm.invsrat = pa.ratio(pyrat.atm.q, atm.ibulk)
  if atm.molmodel != []:
      nabund = len(atm.molmodel)
      # TBD: set ifree to something when model is equil
      atm.ifree = [
          spec.index(mol)
          for mol,model in zip(atm.molfree,atm.molmodel)
          if model != 'equil'
      ]
      # Abundance free-parameter names:
      mol_pnames = []
      mol_tex_names = []
      for mol in atm.molfree:
          if mol == 'metal':
              mol_pnames.append('[Z/H]')
              mol_tex_names.append('[Z/H]')
          elif '_' in mol:
              mol_pnames.append(mol.replace('_','/'))
              mol_tex_names.append(mol.replace('_','/'))
          else:
              mol_pnames.append(f'log({mol})')
              mol_tex_names.append(fr'$\log_{{10}}(X_{{\rm {mol}}})$')

  else:
      nabund = 0
      mol_pnames, mol_tex_names = [], []
      atm.ibulk = None

  # Read stellar spectrum model (starspec, kurucz, or blackbody(tstar)):
  if phy.starspec is not None:
      starwn, starflux = io.read_spectrum(phy.starspec)
  elif phy.kurucz is not None:
      if phy.tstar is None:
          log.error(
              'Undefined stellar temperature (tstar), required for Kurucz model'
          )
      if phy.gstar is None:
          log.error('Undefined stellar gravity (gstar), required for '
                    'Kurucz model.')
      starflux, starwn, kuruczt, kuruczg = ps.read_kurucz(
          phy.kurucz, phy.tstar, np.log10(phy.gstar),
      )
      logg = np.log10(phy.gstar)
      log.msg(
          f'Input stellar params: T={phy.tstar:7.1f} K, log(g)={logg:4.2f}\n'
          f'Best Kurucz match:    T={kuruczt:7.1f} K, log(g)={kuruczg:4.2f}'
      )
  elif phy.marcs is not None:
      pass
  elif phy.phoenix is not None:
      pass
  elif phy.tstar is not None:
      starwn = pyrat.spec.wn
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
          atm.tmodelname, pressure=pyrat.atm.press, nlayers=pyrat.atm.nlayers,
      )
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
      ret.irad = np.arange(ret.nparams, ret.nparams + 1)
      ret.pnames   += [f'Rp ({atm.runits})']
      ret.texnames += [fr'$R_{{\rm planet}}$ ({utex[atm.runits]})']
      ret.nparams += 1
  if 'press' in ret.retflag:
      ret.ipress = np.arange(ret.nparams, ret.nparams + 1)
      ret.pnames   += ['log(P_ref) (bar)']
      ret.texnames += [r'$\log p_{{\rm ref}}$']
      ret.nparams += 1
  if 'mol' in ret.retflag:
      ret.imol = np.arange(ret.nparams, ret.nparams + nabund)
      ret.pnames += mol_pnames
      ret.texnames += mol_tex_names
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
          f'not match\nthe number of required parameters ({ret.nparams}).'
      )
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
