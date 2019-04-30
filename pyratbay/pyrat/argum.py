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

from .  import haze      as hz
from .  import rayleigh  as ray
from .  import alkali    as al

sys.path.append(pc.ROOT + "/pyratbay/atmosphere/")
import MadhuTP


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
      pyrat.mol.molfile = pc.ROOT + "inputs/molecules.dat"

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
  if (pyrat.voigt.Dmin is not None and pyrat.voigt.Dmax is not None
      and pyrat.voigt.Dmax <= pyrat.voigt.Dmin):
      log.error("Dmax ({:g} cm-1) must be > Dmin ({:g} cm-1).".
                format(pyrat.voigt.Dmax, pyrat.voigt.Dmin))

  if (pyrat.voigt.Lmin is not None and pyrat.voigt.Lmax is not None
      and pyrat.voigt.Lmax <= pyrat.voigt.Lmin):
      log.error("Lmax ({:g} cm-1) must be > Lmin ({:g} cm-1).".
                format(pyrat.voigt.Lmax, pyrat.voigt.Lmin))

  if pyrat.runmode == "opacity" or pt.isfile(pyrat.ex.extfile) == 0:
      if pyrat.ex.tmin is None:
          log.error('Undefined lower temperature boundary (tmin) for '
                    'extinction-coefficient grid.')
      if pyrat.ex.tmax is None:
          log.error('Undefined upper temperature boundary (tmax) for '
                    'extinction-coefficient grid.')

  # Check haze models:
  if pyrat.haze.model_names is not None:
      nhpars = 0
      for hmodel in pyrat.haze.model_names:
          ihaze = np.where(hz.hnames == hmodel)[0][0]
          pyrat.haze.model.append(hz.hmodels[ihaze])
          pyrat.haze.nmodels += 1
          nhpars += pyrat.haze.model[-1].npars
      # Process the haze parameters:
      if pyrat.haze.pars is not None:
          if nhpars != len(pyrat.haze.pars):
              log.error('Number of input haze parameters ({:d}) does not '
                        'match the number of required model parameters ({:d}).'.
                        format(len(pyrat.haze.pars), nhpars))
          j = 0
          for i in np.arange(pyrat.haze.nmodels):
              npars = pyrat.haze.model[i].npars
              pyrat.haze.model[i].pars = pyrat.haze.pars[j:j+npars]
              j += npars

  # Check Rayleigh models:
  if pyrat.rayleigh.model_names is not None:
      nrpars = 0
      for rmodel in pyrat.rayleigh.model_names:
          iray = np.where(ray.rnames == rmodel)[0][0]
          pyrat.rayleigh.model.append(ray.rmodels[iray])
          pyrat.rayleigh.nmodels += 1
          nrpars += pyrat.rayleigh.model[-1].npars
      # Process the Rayleigh parameters:
      if pyrat.rayleigh.pars is not None:
          if nrpars != len(pyrat.rayleigh.pars):
              log.error('Number of input Rayleigh parameters ({:d}) does not '
                        'match the number of required model parameters ({:d}).'.
                        format(len(pyrat.rayleigh.pars), nrpars))
          j = 0
          for i in np.arange(pyrat.rayleigh.nmodels):
              npars = pyrat.rayleigh.model[i].npars
              pyrat.rayleigh.model[i].pars = pyrat.rayleigh.pars[j:j+npars]
              j += npars

  # Check alkali arguments:
  if pyrat.alkali.model_names is not None:
      for amodel in pyrat.alkali.model_names:
          ialkali = np.where(al.mnames == amodel)[0][0]
          pyrat.alkali.model.append(al.models[ialkali])
          pyrat.alkali.nmodels += 1

  # Accept ray-path argument:
  if pyrat.runmode in ["spectrum", "mcmc"] and pyrat.od.path is None:
      log.error("Undefined observing geometry (path).  Select between "
                "'transit' or 'eclipse'.")

  # Check system arguments:
  if pyrat.od.path == "transit" and phy.rstar is None:
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
  if obs.filter is not None:
      obs.nfilters = len(obs.filter)
  # Number checks:
  if pyrat.obs.uncert is not None and pyrat.obs.ndata != len(pyrat.obs.uncert):
      log.error("Number of data uncertainty values ({:d}) does not match "
                "the number of data points ({:d}).".
                format(len(pyrat.obs.uncert), pyrat.obs.ndata))
  if (obs.filter is not None and obs.ndata > 0 and obs.ndata != obs.nfilters):
      log.error("Number of filter bands ({:d}) does not match the "
                "number of data points ({:d}).".
                format(obs.nfilters, obs.ndata))

  # Retrieval variables:
  if pyrat.ret.params is not None:
      pyrat.ret.nparams = len(pyrat.ret.params)

  if atm.tmodelname == 'TCEA':
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
      log.warning("Number of requested CPUs ({:d}) is >= than the number "
                  "of available CPUs ({:d}).  Enforced ncpu to {:d}.".
                  format(pyrat.ncpu, mp.cpu_count(), mp.cpu_count()-1))
      pyrat.ncpu = mp.cpu_count() - 1
  log.msg("Check spectrum done.")


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
  log = pyrat.log

  # Setup bulk and variable-abundance species:
  species = pyrat.mol.name
  # Non-overlapping species:
  if atm.bulk is not None and len(np.setdiff1d(atm.bulk, species)) > 0:
      log.error("These bulk species are not present in the atmosphere: {:s}.".
                format(str(np.setdiff1d(atm.bulk, species))))

  if atm.molmodel is not None:
      if atm.molfree is None:
          log.error("molmodel is set, but there are no molfree.")
      if len(atm.molmodel) != len(atm.molfree):
          log.error("There should be one molfree for each molmodel:\n"
              "molmodel: {}\nmolfree: {}".format(atm.molmodel, atm.molfree))
      if len(np.setdiff1d(atm.molfree, species)) > 0:
          log.error("These species are not present in the atmosphere: {:s}.".
                    format(str(np.setdiff1d(atm.molfree, species))))

  # Overlapping species:
  if (atm.bulk is not None and atm.molfree is not None
      and len(np.intersect1d(atm.bulk, atm.molfree)) > 0):
      log.error("These species were marked as both bulk and variable-"
          "abundance: {}.".format(np.intersect1d(atm.bulk, atm.molfree)))

  if pyrat.runmode == "mcmc":
      if ret.retflag is None:
          log.error("Undefined retrieval model flags.  Set the retflag list "
                    "of models selecting from: {:s}.".format(ret.rmodels))
      elif not np.all(np.in1d(ret.retflag, ret.rmodels)):
          log.error("Invalid retrieval model flags in retflag={}.  Available "
                    "options are: {}.".format(ret.retflag, ret.rmodels))
      if atm.bulk is None and "mol" in ret.retflag:
          log.error("Undefined bulk species list (bulk).")
      if atm.molmodel is None and 'mol' in ret.retflag:
          log.error("Species abundances included for retrieval (retflag "
              "contains 'mol') but there are no abundance model (molmodel).")

  # Obtain abundance ratios between the bulk species:
  spec = list(species)
  if atm.bulk is not None:
      atm.ibulk = [spec.index(mol) for mol in atm.bulk]
      atm.bulkratio, atm.invsrat = pa.ratio(pyrat.atm.q, atm.ibulk)
  if atm.molmodel is not None:
      atm.ifree = [spec.index(mol) for mol in atm.molfree]
      nabund = len(atm.ifree)
      # Abundance free-parameter names:
      mpnames   = ["log({:s})".format(mol) for mol in atm.molfree]
      mtexnames = [r"$\log_{{10}}(f_{{\rm {:s}}})$".format(mol)
                   for mol in atm.molfree]
  else:
      nabund = 0
      mpnames, mtexnames = [], []

  # Read stellar spectrum model:
  if phy.starspec is not None:
      starwn, starflux = io.read_spectrum(phy.starspec)
  # Kurucz stellar model:
  elif phy.kurucz is not None:
      if phy.tstar is None:
          log.error('Undefined stellar temperature (tstar), required for '
                    'Kurucz model.')
      if phy.gstar is None:
          log.error('Undefined stellar gravity (gstar), required for '
                    'Kurucz model.')
      starflux, starwn, kuruczt, kuruczg = ps.read_kurucz(phy.kurucz,
          phy.tstar, np.log10(phy.gstar))
      log.msg("Input stellar params: T={:7.1f} K, log(g)={:4.2f}\n"
              "Best Kurucz match:    T={:7.1f} K, log(g)={:4.2f}".
              format(phy.tstar, np.log10(phy.gstar), kuruczt, kuruczg), verb=2)
  # MARCS stellar model:
  elif phy.marcs is not None:
      pass
  # PHOENIX stellar model:
  elif phy.phoenix is not None:
      pass
  # Blackbody stellar model:
  elif phy.tstar is not None:
      starwn   = pyrat.spec.wn
      starflux = ps.bbflux(starwn, phy.tstar)
  else:
      starflux, starwn = None, None

  # Store input stellar spectrum into pyrat:
  phy.starflux  = starflux
  phy.starwn    = starwn
  # Store interpolated stellar spectrum:
  if phy.starflux is not None:
      sinterp = si.interp1d(phy.starwn, phy.starflux)
      pyrat.spec.starflux = sinterp(pyrat.spec.wn)

  # Set observational variables (for given filters and other parameters):
  setfilters(pyrat.obs, pyrat.spec, pyrat.phy)

  # Planet-to-star radius ratio:
  if phy.rplanet is not None and phy.rstar is not None:
      phy.rprs = phy.rplanet/phy.rstar

  # Temperature models and arguments:
  if atm.tmodelname == "TCEA":
      ntemp = 5
      atm.tmodel = pa.temp_TCEA
      atm.targs  = [pyrat.atm.press, phy.rstar, phy.tstar, phy.tint,
                    phy.gplanet, phy.smaxis]
      tpnames   = ["log(kappa)", "log(gamma1)", "log(gamma2)", "alpha", "beta"]
      ttexnames = [r"$\log_{10}(\kappa)$", r"$\log_{10}(\gamma_1)$",
                   r"$\log_{10}(\gamma2)$", r"$\alpha$", r"$\beta$"]
  elif atm.tmodelname == "isothermal":
      ntemp = 1
      atm.tmodel = pa.temp_isothermal
      atm.targs = [pyrat.atm.nlayers]
      tpnames   = ["T (K)"]
      ttexnames = [r"$T\ ({\rm K})$"]
  elif atm.tmodelname == "MadhuNoInv":
      ntemp = 5
      atm.tmodel = MadhuTP.no_inversion
      atm.targs  = [pyrat.atm.press*1e-6]
      tpnames    = ["a1", "a2", "p1", "p3", "T3"]
      ttexnames  = [r"$a_1$", r"$a_2$", r"$p_1$", r"$p_3$", r"$T_3$"]
  elif atm.tmodelname == "MadhuInv":
      ntemp = 6
      atm.tmodel = MadhuTP.inversion
      atm.targs  = [pyrat.atm.press*1e-6]
      tpnames    = ["a1", "a2", "p1", "p2", "p3", "T3"]
      ttexnames  = [r"$a_1$", r"$a_2$", r"$p_1$", r"$p_2$", r"$p_3$", r"$T_3$"]
  else:
      ntemp = 0
      tpnames, ttexnames = [], []

  # Rayleigh models:
  nray     = 0
  rpnames, rtexnames = [], []
  for i in np.arange(pyrat.rayleigh.nmodels):
    rpnames   += pyrat.rayleigh.model[i].pnames
    rtexnames += pyrat.rayleigh.model[i].texnames
    nray += pyrat.rayleigh.model[i].npars

  # Haze models:
  nhaze    = 0
  hpnames, htexnames = [], []
  for i in np.arange(pyrat.haze.nmodels):
      hpnames   += pyrat.haze.model[i].pnames
      htexnames += pyrat.haze.model[i].texnames
      nhaze     += pyrat.haze.model[i].npars

  # Indices to parse the array of fitting parameters:
  if ret.retflag is None:
      ret.retflag = []
  nparams = 0
  ret.pnames, ret.texnames = [], []
  if "pt" in ret.retflag:
      ret.itemp  = np.arange(nparams, nparams + ntemp)
      ret.pnames   += tpnames
      ret.texnames += ttexnames
      nparams += ntemp
  if "rad" in ret.retflag:
      ret.irad   = np.arange(nparams, nparams + 1)  # nrad is always 1
      ret.pnames   += ["Radius (km)"]
      ret.texnames += [r"${\rm Radius\ (km)}$"]
      nparams += 1
  if "mol" in ret.retflag:
      ret.iabund = np.arange(nparams, nparams + nabund)
      ret.pnames   += mpnames
      ret.texnames += mtexnames
      nparams += nabund
  if "ray" in ret.retflag:
      ret.iray   = np.arange(nparams, nparams + nray)
      ret.pnames   += rpnames
      ret.texnames += rtexnames
      nparams += nray
  if "haze" in ret.retflag:
      ret.ihaze  = np.arange(nparams, nparams + nhaze)
      ret.pnames   += hpnames
      ret.texnames += htexnames
      nparams += nhaze
  #if "cloud" in ret.retflag:
  #  ret.icloud  = np.arange(nparams, nparams + ncloud)
  #  ret.pnames   += cpnames
  #  ret.texnames += ctexnames
  #  nparams += ncloud
  if "patchy" in ret.retflag:
      ret.ipatchy = np.arange(nparams, nparams + 1)  # npatchy is always 1
      ret.pnames   += ["f_patchy"]
      ret.texnames += [r"$f_{\rm patchy}$"]
      nparams += 1

  if pyrat.runmode == "mcmc":
      if ret.nparams != nparams:
          log.error("The input number of fitting parameters ({:d}) does not "
                    "match the number of model parameters ({:d}).".
                    format(ret.nparams, nparams))

  # Check for non-retrieval model/parameters:
  if (pyrat.rayleigh.nmodels > 0
      and (pyrat.runmode != "mcmc" or "ray" not in ret.retflag)
      and pyrat.rayleigh.pars is None):
      log.error("Rayleigh parameters (rpars) have not been specified.")
  if (pyrat.haze.nmodels > 0
      and (pyrat.runmode != "mcmc" or "haze" not in ret.retflag)
      and pyrat.haze.pars is None):
      log.error("Haze parameters (hpars) have not been specified.")


def setfilters(obs, spec, phy):
  """
  Set observational variables (pyrat.obs) based on given parameters.
  """
  # Skip if there are no filter bands:
  if obs.filter is None:
      return
  # Load filters:
  bandidx   = []  # Filter wavenumber indices
  starflux  = []  # Interpolated stellar flux
  bandtrans = []  # Normalized interpolated filter transmission
  bandwn    = []  # Band's mean wavenumber
  for i in np.arange(obs.nfilters):
      # Read filter wavenumber and transmission curves:
      filterwn, filtertr = io.read_spectrum(obs.filter[i])
      # Resample the filters into the stellar wavenumber array:
      btrans, bidx = pt.resample(filtertr, filterwn, spec.wn, normalize=True)
      sflux, dummy = pt.resample(phy.starflux, phy.starwn, spec.wn)
      bandidx.append(bidx)
      bandtrans.append(btrans)
      starflux.append(sflux)
      bandwn.append(np.sum(filterwn*filtertr)/np.sum(filtertr))

  # Per-band variables:
  obs.bandidx   = bandidx
  obs.bandtrans = bandtrans
  obs.starflux  = starflux
  obs.bandwn    = np.asarray(bandwn)
  obs.bandflux  = np.zeros(obs.nfilters, np.double)
