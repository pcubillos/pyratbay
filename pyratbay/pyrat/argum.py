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


def checkinputs(pyrat):
  """
  Check that user input arguments make sense.
  """
  # Shortcuts:
  inputs = pyrat.inputs
  phy    = pyrat.phy
  log    = pyrat.log
  spec   = pyrat.spec
  atm    = pyrat.atm

  # Pyrat runmode:
  if inputs.runmode not in pc.rmodes:
      log.error("Invalid runmode ({}). Select from: {:s}.".
                format(inputs.runmode, str(pc.rmodes)))
  pyrat.runmode = inputs.runmode

  # Check that input files exist:
  if inputs.molfile is None:
      inputs.molfile = pc.ROOT + "/inputs/molecules.dat"

  pyrat.atm.atmfile = inputs.get_path('atmfile', 'Atmospheric',    exists=True)
  pyrat.mol.molfile = inputs.get_path('molfile', 'Molecular data', exists=True)
  pyrat.lt.tlifile  = inputs.get_path('tlifile', 'TLI',            exists=True)
  pyrat.cs.files    = inputs.get_path('csfile',  'Cross-section',  exists=True)
  pyrat.ex.extfile  = inputs.get_path('extfile', 'Extinction-coefficient')

  # Check spectrum arguments:
  spec.wnunits = inputs.get_default('wnunits', 'Wavenumber units', 'cm')
  spec.wlunits = inputs.get_default('wlunits', 'Wavelength units', 'um')

  spec.wllow  = inputs.get_param('wllow',  spec.wlunits,
      'Wavelength lower boundary',  gt=0.0)
  spec.wlhigh = inputs.get_param('wlhigh', spec.wlunits,
      'Wavelength higher boundary', gt=0.0)

  # Wavenumber must be taken care differently (take inverse of units):
  if inputs.wnlow is not None:
      if len(inputs.wnlow.split()) == 2:
          wnunits = inputs.wnlow.split()[1]
      else:
          wnunits = spec.wnunits
      wnlow = float(inputs.wnlow.split()[0])
      if wnlow < 0.0:
          log.error("Low wavenumber boundary ({:.2e} {:s}-1) must be >= 0.".
                    format(wnlow, wnunits))
      spec.wnlow = wnlow / pt.u(wnunits)

  if inputs.wnhigh is not None:
      if len(inputs.wnhigh.split()) == 2:
          wnunits = inputs.wnhigh.split()[1]
      else:
          wnunits = spec.wnunits
      wnhigh = float(inputs.wnhigh.split()[0])
      if wnhigh <= 0.0:
          log.error("High wavenumber boundary ({:.2e} {:s}-1) must be > 0.".
                    format(wnhigh, wnunits))
      spec.wnhigh = wnhigh / pt.u(wnunits)

  wnstep = inputs.get_default('wnstep', 'Wavenumber sampling step', '1.0 cm')
  if len(wnstep.split()) == 2:
      wnunits = wnstep.split()[1]
  else:
      wnunits = spec.wnunits
  wnstep = float(wnstep.split()[0])
  if wnstep <= 0:
      log.error("Wavenumber sampling step ({:.2e} {:s}-1) must be be > 0.".
                format(wnstep, wnunits))
  spec.wnstep = wnstep / pt.u(wnunits)

  spec.wnosamp = inputs.get_default('wnosamp',
      'Wavenumber oversampling factor', 2160, ge=1)
  spec.resolution = inputs.get_default('resolution',
      'Spectral resolution', gt=0.0)

  # Check atmospheric layers arguments:
  atm.punits = inputs.get_default('punits', 'Pressure units', 'bar')
  atm.runits = inputs.get_default('runits', 'Distance units', 'km')
  atm.nlayers = inputs.get_default('nlayers',
      'Number of atmospheric layers', gt=0)

  # Pressure boundaries:
  atm.pbottom = inputs.get_param('pbottom', atm.punits,
      'Pressure at bottom of atmosphere', gt=0.0)
  atm.ptop = inputs.get_param('ptop', atm.punits,
      'Pressure at top of atmosphere', gt=0.0)

  # Radius boundaries:
  atm.radlow = inputs.get_param('radlow', atm.runits,
      'Radius at bottom of atmosphere', ge=0.0)
  atm.radhigh = inputs.get_param('radhigh', atm.runits,
      'Radius at top of atmosphere', gt=0.0)
  atm.radstep = inputs.get_param('radstep', atm.runits,
      'Radius sampling step', gt=0.0)

  # System physical parameters:
  phy.rplanet = inputs.get_param('rplanet', atm.runits,
      'Planetary radius', gt=0.0)
  atm.refpressure = inputs.get_param('refpressure', atm.punits,
      'Planetary reference pressure level', gt=0.0)
  phy.gplanet = inputs.get_default('gplanet',
      'Planetary surface gravity (cm s-2)', gt=0.0)
  phy.mplanet = inputs.get_param('mplanet', "gram",
      'Planetary mass', gt=0.0)

  phy.rstar = inputs.get_param('rstar', atm.runits,
      'Stellar radius', gt=0.0)
  phy.gstar = inputs.get_default('gstar',
      'Stellar surface gravity', gt=0.0)
  phy.tstar = inputs.get_default('tstar',
      'Stellar effective temperature (K)', gt=0.0)
  phy.smaxis = inputs.get_param('smaxis', atm.runits,
      'Planetary radius', gt=0.0)
  phy.mstar = inputs.get_param('mstar', 'msun',
      'Stellar mass', gt=0.0)
  pyrat.phy.tint = inputs.get_default('tint',
      'Planetary internal temperature', 100.0, gt=0)

  # Check planetary surface gravity/mass/radius:
  if phy.mplanet is not None:
      atm.hydrom = True  # Use mass value for hydrostatic equilibrium
      if phy.rplanet is None and phy.gplanet is not None:
          phy.rplanet = np.sqrt(pc.G * phy.mplanet / phy.gplanet)
      if phy.rplanet is not None:
          gplanet = pc.G * phy.mplanet / phy.rplanet**2
          if phy.gplanet is None:
              phy.gplanet = gplanet
          elif np.abs(gplanet-phy.gplanet)/phy.gplanet > 0.05:
              log.error("Both mplanet and gplanet were provided, but values "
                  "are inconsistent (>5%): g(mplanet) = {:7.1f} cm s-2 and "
                  "gplanet = {:7.1f} cm s-2.".format(gplanet, phy.gplanet))
  elif phy.gplanet is not None and phy.rplanet is not None:
      phy.mplanet = phy.gplanet * phy.rplanet**2 / pc.G

  # Compute the Hill radius for the planet:
  if (phy.mstar is not None and phy.mplanet is not None
      and phy.smaxis is not None):
      phy.rhill = phy.smaxis * (phy.mplanet/(3*phy.mstar))**(1.0/3.0)

  # Check Voigt-profile arguments:
  pyrat.voigt.extent = inputs.get_default('vextent',
      'Voigt profile extent', 20.0, ge=1.0)

  # Doppler width:
  pyrat.voigt.nDop = inputs.get_default('nDop',
      'Number of Doppler-width samples', 40, ge=1)
  pyrat.voigt.Dmin = inputs.get_default('Dmin',
      'Minimum Doppler HWHM (cm-1)', gt=0.0)
  pyrat.voigt.Dmax = inputs.get_default('Dmax',
      'Maximum Doppler HWHM (cm-1)', gt=0.0)

  if (pyrat.voigt.Dmin is not None and pyrat.voigt.Dmax is not None
      and pyrat.voigt.Dmax <= pyrat.voigt.Dmin):
      log.error("Dmax ({:g} cm-1) must be > Dmin ({:g} cm-1).".
                format(pyrat.voigt.Dmax, pyrat.voigt.Dmin))

  # Lorentz width:
  pyrat.voigt.nLor = inputs.get_default('nLor',
      'Number of Lorentz-width samples', 40, ge=1)
  pyrat.voigt.Lmin = inputs.get_default('Lmin',
      'Minimum Lorentz HWHM (cm-1)', gt=0.0)
  pyrat.voigt.Lmax = inputs.get_default('Lmax',
      'Maximum Lorentz HWHM (cm-1)', gt=0.0)

  if (pyrat.voigt.Lmin is not None and pyrat.voigt.Lmax is not None
      and pyrat.voigt.Lmax <= pyrat.voigt.Lmin):
      log.error("Lmax ({:g} cm-1) must be > Lmin ({:g} cm-1).".
                format(pyrat.voigt.Lmax, pyrat.voigt.Lmin))

  pyrat.voigt.DLratio = inputs.get_default('DLratio',
      'Doppler/Lorentz-width ratio threshold', 0.1, gt=0)

  # Check extinction-coefficient arguments:
  pyrat.ex.ethresh = inputs.get_default('ethresh',
      'Extinction-cofficient threshold', 1e-15, gt=0.0)

  # Require tmin, tmax:
  if (pyrat.runmode == "opacity"
      or (pyrat.runmode in ["spectrum", "mcmc"]
          and pyrat.ex.extfile is not None
          and not os.path.isfile(pyrat.ex.extfile))):
      if inputs.tmin is None:
          log.error("Undefined lower boundary (tmin) of temperature grid for "
                    "extinction-coefficient grid.")
      if inputs.tmax is None:
          log.error("Undefined upper boundary (tmax) of temperature grid for "
                    "extinction-coefficient grid.")

  if inputs.tmin is not None:
      pyrat.ex.tmin = inputs.get_param('tmin', 'kelvin',
          'Minimum temperature of opacity grid', gt=0.0)
  if inputs.tmax is not None:
      pyrat.ex.tmax = inputs.get_param('tmax', 'kelvin',
          'Maximum temperature of opacity grid', gt=pyrat.ex.tmin)
      pyrat.ex.tstep = inputs.get_default('tstep',
          "Opacity grid's temperature sampling step in K", 100.0, gt=0.0)

  # Check haze models:
  if inputs.hazes is not None:
      nhpars = 0
      for hmodel in inputs.hazes:
          if hmodel not in hz.hnames:
              log.error("Haze model '{:s}' is not in the list of available "
                        "models:\n{:s}".format(hmodel, hz.hnames))
          else:
              ihaze = np.where(hz.hnames == hmodel)[0][0]
              pyrat.haze.model.append(hz.hmodels[ihaze])
              pyrat.haze.nmodels += 1
              nhpars += pyrat.haze.model[-1].npars
      # Process the haze parameters:
      pyrat.haze.pars = inputs.hpars
      if pyrat.haze.pars is not None:
          if nhpars != len(pyrat.haze.pars):
              log.error("Number of input haze parameters ({:d}) does not "
                        "match the number of required model parameters ({:d}).".
                        format(len(pyrat.haze.pars), nhpars))
          j = 0
          for i in np.arange(pyrat.haze.nmodels):
              npars = pyrat.haze.model[i].npars
              pyrat.haze.model[i].pars = pyrat.haze.pars[j:j+npars]
              j += npars

  pyrat.haze.fpatchy = inputs.get_default('fpatchy',
      'Patchy-cloud fraction', ge=0.0, le=1.0)

  # Check Rayleigh models:
  if inputs.rayleigh is not None:
      nrpars = 0
      for rmodel in inputs.rayleigh:
          if rmodel not in ray.rnames:
              log.error("Rayleigh model '{:s}' is not in the list of available "
                        "models:\n{:s}".format(rmodel, ray.rnames))
          j = np.where(ray.rnames == rmodel)[0][0]
          pyrat.rayleigh.model.append(ray.rmodels[j])
          pyrat.rayleigh.nmodels += 1
          nrpars += pyrat.rayleigh.model[-1].npars
      # Process the Rayleigh parameters:
      pyrat.rayleigh.pars = inputs.rpars
      if pyrat.rayleigh.pars is not None:
          if nrpars != len(pyrat.rayleigh.pars):
              log.error("Number of input Rayleigh parameters ({:d}) does not "
                        "match the number of required model parameters ({:d}).".
                        format(len(pyrat.rayleigh.pars), nrpars))
          j = 0
          for i in np.arange(pyrat.rayleigh.nmodels):
              npars = pyrat.rayleigh.model[i].npars
              pyrat.rayleigh.model[i].pars = pyrat.rayleigh.pars[j:j+npars]
              j += npars

  # Check alkali arguments:
  if inputs.alkali is not None:
      for amodel in inputs.alkali:
          if amodel not in al.mnames:
              log.error("Alkali model '{:s}' is not in the list of available "
                        "models:\n{:s}.".format(amodel, al.mnames))
          ialkali = np.where(al.mnames == amodel)[0][0]
          pyrat.alkali.model.append(al.models[ialkali])
          pyrat.alkali.nmodels += 1

  # Check optical-depth arguments:
  pyrat.od.maxdepth = inputs.get_default('maxdepth',
      'Maximum optical-depth', 10.0, ge=0.0)

  # Accept ray-path argument:
  pyrat.od.path  = inputs.path
  if pyrat.runmode in ["spectrum", "mcmc"]:
      if pyrat.od.path is None:
          log.error("Undefined observing geometry (path).  Select between "
                    "'transit' or 'eclipse'.")
      elif pyrat.od.path not in ['transit', 'eclipse']:
          log.error("Unknown observing geometry (path = {:s}).  Select between "
                    "'transit' or 'eclipse'.".format(pyrat.od.path))

  # Accept output files:
  spec.outspec = inputs.get_default('outspec',
      'Output spectrum filename', 'spectrum.dat')

  # Check system arguments:
  if pyrat.od.path == "transit" and pyrat.phy.rstar is None:
      log.error("Undefined stellar radius (rstar), required for transmission "
                "spectrum calculation.")
  # Stellar-spectrum models:
  phy.starspec = inputs.get_path('starspec', 'Stellar spectrum', exists=True)
  phy.kurucz   = inputs.get_path('kurucz',   'Kurucz model',     exists=True)
  phy.marcs    = inputs.get_path('marcs',    'MARCS model',      exists=True)
  phy.phoenix  = inputs.get_path('phoenix',  'PHOENIX model',    exists=True)

  # Check raygrid:
  raygrid = inputs.get_default('raygrid',
      'Emission raygrid (deg)', np.array([0, 20, 40, 60, 80.]))
  if raygrid[0] != 0:
      log.error("First angle in raygrid must be 0.0 (normal to surface).")
  if np.any(raygrid < 0) or np.any(raygrid > 90):
      log.error("raygrid angles must lie between 0 and 90 deg.")
  if np.any(np.ediff1d(raygrid) <= 0):
      log.error("raygrid angles must be monotonically increasing.")
  # Store raygrid values in radians:
  spec.raygrid = raygrid * sc.degree

  # Gauss quadrature integration variables:
  spec.quadrature = inputs.get_default('quadrature',
      'Number of Gaussian-quadrature points', ge=1)
  if inputs.quadrature is not None:
      qnodes, qweights = ss.p_roots(inputs.quadrature)
      spec.qnodes   = 0.5*(qnodes + 1.0)
      spec.qweights = 0.5 * qweights

  # Observational parameters:
  pyrat.obs.data   = inputs.data
  pyrat.obs.uncert = inputs.uncert
  pyrat.obs.filter = inputs.filter
  # Number of datapoints and filters:
  if inputs.data is not None:
      pyrat.obs.ndata = len(inputs.data)
  if inputs.filter is not None:
      pyrat.obs.nfilters = len(inputs.filter)
  # Number checks:
  if pyrat.obs.uncert is not None and pyrat.obs.ndata != len(pyrat.obs.uncert):
      log.error("Number of data uncertainty values ({:d}) does not match "
                "the number of data points ({:d}).".
                format(len(pyrat.obs.uncert), pyrat.obs.ndata))
  if pyrat.obs.filter is not None:
      for f in pyrat.obs.filter:
          if not os.path.isfile(f):
              log.error("Filter file: '{:s}' does not exist.".format(f))
      if pyrat.obs.ndata > 0  and  pyrat.obs.ndata != pyrat.obs.nfilters:
          log.error("Number of filter bands ({:d}) does not match the "
                    "number of data points ({:d}).".
                    format(pyrat.obs.nfilters, pyrat.obs.ndata))

  # Retrieval variables:
  # Accept species lists, check after we load the atmospheric model:
  pyrat.ret.retflag  = inputs.retflag
  pyrat.ret.qcap = inputs.get_default('qcap',
      'Metals abundance cap', 0.25, gt=0, lt=1.0)
  pyrat.ret.params = inputs.params
  if pyrat.ret.params is not None:
      pyrat.ret.nparams = len(pyrat.ret.params)
  pyrat.ret.stepsize = inputs.stepsize # FINDME checks
  pyrat.ret.tlow  = inputs.get_default('tlow',
      'Retrieval low-temperature (K) bound', 0)
  pyrat.ret.thigh = inputs.get_default('thigh',
      'Retrieval high-temperature (K) bound', np.inf)

  # Purely-MCMC variables:
  pyrat.ret.mcmcfile = inputs.mcmcfile
  pyrat.ret.walk     = inputs.walk
  pyrat.ret.pmin     = inputs.pmin
  pyrat.ret.pmax     = inputs.pmax
  pyrat.ret.nsamples = inputs.get_default('nsamples',
      'Number of MCMC samples', gt=0)
  pyrat.ret.burnin   = inputs.get_default('burnin',
      'Number of burn-in samples per chain', gt=0)
  pyrat.ret.thinning = inputs.get_default('thinning',
      'MCMC posterior thinning', 1)
  pyrat.ret.nchains  = inputs.get_default('nchains',
      'Number of MCMC parallel chains', ge=1)
  pyrat.ret.grbreak  = inputs.get_default('grbreak',
      'Gelman-Rubin convergence criteria', 0.0, ge=0)
  pyrat.ret.grnmin   = inputs.get_default('grnmin',
      'Gelman-Rubin convergence fraction', 0.5, gt=0.0)
  pyrat.ret.prior    = inputs.prior
  pyrat.ret.priorlow = inputs.priorlow
  pyrat.ret.priorup  = inputs.priorup
  # TBD: implement checks

  # Atmospheric model:
  atm.molmodel = inputs.molmodel
  atm.molfree  = inputs.molfree
  atm.molpars  = inputs.molpars
  atm.bulk     = inputs.bulk
  if inputs.tmodel is not None and inputs.tmodel not in \
          ["TCEA", "isothermal", "MadhuInv", "MadhuNoInv"]:
      log.error("Invalid temperature model '{:s}'.  Select from: "
          "TCEA, MadhuInv, MadhuNoInv or isothermal.".format(inputs.tmodel))
  atm.tmodelname = inputs.tmodel
  atm.tpars = inputs.tpars

  if atm.tmodelname == "TCEA":
      if pyrat.phy.rstar is None:
          log.error("Undefined stellar radius (rstar), required for "
                    "temperature model.")
      if pyrat.phy.tstar is None:
          log.error("Undefined stellar temperature (tstar), required for "
                    "temperature model.")
      if pyrat.phy.smaxis is None:
          log.error("Undefined orbital semi-major axis (smaxis), required for "
                    "temperature model.")
      if pyrat.phy.gplanet is None:
          log.error("Undefined planetary surface gravity (gplanet), required "
                    "for temperature model.")

  pyrat.ncpu = inputs.get_default('ncpu', 'Number of processors', 1, ge=1)
  if pyrat.ncpu >= mp.cpu_count():
      log.warning("Number of requested CPUs ({:d}) is >= than the number "
                  "of available CPUs ({:d}).  Enforced ncpu to {:d}.".
                  format(pyrat.ncpu, mp.cpu_count(), mp.cpu_count()-1))
      pyrat.ncpu = mp.cpu_count() - 1
  log.msg("Check inputs done.")


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
          "abundance: {:s}.".format(np.intersect1d(atm.bulk, atm.molfree)))

  if pyrat.runmode == "mcmc":
      if ret.retflag is None:
          log.error("Unspecified retrieval model flags.  Set the retflag list "
                    "of models selecting from: {:s}.".format(ret.rmodels))
      elif not np.all(np.in1d(ret.retflag, ret.rmodels)):
          log.error("Invalid retrieval model flags in retflag={}.  Available "
                    "options are: {}.".format(ret.retflag, ret.rmodels))
      if atm.bulk is None and "mol" in ret.retflag:
          log.error("Undefined bulk species list (bulk).")
      if atm.molmodel is None and "mol" in ret.retflag:
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
          log.error("Undefined stellar temperature (tstar), required for "
                    "Kurucz model.")
      if phy.gstar is None:
          log.error("Undefined stellar gravity (gstar), required for "
                    "Kurucz model.")
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
