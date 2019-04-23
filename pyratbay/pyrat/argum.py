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
  if inputs.atmfile is None:
      log.error("Undefined atmospheric file (atmfile).")
  elif not os.path.isfile(inputs.atmfile):
      log.error("Atmospheric file: '{:s}' does not exist.".
                format(inputs.atmfile))
  atm.atmfile = inputs.atmfile

  if inputs.tlifile is not None:
      for tli in inputs.tlifile:
          if not os.path.isfile(tli):
              log.error("TLI file: '{:s}' does not exist.".format(tli))
  pyrat.lt.tlifile = pyrat.inputs.tlifile

  if inputs.csfile is not None:
      for cs in pyrat.inputs.csfile:
          if not os.path.isfile(cs):
              log.error("Cross-section file: '{:s}' does not exist.".format(cs))
  pyrat.cs.files = pyrat.inputs.csfile

  if inputs.molfile is None: # Set default
      inputs.molfile = os.path.realpath(pc.ROOT + "/inputs/molecules.dat")
  if not os.path.isfile(inputs.molfile):
      log.error("Molecular-data file: '{:s}' does not exist.".
                format(inputs.molfile))
  pyrat.mol.molfile = os.path.realpath(inputs.molfile)

  if inputs.extfile is not None:
    if not os.path.exists(os.path.realpath(os.path.dirname(inputs.extfile))):
      log.error("Directory for extinction-coefficient file '{:s}' does "
                "not exist.".format(inputs.extfile))
    pyrat.ex.extfile = os.path.realpath(inputs.extfile)

  # Check spectrum arguments:
  spec.wnunits = pt.defaultp(inputs.wnunits, "cm",
      "wnunits input variable defaulted to '{:s}'.", log)
  spec.wlunits = pt.defaultp(inputs.wlunits, "um",
      "wlunits input variable defaulted to '{:s}'.", log)

  spec.wllow = pt.getparam('wllow', inputs.wllow, spec.wlunits, log)
  isgreater(spec.wllow, "um", 0, False,
      "Low wavelength boundary ({:.2e} um) must be >= 0.", log)

  spec.wlhigh = pt.getparam('wlhigh', inputs.wlhigh, spec.wlunits, log)
  isgreater(spec.wlhigh, "um", 0, True,
      "High wavelength boundary ({:.2e} um) must be >= 0.", log)

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

  wnstep = pt.defaultp(inputs.wnstep, "1.0 cm",
      "Input wavenumber sampling step (wnstep) defaulted to '{:s}'.", log)
  if len(wnstep.split()) == 2:
      wnunits = wnstep.split()[1]
  else:
      wnunits = spec.wnunits
  wnstep = float(wnstep.split()[0])
  if wnstep <= 0:
      log.error("Wavenumber sampling step ({:.2e} {:s}-1) must be be > 0.".
                format(wnstep, wnunits))
  spec.wnstep = wnstep / pt.u(wnunits)

  spec.wnosamp = pt.defaultp(inputs.wnosamp, 2160,
      "Input wavenumber oversampling factor (wnosamp) defaulted to {:d}.", log)
  isgreater(spec.wnosamp, "none", 1, False,
      "Wavenumber oversampling factor ({:d}) must be >= 1.", log)

  spec.resolution = inputs.resolution

  # Check atmospheric layers arguments:
  atm.punits = pt.defaultp(inputs.punits, "bar",
     "Input pressure units (punits) defaulted to '{:s}'.", log)
  atm.runits = pt.defaultp(inputs.runits, "km",
     "Input radius units (punits) defaulted to '{:s}'.", log)

  # Pressure boundaries:
  atm.pbottom = pt.getparam('pbottom', inputs.pbottom, atm.punits, log)
  isgreater(atm.pbottom, "bar", 0, True,
            "High atm pressure boundary ({:.2e} bar) must be > 0.0", log)
  atm.ptop = pt.getparam('ptop', inputs.ptop,  atm.punits, log)
  isgreater(atm.ptop, "bar",  0, True,
            "Low atm pressure boundary ({:.2e} bar) must be > 0.0", log)
  # Radius boundaries:
  atm.radlow = pt.getparam('radlow', inputs.radlow,  atm.runits, log)
  isgreater(atm.radlow, "cm", 0, False,
            "Low atm radius boundary ({:.2e} cm) must be >= 0.0", log)
  atm.radhigh = pt.getparam('radhigh', inputs.radhigh, atm.runits, log)
  isgreater(atm.radhigh, "cm", 0, True,
            "High atm radius boundary ({:.2e} cm) must be > 0.0", log)
  atm.radstep = pt.getparam('radstep', inputs.radstep, atm.runits, log)
  isgreater(atm.radstep, "cm", 0, True,
            "Radius step size ({:.2f} cm) must be > 0.", log)
  # System physical parameters:
  phy.rplanet = pt.getparam('rplanet', inputs.rplanet, atm.runits, log)
  isgreater(phy.rplanet, "cm",   0, True,
            "Planetary radius ({:.3e} cm) must be > 0.", log)

  atm.refpressure = pt.getparam('refpressure', inputs.refpressure, atm.punits,
                                log)
  isgreater(atm.refpressure, "bar", 0, True,
      "Planetary reference pressure level ({:8g} bar) must be > 0.", log)

  phy.gplanet = pt.getparam('gplanet', inputs.gplanet, "none", log)
  isgreater(phy.gplanet, "none", 0, True,
            "Planetary surface gravity ({:.2f} cm s-2) must be > 0.", log)

  phy.mplanet = pt.getparam('mplanet', inputs.mplanet, "gram", log)
  isgreater(phy.mplanet, "mearth", 0, True,
            "Planetary mass ({:.2f} Mearth) must be > 0.", log)

  # Check planetary surface gravity:
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

  phy.rstar = pt.getparam('rstar', inputs.rstar, atm.runits, log)
  isgreater(pyrat.phy.rstar, "cm",   0, True,
      "Stellar radius ({:.3e} cm) must be > 0.", log)

  pyrat.phy.gstar = pt.getparam('gstar', inputs.gstar, "none", log)
  isgreater(pyrat.phy.gstar, "none", 0, True,
      "Stellar surface gravity ({:.2f} cm s-2) must be > 0.", log)

  pyrat.phy.tstar = pt.getparam('tstar', inputs.tstar, "none", log)
  isgreater(pyrat.phy.tstar, "none", 0, True,
      "Stellar effective temperature ({:.1f} K) must be > 0.", log)

  pyrat.phy.smaxis = pt.getparam('smaxis', inputs.smaxis, atm.runits, log)
  isgreater(pyrat.phy.smaxis, "cm",   0, True,
      "Planetary radius ({:.3e} cm) must be > 0.", log)

  phy.mstar = pt.getparam('mstar', inputs.mstar, "gram", log)
  isgreater(phy.mstar, "msun", 0, True,
      "Stellar mass ({:.2f} Msun) must be > 0.", log)

  pyrat.phy.tint = pt.defaultp(inputs.tint, 100.0,
      "Planetary internal temperature (tint) defaulted to {:.1f} K.", log)
  isgreater(phy.tint, "none", 0, True,
      "Planetary internal temperature ({:.1f} K) must be > 0.", log)

  # Compute the Hill radius for the planet:
  if (phy.mstar is not None and phy.mplanet is not None and
      phy.smaxis is not None):
      phy.rhill = phy.smaxis * (phy.mplanet/(3*phy.mstar))**(1.0/3.0)

  atm.nlayers = pt.getparam('nlayers', inputs.nlayers, 'none', log,integer=True)
  isgreater(atm.nlayers, "none", 0, True,
      "The number of atmospheric layers ({:d}) must be > 0.", log)

  # Check Voigt-profile arguments:
  pyrat.voigt.extent = pt.defaultp(inputs.vextent, 20.0,
      "Input Voigt extent (vextent) defaulted to {:g}.", log)
  isgreater(pyrat.voigt.extent, "none", 1, False,
      "Voigt extent ({:g}) must be >= 1.0", log)

  # Doppler width:
  pyrat.voigt.nDop = pt.defaultp(inputs.nDop, 40,
      "Number of Doppler-width samples (nDop) defaulted to {:d}.", log)
  isgreater(pyrat.voigt.nDop, "none", 1, False,
      "The number of Doppler-width samples ({:d}) must be >= 1", log)

  pyrat.voigt.Dmin = pt.getparam('Dmin', inputs.Dmin, "none", log)
  isgreater(pyrat.voigt.Dmin, "none", 0, True,
      "Dmin ({:g} cm-1) must be > 0.", log)

  pyrat.voigt.Dmax = pt.getparam('Dmax', inputs.Dmax, "none", log)
  isgreater(pyrat.voigt.Dmax, "none", 0, True,
      "Dmax ({:g} cm-1) must be > 0.", log)

  if (pyrat.voigt.Dmin is not None and pyrat.voigt.Dmax is not None and
      pyrat.voigt.Dmax <= pyrat.voigt.Dmin):
      log.error("Dmax ({:g} cm-1) must be > Dmin ({:g} cm-1).".
                format(pyrat.voigt.Dmax, pyrat.voigt.Dmin))

  # Lorentz width:
  pyrat.voigt.nLor = pt.defaultp(inputs.nLor, 40,
      "Number of Lorentz-width samples (nLor) defaulted to {:d}.", log)
  isgreater(pyrat.voigt.nLor, "none", 1, False,
      "The number of Lorentz-width samples ({:d}) must be >= 1", log)

  pyrat.voigt.Lmin = pt.getparam('Lmin', inputs.Lmin, "none", log)
  isgreater(pyrat.voigt.Lmin, "none", 0, True,
      "Lmin ({:g} cm-1) must be > 0.", log)

  pyrat.voigt.Lmax = pt.getparam('Lmax', inputs.Lmax, "none", log)
  isgreater(pyrat.voigt.Lmax, "none", 0, True,
      "Lmax ({:g} cm-1) must be > 0.", log)

  if (pyrat.voigt.Lmin is not None and pyrat.voigt.Lmax is not None and
      pyrat.voigt.Lmax <= pyrat.voigt.Lmin):
      log.error("Lmax ({:g} cm-1) must be > Lmin ({:g} cm-1).".
                format(pyrat.voigt.Lmax, pyrat.voigt.Lmin))

  pyrat.voigt.DLratio = pt.defaultp(inputs.DLratio, 0.1,
      "Doppler/Lorentz-width ratio threshold (DLratio) defaulted to {:g}.", log)
  isgreater(pyrat.voigt.DLratio, "none", 0, True,
      "Doppler/Lorentz-width ratio threshold ({:g}) must be > 0.", log)

  # Check extinction-coefficient arguments:
  pyrat.ex.ethresh = pt.defaultp(inputs.ethresh, 1e-15,
      "ethresh defaulted to {}.", log)

  isgreater(pyrat.ex.ethresh, "none", 0, True,
      "Extinction-coefficient threshold ({:g}) must be positive.", log)
  # Require tmin, tmax:
  if (pyrat.runmode == "opacity" or
      (pyrat.runmode in ["spectrum", "mcmc"] and
       pyrat.ex.extfile is not None and
       not os.path.isfile(pyrat.ex.extfile))):
    if inputs.tmin is None:
      log.error("Undefined lower boundary (tmin) of temperature grid for "
               "extinction-coefficient grid.")
    if inputs.tmax is None:
      log.error("Undefined upper boundary (tmax) of temperature grid for "
               "extinction-coefficient grid.")

  if inputs.tmin is not None:
    pyrat.ex.tmin = pt.getparam('tmin', inputs.tmin, "kelvin", log)
    isgreater(pyrat.ex.tmin,  "kelvin", 0, True,
          "Minimum temperature sample ({:g} K) must be positive.", log)
  if inputs.tmax is not None:
    pyrat.ex.tmax  = pt.getparam('tmax', inputs.tmax, "kelvin", log)
    isgreater(pyrat.ex.tmax,  "kelvin", 0, True,
          "Maximum temperature sample ({:g} K) must be positive.", log)

    pyrat.ex.tstep = pt.defaultp(inputs.tstep, 100,
      "Extinction-coefficient grid's temperature sampling interval (tstep) "
      "defaulted to {:g} K.", log)

    isgreater(pyrat.ex.tstep, "kelvin", 0, True,
      "Temperature sample step interval ({:g} K) must be positive.", log)

    if pyrat.ex.tmax <= pyrat.ex.tmin:
      log.error("Extinction-coefficient grid's maximum temperature ({:g} K) "
               "must be > minimum temperature ({:g} K).".
                format(pyrat.ex.tmax, pyrat.ex.tmin))

  # Check haze models:
  if inputs.hazes is not None:
    nhpars = 0
    for hmodel in inputs.hazes:
      if hmodel not in hz.hnames:
        log.error("Haze model '{:s}' is not in the list of available models:"
                 "\n{:s}".format(hmodel, hz.hnames))
      else:
        ihaze = np.where(hz.hnames == hmodel)[0][0]
        pyrat.haze.model.append(hz.hmodels[ihaze])
        pyrat.haze.nmodels += 1
        nhpars += pyrat.haze.model[-1].npars
    # Process the haze parameters:
    pyrat.haze.pars = inputs.hpars
    if pyrat.haze.pars is not None:
      if nhpars != len(pyrat.haze.pars):
        log.error("The number of input haze parameters ({:d}) does not match "
                 "the number of required haze parameters ({:d}).".
                 format(len(pyrat.haze.pars), nhpars))
      j = 0
      for i in np.arange(pyrat.haze.nmodels):
        npars = pyrat.haze.model[i].npars
        pyrat.haze.model[i].pars = pyrat.haze.pars[j:j+npars]
        j += npars

  if inputs.fpatchy is not None:
    if inputs.fpatchy < 0 or inputs.fpatchy > 1:
      log.error("Invalid patchy-cloud fraction ({:g}).  fpatchy must be "
               "in the range 0--1.".format(inputs.fpatchy))
    pyrat.haze.fpatchy = inputs.fpatchy

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
        log.error("The number of input Rayleigh parameters ({:d}) does not "
                 "match the number of required parameters ({:d}).".
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
        log.error("Alkali model '{:s}' is not in the list of available models:"
                 "\n{:s}.".format(amodel, al.mnames))
      ialkali = np.where(al.mnames == amodel)[0][0]
      pyrat.alkali.model.append(al.models[ialkali])
      pyrat.alkali.nmodels += 1

  # Check optical-depth arguments:
  pyrat.od.maxdepth = pt.defaultp(inputs.maxdepth, 10.0,
   "Maximum optical-depth (maxdepth) defaulted to {:g}.", log)
  isgreater(pyrat.od.maxdepth, "none", 0, False,
            "Maximum optical-depth limit ({:g}) must be >= 0.0", log)

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
  spec.outspec = pt.defaultp(inputs.outspec, 'outspec.dat',
      "Output spectrum filename (outspec) defaulted to '{:s}'.", log)

  # Check system arguments:
  if pyrat.od.path == "transit" and pyrat.phy.rstar is None:
    log.error("Undefined stellar radius (rstar), required for transmission "
             "spectrum calculation.")
  # Stellar-spectrum models:
  pyrat.phy.starspec = inputs.starspec
  pyrat.phy.kurucz   = inputs.kurucz
  pyrat.phy.marcs    = inputs.marcs
  pyrat.phy.phoenix  = inputs.phoenix
  if inputs.starspec is not None and not os.path.isfile(inputs.starspec):
    log.error("Stellar-spectrum model file: '{:s}' does not exist.".
             format(inputs.starspec))
  if inputs.kurucz  is not None and not os.path.isfile(inputs.kurucz):
    log.error("Stellar Kurucz model file: '{:s}' does not exist.".
             format(inputs.kurucz))
  if inputs.marcs   is not None and not os.path.isfile(inputs.marcs):
    log.error("Stellar MARCS model file: '{:s}' does not exist.".
             format(inputs.marcs))
  if inputs.phoenix is not None and not os.path.isfile(inputs.phoenix):
    log.error("Stellar PHOENIX model file: '{:s}' does not exist.".
             format(inputs.phoenix))

  # Check raygrid:
  if inputs.raygrid is None:
    raygrid = pt.defaultp(inputs.raygrid, np.array([0, 20, 40, 60, 80.]),
        "Defaulted emission raygrid to {}.", log)
  else:
    raygrid = inputs.raygrid
    if raygrid[0] != 0:
      log.error("First angle in raygrid must be 0.0 (normal to surface).")
    if np.any(raygrid < 0) or np.any(raygrid > 90):
      log.error("raygrid angles must lie between 0 and 90 deg.")
    if np.any(np.ediff1d(raygrid) <= 0):
      log.error("raygrid angles must be monotonically increasing.")
  # Store raygrid values in radians:
  spec.raygrid = raygrid * sc.degree

  # Gauss quadrature integration variables:
  spec.quadrature = inputs.quadrature
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
    log.error("The number of data uncertainty values ({:d}) does not match "
       "the number of data points ({:d}).".
        format(len(pyrat.obs.uncert), pyrat.obs.ndata))
  if pyrat.obs.filter is not None:
    for f in pyrat.obs.filter:
      if not os.path.isfile(f):
        log.error("Filter file: '{:s}' does not exist.".format(f))
    if pyrat.obs.ndata > 0  and  pyrat.obs.ndata != pyrat.obs.nfilters:
      log.error("The number of filter bands ({:d}) does not match the number "
          "of data points ({:d}).".format(pyrat.obs.nfilters, pyrat.obs.ndata))

  # Retrieval variables:
  # Accept species lists, check after we load the atmospheric model:
  pyrat.ret.retflag  = inputs.retflag
  pyrat.ret.qcap     = inputs.qcap
  pyrat.ret.qcap = pt.defaultp(inputs.qcap, 0.25, "qcap defaulted to '{}'.", log)
  pyrat.ret.params   = inputs.params
  if pyrat.ret.params is not None:
      pyrat.ret.nparams = len(pyrat.ret.params)
  pyrat.ret.stepsize = inputs.stepsize # FINDME checks
  pyrat.ret.tlow  = pt.defaultp(inputs.tlow, 0.0, "tlow defaulted to {}.", log)
  pyrat.ret.thigh = pt.defaultp(inputs.thigh, np.inf,
                                "thigh defaulted to {}.", log)

  # Purely-MCMC variables:
  pyrat.ret.mcmcfile = inputs.mcmcfile
  pyrat.ret.walk     = inputs.walk
  pyrat.ret.pmin     = inputs.pmin
  pyrat.ret.pmax     = inputs.pmax
  pyrat.ret.nsamples = inputs.nsamples
  pyrat.ret.burnin   = inputs.burnin
  pyrat.ret.thinning = pt.defaultp(inputs.thinning, 1,
      "thinning defaulted to {}.", log)
  pyrat.ret.nchains  = inputs.nchains
  pyrat.ret.grbreak  = pt.defaultp(inputs.grbreak, 0.0,
      "grbreak defaulted to {}.", log)
  pyrat.ret.grnmin   = pt.defaultp(inputs.grnmin, 0.5,
      "grnmin defaulted to {}.", log)
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

  if np.abs(pyrat.ret.qcap-0.5) > 0.5:
    log.error("Trace abundances cap (qcap={:.3f}) must lie in the range "
             "between 0.0 and 1.0.".format(pyrat.ret.qcap))
  if atm.tmodelname == "TCEA":
    if pyrat.phy.rstar is None:
      log.error("Undefined stellar radius (rstar), required for temperature "
               "model.")
    if pyrat.phy.tstar is None:
      log.error("Undefined stellar temperature (tstar), required for "
               "temperature model.")
    if pyrat.phy.smaxis is None:
      log.error("Undefined orbital semi-major axis (smaxis), required for "
               "temperature model.")
    if pyrat.phy.gplanet is None:
      log.error("Undefined planetary surface gravity (gplanet), required for "
               "temperature model.")

  # Number of processors:
  if inputs.ncpu is None:
      pyrat.ncpu = 1
  else:
      pyrat.ncpu = pt.getparam('ncpu', inputs.ncpu, "none", log, integer=True)
  isgreater(pyrat.ncpu, "none", 1, False,
            "The number of processors ({:d}) must be >= 1.", log)
  if pyrat.ncpu >= mp.cpu_count():
    log.warning("The number of requested CPUs ({:d}) is >= than the number "
       "of available CPUs ({:d}).  Enforced ncpu to {:d}.".
       format(pyrat.ncpu, mp.cpu_count(), mp.cpu_count()-1))
    pyrat.ncpu = mp.cpu_count() - 1
  log.msg("Check inputs done.")


def isgreater(value, units, thresh, equal=False, text="", log=None):
  """
  Check that value (if not None) is greater than thresh.
  Throw error if not.

  Parameters
  ----------
  value: Scalar
    The value being tested.
  units: String
    The units of the value.
  thresh: Scalar
    Threshold against which value is being compared.
  equal: Boolean
    If True, strictly require greater than.
  text: String
    Text to show if condition is not satisfied.
  log: File
    Pyrat screen-output log file.

  Returns
  -------
  The value in the pyrat units.
  """
  # Set comparison command:
  if equal:
    compare = np.less_equal
  else:
    compare = np.less

  # Check value:
  if value is None:
    return

  if compare(value, thresh):
    log.error(text.format(value/pt.u(units)), tracklev=-3)


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
  if atm.bulk is not None  and  len(np.setdiff1d(atm.bulk, species)) > 0:
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
  if (atm.bulk is not None  and  atm.molfree is not None  and
      len(np.intersect1d(atm.bulk, atm.molfree)) > 0):
    log.error("These species were marked as both bulk and variable-abundance: "
             "{:s}.".format(np.intersect1d(atm.bulk, atm.molfree)))

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
      log.error("Species abundances included for retrieval (retflag contains "
               "'mol') but there are no abundance model (molmodel).")

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
      log.error("Undefined stellar temperature (tstar), required for Kurucz "
               "model.")
    if phy.gstar is None:
      log.error("Undefined stellar gravity (gstar), required for Kurucz model.")
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
    nhaze += pyrat.haze.model[i].npars

  # Indices to parse the array of fitting parameters:
  if ret.retflag is None:
    ret.retflag = []
  nparams = 0
  ret.pnames, ret.texnames = [], []
  if "pt" in ret.retflag:
    ret.itemp  = np.arange(nparams, nparams + ntemp)
    ret.pnames    += tpnames
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
    ret.pnames    += rpnames
    ret.texnames += rtexnames
    nparams += nray
  if "haze" in ret.retflag:
    ret.ihaze  = np.arange(nparams, nparams + nhaze)
    ret.pnames   += hpnames
    ret.texnames += htexnames
    nparams += nhaze
  #if "cloud" in ret.retflag:
  #  ret.icloud  = np.arange(nparams, nparams + ncloud)
  #  ret.pnames    += cpnames
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
  if (pyrat.rayleigh.nmodels > 0 and
      (pyrat.runmode != "mcmc" or "ray" not in ret.retflag)):
    if pyrat.rayleigh.pars is None:
      log.error("Rayleigh parameters (rpars) have not been specified.")
  if (pyrat.haze.nmodels > 0 and
      (pyrat.runmode != "mcmc" or "haze" not in ret.retflag)):
    if pyrat.haze.pars is None:
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
