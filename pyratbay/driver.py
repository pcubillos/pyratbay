# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["run"]

import os
import sys

import numpy as np

from . import tools      as pt
from . import constants  as pc
from . import lineread   as lr
from . import plots      as pp
from . import atmosphere as pa
from . import io         as io
from .pyrat import Pyrat

sys.path.append(pc.ROOT + "modules/MCcubed/")
import MCcubed as mc3


@pt.ignore_system_exit
def run(cfile, init=False):
  """
  Pyrat Bay initialization driver.

  Parameters
  ----------
  cfile: String
      A Pyrat Bay configuration file.
  init: Bool
      If True, only initialize a Pyrat object (no spectra calculation).
      This is useful when computing spectra interactively.
  """
  pyrat = Pyrat(cfile)
  log = pyrat.log
  phy = pyrat.phy
  atm = pyrat.atm
  ret = pyrat.ret
  inputs = pyrat.inputs

  # Call lineread package:
  if pyrat.runmode == 'tli':
      if pyrat.lt.tlifile is None:
          log.error('Undefined TLI file (tlifile).')
      lr.makeTLI(inputs.dblist, inputs.pflist, inputs.dbtype,
                 pyrat.lt.tlifile[0], pyrat.spec.wllow, pyrat.spec.wlhigh,
                 pyrat.spec.wlunits, log)
      return


  # Get gplanet from mplanet and rplanet if necessary:
  mplanet = phy.mplanet is not None
  gplanet = phy.gplanet is not None
  rplanet = phy.rplanet is not None

  # Check planetary surface gravity/mass/radius:
  if mplanet and rplanet and gplanet:
      gplanet = pc.G * phy.mplanet / phy.rplanet**2
      if np.abs(gplanet-phy.gplanet)/phy.gplanet > 0.05:
          log.error("All mplanet, rplanet, and gplanet were provided, but "
              "values are inconsistent (>5%): g(M,R) = {:7.1f} cm s-2 and "
              "gplanet = {:7.1f} cm s-2.".format(gplanet, phy.gplanet))
  elif not mplanet and rplanet and gplanet:
      phy.mplanet = phy.gplanet * phy.rplanet**2 / pc.G
  elif mplanet and not rplanet and gplanet:
      phy.rplanet = np.sqrt(pc.G * phy.mplanet / phy.gplanet)
  elif mplanet and rplanet and not gplanet:
      phy.gplanet = pc.G * phy.mplanet / phy.rplanet**2


  # Compute pressure-temperature profile:
  if pyrat.runmode in ['pt', 'atmosphere'] or pt.isfile(atm.atmfile) != 1:
      if pt.isfile(atm.ptfile) == 1:
          log.msg("\nReading pressure-temperature file: '{:s}'.".
                  format(atm.ptfile))
          pressure, temperature = io.read_pt(atm.ptfile)
      else:
          check_pressure(pyrat)
          pressure = pa.pressure(atm.ptop, atm.pbottom, atm.nlayers,
              'barye', log)
          check_temp(pyrat)
          temperature = pa.temperature(atm.tmodelname, pressure,
               phy.rstar, phy.tstar, phy.tint, phy.gplanet, phy.smaxis,
               atm.runits, atm.nlayers, log, atm.tpars)

  # Return temperature-pressure if requested:
  if pyrat.runmode == 'pt':
      return pressure, temperature


  # Compute atmospheric abundances:
  if pyrat.runmode == 'atmosphere' or pt.isfile(atm.atmfile) != 1:
      check_atm(pyrat)
      abundances = pa.abundances(atm.atmfile, pressure, temperature,
          inputs.species, inputs.elements, inputs.uniform, atm.punits,
          inputs.xsolar, inputs.solar, log)

  # Return atmospheric model if requested:
  if pyrat.runmode == 'atmosphere':
      return pressure, temperature, abundances


  # Check status of extinction-coefficient file if necessary:
  if pyrat.runmode != 'spectrum' and pt.isfile(pyrat.ex.extfile) == -1:
      log.error("Undefined extinction-coefficient file (extfile).")

  # Force to re-calculate extinction-coefficient file if requested:
  if pyrat.runmode == 'opacity' and pt.isfile(pyrat.ex.extfile) == 1:
      os.remove(pyrat.ex.extfile)

  if pyrat.runmode == 'mcmc' and pyrat.ret.mcmcfile is None:
      log.error('Undefined MCMC file (mcmcfile).')

  # Initialize pyrat object:
  pyrat.setup_spectrum()
  #if pyrat.inputs.resume: # Bypass writting all of the initialization log:
  #    pyrat = Pyrat(args, log=None)
  #    pyrat.log = log
  #else:
  #    pyrat = Pyrat(args, log=log)

  # Stop and return if requested:
  if init or pyrat.runmode == 'opacity':
      return pyrat

  # Compute spectrum and return pyrat object if requested:
  if pyrat.runmode == "spectrum":
      pyrat.run()
      return pyrat


  muted_log = mc3.utils.Log(None, verb=0, width=80)
  pyrat.log = muted_log      # Mute logging in PB, but not in MC3
  pyrat.spec.outspec = None  # Avoid writing spectrum file during MCMC
  retmodel = False  # Return only the band-integrated spectrum
  # Basename of the output files (no path, no extension):
  outfile = os.path.splitext(os.path.basename(pyrat.ret.mcmcfile))[0]

  # Run MCMC:
  mc3_out = mc3.mcmc(data=pyrat.obs.data, uncert=pyrat.obs.uncert,
      func=pyrat.eval, indparams=[retmodel], params=ret.params,
      pmin=ret.pmin, pmax=ret.pmax, stepsize=ret.stepsize,
      prior=ret.prior, priorlow=ret.priorlow, priorup=ret.priorup,
      walk=ret.walk, nsamples=ret.nsamples,
      nchains=ret.nchains, burnin=ret.burnin, thinning=ret.thinning,
      grtest=True, grbreak=ret.grbreak, grnmin=ret.grnmin,
      hsize=10, kickoff='normal', log=log, nproc=pyrat.ncpu,
      plots=True, showbp=False,
      pnames=ret.pnames, texnames=ret.texnames,
      resume=inputs.resume, savefile=ret.mcmcfile)

  if mc3_out is None:
      log.error("Error in MC3.")

  bestp, CRlo, CRhi, stdp, posterior, Zchain = mc3_out
  pyrat.ret.posterior = posterior
  pyrat.ret.bestp     = bestp

  # Best-fitting model:
  pyrat.spec.outspec = "{:s}_bestfit_spectrum.dat".format(outfile)
  pyrat.ret.spec_best, pyrat.ret.bestbandflux = pyrat.eval(bestp)

  header  = "# MCMC best-fitting atmospheric model.\n\n"
  bestatm = "{:s}_bestfit_atmosphere.atm".format(outfile)
  pa.writeatm(bestatm, pyrat.atm.press, pyrat.atm.temp,
              pyrat.mol.name, pyrat.atm.q, pyrat.atm.punits,
              header, radius=pyrat.atm.radius, runits='km')

  pyrat.plot_spectrum(spec='best',
      filename='{:s}_bestfit_spectrum.png'.format(outfile))

  if pyrat.atm.tmodelname in ['tcea', 'madhu_inv', 'madhu_noinv']:
      pyrat.plot_posterior_pt('{:s}_posterior_PT_profile.png'.format(outfile))

  if pyrat.od.path == "eclipse":
      cf  = pt.cf(pyrat.od.depth, pyrat.atm.press, pyrat.od.B)
      bcf = pt.bandcf(cf, pyrat.obs.bandtrans, pyrat.spec.wn, pyrat.obs.bandidx)
  elif pyrat.od.path == "transit":
      transmittance = pt.transmittance(pyrat.od.depth, pyrat.od.ideep)
      bcf = pt.bandcf(transmittance, pyrat.obs.bandtrans, pyrat.spec.wn,
                      pyrat.obs.bandidx)
  pp.cf(bcf, 1.0/(pyrat.obs.bandwn*pc.um), pyrat.od.path,
        pyrat.atm.press, pyrat.atm.radius,
        pyrat.atm.rtop, filename="{:s}_bestfit_cf.png".format(outfile))

  pyrat.log = log  # Un-mute
  log.msg("\nOutput MCMC posterior results, log, bestfit atmosphere, "
          "and spectrum:\n'{:s}.npz',\n'{:s}',\n'{:s}',\n'{:s}'.\n\n".
          format(outfile, os.path.basename(inputs.logfile), bestatm,
                 pyrat.spec.outspec))
  return pyrat


def check_pressure(pyrat):
  """
  Check the input arguments to calculate the pressure profile.
  """
  if pyrat.atm.nlayers is None:
      pyrat.log.error("Undefined number of atmospheric layers (nlayers).")
  if pyrat.atm.ptop is None:
      pyrat.log.error("Undefined atmospheric top pressure (ptop)")
  if pyrat.atm.pbottom is None:
      pyrat.log.error("Undefined atmospheric bottom pressure (pbottom)")


def check_temp(pyrat):
  """
  Check the input arguments to calculate the temperature profile.
  """
  log = pyrat.log
  atm = pyrat.atm
  if atm.tmodelname is None:
      log.error("Undefined temperature model (tmodel).")
  if atm.tpars is None:
      log.error("Undefined temperature-model parameters (tpars).")

  if atm.tmodelname == 'isothermal':
      if len(atm.tpars) != 1:
          log.error("Wrong number of parameters ({:d}) for the isothermal "
                    "temperature model (1).".format(len(atm.tpars)))

  elif atm.tmodelname == 'tcea':
      if len(atm.tpars) != 5:
          log.error("Wrong number of parameters ({:d}) for the tcea "
                    "temperature model (5).".format(len(atm.tpars)))
      if pyrat.phy.rstar is None:
          log.error("Undefined stellar radius (rstar).")
      if pyrat.phy.tstar is None:
          log.error("Undefined stellar temperature (tstar).")
      if pyrat.phy.smaxis is None:
          log.error("Undefined orbital semi-major axis (smaxis).")
      if pyrat.phy.gplanet is None:
          log.error("Undefined planetary surface gravity, set either "
                    "gplanet or mplanet and rplanet.")


def check_atm(pyrat):
  """
  Check the input arguments to calculate the atmospheric model.
  """
  atm = pyrat.atm
  if atm.atmfile is None:
      pyrat.log.error("Undefined atmospheric file (atmfile).")
  if pyrat.inputs.species is None:
      pyrat.log.error("Undefined atmospheric species list (species).")
  # Uniform-abundances profile:
  if pyrat.inputs.uniform is not None:
      if len(pyrat.inputs.uniform) != len(pyrat.inputs.species):
          pyrat.log.error("Number of uniform abundances ({:d}) does not "
                          "match the number of species ({:d}).".
                          format(len(pyrat.inputs.uniform), len(pyrat.inputs.species)))
      return
  # TEA abundances:
  if pyrat.inputs.elements is None:
      pyrat.log.error("Undefined atmospheric atomic composition (elements).")
  pyrat.inputs.solar = pyrat.inputs.get_default('solar', 'Solar-abundance file',
      pc.ROOT+'inputs/AsplundEtal2009.txt')
  pyrat.inputs.atomicfile = pyrat.inputs.get_default('atomicfile',
      'Atomic-composition file', './atomic.tea')
  pyrat.inputs.patm = pyrat.inputs.get_default('patm',
      'Pre-atmospheric file', './preatm.tea')
  pyrat.inputs.xsolar = pyrat.inputs.get_default('xsolar',
      'Solar-metallicity scaling factor', 1.0, gt=0.0, wflag=True)
