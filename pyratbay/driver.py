# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["run"]

import os
import sys

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
  Pyrat Bay (Python Radiative Transfer in a Bayesian framework)
  initialization driver.

  Parameters
  ----------
  cfile: String
      A Pyrat Bay configuration file.
  init: Bool
      If True, only initialize a Pyrat object (no spectra calculation).
      This is useful when computing spectra interactively.
  """
  # Parse command line arguments:
  args, log = pt.parse(cfile)

  # Check run mode:
  if args.runmode not in pc.rmodes:
      log.error("Invalid runmode ({}). Select from: {:s}.".
                format(args.runmode, str(pc.rmodes)))

  # Call lineread package:
  if args.runmode == "tli":
      if args.tlifile is None:
          log.error('No output TLI file specified.')
      if args.wlunits is None:
          args.wlunits = 'um'
      lr.makeTLI(args.dblist, args.pflist, args.dbtype, args.tlifile[0],
                 args.wllow, args.wlhigh, args.wlunits, log)
      return


  # Get gplanet from mplanet and rplanet if necessary:
  if (args.gplanet is None and args.rplanet is not None and
      args.mplanet is not None):
      mplanet = args.get_param('mplanet', None,        'Planet mass',   gt=0.0)
      rplanet = args.get_param('rplanet', args.runits, 'Planet radius', gt=0.0)
      args.gplanet = pc.G * mplanet / rplanet**2

  # Compute pressure-temperature profile:
  if args.runmode in ["pt", "atmosphere"] or pt.isfile(args.atmfile) != 1:
      # Check if PT file is provided:
      if args.ptfile is None:
          check_pressure(args, log)
          pressure = pa.pressure(args.ptop, args.pbottom, args.nlayers,
                                 args.punits, log)
          check_temp(args, log)
          temperature = pa.temperature(args.tmodel, pressure,
               args.rstar, args.tstar, args.tint, args.gplanet, args.smaxis,
               args.runits, args.nlayers, log, args.tpars)
      # If PT file is provided, read it:
      elif os.path.isfile(args.ptfile):
          log.msg("\nReading pressure-temperature file: '{:s}'.".
                  format(args.ptfile))
          pressure, temperature = io.read_pt(args.ptfile)

  # Return temperature-pressure if requested:
  if args.runmode == "pt":
      return pressure, temperature


  # Compute atmospheric abundances:
  if args.runmode == "atmosphere" or pt.isfile(args.atmfile) != 1:
      check_atm(args, log)
      xsolar = args.get_default('xsolar', 'Solar-metallicity scale factor')
      abundances = pa.abundances(args.atmfile, pressure, temperature,
          args.species, args.elements, args.uniform, args.punits, xsolar,
          args.solar, log)

  # Return atmospheric model if requested:
  if args.runmode == "atmosphere":
      return pressure, temperature, abundances

  # Check status of extinction-coefficient file if necessary:
  if args.runmode != "spectrum" and pt.isfile(args.extfile) == -1:
      log.error("Unspecified extinction-coefficient file (extfile).")

  # Force to re-calculate extinction-coefficient file if requested:
  if args.runmode == "opacity" and pt.isfile(args.extfile) == 1:
      os.remove(args.extfile)

  if args.runmode == "mcmc" and args.mcmcfile is None:
      log.error('No MCMC file specified.')

  # Initialize pyrat object:
  if args.resume: # Bypass writting all of the initialization log:
      pyrat = Pyrat(args, log=None)
      pyrat.log = log
  else:
      pyrat = Pyrat(args, log=log)

  # Stop and return if requested:
  if init or pyrat.runmode == 'opacity':
      return pyrat

  # Compute spectrum and return pyrat object if requested:
  if args.runmode == "spectrum":
      pyrat.run()
      return pyrat

  # Retrieval checks:
  if pyrat.od.path == "eclipse" and pyrat.phy.starflux is None:
      log.error("Unspecified stellar flux model.")
  if pyrat.od.path == "eclipse" and pyrat.phy.rprs is None:
      log.error("Undefined Rp/Rs.")
  if pyrat.obs.filter is None:
      log.error("Undefined transmission filters.")
  if pyrat.obs.data is None:
      log.error("Undefined data.")
  if pyrat.obs.uncert is None:
      log.error("Undefined data uncertainties.")

  muted_log = mc3.utils.Log(None, verb=0, width=80)
  pyrat.log = muted_log    # Mute logging in PB, but not in MC3
  pyrat.spec.outspec = None  # Avoid writing spectrum file during MCMC
  retmodel = False  # Return only the band-integrated spectrum
  # Basename of the output files (no path, no extension):
  outfile = os.path.splitext(os.path.basename(pyrat.ret.mcmcfile))[0]
  ret = pyrat.ret

  # Run MCMC:
  mc3_out = mc3.mcmc(data=args.data, uncert=args.uncert,
      func=pyrat.eval, indparams=[retmodel], params=ret.params,
      pmin=ret.pmin, pmax=ret.pmax, stepsize=ret.stepsize,
      prior=ret.prior, priorlow=ret.priorlow, priorup=ret.priorup,
      walk=ret.walk, nsamples=ret.nsamples,
      nchains=ret.nchains, burnin=ret.burnin, thinning=ret.thinning,
      grtest=True, grbreak=ret.grbreak, grnmin=ret.grnmin,
      hsize=10, kickoff='normal', log=log, nproc=pyrat.ncpu,
      plots=True, showbp=False,
      pnames=ret.pnames, texnames=ret.texnames,
      resume=args.resume, savefile=ret.mcmcfile)

  if mc3_out is None:
      log.error("Error in MC3.")
  else:
      bestp, CRlo, CRhi, stdp, posterior, Zchain = mc3_out

  # Best-fitting model:
  pyrat.spec.outspec = "{:s}_bestfit_spectrum.dat".format(outfile)
  dummy = pyrat.eval(bestp, retmodel=False)

  # Best-fit atmfile header:
  header = "# MCMC best-fitting atmospheric model.\n\n"
  # Write best-fit atmfile:
  bestatm = "{:s}_bestfit_atmosphere.atm".format(outfile)
  pa.writeatm(bestatm, pyrat.atm.press, pyrat.atm.temp,
              pyrat.mol.name, pyrat.atm.q, pyrat.atm.punits,
              header, radius=pyrat.atm.radius, runits='km')

  # Best-fitting spectrum:
  pp.spectrum(pyrat=pyrat, logxticks=args.logxticks, yran=args.yran,
              filename="{:s}_bestfit_spectrum.png".format(outfile))
  # Posterior PT profiles:
  if pyrat.atm.tmodelname in ["TCEA", "MadhuInv", "MadhuNoInv"]:
      pp.PT(posterior, besttpars=bestp[pyrat.ret.itemp], pyrat=pyrat,
            filename="{:s}_PT_posterior_profile.png".format(outfile))
  # Contribution or transmittance functions:
  if   pyrat.od.path == "eclipse":
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
          format(outfile, os.path.basename(args.logfile), bestatm,
                 pyrat.spec.outspec))
  log.close()
  return pyrat, bestp


def check_pressure(args, log):
  """
  Check the input arguments to calculate the pressure profile.
  """
  args.punits = args.get_default('punits', 'Pressure units', 'bar')
  if args.nlayers is None:
      log.error("Undefined number of atmospheric layers (nlayers).")
  if args.ptop is None:
      log.error("Undefined atmospheric top pressure (ptop)")
  if args.pbottom is None:
      log.error("Undefined atmospheric bottom pressure (pbottom)")
  args.get_default('nlayers', 'Number of atmospheric layers', gt=0)


def check_temp(args, log):
  """
  Check the input arguments to calculate the temperature profile.
  """
  if args.tmodel is None:
      log.error("Undefined temperature model (tmodel).")
  if args.tpars is None:
      log.error("Undefined temperature-model parameters (tpars).")

  if args.tmodel == "TCEA":
      if len(args.tpars) != 5:
          log.error("Wrong number of parameters ({:d}) for the TCEA "
                    "temperature model (5).".format(len(args.tpars)))
      if args.rstar is None:
          log.error("Undefined stellar radius (rstar).")
      if args.tstar is None:
          log.error("Undefined stellar temperature (tstar).")
      if args.smaxis is None:
          log.error("Undefined orbital semi-major axis (smaxis).")
      if (args.gplanet is None and
          (args.mplanet is None or args.rplanet is None)):
          log.error("Undefined planetary surface gravity, set either "
                    "gplanet or mplanet and rplanet.")
      args.tint = args.get_default('tint',
          'Planetary internal temperature (K)', 100.0)
      args.runits = args.get_default('runits', 'Distance units', 'cm')

  elif args.tmodel == "isothermal":
      if len(args.tpars) != 1:
          log.error("Wrong number of parameters ({:d}) for the isothermal "
                    "temperature model (1).".format(len(args.tpars)))


def check_atm(args, log):
  """
  Check the input arguments to calculate the atmospheric model.
  """
  if args.atmfile is None:
      log.error("Undefined atmospheric file (atmfile).")
  if args.species is None:
      log.error("Undefined atmospheric species list (species).")
  args.punits = args.get_default('punits', 'Pressure units', 'bar')

  # Uniform-abundances profile:
  if args.uniform is not None:
      if len(args.uniform) != len(args.species):
          log.error("Number of uniform abundances ({:d}) does not match the "
                    "number of species ({:d}).".
                    format(len(args.uniform), len(args.species)))
      return
  else:  # TEA abundances:
      if args.elements is None:
          log.error("Undefined atmospheric atomic-composition list (elements).")
      args.solar = args.get_default('solar', 'Solar-abundance file',
          pc.ROOT+'inputs/AsplundEtal2009.txt')

      args.atomicfile = args.get_default('atomicfile',
          'Atomic-composition file', './atomic.tea')
      args.patm = args.get_default('patm',
          'Pre-atmospheric file', './preatm.tea')
      args.xsolar = args.get_default('xsolar',
          'Solar-metallicity scaling factor', 1.0)
