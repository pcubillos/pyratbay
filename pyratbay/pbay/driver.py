# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt

from .. import tools      as pt
from .. import constants  as pc
from .. import plots      as pp
from .. import pyrat      as py
from .. import lineread   as lr
from .. import atmosphere as pa

from .  import argum     as ar
from .  import pyratfit  as pf

rootdir = os.path.realpath(os.path.dirname(__file__) + "/../../")
MC3dir = rootdir + "/modules/MCcubed/"
sys.path.append(MC3dir)
import MCcubed as mc3

__all__ = ["run"]


def run(argv, main=False):
  """
  Pyrat Bay (Python Radiative Transfer in a Bayesian framework)
  initialization driver.

  Parameters
  ----------
  argv: List or string
     If called from the shell, the list of command line arguments; if
     called from the Python interpreter, the configuration-file name.
  main: Bool
     Flag to indicate if Pyrat was called from the shell (True) or from
     the Python interpreter.
  """

  # Put everything into a try--except to catch the sys.exit() Traceback.
  try:
    # Setup the command-line-arguments input:
    if main is False:
      sys.argv = ['pbay.py', '-c', argv]

    # Warnings log:
    wlog = []

    # Setup time tracker:
    timestamps = []
    timestamps.append(time.time())

    # Parse command line arguments:
    args, log = ar.parse(wlog)
    timestamps.append(time.time())

    # Check run mode:
    if args.runmode not in pc.rmodes:
      pt.error("Invalid runmode ({:s}). Select from: {:s}.".
                format(args.runmode, str(pc.rmodes)))

    # Call lineread package:
    if args.runmode == "tli":
      parser = lr.parser()
      lr.makeTLI(parser.dblist,  parser.pflist, parser.dbtype,
                 parser.outfile, parser.iwl, parser.fwl, args.verb)
      return

    # Get gplanet from mplanet and rplanet if necessary:
    if (args.gplanet is None and args.rplanet is not None and
        args.mplanet is not None):
      args.gplanet = (pc.G * pt.getparam(args.mplanet, "gram") /
                      pt.getparam(args.rplanet, args.radunits)**2)

    # Compute pressure-temperature profile:
    if args.runmode in ["pt", "atmosphere"] or pt.isfile(args.atmfile) != 1:
      # Check if PT file is provided:
      if args.ptfile is None:
        ar.checkpressure(args, log, wlog)  # Check pressure inputs
        pressure = pa.pressure(args.ptop, args.pbottom, args.nlayers,
                               args.punits, log)
        ar.checktemp(args, log, wlog)      # Check temperature inputs
        temperature = pa.temperature(args.tmodel, pressure,
           args.rstar, args.tstar, args.tint, args.gplanet, args.smaxis,
           args.radunits, args.nlayers, log, args.tparams)

      # If PT file is provided, read it:
      elif os.path.isfile(args.ptfile):
        pt.msg(args.verb-3, "\nReading pressure-temperature file:"
               " '{:s}'.".format(args.ptfile), log)
        pressure, temperature = pa.read_ptfile(args.ptfile)

    # Return temperature-pressure if requested:
    if args.runmode == "pt":
      plt.figure(1)
      plt.semilogy(temperature, pressure)
      plt.ylim(np.max(pressure), np.min(pressure))
      plt.xlabel("Temperature  (K)")
      plt.ylabel("Pressure  (barye)")
      plt.savefig("tmp_pt.pdf")
      return pressure, temperature

    # Compute or read atmospheric abundances:
    if args.runmode == "atmosphere" or pt.isfile(args.atmfile) != 1:
      ar.checkatm(args, log, wlog)
      xsolar = pt.getparam(args.xsolar, "none")
      pa.calcatm(args.atmfile, pressure, temperature, args.species,
                 args.punits, args.solar, xsolar, args.uniform,
                 args.elements, 1, log, int(args.verb>0))

    # Return atmospheric model if requested:
    if args.runmode == "atmosphere":
      species, pressure, temperature, abundances = pa.readatm(args.atmfile)
      return pressure, temperature, abundances

    # Check status of extinction-coefficient file if necessary:
    if args.runmode != "spectrum" and pt.isfile(args.extfile) == -1:
      pt.error("Unspecified extinction-coefficient file (extfile).", log)

    # Force to re-calculate extinction-coefficient file if requested:
    if args.runmode == "opacity" and pt.isfile(args.extfile):
      os.remove(args.extfile)

    # Initialize pyrat object:
    if args.resume: # Bypass writting all of the initialization log:
      nolog = open("deleteme.log", "w")
      pyrat = py.init(args.cfile, log=nolog)
      nolog.close()
      os.remove("deleteme.log")
      pyrat.log = log
    else:
      pyrat = py.init(args.cfile, log=log)

    # Compute spectrum and return pyrat object if requested:
    if args.runmode == "spectrum":
      pyrat = py.run(pyrat)
      return pyrat

    # End if necessary:
    if args.runmode == "opacity":
      return pyrat

    # Parse retrieval info into the Pyrat object:
    pf.init(pyrat, args, log)
    verb = pyrat.verb     # Mute pyrat
    pyrat.verb = 0        # Mute pyrat
    pyrat.outspec = None  # Avoid writing spectrum file during MCMC

    # Basename of the output files:
    outfile = os.path.splitext(os.path.basename(log.name))[0]
    # Run MCMC:
    freeze   = True  # Freeze abundances evoer iterations
    retmodel = False # Return only the band-integrated spectrum
    mc3_out = mc3.mcmc(data=args.data, uncert=args.uncert,
           func=pf.fit, indparams=[pyrat,freeze,retmodel], params=args.params,
           pmin=args.pmin, pmax=args.pmax, stepsize=args.stepsize,
           prior=args.prior, priorlow=args.priorlow, priorup=args.priorup,
           walk=args.walk, nsamples=args.nsamples, nchains=args.nchains,
           burnin=args.burnin, thinning=args.thinning,
           grtest=True, grbreak=args.grbreak, grnmin=args.grnmin,
           hsize=10, kickoff='normal', log=log, nproc=args.nproc,
           plots=True, parname=pyrat.ret.parname, showbp=False,
           resume=args.resume, savefile="{:s}.npz".format(outfile))

    if mc3_out is None:
      pt.error("Error in MC3.", pyrat.log)
    else:
      bestp, CRlo, CRhi, stdp, posterior, Zchain = mc3_out

    # Best-fitting model:
    pyrat.outspec = "{:s}_bestfit_spectrum.dat".format(outfile)
    bestbandflux = pf.fit(bestp, pyrat, retmodel=False)

    # Best-fit atmfile header:
    header = "# MCMC best-fitting atmospheric model.\n\n"
    # Write best-fit atmfile:
    bestatm = "{:s}_bestfit_atmosphere.atm".format(outfile)
    pa.writeatm(bestatm, pyrat.atm.press, pyrat.atm.temp,
                pyrat.mol.name, pyrat.atm.q, pyrat.atm.punits,
                header, radius=pyrat.atm.radius, runits='km')

    pyrat.verb = verb  # Un-mute

    # Best-fitting spectrum:
    pp.spectrum(pyrat=pyrat, logxticks=args.logxticks, yran=args.yran,
                filename="{:s}_bestfit_spectrum.png".format(outfile))
    # Posterior PT profiles:
    if pyrat.ret.tmodelname in ["TCEA", "MadhuInv", "MadhuNoInv"]:
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

    pt.msg(pyrat.verb-3, "\nOutput MCMC posterior results, log, bestfit "
      "atmosphere, and spectrum:\n'{:s}.npz',\n'{:s}',\n'{:s}',\n'{:s}'.\n\n".
      format(outfile, os.path.basename(args.logfile), bestatm,
             pyrat.outspec), pyrat.log, 0)
    log.close()
    return pyrat, bestp

  # Avoid printing to screeen the System-Exit Traceback error:
  except SystemExit:
    return None
