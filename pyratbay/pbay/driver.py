# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import time
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt

from .. import tools      as pt
from .. import constants  as pc
from .. import plots      as pp
from .. import pyrat      as py
from .. import lineread   as lr
from .. import atmosphere as atm

from .  import argum     as ar
from .  import makecfg   as mc
from .  import pyratfit  as pf

rootdir = os.path.realpath(os.path.dirname(__file__) + "/../../")
TEAdir = rootdir + "/modules/TEA/"
MC3dir = rootdir + "/modules/MCcubed/"

sys.path.append(MC3dir)
import MCcubed as mc3


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
  if args.runmode not in ['tli','pt','atmosphere','opacity','spectrum','mcmc']:
    pt.error("Invalid runmode ({:s}). Select from: 'tli', 'pt', 'atmosphere', "
             "'opacity', 'spectrum', 'mcmc'.".format(args.runmode))

  # Call lineread package:
  if args.runmode == "tli":
    parser = lr.parser()
    lr.makeTLI(parser.dblist,  parser.pflist, parser.dbtype,
               parser.outfile, parser.iwl, parser.fwl, parser.verb)
    return

  # Get gplanet from mplanet and rplanet if necessary:
  if (args.gplanet is None and args.rplanet is not None and
      args.mplanet is not None):
    args.gplanet = (pc.G * pt.getparam(args.mplanet, "gram") /
                    pt.getparam(args.rplanet, args.radunits)**2)

  # Compute pressure-temperature profile:
  if args.runmode in ["pt", "atmosphere"] or pt.isfile(args.atmfile) != 1:
    ar.checkpressure(args, log, wlog)  # Check pressure inputs
    pressure = atm.pressure(args.ptop, args.pbottom, args.nlayers,
                           args.punits, log)
    ar.checktemp(args, log, wlog)      # Check temperature inputs
    temperature = atm.temperature(args.tmodel, pressure,
       args.rstar, args.tstar, args.tint, args.gplanet, args.smaxis,
       args.radunits, args.nlayers, log, args.tparams)

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
    calcatm(args, pressure, temperature, log, wlog)

  # Return atmospheric model if requested:
  if args.runmode == "atmosphere":
    species, pressure, temperature, abundances = atm.readatm(args.atmfile)
    return pressure, temperature, abundances

  # Check status of extinction-coefficient file if necessary:
  if args.runmode != "spectrum" and pt.isfile(args.extfile) == -1:
    pt.error("Unspecified extinction-coefficient file (extfile).", log)

  # Force to re-calculate extinction-coefficient file if requested:
  if args.runmode == "opacity" and pt.isfile(args.extfile):
    os.remove(args.extfile)
    pass

  # Initialize pyrat object:
  pyrat = py.init(args.cfile)

  # Compute spectrum and return pyrat object if requested:
  if args.runmode == "spectrum":
    pyrat = py.run(pyrat)
    return pyrat

  # End if necessary:
  if args.runmode == "opacity":
    return

  # Parse retrieval into the Pyrat object:
  pf.init(pyrat, args, log)
  pyrat.verb = 0  # Mute pyrat

  # Run MCMC:
  bestp, uncertp, posterior, Zchain = mc3.mcmc(data=args.data,
         uncert=args.uncert,
         func=pf.fit, indparams=[pyrat, True], params=args.params,
         pmin=args.pmin, pmax=args.pmax, stepsize=args.stepsize,
         walk=args.walk, nsamples=args.nsamples, nchains=args.nchains,
         burnin=args.burnin, thinning=args.thinning, grtest=True,
         hsize=10, kickoff='normal', log='MCMC.log',
         plots=True, savefile="test.npz")

  # Best-fitting model:
  bestbandflux = pf.fit(bestp, pyrat)

  # Best-fitting spectrum:
  pp.spectrum(pyrat=pyrat, filename="bestfit_spectrum.png")
  # Posterior PT profiles:
  if pyrat.ret.tmodelname == "TCEA":
    pp.TCEA(posterior, besttpars=bestp[pyrat.ret.itemp], pyrat=pyrat)
  # Contribution functions:
  cf  = pt.cf(pyrat.od.depth, pyrat.atm.press, pyrat.od.B)
  bcf = pt.bandcf(cf, pyrat.obs.bandtrans, pyrat.obs.bandidx)
  pp.cf(bcf, 1.0/(pyrat.obs.bandwn*pc.um), pyrat.atm.press,
        filename="bestfit_cf.png")


  log.close()
  return pyrat, bestp


def calcatm(args, pressure, temperature, log, wlog):
  """
  Compute atmospheric abundaces for given pressure, temperature profile:
  """
  ar.checkatm(args, log, wlog)

  # Uniform-abundances profile:
  if args.uniform is not None:
    atm.uniform(args.atmfile, pressure, temperature, args.species,
               args.uniform, args.punits)
    pt.msg(1, "\nProduced uniform-abundances atmospheric file: '{:s}'.".
              format(args.atmfile), log)
  # TEA abundances:
  else:
    pt.msg(1, "\nRun TEA to compute thermochemical-equilibrium "
              "abundances.", log)
    xsolar = pt.getparam(args.xsolar, "none")
    swap   = None
    atm.makeatomic(args.solar, args.atomicfile, xsolar, swap)
    # Pre-atmospheric file:
    atm.makepreatm(pressure/pt.u(args.punits), temperature, args.atomicfile,
                  args.elements, args.species, args.patm)
    # Run TEA:
    mc.makeTEA(abun_file=args.atomicfile)
    proc = subprocess.Popen([TEAdir + "tea/runatm.py", args.patm, "TEA"])
    proc.communicate()
    # Reformat the TEA output into the pyrat format:
    atm.TEA2pyrat("./TEA/TEA/results/TEA.tea", args.atmfile)
    shutil.rmtree("TEA")
    pt.msg(1, "Produced TEA atmospheric file '{:s}'.".format(args.atmfile), log)

