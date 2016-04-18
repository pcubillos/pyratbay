#!/usr/bin/env python

# FINDME a LICENSE

import os
import sys
import time
import shutil
import subprocess
import numpy as np

from .. import tools     as pt
from .. import constants as pc
from .  import argum     as ar
from .  import makeatm   as ma
from .  import makecfg   as mc

rootdir = os.path.realpath(os.path.dirname(__file__) + "/../../")
TEAdir = rootdir + "/modules/TEA/"

sys.path.append(rootdir + '/pyratbay/lib/')
import pt as PT


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

  ar.checkinputs(args, log, wlog)

  # Unpack pressure input variables:
  punits  = args.punits
  ptop    = pt.getparam(args.ptop,    args.punits)
  pbottom = pt.getparam(args.pbottom, args.punits)
  if ptop >= pbottom:
    pt.error("Bottom-layer pressure ({:.2e} bar) must be higher than the"
      "top-layer pressure ({:.2e} bar).".format(pbottom/pt.u("bar"),
                                                ptop/pt.u("bar")))
  nlayers = args.nlayers

  # Create pressure array in barye (CGS) units:
  pressure = np.logspace(np.log10(ptop), np.log10(pbottom), nlayers)
  pt.msg(1, "Creating {:d}-layer atmospheric model between {:.1e} "
      "and {:.1e} bar.".format(nlayers, ptop/pc.bar, pbottom/pc.bar), log)

  # Compute the temperature profile:
  if args.tmodel == "TCEA":
    rstar   = pt.getparam(args.rstar, args.radunits)
    tstar   = pt.getparam(args.tstar, "kelvin")
    tint    = pt.getparam(args.tint,  "kelvin")
    pgrav   = pt.getparam(args.pgrav, "none")
    smaxis  = pt.getparam(args.smaxis, args.radunits)
    tparams = args.tparams
    temp    = PT.TCEA(tparams, pressure, rstar, tstar, tint, smaxis, pgrav)
    pt.msg(1, "\nComputed TCEA temperature model.", log)
  elif args.tmodel == "isothermal":
    temp = np.tile(tparams[0], nlayers)
    pt.msg(1, "\nComputed isothermal temperature model.", log)

  # Make the atmospheric file:
  atmfile  = args.atmfile
  uniform  = args.uniform
  species  = args.species
  elements = args.elements
  solar    = args.solar
  atomicfile = args.atomicfile
  xsolar   = pt.getparam(args.xsolar, "none")
  patm     = args.patm

  # Uniform-abundances profile:
  if uniform is not None:
    if len(uniform) != len(species):
      pt.error("Number of uniform abundances ({:d}) does not match the "
        "number of species ({:d}).".format(len(uniform), len(species)), log)
    ma.uniform(atmfile, pressure, temp, species, uniform, punits)
    pt.msg(1, "\nProduced uniform-abundances atmospheric file: '{:s}'.".
              format(atmfile), log)
  # TEA abundances:
  else:
    pt.msg(1, "\nRun TEA to compute thermochemical-equilibrium "
              "abundances.", log)
    swap   = None
    ma.makeatomic(solar, atomicfile, xsolar, swap)
    # Pre-atmospheric file:
    ma.makepreatm(pressure/pt.u(punits), temp, atomicfile,
                  elements, species, patm)
    # Run TEA:
    mc.makeTEA(abun_file=atomicfile)
    proc = subprocess.Popen([TEAdir + "tea/runatm.py", patm, "TEA"])
    proc.communicate()
    # Reformat the TEA output into the pyrat format:
    ma.TEA2pyrat("./TEA/TEA/results/TEA.tea", atmfile)
    shutil.rmtree("TEA")
    pt.msg(1, "Produced TEA atmospheric file: '{:s}'.".format(atmfile), log)

  # Compute an opacity grid:
  pass

  # Full Pyrat Bay run:
  if False:
    # FINDME: So in principle, I could do the initialization here ...
    pass
    # Run MCMC:
    posterior, bestp = mc3.mcmc(data=data, func=func, indparams=indparams,
                    params=params,
                    numit=numit, nchains=nchains, walk=walk, grtest=grtest,
                    leastsq=leastsq, chisqscale=chisqscale,
                    burnin=burnin, plots=plots, savefile=savefile,
                    savemodel=savemodel)

  # Post processing:
  pass

  log.close()
