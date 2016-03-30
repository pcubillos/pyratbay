#!/usr/bin/env python

# FINDME a LICENSE

import os
import sys
import time
import shutil
import subprocess
import numpy as np

from .. import tools   as pt
from .  import argum   as ar
from .  import makeatm as ma
from .  import makecfg as mc

rootdir = os.path.realpath(os.path.dirname(__file__) + "/../../")
TEAdir = rootdir + "/modules/TEA/"

sys.path.append(rootdir + '/pyratbay/lib/')
import pt as PT


def run(cfile, flag, main=False):
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
    sys.argv = ['pbay.py', '-c', cfile]

  # Warnings log:
  wlog = []

  # Setup time tracker:
  timestamps = []
  timestamps.append(time.time())

  # Parse command line arguments:
  args, log = ar.parse(wlog)
  timestamps.append(time.time())

  # Compute an atmospheric file:

  # Unpack pressure input variables:
  punits  = args.punits
  ptop    = pt.getparam(args.ptop,    args.punits)
  pbottom = pt.getparam(args.pbottom, args.punits)
  nlayers = args.nlayers

  # Create pressure array in barye (CGS) units:
  pressure = np.logspace(np.log10(ptop), np.log10(pbottom), nlayers)

  # Compute the temperature profile:
  if args.tmodel == "TCEA":
    rstar   = pt.getparam(args.rstar, args.radunits)
    tstar   = pt.getparam(args.tstar, "kelvin")
    tint    = pt.getparam(args.tint,  "kelvin")
    gplanet = pt.getparam(args.surfgravity, "none")
    smaxis  = pt.getparam(args.smaxis, args.radunits)
    tparams = args.tparams
    temp    = PT.TCEA(tparams, pressure, rstar, tstar, tint, smaxis, gplanet)

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
    ma.uniform(atmfile, pressure, temp, species, uniform, punits)
  # TEA abundances:
  else:
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

  # Compute an opacity grid:
  pass

  # Full Pyrat Bay run:
  pass

  # Post processing:
  pass

