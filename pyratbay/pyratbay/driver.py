import os
import sys
import time
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt

from .. import tools     as pt
from .. import constants as pc
from .. import pyrat     as pyrat
from .. import lineread  as lr

from .  import argum     as ar
from .  import makeatm   as ma
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

  # Call lineread package:
  if args.runmode == "tli":
    parser = lr.parser()
    lr.makeTLI(parser.dblist,  parser.pflist, parser.dbtype,
               parser.outfile, parser.iwl, parser.fwl, parser.verb)
    return

  # Compute pressure-temperature profile:
  if args.runmode in ["pt", "atmosphere"] or pt.isfile(args.atmfile) != 1:
    pressure    = calcp(args, log, wlog)
    temperature = ma.calct(args, pressure, log, wlog)

  # Return temperature-pressure if requested:
  if args.runmode == "pt":
    plt.figure(1)
    plt.semilogy(temperature, pressure)
    plt.ylim(np.max(pressure), np.min(pressure))
    plt.savefig("tmp_pt.pdf")
    return pressure, temperature

  # Compute or read atmospheric abundances:
  if args.runmode == "atmosphere" or pt.isfile(args.atmfile) != 1:
    calcatm(args, pressure, temperature, log, wlog)

  # Return atmospheric model if requested:
  if args.runmode == "atmosphere":
    species, pressure, temperature, abundances = ma.readatm(args.atmfile)
    return pressure, temperature, abundances

  # Check status of extinction-coefficient file if necessary:
  if args.runmode != "spectrum" and pt.isfile(args.extfile) == -1:
    pt.error("Unspecified extinction-coefficient file (extfile).", log)

  # Force to re-calculate extinction-coefficient file if requested:
  if args.runmode == "opacity":
    # FINDME: os.remove(args.extfile)
    pass

  # Initialize pyrat object:
  py = pyrat.init(args.cfile)

  # Compute spectrum and return pyrat object if requested:
  if args.runmode == "spectrum":
    py = pyrat.run(py)
    return py

  # End if necessary:
  if args.runmode == "opacity":
    return

  # Parse retrieval into the Pyrat object:
  pf.init(py, args, log)
  return py, args, pf.fit


  # Full Pyrat Bay run:
  if False:
    # Run MCMC:
    bestp, uncertp, posterior, Zchain = mc3.mcmc(data=args.data,
         uncert=args.uncert,
         func=pf.fit, indparams=indparams, params=params,
         nsamples=nsamples, nchains=nchains, walk="snooker", grtest=grtest,
         burnin=burnin, plots=plots, savefile=savefile)

  # Post processing:
  pass

  #log.close()


def calcp(args, log, wlog):
  """
  Calculate a pressure profile.

  Parameters
  ----------
  args:
  log:
  wlog:

  Returns
  -------
  pressure: 1D float ndarray
    Atmospheric pressure profile in Barye units.
  """
  # Check pressure inputs:
  ar.checkpressure(args, log, wlog)
  # Unpack pressure input variables:
  punits  = args.punits
  ptop    = pt.getparam(args.ptop,    args.punits)
  pbottom = pt.getparam(args.pbottom, args.punits)
  if ptop >= pbottom:
    pt.error("Bottom-layer pressure ({:.2e} bar) must be higher than the"
      "top-layer pressure ({:.2e} bar).".format(pbottom/pt.u("bar"),
                                                ptop/pt.u("bar")))

  # Create pressure array in barye (CGS) units:
  pressure = np.logspace(np.log10(ptop), np.log10(pbottom), args.nlayers)
  pt.msg(1, "Creating {:d}-layer atmospheric model between {:.1e} "
      "and {:.1e} bar.".format(args.nlayers, ptop/pc.bar, pbottom/pc.bar), log)
  return pressure


def calcatm(args, pressure, temperature, log, wlog):
  """
  Compute atmospheric abundaces for given pressure, temperature profile:
  """
  ar.checkatm(args, log, wlog)

  # Uniform-abundances profile:
  if args.uniform is not None:
    ma.uniform(args.atmfile, pressure, temperature, args.species,
               args.uniform, args.punits)
    pt.msg(1, "\nProduced uniform-abundances atmospheric file: '{:s}'.".
              format(args.atmfile), log)
  # TEA abundances:
  else:
    pt.msg(1, "\nRun TEA to compute thermochemical-equilibrium "
              "abundances.", log)
    xsolar = pt.getparam(args.xsolar, "none")
    swap   = None
    ma.makeatomic(args.solar, args.atomicfile, xsolar, swap)
    # Pre-atmospheric file:
    ma.makepreatm(pressure/pt.u(args.punits), temperature, args.atomicfile,
                  args.elements, args.species, args.patm)
    # Run TEA:
    mc.makeTEA(abun_file=args.atomicfile)
    proc = subprocess.Popen([TEAdir + "tea/runatm.py", args.patm, "TEA"])
    proc.communicate()
    # Reformat the TEA output into the pyrat format:
    ma.TEA2pyrat("./TEA/TEA/results/TEA.tea", args.atmfile)
    shutil.rmtree("TEA")
    pt.msg(1, "Produced TEA atmospheric file '{:s}'.".format(args.atmfile), log)

