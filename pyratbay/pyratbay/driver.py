import os
import sys
import time
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt

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

  #ar.checkinputs(args, log, wlog)

  # Call lineread package:
  if args.runmode == "tli":
    return

  # Compute pressure-temperature profile:
  if args.runmode in ["pt", "atmosphere"] or pt.isfile(args.atmfile) != 1:
    pressure    = calcp(args, log, wlog)
    temperature = calct(args, pressure, log, wlog)

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

  #log.close()


def calcp(args, log, wlog):
  """
  Calculate the pressure profile.
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


def calct(args, pressure, log, wlog):
  """
  Calculate the temperature profile.
  """
  ar.checktemp(args, log, wlog)
  if args.tmodel == "TCEA":
    rstar   = pt.getparam(args.rstar, args.radunits)
    tstar   = pt.getparam(args.tstar, "kelvin")
    tint    = pt.getparam(args.tint,  "kelvin")
    gplanet = pt.getparam(args.gplanet, "none")
    smaxis  = pt.getparam(args.smaxis, args.radunits)
    temperature = PT.TCEA(args.tparams, pressure,
                          rstar, tstar, tint, smaxis, gplanet)
    pt.msg(1, "\nComputed TCEA temperature model.", log)
  elif args.tmodel == "isothermal":
    temperature = np.tile(args.tparams[0], args.nlayers)
    pt.msg(1, "\nComputed isothermal temperature model.", log)

  return temperature


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
    pt.msg(1, "Produced TEA atmospheric file: '{:s}'.".format(args.atmfile), log)

