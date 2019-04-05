# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

from .. import tools     as pt
from .. import constants as pc


def checkpressure(args, log):
  """
  Check the input arguments to calculate the pressure profile.
  """
  args.punits = pt.defaultp(args.punits, "bar",
     "punits input variable defaulted to '{:s}'.", log)
  args.nlayers = pt.defaultp(args.nlayers, 100,
     "Number of atmospheric-model layers defaulted to {:d}.", log)
  args.ptop = pt.defaultp(args.ptop, "1e-8 bar",
     "Atmospheric-model top-pressure boundary defaulted to {:s}.", log)
  args.pbottom = pt.defaultp(args.pbottom, "100 bar",
     "Atmospheric-model bottom-pressure boundary defaulted to {:s}.", log)


def checktemp(args, log):
  """
  Check the input arguments to calculate the temperature profile.
  """
  if args.tmodel is None:
    log.error("Undefined temperature model (tmodel).")
  if args.tpars is None:
    log.error("Undefined temperature-model parameters (tpars).")

  if args.tmodel == "TCEA":
    if len(args.tpars) != 5:
      log.error("Wrong number of parameters ({:d}) for the TCEA temperature "
                "model (5).".format(len(args.tpars)))
    if args.rstar is None:
      log.error("Undefined stellar radius (rstar).")
    if args.tstar is None:
      log.error("Undefined stellar temperature (tstar).")
    if args.smaxis is None:
      log.error("Undefined orbital semi-major axis (smaxis).")
    if (args.gplanet is None and
        (args.mplanet is None or args.rplanet is None)):
      log.error("Undefined planetary surface gravity, set either gplanet or "
                "mplanet and rplanet.")
    args.tint = pt.defaultp(args.tint, "100.0",
        "Planetary internal temperature defaulted to {:s} K.", log)
    args.runits = pt.defaultp(args.runits, "cm",
        "runits input variable defaulted to '{:s}'.", log)

  elif args.tmodel == "isothermal":
    if len(args.tpars) != 1:
      log.error("Wrong number of parameters ({:d}) for the isothermal "
                "temperature model (1).".format(len(args.tpars)))


def checkatm(args, log):
  """
  Check the input arguments to calculate the atmospheric model.
  """
  if args.atmfile is None:
    log.error("Undefined atmospheric file (atmfile).")
  if args.species is None:
    log.error("Undefined atmospheric species list (species).")
  args.punits = pt.defaultp(args.punits, "bar",
     "punits input variable defaulted to '{:s}'.", log)

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
    args.solar = pt.defaultp(args.solar, pc.ROOT+"inputs/AsplundEtal2009.txt",
      "Solar-abundances file defaulted to '{:s}'.", log)
    args.atomicfile = pt.defaultp(args.atomicfile, "./atomic.tea",
      "Atomic-composition file defaulted to '{:s}'.", log)
    args.patm = pt.defaultp(args.patm, "./preatm.tea",
      "Pre-atmospheric file defaulted to '{:s}'.", log)
    args.xsolar = pt.defaultp(args.xsolar, 1.0,
      "Solar-metallicity scaling factor defaulted to {:.2f}.", log)


def checkinputs(args, log):
  """
  Check that the input values (args) make sense.
  """
  # Stellar model:
  if args.starspec is not None:
    # Check file exists
    pass
  elif args.kurucz is not None:
    # Check file exists
    # Check gstar exists
    pass
  else:
    #log.error("Stellar spectrum model was not specified.")
    pass
