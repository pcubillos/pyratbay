# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["init"]


def init(pyrat, args, log):
  """
  Initialize variables that will be used in the atmospheric retrieval.
  """
  # Check stellar spectrum model:
  if pyrat.od.path == "eclipse" and pyrat.phy.starflux is None:
    log.error("Unspecified stellar flux model.")

  # Check filter files and data:
  if pyrat.obs.filter is None:
    log.error("Undefined waveband transmission filters (filter).")
  if pyrat.obs.data is None:
    log.error("Undefined data.")
  if pyrat.obs.uncert is None:
    log.error("Undefined data uncertainties.")

  # Check rprs:
  if pyrat.od.path == "eclipse" and pyrat.phy.rprs is None:
    log.error("Undefined Rp/Rs.")

  # Boundaries for temperature profile:
  pyrat.ret.tlow  = args.tlow
  pyrat.ret.thigh = args.thigh
