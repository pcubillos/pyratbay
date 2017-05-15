# Copyright (c) 2016-2017 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import time
import numpy as np

from .. import tools as pt

from . import argum      as ar
from . import makesample as ms
from . import readatm    as ra
from . import readlinedb as rl
from . import voigt      as v
from . import extinction as ex
from . import crosssec   as cs
from . import haze       as hz
from . import rayleigh   as ray
from . import alkali     as al
from . import optdepth   as od
from . import spectrum   as sp

from .objects import Pyrat



def init(argv, main=False, log=None):
  """
  PyRaT (Python Radiative Transfer) initialization driver (or should
  I say shipmaster?).

  Parameters
  ----------
  argv: List or string
     If called from the shell, the list of command line arguments; if
     called from the Python interpreter, the configuration-file name.
  main: Bool
     Flag to indicate if Pyrat was called from the shell (True) or from
     the Python interpreter.

  Returns
  -------
  cavendish: Pyrat instance
     The Pyrat object.
  """
  # Setup the command-line-arguments input:
  if main is False:
    sys.argv = ['pyrat.py', '-c', argv]

  # Setup time tracker:
  timestamps = []
  timestamps.append(time.time())

  # Initialize a pyrat object:
  pyrat = Pyrat()
  timestamps.append(time.time())

  # Parse command line arguments into pyrat:
  ar.parse(pyrat, log)
  timestamps.append(time.time())

  # Check that user input arguments make sense:
  ar.checkinputs(pyrat)
  timestamps.append(time.time())

  # Initialize wavenumber sampling:
  ms.makewavenumber(pyrat)
  timestamps.append(time.time())

  # Read the atmospheric file:
  ra.readatm(pyrat)
  timestamps.append(time.time())

  # Read line database:
  rl.readlinedb(pyrat)
  timestamps.append(time.time())

  # Make radius sampling:
  ms.makeradius(pyrat)
  timestamps.append(time.time())

  # Set up observational/retrieval parameters:
  ar.setup(pyrat)

  # Extinction Voigt grid:
  v.voigt(pyrat)
  # Alkali Voigt grid:
  al.init(pyrat)
  timestamps.append(time.time())

  # Calculate extinction-coefficient table:
  ex.exttable(pyrat)
  timestamps.append(time.time())
  if pyrat.runmode == "opacity":  # Early way out
    return pyrat

  # Read CIA files:
  cs.read(pyrat)
  timestamps.append(time.time())

  pyrat.timestamps = list(np.ediff1d(timestamps))
  return pyrat


def run(pyrat, inputs=None):
  """
  PyRaT driver to calculate the spectrum.

  Parameters:
  -----------
  pyrat: A Pyrat instance
  inputs: list
     A list containing a 1D float array of temperatures, and a 2D array of
     the species abundances.
  """
  timestamps = []
  timestamps.append(time.time())

  # Re-calculate atmospheric properties if required:
  if inputs is not None:
    status = ra.reloadatm(pyrat, *inputs)
    if status == 0:
      return

  # Interpolate CIA absorption:
  cs.interpolate(pyrat)
  timestamps.append(time.time())

  # Calculate the haze absorption:
  hz.absorption(pyrat)
  ray.absorption(pyrat)
  timestamps.append(time.time())

  # Calculate the alkali absorption:
  al.absorption(pyrat)
  timestamps.append(time.time())

  # Calculate the optical depth:
  od.opticaldepth(pyrat)
  timestamps.append(time.time())

  # Calculate the spectrum:
  sp.spectrum(pyrat)
  timestamps.append(time.time())

  pyrat.timestamps += list(np.ediff1d(timestamps))
  dtime = pyrat.timestamps[0:9] + pyrat.timestamps[-6:]
  pt.msg(pyrat.verb-4, "\nTimestamps:\n"
        " Init:     {:10.6f}\n Parse:    {:10.6f}\n Inputs:   {:10.6f}\n"
        " Wnumber:  {:10.6f}\n Atmosph:  {:10.6f}\n TLI:      {:10.6f}\n"
        " Layers:   {:10.6f}\n Voigt:    {:10.6f}\n Extinct:  {:10.6f}\n"
        " CIA read: {:10.6f}\n CIA intp: {:10.6f}\n Haze:     {:10.6f}\n"
        " Alkali:   {:10.6f}\n O.Depth:  {:10.6f}\n Spectrum: {:10.6f}".
        format(*dtime), pyrat.log)

  if len(pyrat.wlog) > 0:
    # Write all warnings to file:
    wpath, wfile = os.path.split(pyrat.logfile)
    wfile = "{:s}/warnings_{:s}".format(wpath, wfile)
    with open(wfile, "w") as warns:
      warns.write("Warnings log:\n\n{:s}\n".format(pt.sep))
      warns.write("\n\n{:s}\n".format(pt.sep).join(pyrat.wlog))
    # Report it:
    pyrat.warning("There were {:d} warnings raised.  See '{:s}'.".
                   format(len(pyrat.wlog), wfile))
    pyrat.log.flush()
  return pyrat

