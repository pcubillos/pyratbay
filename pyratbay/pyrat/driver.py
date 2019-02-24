# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import sys
import time
import numpy as np

from . import argum      as ar
from . import makesample as ms
from . import readatm    as ra
from . import readlinedb as rl
from . import voigt      as v
from . import extinction as ex
from . import crosssec   as cs
from . import alkali     as al

from .objects import Pyrat


def init(argv, main=False, log=None):
  """
  Pyrat (Python Radiative Transfer) initialization driver.

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
  pyrat: Pyrat instance
     A Pyrat object.
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

  # Make atmospheric profiles (pressure, radius, temperature, abundances):
  ms.make_atmprofiles(pyrat)
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

  # Read CIA files:
  cs.read(pyrat)
  timestamps.append(time.time())

  pyrat.timestamps = list(np.ediff1d(timestamps))
  return pyrat
