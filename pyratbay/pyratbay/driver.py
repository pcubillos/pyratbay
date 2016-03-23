#!/usr/bin/env python

# FINDME a LICENSE

import os
import sys
import time
import numpy as np

from .. import tools as pt
from .. import VERSION   as ver

from . import argum      as ar

#from . import makesample as ms
#from . import readatm    as ra
#from . import readlinedb as rl
#from . import voigt      as v
#from . import extinction as ex
#from . import crosssec   as cs
#from . import haze       as hz
#from . import alkali     as al
#from . import optdepth   as od
#from . import spectrum   as sp
#
#from .objects import Pyrat


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

  Returns
  -------
  cavendish: Pyrat instance
     The Pyrat object.
  """

  pt.msg(1, "{:s}\n"
            "  Python Radiative Transfer in a Bayesian framework (Pyrat Bay).\n"
            "  Version {:d}.{:d}.{:d}.\n"
            "  Copyright (c) 2016 Patricio Cubillos and collaborators.\n"
            "  Pyrat Bay is open-source software under the RR license.\n"
            "{:s}\n\n".format(pt.sep, ver.PBAY_VER, ver.PBAY_MIN,
                                      ver.PBAY_REV, pt.sep), None)#pb.log)

  # Setup the command-line-arguments input:
  if main is False:
    sys.argv = ['pbay.py', '-c', cfile]

  # Setup time tracker:
  timestamps = []
  timestamps.append(time.time())

  # Parse command line arguments:
  args = ar.parse()
  timestamps.append(time.time())

  # Compute an atmospheric file:
  pass

  # Compute an opacity grid:
  pass

  # Full Pyrat Bay run:
  pass

  # Post processing:
  pass

