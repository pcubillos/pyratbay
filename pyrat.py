#!/usr/bin/env python
import sys, os
import time
import numpy as np

# Main dir (where pyrat.py is located):
maindir = os.path.dirname(os.path.realpath(__file__))
# Add paths to Python a C source folders:
sys.path.append(maindir + '/src_Py/')
sys.path.append(maindir + '/src_C/lib/')

import argum      as ar
import makesample as ms
import readatm    as ra
import readlinedb as rl
import voigt      as v
import extinction as ex
import crosssec   as cs
import optdepth   as od
import spectrum   as sp

from objects import Pyrat


def init(argv, main=False):
  """
  PyRaT (Python Radiative Transfer) initialization driver.

  Parameters:
  -----------
  argv: String
     If ran from the Shell: the list of command line arguments.  If ran
     from the Python  Interpreter: the name of the configuration file.
  main: Bool
     Flag to indicate if Pyrat was called from the Shell (True) or from
     the Python Interpreter.

  Returns:
  --------
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
  ar.parse(pyrat)
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

  # Extinction gridding:
  v.voigt(pyrat)
  timestamps.append(time.time())

  # Read CIA files:
  cs.read(pyrat)
  timestamps.append(time.time())

  # Calculate extinction-coefficient table:
  ex.exttable(pyrat)
  timestamps.append(time.time())

  pyrat.timestamps = timestamps
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
  timestamps = pyrat.timestamps
  # Re-calculate atmospheric properties if required:
  if inputs is not None:
    ms.reloadatm(pyrat, *inputs)

  # Interpolate CIA absorption:
  cs.interpolate(pyrat)
  timestamps.append(time.time())

  # Calculate the optical depth:
  od.opticaldepth(pyrat)
  timestamps.append(time.time())

  # Calculate the modulation (transit) or emission (eclipse) spectrum:
  sp.spectrum(pyrat)
  timestamps.append(time.time())

  print("\nTimestamps:\n"
        " Init:     {:10.6f}\n Parse:    {:10.6f}\n Inputs:   {:10.6f}\n"
        " Wnumber:  {:10.6f}\n Atmosph:  {:10.6f}\n TLI:      {:10.6f}\n"
        " Layers:   {:10.6f}\n Voigt:    {:10.6f}\n CIA read: {:10.6f}\n"
        " Extinct:  {:10.6f}\n CIA intp: {:10.6f}\n O.Depth:  {:10.6f}\n"
        " Spectrum: {:10.6f}".format(*np.ediff1d(timestamps)))
  return pyrat


if __name__ == "__main__":
  pyrat = init(sys.argv, True)
  run(pyrat)

