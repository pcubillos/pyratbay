#!/usr/bin/env python
import sys, os
import time
import numpy as np

# Main dir (where pyrat.py is located):
maindir = os.path.dirname(os.path.realpath(__file__))
# Add paths to Python a C folders:
sys.path.append(maindir + '/src_Py/')
sys.path.append(maindir + '/src_C/lib/')

import argum      as ar
import makesample as ms
import readatm    as ra
import readlinedb as rl
import extinction as ex
import cia        as cia
import optdepth   as od
import spectrum   as sp

from objects import pyrat

def main(argv):
  """
  PyRaT driver.
  """

  timestamps = []
  timestamps.append(time.time())
  # Initialize a pyrat object:
  cavendish = pyrat()
  timestamps.append(time.time())

  # Parse command line arguments into pyrat:
  ar.parse(cavendish)
  timestamps.append(time.time())

  # Check that user input arguments make sense:
  ar.checkinputs(cavendish)
  timestamps.append(time.time())

  # Initialize wavenumber sampling:
  ms.makewavenumber(cavendish)
  timestamps.append(time.time())

  # Read the atmospheric file:
  ra.readatm(cavendish)
  timestamps.append(time.time())

  # Read line database:
  rl.readlinedb(cavendish)
  timestamps.append(time.time())

  # Make radius sampling:
  ms.makeradius(cavendish)
  timestamps.append(time.time())

  # Extinction gridding:
  ex.voigt(cavendish)
  timestamps.append(time.time())

  # Read CIA files:
  cia.read(cavendish)
  timestamps.append(time.time())

  # Calculate extinction-coefficient table:
  ex.exttable(cavendish)
  timestamps.append(time.time())

  # ::: Pyrat-Bay loop boundary :::
  # Interpolate CIA absorption:
  cia.interpolate(cavendish)
  timestamps.append(time.time())

  # Calculate the optical depth:
  od.opticaldepth(cavendish)
  timestamps.append(time.time())

  # Calculate the modulation (transit) or emission (hem) spectrum:
  sp.spectrum(cavendish)
  timestamps.append(time.time())

  print("Next one!")
  print("Timestamps:\n"
        " Init:    {:10.6f}\n Parse:   {:10.6f}\n Inputs:  {:10.6f}\n"
        " Wnumber: {:10.6f}\n Atm:     {:10.6f}\n TLI:     {:10.6f}\n"
        " Layers:  {:10.6f}\n Voigt:   {:10.6f}\n CIAr:    {:10.6f}\n"
        " Ext:     {:10.6f}\n CIAi:    {:10.6f}\n Depth:   {:10.6f}\n"
        " Spect:   {:10.6f}".format(*np.ediff1d(timestamps)))
  return 1

if __name__ == "__main__":
  main(sys.argv)
