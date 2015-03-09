#!/usr/bin/env python
import sys, os

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

from objects import pyrat

def main(argv):
  """
  Doc me!

  Modification History:
  ---------------------
  2014-04-26  patricio  Initial implementation.
  2014-08-17  patricio  Added ex.voigt step.
  2014-08-31  patricio  Added cia steps.
  """

  # Initialize a pyrat object:
  cavendish = pyrat()

  # Parse command line arguments into pyrat:
  ar.parse(cavendish)

  # Check that user input arguments make sense:
  ar.checkinputs(cavendish)

  # Initialize wavenumber sampling:
  ms.makewavenumber(cavendish)

  # Read the atmospheric file:
  ra.readatm(cavendish)

  # Read line database:
  rl.readlinedb(cavendish)

  # Make radius sampling:
  ms.makeradius(cavendish)

  # Extinction gridding:
  ex.voigt(cavendish)

  # Read CIA files:
  cia.read(cavendish)

  # Calculate extinction-coefficient table:
  ex.exttable(cavendish)

  # :: Pyrat-Bay loop boundary ::

  # Interpolate CIA absorption:
  cia.interpolate(cavendish)

  # Calculate the optical depth:
  od.opticaldepth(cavendish)

  print("Next one!")
  return 1

if __name__ == "__main__":
  main(sys.argv)
