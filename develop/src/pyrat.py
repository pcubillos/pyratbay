#!/usr/bin/env python
import sys
import argum      as ar
import makesample as ms
import readatm    as ra
import readlinedb as rl
import extinction as ex

from objects import pyrat

def main(argv):
  """
  Doc me!

  Modification History:
  ---------------------
  2014-04-26  patricio  Initial implementation.
  2014-08-17  patricio  Added ex.voigt step.
  """

  # Initialyze a pyrat object:
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

  print("Next one!")

if __name__ == "__main__":
  main(sys.argv)
