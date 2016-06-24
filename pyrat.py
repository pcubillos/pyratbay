#!/usr/bin/env python

# FINDME a LICENSE

import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import pyratbay as pb


def main():
  """
  PyRaT (Python Radiative Transfer) initialization driver.
  """

  # Initialize the Pyrat Object:
  pyrat = pb.pyrat.init(sys.argv, True)
  # Compute the spectrum:
  pb.pyrat.run(pyrat)


if __name__ == "__main__":
  main()

