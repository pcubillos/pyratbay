#!/usr/bin/env python

# FINDME a LICENSE

import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import pyratbay as pbay


def main():
  """
  PyRaT (Python Radiative Transfer) initialization shipmaster.
  """

  # Initialize the Pyrat Object:
  pyrat = pbay.pyrat.init(sys.argv, True)
  # Compute the spectrum:
  pbay.pyrat.run(pyrat)


if __name__ == "__main__":
  main()

