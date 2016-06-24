#!/usr/bin/env python

# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

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

