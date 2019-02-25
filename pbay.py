#! /usr/bin/env python

# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import pyratbay as pb


def main():
    """
    Pyrat Bay: Python Radiative Transfer in a Bayesian framework

    This code was initially developed in parallel with the Bayesian
    Atmospheric Radiative Transfer (BART) code, developed at UCF:
    https://github.com/exosports/BART
    """
    pb.pbay.run(sys.argv, True)


if __name__ == "__main__":
    main()
