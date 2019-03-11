#! /usr/bin/env python

# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import argparse

import pyratbay as pb


def main():
    """
    Pyrat Bay: Python Radiative Transfer in a Bayesian framework

    This code was initially developed in parallel with the Bayesian
    Atmospheric Radiative Transfer (BART) code, developed at UCF:
    https://github.com/exosports/BART
    """
    # Parse configuration file:
    parser = argparse.ArgumentParser(description=__doc__, add_help=True,
                  formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--version", action="version",
                       help="Show Pyrat Bay's version.",
                       version='Pyrat Bay version {:s}.'.format(pb.__version__))
    parser.add_argument("-c", dest='cfile', default=None,
                       help="Configuration file.")

    # Parse command-line args:
    args, unknown = parser.parse_known_args()

    # Make calls:
    if args.cfile is not None:
        pb.pbay.run(args.cfile)


if __name__ == "__main__":
    main()
