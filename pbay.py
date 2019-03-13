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

    Examples
    --------
    # Run Pyrat Bay:
    python pbay.py -c config.cfg

    # Re-format partition function files
    python pbay.py -pf exomol 14N-1H3__BYTe.pf 15N-1H3__BYTe-15.pf
    python pbay.py -pf kurucz h2opartfn.dat
    """
    # Parse configuration file:
    parser = argparse.ArgumentParser(description=__doc__, add_help=True,
                  formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--version", action="version",
                       help="Show Pyrat Bay's version.",
                       version='Pyrat Bay version {:s}.'.format(pb.__version__))
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-c", dest='cfile', default=None,
                       help="Configuration file.")
    group.add_argument("-pf", dest='pf', default=None, nargs='+',
                       help="Format a partition-function file.")
    group.add_argument("-cs", dest='cs', default=None, nargs='+',
                       help="Format a cross-section file.")
    # Parse command-line args:
    args, unknown = parser.parse_known_args()

    # Partition-function reformatting:
    if args.pf is not None:
        if args.pf[0] == 'exomol':
            pb.tools.pf_exomol(args.pf[1:])
        elif args.pf[0] == 'kurucz':
            pb.tools.pf_kurucz(args.pf[1])
    # Cross-section reformatting:
    elif args.cs is not None:
        pass
    # Pyrat-Bay run:
    elif args.cfile is not None:
        pb.pbay.run(args.cfile)


if __name__ == "__main__":
    main()
