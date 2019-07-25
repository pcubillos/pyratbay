#! /usr/bin/env python

# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import argparse

import pyratbay as pb


def main():
    """
    Pyrat Bay: Python Radiative Transfer in a Bayesian framework

    Pyrat Bay is a modification of the GNU GPL v2-licensed transit code
    (developed by Patricio Rojo); therefore, Pyrat Bay (v1.0+) is also
    released under the GNU GLP v2 license.

    Examples
    --------
    # Run Pyrat Bay:
    python pbay.py -c config.cfg

    # Re-format partition-function files
    python pbay.py -pf exomol 14N-1H3__BYTe.pf 15N-1H3__BYTe-15.pf
    python pbay.py -pf kurucz h2opartfn.dat
    python pbay.py -pf tips H2O

    # Re-format cross-section files:
    python pbay.py -cs hitran H2-H2_2011.cia 2 10
    python pbay.py -cs borysow ciah2he_dh_quantmech H2 He
    """
    # Parse configuration file:
    parser = argparse.ArgumentParser(description=__doc__, add_help=True,
                  formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--version", action="version",
                       help="Show Pyrat Bay's version.",
                       version='Pyrat Bay version {:s}.'.format(pb.__version__))
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-c", dest='cfile', default=None,
                       help="Run Pyrat Bay for given configuration file.")
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
        elif args.pf[0] == 'tips':
            print('TIPS partition function is TBD.')
        else:
            print('Invalid partition-function type.')

    # Cross-section reformatting:
    elif args.cs is not None:
        if args.cs[0] == 'hitran':
            pb.tools.cia_hitran(args.cs[1], int(args.cs[2]), int(args.cs[3]))
        elif args.cs[0] == 'borysow':
            pb.tools.cia_borysow(*args.cs[1:])
        else:
            print('Invalid cross-section type.')

    # Pyrat-Bay run:
    elif args.cfile is not None:
        pb.run(args.cfile)


if __name__ == "__main__":
    main()
