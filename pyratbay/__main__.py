#! /usr/bin/env python

# Copyright (c) 2021-2024 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import argparse
import matplotlib

import pyratbay as pb


def main():
    """
    Pyrat Bay: Python Radiative Transfer in a Bayesian framework

    Pyrat Bay is open-source software released under the GNU GPL v2
    license.  The code is compatible with Python >= 3.6
    Find the source code at https://github.com/pcubillos/pyratbay
    Find the documentation at https://pyratbay.readthedocs.io/en/latest
    Part of Pyrat Bay is based on the Transit code written by Patricio Rojo.

    Examples
    --------
    # Run Pyrat Bay:
    pbay -c config.cfg

    # Re-format partition-function files
    pbay -pf exomol 14N-1H3__BYTe.pf 15N-1H3__BYTe-15.pf
    pbay -pf kurucz h2opartfn.dat
    pbay -pf tips H2O

    # Re-format cross-section files:
    pbay -cs hitran H2-H2_2011.cia 2 10
    pbay -cs borysow ciah2he_dh_quantmech H2 He
    """
    matplotlib.use('Agg')
    # Parse configuration file:
    parser = argparse.ArgumentParser(
        description=__doc__,
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        help="Show Pyrat Bay's version.",
        version=f'Pyrat Bay version {pb.__version__}',
    )
    parser.add_argument(
        '--root',
        action='version',
        help="Show Pyrat Bay's ROOT directory.",
        version=pb.constants.ROOT,
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '-c', dest='cfile', default=None,
        help='Run Pyrat Bay for given configuration file.',
    )
    group.add_argument(
        '--post', dest='post', default=None,
        help='Post-processing posterior data after a retrieval run.',
    )
    group.add_argument(
        '-pf', dest='pf', default=None, nargs='+',
        help='Format a partition-function file.',
    )
    group.add_argument(
        '-cs', dest='cs', default=None, nargs='+',
        help='Format a cross-section file.',
    )

    parser.add_argument(
        '-suf', dest='suffix', default=None,
        help='Suffix for post-processed file.',
    )

    # Parse command-line args:
    args, unknown = parser.parse_known_args()

    # Partition-function reformatting:
    if args.pf is not None:
        outfile = 'default'
        if args.pf[0] == 'exomol':
            pb.opacity.partitions.exomol(args.pf[1:], outfile=outfile)
        elif args.pf[0] == 'kurucz':
            pb.opacity.partitions.kurucz(args.pf[1], outfile=outfile)
        elif args.pf[0] == 'tips':
            dbtype = args.pf[2] if len(args.pf) > 2 else 'as_tips'
            pb.opacity.partitions.tips(
                args.pf[1], outfile=outfile, db_type=dbtype)
        else:
            print('Invalid partition-function type.')

    # Cross-section reformatting:
    elif args.cs is not None:
        if args.cs[0] == 'hitran':
            if len(args.cs) < 3:
                t_step = 1
            else:
                t_step = int(args.cs[2])
            if len(args.cs) < 4:
                wn_step = 1
            else:
                wn_step = int(args.cs[3])
            pb.tools.cia_hitran(args.cs[1], t_step, wn_step)
        elif args.cs[0] == 'borysow':
            pb.tools.cia_borysow(*args.cs[1:])
        else:
            print('Invalid cross-section type.')

    # Pyrat-Bay run:
    elif args.cfile is not None:
        pb.run(args.cfile)

    # Post processing:
    elif args.post is not None:
        suffix = '' if args.suffix is None else args.suffix
        pb.tools.posterior_post_processing(cfg_file=args.post, suffix=suffix)


if __name__ == '__main__':
    matplotlib.pyplot.ioff()
    main()
