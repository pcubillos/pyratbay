#!/usr/bin/env python

# Copyright (c) 2016 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).


import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import pyratbay as pb


def main():
  """
  Lineread shell-executable driver.

  Usage
  -----
  Execute from the shell:
  ./lineread.py [--option <args>]

  To display the list of command-line arguments execute:
  ./lineread.py --help
  """

  # Parse the command-line arguments:
  parser = pb.lineread.parser()

  # Call the Lineread driver:
  pb.lineread.makeTLI(parser.dblist,  parser.pflist, parser.dbtype,
                      parser.outfile, parser.iwl, parser.fwl,
                      parser.verb)


if __name__ == "__main__":
  main()
