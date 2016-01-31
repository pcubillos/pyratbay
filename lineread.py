#!/usr/bin/env python

# FINDME a LICENSE

import sys, os
import ConfigParser, argparse
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
#import pyratbay.lineread as lr
import pyratbay as pbay

import numpy as np


def main():
  """
  Usage:
  ------
  Execute from the shell:
  ./lineread.py [--option <args>]

  To display the list of command-line arguments execute:
  ./lineread.py --help
  """

  # Parse the command-line arguments:
  parser = pbay.lineread.parser()

  print("main holi")
  #print(dir(parser))
  #print(parser.fwav)
  #print(parser.iwav)
  print(parser)

  # Call the Lineread driver:
  pbay.lineread.makeTLI(parser.dblist,  parser.pflist, parser.dbtype,
                        parser.outfile, parser.iwl, parser.fwl,
                        parser.verb)


if __name__ == "__main__":
  main()
