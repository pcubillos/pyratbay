#!/usr/bin/env python

# FINDME a LICENSE

import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import pyratbay as pbay


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
  parser = pbay.lineread.parser()

  # Call the Lineread driver:
  pbay.lineread.makeTLI(parser.dblist,  parser.pflist, parser.dbtype,
                        parser.outfile, parser.iwl, parser.fwl,
                        parser.verb)


if __name__ == "__main__":
  main()
