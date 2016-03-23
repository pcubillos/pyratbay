import sys, os
import ConfigParser, argparse
import numpy as np

from .. import tools as pt

def parse():
  """
  Read the command line arguments.
  """

  # Parser to process a configuration file:
  cparser = argparse.ArgumentParser(description=__doc__, add_help=False,
                           formatter_class=argparse.RawDescriptionHelpFormatter)
  # Add config file option:
  cparser.add_argument("-c", "--cfile",
                       help="Configuration filename (string).", metavar="FILE")
  # remaining_argv contains all other command-line-arguments:
  args, remaining_argv = cparser.parse_known_args()

  if args.cfile is None:
    pt.error("No configuration file specified.")
  elif not os.path.isfile(args.cfile):
    pt.error("Configuration file '{:s}' does not exist.".format(args.cfile))

  config = ConfigParser.SafeConfigParser()
  config.read([args.cfile])
  if "pbay" not in config.sections():
    pt.error("Invalid configuration file: '{:s}'.".format(args.cfile))

  defaults = dict(config.items("pbay"))

  # Inherit options from cparser:
  parser = argparse.ArgumentParser(parents=[cparser])

  parser.set_defaults(**defaults)
  args = parser.parse_args(remaining_argv)

  return args
