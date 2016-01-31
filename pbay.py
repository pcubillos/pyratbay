#! /usr/bin/env python

import sys
import os
#import re
#import shutil
#import time
#import subprocess
#import argparse
#import ConfigParser 
#import numpy as np

# Directory of BART.py file:
maindir = os.path.dirname(os.path.realpath(__file__))

# Add path to submodules and import:
sys.path.append(maindir + '/src_Py/')
import ptools as pt

sys.path.append(maindir + '/src_C/lib/')

#TEAdir     = maindir + "/modules/TEA/"
#MC3dir     = maindir + "/modules/MCcubed/src"

#import makeP     as mp
#import InitialPT as ipt
#import makeatm   as mat
#import makecfg   as mc
#import bestFit   as bf

#sys.path.append(MC3dir)
#import mcutils   as mu

def main():
  """
  Pyrat Bay: Python Radiative Transfer in a Bayesian framework

  One function to run them all.

  Main developer:
  ---------------
  Patricio Cubillos  pcubillos@fulbrightmail.org

  Notes:
  ------
  This code is based on the Bayesian Atmospheric Radiative Transfer (BART)
  code, developed at UCF:  https://github.com/joeharr4/BART
  """

  pt.msg(1,
    "\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
    "\n  Python Radiative Transfer in a Bayesian framework (Pyrat Bay)"
    "\n  -------------------------------------------------------------"
  #"\n\nA code to infer planetary atmospheric properties based on observed  "
  #  "\nspectroscopic information."
  "\n\n  Copyright (C) 2016  Patricio Cubillos.  All rights reserved."
    "\n                      Space Research Institute "
    "\n                      (Institut fuer Weltraumforschung, IWF)."
    "\n                      Graz, Austria."
  "\n\n  Contact:  Patricio Cubillos  pyratbay@gmail.com"
    "\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")

  # Parse arguments:
  # Initialization:
  # MCMC setup:
  # MCMC run:
  # Post-MCMC:

if __name__ == "__main__":
  main()
