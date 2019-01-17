# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

"""
Waveband-integrated emission utilities.

This set of routines read waveband filters and compute band-integrated
fluxes over the filter transmission curve.
"""

__all__ = ["Bwn", "Bwn2D"]

import os
import sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../lib')
from _blackbody import Bwn, Bwn2D


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__ ):
        del locals()[varname]
del(varname)
