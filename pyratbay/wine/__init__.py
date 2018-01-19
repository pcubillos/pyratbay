# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

"""
Waveband-integrated emission utilities.

This set of routines read waveband filters and compute band-integrated
fluxes over the filter transmission curve.
"""

from .wine import *

__all__ = wine.__all__

# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__ ):
        del locals()[varname]
del(varname)
