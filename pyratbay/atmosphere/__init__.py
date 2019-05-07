# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

"""
Atmospheric modeling functions.
"""

from .qscaling   import *
from .atmosphere import *
from . import tmodels

__all__ = qscaling.__all__ + atmosphere.__all__ + ['tmodels']

# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__ ):
        del locals()[varname]
del(varname)
