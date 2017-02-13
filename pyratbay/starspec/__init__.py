# Copyright (c) 2016-2017 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

"""
Stellar-spectrum models.
"""

__all__ = ["readpyrat", "readkurucz", "kunpack"]

from .readpyrat import *
from .kurucz    import *
#from .marcs     import *
#from .phoenix   import *

# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__ ):
        del locals()[varname]
del(varname)
