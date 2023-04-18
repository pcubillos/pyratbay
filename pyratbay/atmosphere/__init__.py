# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

"""
Atmospheric modeling functions.
"""

from .atmosphere import *
from .vmr_scaling import *
from . import tmodels
from . import clouds

__all__ = (
    atmosphere.__all__
  + vmr_scaling.__all__
  + ['tmodels']
  + ['clouds']
)


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__):
        del locals()[varname]
del(varname)
