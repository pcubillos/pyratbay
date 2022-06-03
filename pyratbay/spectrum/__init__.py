# Copyright (c) 2021-2022 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

from .blackbody import *
from .kurucz import *
#from .marcs import *
#from .phoenix import *
from .spec_tools import *
from .contribution_funcs import *
from .convection import *
from .radiative_transfer import *

__all__ = (
    []  # empty placeholder
    + blackbody.__all__
    + kurucz.__all__
    + spec_tools.__all__
    + contribution_funcs.__all__
    + convection.__all__
    + radiative_transfer.__all__
)


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for var in dir():
    if not ((var.startswith('__') and var.endswith('__')) or var in __all__):
        del locals()[var]
del(var)
