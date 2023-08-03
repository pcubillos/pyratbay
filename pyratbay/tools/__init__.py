# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

from .data import *
from .mpi_tools import *
from .parser import *
from .retrieval_tools import *
from .tools import *

__all__ = (
    data.__all__
    + mpi_tools.__all__
    + parser.__all__
    + retrieval_tools.__all__
    + tools.__all__
)

# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__):
        del locals()[varname]
del(varname)
