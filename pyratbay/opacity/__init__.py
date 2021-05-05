# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'make_tli',
    'linelist',
    'partitions',
    'broadening',
    ]

from .lread import *
from . import linelist
from . import partitions
from . import broadening

from .. import VERSION as ver


# Lineread version:
__version__ = f'{ver.LR_VER}.{ver.LR_MIN}.{ver.LR_REV}'


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__):
        del locals()[varname]
del(varname)
