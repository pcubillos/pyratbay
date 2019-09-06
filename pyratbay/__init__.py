# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = [
    'constants',
    'io',
    'tools',
    'blackbody',
    'broadening',
    'lineread',
    'plots',
    'starspec',
    'atmosphere',
    'Pyrat',
    'run',
    ]

from . import constants
from . import io
from . import tools
from . import blackbody
from . import broadening
from . import lineread
from . import plots
from . import starspec
from . import atmosphere
from .pyrat import Pyrat
from .driver import run

# Pyrat Bay version:
from .VERSION import __version__


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__):
        del locals()[varname]
del(varname)
