# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = [
    "constants",
    "io",
    "blackbody",
    "broadening",
    "lineread",
    "plots",
    "starspec",
    "tools",
    "atmosphere",
    "Pyrat",
    "init",
    "run",
]

from . import constants
from . import io
from . import blackbody
from . import broadening
from . import lineread
from . import plots
from . import starspec
from . import tools
from . import atmosphere
from .pyrat import Pyrat
from .driver import init, run

from . import VERSION as ver


# Pyrat-Bay version:
__version__ = "{:d}.{:d}.{:d}".format(ver.PBAY_VER, ver.PBAY_MIN, ver.PBAY_REV)


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__ ):
        del locals()[varname]
del(varname)
