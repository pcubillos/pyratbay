# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = [
    "Pyrat",
    "constants",
    "io",
    "blackbody",
    "broadening",
    "lineread",
    "wine",
    "plots",
    "starspec",
    "tools",
    "atmosphere",
    "pbay",
]

# No pre-requirements:
from . import constants
from . import io
from . import blackbody
# Require constants:
from . import broadening
from . import lineread
from . import wine
# Require constants, wine:
from . import plots
# Require constants, blackbody:
from . import starspec
# Require constants, io, starspec, (pbay):
from . import tools
# Require constants, tools:
from . import atmosphere
# Require many things:
from .pyrat import Pyrat
# Require even more things:
from . import pbay

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
