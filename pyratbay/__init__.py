# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["pbay", "lineread", "pyrat", "constants", "tools", "wine",
           "blackbody", "broadening", "starspec", "atmosphere", "plots"]

from . import pbay
from . import pyrat
from . import lineread
from . import constants
from . import tools
from . import wine
from . import blackbody
from . import broadening
from . import starspec
from . import atmosphere
from . import plots

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
