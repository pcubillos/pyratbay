# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["makeTLI", "parser", "db"]

from .lread import makeTLI, parser
from . import db

from .. import VERSION as ver

# Lineread version:
__version__ = "{:d}.{:d}.{:d}".format(ver.LR_VER, ver.LR_MIN, ver.LR_REV)


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__ ):
        del locals()[varname]
del(varname)
