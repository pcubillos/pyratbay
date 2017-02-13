# Copyright (c) 2016-2017 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["init", "run"]

from .driver import init, run

from .. import VERSION as ver

# PyRaT version:
__version__ = "{:d}.{:d}.{:d}".format(ver.PYRAT_VER, ver.PYRAT_MIN,
                                      ver.PYRAT_REV)


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__ ):
        del locals()[varname]
del(varname)
