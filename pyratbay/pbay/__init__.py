# Copyright (c) 2016-2018 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["run", "mc", "ar", "pf"]

from .driver import run

from .  import makecfg as mc
from .  import pyratfit as pf
from .. import VERSION as ver

# Temporary:
from . import argum as ar

# Pyrat Bay version:
__version__ = "{:d}.{:d}.{:d}".format(ver.PBAY_VER, ver.PBAY_MIN,
                                      ver.PBAY_REV)


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__ ):
        del locals()[varname]
del(varname)
