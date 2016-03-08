__all__ = ["pyratbay", "lineread", "pyrat", "constants", "tools"]

from . import pyratbay
from . import pyrat
from . import lineread
from . import constants
from . import tools

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
