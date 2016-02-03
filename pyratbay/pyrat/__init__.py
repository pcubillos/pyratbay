__all__ = ["init", "run"]

from .shipmaster import init, run

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
