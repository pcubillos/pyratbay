# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

"""Cloud models."""

from .gray import *

__all__ = gray.__all__ + ['get_model']

def get_model(name):
    """Get a cloud model by its name."""
    for model in __all__:
        if model.lower() == name:
            return eval('{:s}()'.format(model))


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__ ):
        del locals()[varname]
del(varname)
