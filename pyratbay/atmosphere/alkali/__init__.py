# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

"""Alkali resonant-line models."""

from .alkali import *

__all__ = alkali.__all__ + ['get_model']

def get_model(name):
    """Get an alkali model by its name."""
    if name == 'sodium_vdw':
        return SodiumVdW()
    if name == 'potassium_vdw':
        return PotassiumVdW()


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__ ):
        del locals()[varname]
del(varname)
