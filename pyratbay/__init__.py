# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'constants',
    'io',
    'tools',
    'opacity',
    'plots',
    'spectrum',
    'atmosphere',
    'Pyrat',
    'run',
    ]

from . import constants
from . import io
from . import tools
from . import opacity
from . import plots
from . import spectrum
from . import atmosphere
from .pyrat import Pyrat
from .driver import run

from .pyrat import read_atm as _ra
from .pyrat import crosssec as _cs
from .pyrat import clouds   as _cl
from .pyrat import rayleigh as _ray
from .pyrat import alkali   as _al
from .pyrat import optical_depth as _od
from .pyrat import spectrum as _sp
__all__ += ['_ra', '_cs', '_cl', '_ray', '_al', '_od', '_sp']

from .VERSION import __version__


# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__):
        del locals()[varname]
del(varname)
