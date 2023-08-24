# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

__all__ = [
    'constants',
    'io',
    'tools',
    'opacity',
    'plots',
    'spectrum',
    'atmosphere',
    'Pyrat',
    'Atmosphere',
    'run',
]

import warnings

from . import constants
from . import io
from . import tools
from . import opacity
from . import plots
from . import spectrum
from . import atmosphere
from .pyrat import Pyrat
from .pyrat.atmosphere import Atmosphere
from .driver import run

from .pyrat import optical_depth as _od
from .pyrat import spectrum as _sp
__all__ += ['_od', '_sp']

from .version import __version__


# Display deprecation warnings when running in the interpreter:
warnings.filterwarnings(
    'default',
    category=DeprecationWarning,
    module=fr'^{__name__}\.',
)

# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__):
        del locals()[varname]
del(varname)
