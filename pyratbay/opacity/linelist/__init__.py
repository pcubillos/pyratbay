# Copyright (c) 2021-2025 Cubillos & Blecic
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

from .hitran import Hitran
from .exomol import Exomol
from .repack import Repack
from .pands import Pands
from .tioschwenke import Tioschwenke
from .voplez import Voplez
from .vald import Vald

__all__ = [
    'Hitran',
    'Exomol',
    'Repack',
    'Pands',
    'Tioschwenke',
    'Voplez',
    'Vald',
    ]

# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __all__):
        del locals()[varname]
del(varname)
