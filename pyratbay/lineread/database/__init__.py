# Copyright (c) 2016-2020 Patricio Cubillos.
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE).

__all__ = ["hitran", "pands", "tioschwenke", "voplez", "vald", "exomol",
           "repack"]

from .hitran      import hitran
from .pands       import pands
from .tioschwenke import tioschwenke
from .voplez      import voplez
from .vald        import vald
from .exomol      import exomol
from .repack      import repack

# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
  if not ((varname.startswith('__') and varname.endswith('__')) or
          varname in __all__ ):
    del locals()[varname]
del(varname)
