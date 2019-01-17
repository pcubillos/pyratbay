# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

__all__ = ["bbflux"]

import os
import sys
import numpy as np

from .. import blackbody as bb


def bbflux(wn, Teff):
  """
  Compute the emission flux of a blackbody at temperature Teff
  in wavenumber space.

  Parameters
  ----------
  wn: 1D float ndarray
     Wavenumber array where to evaluate the flux (cm-1).
  Teff: Float
     The effective temperature (Kelvin).

  Return
  ------
  flux: 1D float ndarray
     blackbody flux (erg s-1 cm-2 cm) at wn.

  Note
  ----
  This function adopts the spherical-cow approximation.
  """
  return np.pi * bb.Bwn(wn, Teff)
