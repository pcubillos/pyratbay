# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'blackbody_wn',
    'blackbody_wn_2D',
    'bbflux',
    ]

from numbers import Integral
from collections import Iterable

import numpy as np

from ..lib._blackbody import blackbody_wn, blackbody_wn_2D


def bbflux(wn, teff):
    r"""
    Compute the emission flux of a blackbody at temperature Teff
    in wavenumber space.

    Parameters
    ----------
    wn: 1D float iterable
       Wavenumber array where to evaluate the flux (cm-1).
    teff: Float
       The effective temperature (Kelvin).

    Return
    ------
    flux: 1D float ndarray
       blackbody flux (erg s-1 cm-2 cm) at wn.

    Examples
    --------
    >>> import pyratbay.spectrum as ps
    >>> import pyratbay.constants as pc
    >>> import numpy as np
    >>> tsun = 5772.0
    >>> wn = np.logspace(-1, 5, 30000)
    >>> flux = ps.bbflux(wn, tsun)
    >>> # Solar constant:
    >>> s = np.trapz(flux, wn) * (pc.rsun/pc.au)**2
    >>> print("Solar constant (Teff={:.0f}K): S = {:.1f} W m-2\n"
    >>>       "Wien's displacement law: wn(flux_max) = {:.1f} cm-1\n"
    >>>       "             5.879E10 Hz/K * Teff / c = {:.1f} cm-1".
    >>>       format(tsun, s*1e-3, wn[np.argmax(flux)], 5.879e10*tsun/pc.c))
    Solar constant (Teff=5772K): S = 1361.2 W m-2
    Wien's displacement law: wn(flux_max) = 11318.0 cm-1
                 5.879E10 Hz/K * Teff / c = 11319.0 cm-1
    """
    if not isinstance(wn, Iterable):
        raise ValueError('Input wn must be an iterable.')

    if isinstance(wn, (list, tuple)):
        wn = np.array(wn, np.double)
    if isinstance(wn[0], Integral):
        return np.pi * blackbody_wn(np.array(wn, np.double), float(teff))
    return np.pi * blackbody_wn(wn, float(teff))
