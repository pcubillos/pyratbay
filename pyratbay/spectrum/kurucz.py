# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'read_kurucz',
    ]

import numpy as np

from .. import constants as pc


def read_kurucz(filename, temp=None, logg=None):
    """
    Extract stellar flux models from a Kurucz file.
    Kurucz model files can be found at http://kurucz.harvard.edu/grids.html

    Parameters
    ----------
    filename: String
        Name of a Kurucz model file.
    temp: Float
        Requested surface temperature for the Kurucz model.
        If temp and logg are not None, return the model with the closest
        surface temperature and gravity.
    logg: Float
        Requested log10 of the surface gravity for the Kurucz model
        (where g is in cgs units).

    Returns
    -------
    flux: 1D or 2D float ndarray
        If temp and logg are not None, a 1D array with the kurucz surface
        flux per unit wavenumber (erg s-1 cm-2 cm) of the closest model to
        the input temperature and gravity.
        Else, a 2D array with all kurucz models in file, of shape
        [nmodels, nwave].
    wavenumber: 1D ndarray
        Wavenumber sampling of the flux models (in cm-1 units).
    ktemp: Scalar or 1D float ndarray
        Surface temperature of the output models (in Kelvin degrees).
    klogg: Scalar or 1D float ndarray
        log10 of the stellar surface gravity of the output models (in cm s-2).
    continuum: 2D ndarray
        The models' fluxes with no line absorption.  Same units and
        shape of flux. Returned only if temp and logg are None.

    Examples
    --------
    >>> import pyratbay.spectrum as ps
    >>> import pyratbay.constants as pc
    >>> import numpy as np
    >>> # Download a Kurucz stellar model file from:
    >>> # http://kurucz.harvard.edu/grids/gridp00odfnew/fp00k0odfnew.pck
    >>> # Read a single model from the kurucz file:
    >>> kfile = 'fp00k0odfnew.pck'
    >>> tsun = 5770.0  # Sun's surface temperature
    >>> gsun = 4.44    # Sun's surface gravity (log)
    >>> flux, wn, ktemp, klogg = ps.read_kurucz(kfile, tsun, gsun)
    >>> # Compute brightness at 1 AU from a 1 Rsun radius star:
    >>> s = np.trapz(flux, wn) * (pc.rsun/pc.au)**2
    >>> print("Solar constant [T={:.0f} K, logg={:.1f}]:  S = {:.1f} W m-2".
    >>>       format(ktemp, klogg, s * 1e-3))
    Solar constant [T=5750 K, logg=4.5]:  S = 1340.0 W m-2
    >>> # Pretty close to the solar constant: ~1361 W m-2

    >>> # Read the whole set of models in file:
    >>> # (in this case, ktemp and klogg are 1D arrays)
    >>> fluxes, wn, ktemp, klogg, continua = ps.read_kurucz(kfile)
    """
    # Read file into memory:
    with open(filename, 'r') as f:
        lines = f.readlines()

    iheaders = [i for i,line in enumerate(lines) if line.startswith('TEFF')]
    headers = [lines[i].strip() for i in iheaders]
    ktemp = np.array([line[ 5:12] for line in headers], np.double)
    klogg = np.array([line[22:29] for line in headers], np.double)

    # Get wavelength array (in nm):
    i = 0
    while lines[i].strip() != 'END':
        i += 1
    wl_start = i + 1
    wl_end = iheaders[0]
    wavelength = np.array(''.join(lines[wl_start:wl_end]).split(), np.double)
    wavenumber = 1.0/(wavelength*pc.nm)
    # Sort by increasing wavenumber:
    wavenumber = np.flip(wavenumber, axis=0)

    nmodels = len(headers)
    nwave = len(wavenumber)
    nlines = (iheaders[1] - iheaders[0] - 1) // 2
    vsize = 10

    if temp is not None and logg is not None:
        tmodel = ktemp[np.argmin(np.abs(ktemp-temp))]
        gmodel = klogg[np.argmin(np.abs(klogg-logg))]
        imodels = np.where((ktemp == tmodel) & (klogg == gmodel))[0]
    else:
        imodels = range(nmodels)

    # Read intensity per unit frequency (erg s-1 cm-2 Hz-1 ster-1):
    intensity = np.zeros((nmodels, nwave), np.double)
    continuum = np.zeros((nmodels, nwave), np.double)
    for k,i in enumerate(imodels):
        istart = iheaders[i] + 1
        data = ''.join(lines[istart:istart+nlines]).replace('\n','')
        intensity[k] = [data[j*vsize:(j+1)*vsize] for j in range(nwave)]

        data = ''.join(lines[istart+nlines:istart+2*nlines]).replace('\n','')
        continuum[k] = [data[j*vsize:(j+1)*vsize] for j in range(nwave)]

    # Convert intensity per unit frequency to surface flux per unit
    # wavenumber (erg s-1 cm-2 cm):
    flux      = np.flip(intensity, axis=1) * 4.0*np.pi * pc.c
    continuum = np.flip(continuum, axis=1) * 4.0*np.pi * pc.c

    if temp is not None and logg is not None:
        return flux[0], wavenumber, tmodel, gmodel

    return flux, wavenumber, ktemp, klogg, continuum
