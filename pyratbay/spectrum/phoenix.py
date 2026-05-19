# Copyright (c) 2021-2026 Cubillos & Blecic
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'list_phoenix_files',
    'fetch_phoenix',
    'read_phoenix',
]

import os

import numpy as np
import requests
import h5py

from .. import constants as pc


# https://www.fdr.uni-hamburg.de/record/18108
FDRV1 = 'https://www.fdr.uni-hamburg.de/record/16738/files/'
FDRV3 = 'https://www.fdr.uni-hamburg.de/record/17670/files/'


def list_phoenix_files(teff=None, logg=None, metal=None):
    """
    Get the list of PHOENIX NewEra models with closest parameters to
    input values.  Leave input as None to get all models.

    Parameters
    ----------
    teff: Float
        Effective temperature (K).
    logg: Float
        log of surface gravity.
    metal: Float
        Metallicity [M/H].

    Returns
    -------
    sed_files: List of strings
        PHOENIX SED file names.

    Examples
    --------
    >>> import pyratbay.spectrum as ps
    >>>
    >>> # List all temperature models with given log_g and metallicity:
    >>> teff = None
    >>> logg = 4.57
    >>> metal = 0.35
    >>> sed_models = ps.list_phoenix_files(teff, logg, metal)
    """
    # Fetch PHOENIX SED list if needed:
    sed_list_file = f'{pc.ROOT}pyratbay/data/phoenix_seds.dat'
    if not os.path.exists(sed_list_file):
        url = 'https://www.fdr.uni-hamburg.de/record/18108/files/list_of_available_NewEraV3_models.txt'
        response = requests.get(url, stream=True)
        if not response.ok:
            raise ValueError('HTTP request failed')
        with open(sed_list_file, mode="wb") as file:
            for chunk in response.iter_content(chunk_size=10*1024):
                file.write(chunk)


    files = np.loadtxt(sed_list_file, dtype=str, unpack=True, skiprows=1)[1]
    t_effs = np.array([file[3:8] for file in files], float)
    log_gs = np.array([file[9:13] for file in files], float)
    metals = np.array([file[13:17] for file in files], float)

    mask = np.ones(len(files), bool)

    # alpha!=0 not enabled at the moment
    alpha = None
    if alpha is None:
        has_alpha = np.array(['alpha' in file for file in files])
        mask &= ~has_alpha

    all_t = np.unique(t_effs[mask])
    if teff is not None:
        nearest_t = all_t[np.argmin(np.abs(teff-all_t))]
        t_mask = t_effs==nearest_t
        mask &= t_mask

    all_g = np.unique(log_gs[mask])
    if logg is not None:
        nearest_g = all_g[np.argmin(np.abs(logg-all_g))]
        g_mask = log_gs==nearest_g
        mask &= g_mask

    all_z = np.unique(metals[mask])
    if metal is not None:
        nearest_z = all_z[np.argmin(np.abs(metal-all_z))]
        z_mask = metals==nearest_z
        mask &= z_mask

    return sorted(files[mask].tolist())


def _download_file(url, filename):
    response = requests.get(url, stream=True)
    if not response.ok:
        raise ValueError('HTTP request failed')

    with open(filename, 'wb') as file:
        for chunk in response.iter_content(chunk_size=10*1024):
            file.write(chunk)


def fetch_phoenix(teff, logg, metal, folder='.'):
    """
    Fetch PHOENIX New-Era stellar SED models (url request) that matches the
    closest to the input metallicity, effective temperature, and log(g).

    Parameters
    ----------
    teff: Float
        Effective temperature (K). Set to None to fetch all teff models.
    logg: Float
        log of surface gravity. Set to None to fetch all logg models.
    metal: Float
        Metallicity [M/H]. Set to None to fetch all metallicity models.
    folder: String
        Folder where to store the PHOENIX files.

    Examples
    --------
    >>> import pyratbay.spectrum as ps
    >>>
    >>> # Download all models with fixed log_g and metallicity:
    >>> teff = None
    >>> logg = 4.57
    >>> metal = 0.35
    >>> folder = 'phoenix/'
    >>> ps.fetch_phoenix(metal, teff, logg, folder=folder)
    """
    files = list_phoenix_files(teff, logg, metal)

    for file in files:
        temp = int(file[3:8])
        prefix = FDRV1 if temp < 5000 else FDRV3
        url = f'{prefix}{file}?download=1'
        path = f'{folder}/{file}'
        _download_file(url, path)


def read_phoenix(filename, resolution='lo'):
    """
    Read a PHOENIX New-Era SED model.

    Parameters
    ----------
    filename: str
        Path to PHOENIX SED file.
    resolution: str
        Extract low ('lo', R~10K) or high ('hi', R~1M) resolution spectra.

    Returns
    -------
    wl: 1D float array
        (vacuum) wavelength in microns.
    flux: 1D float array
        Flux in erg s-1 cm-2 cm.

    Examples
    --------
    >>> import pyratbay.spectrum as ps
    >>> filename = 'lte05400-4.50+0.5.PHOENIX-NewEra-ACES-COND-2023.HSR.h5'
    >>> wl, flux = ps.read_phoenix(filename)
    """
    file = h5py.File(filename, 'r')
    if resolution == 'hi':
        wl_key = '/PHOENIX_SPECTRUM/wl'
        f_key = '/PHOENIX_SPECTRUM/flux'
    else:
        wl_key = '/PHOENIX_SPECTRUM_LSR/wl'
        f_key = '/PHOENIX_SPECTRUM_LSR/fl'

    # Read spectrum (wl in microns, flux in erg s-1 cm-2 cm):
    wl = file[wl_key][()] * pc.A / pc.um
    flux = 10.0**file[f_key][()] * (wl*pc.um)**2.0

    return wl, flux
