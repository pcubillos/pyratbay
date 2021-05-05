# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import sys
import itertools
import pytest
import configparser
import pathlib
import tempfile


def pytest_collection_modifyitems(items):
    """Sort tests with a 'sort' mark by order keyword."""
    order = [item.get_closest_marker('sort').kwargs['order']
             if item.get_closest_marker('sort') is not None
             else -1
             for item in items]

    last = itertools.count(1 + max(order) if order else 0)
    order = {item:val if val >= 0 else next(last)
             for item,val in zip(items,order)}
    items[:] = sorted(order, key=order.get)


def make_config(path, cfile, reset={}, remove=[]):
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read([cfile])
    for var, val in reset.items():
        config.set('pyrat', var, val)
    for var in remove:
        config.remove_option('pyrat', var)
    cfg_file = str(path / 'test.cfg')
    with open(cfg_file, 'w') as cfg:
        config.write(cfg)
    return cfg_file


@pytest.fixture
def undefined():
    data = {
        'nlayers': 'Undefined number of atmospheric layers (nlayers).',
        'ptop':    'Undefined atmospheric top pressure (ptop)',
        'pbottom': 'Undefined atmospheric bottom pressure (pbottom)',
        'tmodel':  'Undefined temperature model (tmodel)',
        'tpars':   'Undefined temperature-model parameters (tpars)',
        'rstar':   'Undefined stellar radius (rstar)',
        'tstar':   'Undefined stellar temperature (tstar)',
        'gstar':   'Undefined stellar gravity (gstar)',
        'smaxis':  'Undefined orbital semi-major axis (smaxis)',
        'mplanet': 'Undefined planetary surface gravity, set either '
                   'gplanet or mplanet and\nrplanet.',
        'rplanet': 'Undefined planetary surface gravity, set either '
                   'gplanet or mplanet and\nrplanet.',
        'gplanet': 'Undefined planetary surface gravity (gplanet)',
        'atmfile': 'Undefined atmospheric file (atmfile).',
        'species': 'Undefined atmospheric species list (species).',
        'elements': 'Undefined elemental composition list (elements) for '
                    'tea chemistry model.',
        'uniform': 'Undefined list of uniform volume mixing ratios (uniform) '
                   'for uniform\nchemistry model.',
        'refpressure': 'Cannot compute hydrostatic-equilibrium radius profile.'
                       '  Undefined reference\npressure level (refpressure).',
    }
    return data

@pytest.fixture
def undefined_spec():
    data = {
        'wllow':   'High wavenumber boundary is undefined.  Either set '
                   'wnhigh or wllow.',
        'wlhigh':  'Low wavenumber boundary is undefined.  Either set '
                   'wnlow or wlhigh.',
        'wnstep':  'Undefined wavenumber sampling step size (wnstep).',
        'wnosamp': 'Undefined wavenumber oversampling factor (wnosamp).',
        'rt_path': 'Undefined radiative-transfer observing geometry (rt_path).'
                   '  Select from',
        'specfile': 'Undefined output spectrum file (specfile).',
        'tlifile': 'TLI file (tlifile) does not exist',
         # Transmission
        'rstar': 'Undefined stellar radius (rstar), required for '
                 'transmission calculation.',
        }
    return data


@pytest.fixture
def undefined_opacity():
    data = {
        'tmin': 'Undefined lower temperature boundary (tmin) for '
                'extinction-coefficient grid.',
        'tmax': 'Undefined upper temperature boundary (tmax) for '
                'extinction-coefficient grid.',
        'tstep': 'Undefined temperature sampling step (tstep) for '
                 'extinction-coefficient grid.',
        'tlifile': 'Requested extinction-coefficient table, but there are '
                   'no input TLI files.',
    }
    return data


@pytest.fixture
def undefined_mcmc():
    data = {
        'retflag':"Undefined retrieval model flags.  Select from ['temp', "
                  "'rad', 'mol', 'ray',\n'cloud', 'patchy', 'mass'].",
        'params': 'Undefined retrieval fitting parameters (params).',
        'data':   'Undefined transit/eclipse data (data).',
        'uncert': 'Undefined data uncertainties (uncert).',
        'filters': 'Undefined transmission filters (filters).',
        'sampler': 'Undefined retrieval algorithm (sampler).  Select '
                   'from [snooker].',
        'nsamples': 'Undefined number of retrieval samples (nsamples).',
        'burnin':   'Undefined number of retrieval burn-in samples (burnin).',
        'nchains':  'Undefined number of retrieval parallel chains (nchains).',
        'rstar':    'Undefined radius ratio (need rplanet and rstar).',
        'tmodel':   'Requested temp in retflag, but there is no tmodel.',
        'rayleigh': 'Requested ray in retflag, but there are no rayleigh '
                    'models.',
        'clouds': 'Requested cloud in retflag, but there are no cloud models.',
        'molmodel': "Requested mol in retflag, but there is no 'molmodel'.",
        'bulk': 'Requested mol in retflag, but there are no bulk species.',
    }
    return data


@pytest.fixture
def invalid_raygrid():
    data = {
        '10 60 90': 'First angle in raygrid must be 0.0 (normal to surface).',
        '0 30 60 100': 'raygrid angles must lie between 0 and 90 deg.',
        '0 30 90 60': 'raygrid angles must be monotonically increasing.',
    }
    return data

@pytest.fixture
def invalid():
    data = {
        'runmode': 'Invalid running mode (runmode): invalid. Select from',
        'rayleigh':'Invalid Rayleigh model (rayleigh): invalid. Select from',
        'clouds':  'Invalid cloud model (clouds): invalid. Select from',
        'alkali':  'Invalid alkali model (alkali): invalid. Select from',
        'rt_path': 'Invalid radiative-transfer observing geometry (rt_path):'
                   ' invalid. Select\nfrom',
        'tmodel':  'Invalid temperature model (tmodel): invalid. Select from',
        'molmodel': 'Invalid molecular-abundance model (molmodel): invalid. '
                    'Select from',
        'retflag': 'Invalid retrieval flag (retflag): invalid. Select from'
    }
    return data

@pytest.fixture
def invalid_file():
    data = {
        'atmfile':  'Atmospheric file (atmfile) does not exist',
        'tlifile':  'TLI file (tlifile) does not exist',
        'specfile': 'Spectrum file (specfile) does not exist',
        'mcmcfile': 'MCMC file (mcmcfile) does not exist',
        'extfile':  'Extinction-coefficient file (extfile) does not exist',
        'ptfile':   'Pressure-temperature file (ptfile) does not exist',
        'logfile':  'Log file (logfile) does not exist',
        'starspec': 'Stellar spectrum file (starspec) does not exist',
        'kurucz':   'Kurucz model file (kurucz) does not exist',
        'marcs':    'MARCS model file (marcs) does not exist',
        'phoenix':  'PHOENIX model file (phoenix) does not exist',
        'filters':  'Filter pass-bands file (filters) does not exist',
        'dblist':   'Opacity database file (dblist) does not exist',
        'molfile':  'Molecular data file (molfile) does not exist',
        'csfile':   'Cross-section file (csfile) does not exist',
    }
    return data


@pytest.fixture
def invalid_path():
    data = {
        'atmfile':  'Folder for Atmospheric file (atmfile) does not exist',
        'tlifile':  'Folder for TLI file (tlifile) does not exist',
        'specfile': 'Folder for Spectrum file (specfile) does not exist',
        'mcmcfile': 'Folder for MCMC file (mcmcfile) does not exist',
        'extfile':  'Folder for Extinction-coefficient file (extfile) does '
                    'not exist',
        'ptfile':   'Folder for Pressure-temperature file (ptfile) does '
                    'not exist',
        'logfile':  'Folder for Log file (logfile) does not exist',
    }
    return data


@pytest.fixture
def invalid_temp():
    data = {
        'csfile':  'One or more input temperature values lies out of '
                   'the cross-section\ntemperature boundaries '
                   '(K): [  60.0, 3000.0].',
        'tlifile': 'One or more input temperature values lies out of '
                   'the line-transition\ntemperature boundaries '
                   '(K): [   1.0, 5000.0]',
        'extfile': 'One or more input temperature values lies out of the '
                   'tabulated\nextinction-coefficient temperature boundaries '
                   '(K): [ 300.0, 3000.0].',
    }
    return data


