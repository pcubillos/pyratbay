import sys
import pytest
if sys.version_info.major == 3:
    import configparser
else:
    import ConfigParser as configparser


def make_config(path, cfile, reset={}, remove=[]):
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read([cfile])
    for var in remove:
        config.remove_option('pyrat', var)
    for var, val in reset.items():
        config.set('pyrat', var, val)
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
        'smaxis':  'Undefined orbital semi-major axis (smaxis)',
        'mplanet': 'Undefined planetary surface gravity, set either '
                   'gplanet or mplanet and\nrplanet.',
        'rplanet': 'Undefined planetary surface gravity, set either '
                   'gplanet or mplanet and\nrplanet.',
        'gplanet': 'Undefined planetary surface gravity (gplanet)',
        'atmfile': 'Undefined atmospheric file (atmfile).',
        'species': 'Undefined atmospheric species list (species).',
        'elements': 'Undefined atmospheric atomic composition (elements).',
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
        'path':    "Undefined observing geometry (path).  Select between "
                   "'transit' or 'eclipse'.",
        'outspec': 'Undefined output spectrum file (outspec).',
        'tlifile': 'TLI file (tlifile) does not exist',
        ''
         # Transmission
        'rstar': 'Undefined stellar radius (rstar), required for '
                 'transmission calculation.',
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
        'hazes':   'Invalid aerosol model (hazes): invalid. Select from',
        'alkali':  'Invalid alkali model (alkali): invalid. Select from',
        'path':    'Invalid observing geometry (path): invalid. Select from',
        'tmodel':  'Invalid temperature model (tmodel): invalid. Select from',
        'retflag': 'Invalid retrieval flag (retflag): invalid. Select from'
    }
    return data

@pytest.fixture
def invalid_file():
    data = {
        'atmfile':  'Atmospheric file (atmfile) does not exist',
        'tlifile':  'TLI file (tlifile) does not exist',
        'outspec':  'Output spectrum file (outspec) does not exist',
        'mcmcfile': 'MCMC file (mcmcfile) does not exist',
        'extfile':  'Extinction-coefficient file (extfile) does not exist',
        'ptfile':   'Pressure-temperature file (ptfile) does not exist',
        'logfile':  'Log file (logfile) does not exist',
        'starspec': 'Stellar spectrum file (starspec) does not exist',
        'kurucz':   'Kurucz model file (kurucz) does not exist',
        'marcs':    'MARCS model file (marcs) does not exist',
        'phoenix':  'PHOENIX model file (phoenix) does not exist',
        'filter':   'Filter pass-bands file (filter) does not exist',
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
        'outspec':  'Folder for Output spectrum file (outspec) does not exist',
        'mcmcfile': 'Folder for MCMC file (mcmcfile) does not exist',
        'extfile':  'Folder for Extinction-coefficient file (extfile) does '
                    'not exist',
        'ptfile':   'Folder for Pressure-temperature file (ptfile) does '
                    'not exist',
        'logfile':  'Folder for Log file (logfile) does not exist',
    }
    return data


