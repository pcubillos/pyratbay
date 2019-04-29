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
        'rstar':   'Undefined stellar radius (rstar).',
        'tstar':   'Undefined stellar temperature (tstar).',
        'smaxis':  'Undefined orbital semi-major axis (smaxis).',
        'mplanet': 'Undefined planetary surface gravity, set either '
                   'gplanet or mplanet and\nrplanet.',
        'rplanet': 'Undefined planetary surface gravity, set either '
                   'gplanet or mplanet and\nrplanet.',
#        '':'',
    }

    return data
