import sys
import pytest
if sys.version_info.major == 3:
    import configparser
else:
    import ConfigParser as configparser


def replace(cfg, param, value):
    list_cfg = cfg.split('\n')
    i = 0
    while not list_cfg[i].startswith(param):
        i += 1
    list_cfg[i] = '{:s} = {:s}'.format(param, value)
    return '\n'.join(list_cfg)


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
def configs():
    pt_iso = u"""[pyrat]
runmode = pt
punits  = bar
runits  = km
ptop    = 1e-6 bar
pbottom = 100.0 bar
nlayers = 81
tmodel  = isothermal
tpars   = 1500.0
verb    = 2"""

    pt_tcea = u"""[pyrat]
runmode  = pt
punits   = bar
runits   = km
ptop    = 1e-6 bar
pbottom = 100.0 bar
nlayers = 81
tmodel  = TCEA
tpars   = -1.5 -0.8 -0.8 0.5 1.0
rstar   = 0.756 rsun
tstar   = 5040.0
smaxis  = 0.031 au
mplanet = 1.13 mjup
rplanet = 1.134 rjup
tint    = 100.0
verb    = 2"""

    atm_uniform = u"""[pyrat]
runmode  = atmosphere
punits   = bar
runits   = km
ptop     = 1e-6 bar
pbottom  = 100.0 bar
nlayers  = 81
tmodel   = TCEA
tpars    = -1.5 -0.8 -0.8 0.5 1.0
rstar    = 0.756 rsun
tstar    = 5040.0
smaxis   = 0.031 au
mplanet  = 1.13 mjup
rplanet  = 1.134 rjup
tint     = 100.0
atmfile  = test.atm
species  = H2   He    Na   H2O  CH4  CO   CO2
uniform  = 0.85 0.149 3e-6 4e-4 1e-4 5e-4 1e-7
verb     = 2"""

    atm_tea = u"""[pyrat]
runmode  = atmosphere
punits   = bar
runits   = km
ptop     = 1e-6 bar
pbottom  = 100.0 bar
nlayers  = 81
tmodel   = TCEA
tparams  = -1.5  -0.8  -0.8  0.5  1.0
rstar    = 0.756 rsun
tstar    = 5040.0
smaxis   = 0.031 au
mplanet  = 1.13 mjup
rplanet  = 1.134 rjup
tint     = 100.0
atmfile  = test.atm
elements = H He C N O Na K
xsolar   = 1.0
species  = H2 He Na K H2O CH4 CO CO2 NH3 HCN N2
verb     = 2"""

    data = {
        'pt_iso':  pt_iso,
        'pt_tcea': pt_tcea,
        'atm_uniform': atm_uniform,
        'atm_tea': atm_tea,
    }
    return data


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
