# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import subprocess
import pytest

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay as pb

os.chdir(ROOT+'tests')


pt_iso_cfg = u"""[pyrat]
runmode = pt
punits  = bar
runits  = km
ptop    = 1e-6 bar
pbottom = 100.0 bar
nlayers = 81
tmodel  = isothermal
tpars   = 1500.0
verb    = 2"""

pt_tcea_cfg = u"""[pyrat]
runmode  = pt
punits   = bar
radunits = km
ptop    = 1e-6 bar
pbottom = 100.0 bar
nlayers = 81
verb    = 2
tmodel  = TCEA
tpars   = -1.5 -0.8 -0.8 0.5 1.0
rstar   = 0.756 rsun
tstar   = 5040.0
smaxis  = 0.031 au
mplanet = 1.13 mjup
rplanet = 1.134 rjup
tint    = 100.0"""


def replace(cfg, param, value):
    list_cfg = cfg.split('\n')
    i = 0
    while not list_cfg[i].startswith(param):
        i += 1
    list_cfg[i] = '{:s} = {:s}'.format(param, value)
    return '\n'.join(list_cfg)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check runmode and command_line vs interpreter runs:
@pytest.mark.parametrize('call, runmode',
    [('command_line', 'None'),
     ('command_line', 'invalid'),
     ('interpreter', 'None'),
     ('interpreter', 'invalid')])
def test_run_runmode(capfd, runmode, call):
    config = 'fail_runmode_{:s}.cfg'.format(runmode.lower())
    if call == 'interpreter':
        pyrat = pb.run(config)
        assert pyrat is None
    elif call == 'command_line':
        subprocess.call('../pbay.py -c {:s}'.format(config).split())
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'run', line:",
            "Invalid runmode ({:s}). Select from: ['tli', 'pt', 'atmosphere',".
             format(runmode)]
    for cap in caps:
        assert cap in captured.out


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check output for each case is defined:
def test_run_tli(capfd):
    pyrat = pb.run('fail_tli_tlifile.cfg')
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'run', line:",
            "No output TLI file specified."]
    for cap in caps:
        assert cap in captured.out


@pytest.mark.parametrize('runmode',
    ['opacity', 'mcmc'])
def test_run_opacity_extfile(capfd, runmode):
    pyrat = pb.run('fail_{:s}_extfile.cfg'.format(runmode))
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'run', line:",
            "Unspecified extinction-coefficient file (extfile)."]
    for cap in caps:
        assert cap in captured.out


def test_run_mcmc_mcmcfile(capfd):
    pyrat = pb.run('fail_mcmc_mcmcfile.cfg')
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'run', line:",
            "No MCMC file specified."]
    for cap in caps:
        assert cap in captured.out


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check integer and float data types:
# Same applies to any other integer/float input variables.
@pytest.mark.parametrize('param, value',
    [('nlayers', '10.5'),
     ('nlayers', '10 20'),
     ('nlayers', 'a')])
def test_invalid_integer_type(tmp_path, capfd, param, value):
    cfg = replace(pt_iso_cfg, param, value)
    path = tmp_path / 'pt_fail.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert ("Invalid data type for {:s}, could not convert string to "
            "integer: '{:s}'".format(param, value)) in captured.out


@pytest.mark.parametrize('param, value',
    [('tstar', '100 200'),
     ('tstar', 'a')])
def test_invalid_float_type(tmp_path, capfd, param, value):
    cfg = replace(pt_tcea_cfg, param, value)
    path = tmp_path / 'pt_fail.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert ("Invalid data type for {:s}, could not convert string to "
            "float: '{:s}'".format(param, value)) in captured.out


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# pt runmode fails:
@pytest.mark.parametrize('arg_del, error',
    [('nlayers', 'Undefined number of atmospheric layers (nlayers).'),
     ('ptop',    'Undefined atmospheric top pressure (ptop)'),
     ('pbottom', 'Undefined atmospheric bottom pressure (pbottom)')])
def test_pressure_missing(arg_del, error, capfd, tmp_path):
    cfg = pt_iso_cfg.replace(arg_del, 'dummy')
    path = tmp_path / 'pt_fail.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'driver.py', function: 'check_pressure'" \
           in captured.out
    assert error in captured.out


@pytest.mark.parametrize('value', ['a', '10.0 20.0', '10.0 bar 30.0'])
@pytest.mark.parametrize('param', ['ptop', 'pbottom'])
def test_pressure_invalid_type(capfd, tmp_path, param, value):
    cfg = replace(pt_iso_cfg, param, value)
    path = tmp_path / 'pt_fail.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'atmosphere.py', function: 'pressure'" \
           in captured.out
    assert "Invalid value '{:s}' for parameter {:s}.". \
           format(value, param) in captured.out


@pytest.mark.parametrize('param', ['ptop', 'pbottom'])
def test_negative_pressures(capfd, tmp_path, param):
    cfg = replace(pt_iso_cfg, param, '-10')
    path = tmp_path / 'pt_fail.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'atmosphere.py', function: 'pressure'" \
           in captured.out
    assert "{:s} must be > 0.0".format(param) in captured.out

def test_negative_nlayers(capfd, tmp_path):
    cfg = replace(pt_iso_cfg, 'nlayers', '-10')
    path = tmp_path / 'pt_fail.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'driver.py', function: 'check_pressure'" \
           in captured.out
    assert "Number of atmospheric layers (nlayers) must be > 0" in captured.out


def test_pt_invalid_tmodel(capfd, tmp_path):
    cfg = replace(pt_iso_cfg, 'tmodel', 'invalid')
    path = tmp_path / 'pt_missing.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'atmosphere.py', function: 'temperature'",
            "Invalid input temperature model 'invalid'.  Select from: 'TCEA',",
            "'MadhuInv', 'MadhuNoInv', or 'isothermal'."]
    for cap in caps:
        assert cap in captured.out


@pytest.mark.parametrize('param, error',
    [('tmodel', 'temperature model (tmodel)'),
     ('tpars',  'temperature-model parameters (tpars)')])
def test_pt_temperature_missing(tmp_path, capfd, param, error):
    cfg = pt_iso_cfg.replace(param, 'dummy')
    path = tmp_path / 'pt_missing.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'check_temp', line:",
            "Undefined {:s}.".format(error)]
    for cap in caps:
        assert cap in captured.out


@pytest.mark.parametrize('cfg_file, error',
    [(pt_iso_cfg, 'isothermal temperature model (1).'),
     (pt_tcea_cfg, 'TCEA temperature model (5).')])
def test_pt_tpars(cfg_file, error, capfd, tmp_path):
    cfg = replace(cfg_file, 'tpars', '1000.0 2000.0')
    path = tmp_path / 'pt_fail.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'check_temp', line:",
            "Wrong number of parameters (2) for the {:s}".format(error)]
    for cap in caps:
        assert cap in captured.out


@pytest.mark.parametrize('arg_del, error',
    [('rstar',   'Undefined stellar radius (rstar).'),
     ('tstar',   'Undefined stellar temperature (tstar).'),
     ('smaxis',  'Undefined orbital semi-major axis (smaxis).'),
     ('mplanet', 'Undefined planetary surface gravity, set either '
                 'gplanet or mplanet and\nrplanet.'),
     ('rplanet', 'Undefined planetary surface gravity, set either '
                 'gplanet or mplanet and\nrplanet.')])
def test_tcea_missing(arg_del, error, capfd, tmp_path):
    cfg = pt_tcea_cfg.replace(arg_del, 'dummy')
    path = tmp_path / 'pt_missing.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'driver.py', function: 'check_temp'" \
           in captured.out
    assert error in captured.out

