# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import subprocess
import pytest

from conftest import replace

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay as pb

os.chdir(ROOT+'tests')


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
def test_invalid_integer_type(tmp_path, configs, capfd, param, value):
    cfg = replace(configs['pt_iso'], param, value)
    path = tmp_path / 'test.cfg'
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
def test_invalid_float_type(tmp_path, configs, capfd, param, value):
    cfg = replace(configs['pt_tcea'], param, value)
    path = tmp_path / 'test.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert ("Invalid data type for {:s}, could not convert string to "
            "float: '{:s}'".format(param, value)) in captured.out


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# pt runmode fails:
@pytest.mark.parametrize('missing', ['nlayers', 'ptop', 'pbottom'])
def test_pressure_missing(tmp_path, configs, capfd, undefined, missing):
    cfg = configs['pt_iso'].replace(missing, 'dummy')
    path = tmp_path / 'test.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'driver.py', function: 'check_pressure'" \
           in captured.out
    assert undefined[missing] in captured.out


@pytest.mark.parametrize('value', ['a', '10.0 20.0', '10.0 bar 30.0'])
@pytest.mark.parametrize('param', ['ptop', 'pbottom'])
def test_pressure_invalid_type(tmp_path, configs, capfd, param, value):
    cfg = replace(configs['pt_iso'], param, value)
    path = tmp_path / 'test.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'atmosphere.py', function: 'pressure'" \
           in captured.out
    assert "Invalid value '{:s}' for parameter {:s}.". \
           format(value, param) in captured.out


@pytest.mark.parametrize('param', ['ptop', 'pbottom'])
def test_negative_pressures(tmp_path, configs, capfd, param):
    cfg = replace(configs['pt_iso'], param, '-10')
    path = tmp_path / 'test.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'atmosphere.py', function: 'pressure'" \
           in captured.out
    assert "{:s} must be > 0.0".format(param) in captured.out


def test_negative_nlayers(tmp_path, configs, capfd):
    cfg = replace(configs['pt_iso'], 'nlayers', '-10')
    path = tmp_path / 'test.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'driver.py', function: 'check_pressure'" \
           in captured.out
    assert "Number of atmospheric layers (nlayers) must be > 0" in captured.out


def test_pt_invalid_tmodel(tmp_path, configs, capfd):
    cfg = replace(configs['pt_iso'], 'tmodel', 'invalid')
    path = tmp_path / 'test.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'atmosphere.py', function: 'temperature'",
            "Invalid input temperature model 'invalid'.  Select from: 'TCEA',",
            "'MadhuInv', 'MadhuNoInv', or 'isothermal'."]
    for cap in caps:
        assert cap in captured.out


@pytest.mark.parametrize('missing', ['tmodel', 'tpars'])
def test_pt_temperature_missing(tmp_path, configs, capfd, missing, undefined):
    cfg = configs['pt_iso'].replace(missing, 'dummy')
    path = tmp_path / 'test.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_temp'" \
           in captured.out
    assert undefined[missing] in captured.out


@pytest.mark.parametrize('cfg_file, error',
    [('pt_iso', 'isothermal temperature model (1).'),
     ('pt_tcea', 'TCEA temperature model (5).')])
def test_pt_tpars(cfg_file, configs, error, capfd, tmp_path):
    cfg = replace(configs[cfg_file], 'tpars', '1000.0 2000.0')
    path = tmp_path / 'test.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'check_temp', line:",
            "Wrong number of parameters (2) for the {:s}".format(error)]
    for cap in caps:
        assert cap in captured.out


@pytest.mark.parametrize('missing',
    ['rstar', 'tstar', 'smaxis', 'mplanet', 'rplanet',])
def test_tcea_missing(tmp_path, configs, capfd, missing, undefined):
    cfg = configs['pt_tcea'].replace(missing, 'dummy')
    path = tmp_path / 'test.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'driver.py', function: 'check_temp'" \
           in captured.out
    assert undefined[missing] in captured.out


@pytest.mark.parametrize('missing',
    ['rstar', 'tstar', 'smaxis', 'mplanet', 'rplanet',])
def test_tcea_negatives(tmp_path, configs, capfd, missing):
    cfg = replace(configs['pt_tcea'], missing, '-10.0')
    if missing == 'mplanet':
        cfg = replace(cfg, missing, '-1.0 mjup')
    path = tmp_path / 'test.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "{:s} must be > 0.0".format(missing) in captured.out


def test_tcea_missing_mass_units(tmp_path, configs, capfd):
    cfg = replace(configs['pt_tcea'], 'mplanet', '1.0')
    path = tmp_path / 'test.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Invalid units 'None' for parameter mplanet." in captured.out

