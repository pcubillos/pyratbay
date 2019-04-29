# Copyright (c) 2016-2019 Patricio Cubillos and contributors.
# Pyrat Bay is currently proprietary software (see LICENSE).

import os
import sys
import subprocess
import pytest

from conftest import make_config

ROOT = os.path.realpath(os.path.dirname(__file__) + '/..') + '/'
sys.path.append(ROOT)
import pyratbay as pb

os.chdir(ROOT+'tests')


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check runmode and command_line vs interpreter runs:
@pytest.mark.parametrize('runmode', ['None', 'invalid'])
@pytest.mark.parametrize('call',    ['command_line', 'interpreter'])
def test_run_runmode(tmp_path, capfd, runmode, call):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'runmode':runmode})
    if call == 'interpreter':
        pyrat = pb.run(cfg)
        assert pyrat is None
    elif call == 'command_line':
        subprocess.call('../pbay.py -c {:s}'.format(cfg).split())
    captured = capfd.readouterr()
    caps = ["Error in module: 'parser.py', function: 'parse'",
            "Invalid running mode (runmode): {:s}. Select from: ['tli', 'pt',".format(runmode)]
    for cap in caps:
        assert cap in captured.out


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check output for each case is defined:
def test_run_tli(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/tli_hitran_1.1-1.7um_test.cfg',
        remove=['tlifile'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'run'",
            "Undefined TLI file (tlifile)."]
    for cap in caps:
        assert cap in captured.out


@pytest.mark.parametrize('cfile',
    ['opacity_test.cfg', 'mcmc_transmission_test.cfg'])
def test_run_opacity_extfile(tmp_path, capfd, cfile):
    cfg = make_config(tmp_path, ROOT+'tests/{:s}'.format(cfile),
        remove=['extfile'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'run'",
            "Undefined extinction-coefficient file (extfile)."]
    for cap in caps:
        assert cap in captured.out


def test_run_mcmc_mcmcfile(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/mcmc_transmission_test.cfg',
        remove=['mcmcfile'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'run'",
            "Undefined MCMC file (mcmcfile)."]
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
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={param:value})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert ("Invalid data type for {:s}, could not convert string to "
            "integer: '{:s}'".format(param, value)) in captured.out


@pytest.mark.parametrize('param, value',
    [('tstar', '100 200'),
     ('tstar', 'a')])
def test_invalid_float_type(tmp_path, capfd, param, value):
    cfg = make_config(tmp_path, ROOT+'tests/pt_tcea.cfg',
        reset={param:value})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert ("Invalid data type for {:s}, could not convert string to "
            "float: '{:s}'".format(param, value)) in captured.out


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# pt runmode fails:
@pytest.mark.parametrize('param', ['nlayers', 'ptop', 'pbottom'])
def test_pressure_missing(tmp_path, capfd, undefined, param):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        remove=[param])
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'driver.py', function: 'check_pressure'" \
           in captured.out
    assert undefined[param] in captured.out


@pytest.mark.parametrize('value', ['a', '10.0 20.0', '10.0 bar 30.0'])
@pytest.mark.parametrize('param', ['ptop', 'pbottom'])
def test_pressure_invalid_type(tmp_path, capfd, param, value):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={param:value})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert "Invalid value '{:s}' for parameter {:s}.". \
           format(value, param) in captured.out


@pytest.mark.parametrize('param', ['ptop', 'pbottom'])
def test_negative_pressures(tmp_path, capfd, param):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={param:'-10'})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert "({:s}) must be > 0.0".format(param) in captured.out


def test_negative_nlayers(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={'nlayers':'-10'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert "Number of atmospheric layers (nlayers) must be > 0" in captured.out


def test_pt_invalid_tmodel(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={'tmodel':'invalid'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'parser.py', function: 'parse'",
            "Invalid temperature model (tmodel): invalid. Select from: ['isothermal',",
            "'TCEA', 'MadhuInv', 'MadhuNoInv']."]
    for cap in caps:
        assert cap in captured.out


@pytest.mark.parametrize('param', ['tmodel', 'tpars'])
def test_pt_temperature_missing(tmp_path, capfd, param, undefined):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_temp'" \
           in captured.out
    assert undefined[param] in captured.out


@pytest.mark.parametrize('cfile, error',
    [('pt_isothermal.cfg', 'isothermal temperature model (1).'),
     ('pt_tcea.cfg',       'TCEA temperature model (5).')])
def test_pt_tpars(tmp_path, capfd, cfile, error):
    cfg = make_config(tmp_path, ROOT+'tests/{:s}'.format(cfile),
        reset={'tpars':'100.0 200.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'check_temp'",
            "Wrong number of parameters (2) for the {:s}".format(error)]
    for cap in caps:
        assert cap in captured.out


@pytest.mark.parametrize('param',
    ['rstar', 'tstar', 'smaxis', 'mplanet', 'rplanet',])
def test_tcea_missing(tmp_path, capfd, param, undefined):
    cfg = make_config(tmp_path, ROOT+'tests/pt_tcea.cfg',
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_temp'" \
           in captured.out
    assert undefined[param] in captured.out


@pytest.mark.parametrize('param, value',
    [('rstar', '-1.0 rsun'),
     ('tstar', '-1000'),
     ('smaxis', '-0.01 au'),
     ('mplanet', '-1.0 mjup'),
     ('rplanet', '-1.0 rjup')])
def test_tcea_negatives(tmp_path, capfd, param, value):
    cfg = make_config(tmp_path, ROOT+'tests/pt_tcea.cfg',
        reset={param:value})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "({:s}) must be > 0.0".format(param) in captured.out


def test_tcea_missing_mass_units(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/pt_tcea.cfg',
        reset={'mplanet':'1.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Invalid units 'None' for parameter mplanet." in captured.out

