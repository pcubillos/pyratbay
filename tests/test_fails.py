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


pt_cfg = """[pyrat]
runmode = pt
punits  = bar
runits  = km
ptop    = 1e-6 bar
pbottom = 100.0 bar
nlayers = 81
tmodel  = isothermal
tpars   = 1500.0
verb    = 2"""


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


def test_pt_tmodel(capfd):
    pyrat = pb.run('fail_pt_tmodel.cfg')
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'check_temp', line:",
            "Undefined temperature model (tmodel)."]
    for cap in caps:
        assert cap in captured.out


def test_pt_isothermal_tpars(capfd):
    pyrat = pb.run('fail_isothermal_missing_tpars.cfg')
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'check_temp', line:",
            "Undefined temperature-model parameters (tpars)."]
    for cap in caps:
        assert cap in captured.out


def test_pt_isothermal_tpars(capfd):
    pyrat = pb.run('fail_isothermal_tpars.cfg')
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'check_temp', line:",
            "Wrong number of parameters (2) for the isothermal temperature model (1)."]
    for cap in caps:
        assert cap in captured.out


@pytest.mark.parametrize('arg_del, error',
    [('nlayers', 'Undefined number of atmospheric layers (nlayers).'),
     ('ptop',    'Undefined atmospheric top pressure (ptop)'),
     ('pbottom', 'Undefined atmospheric bottom pressure (pbottom)')])
def test_pressure_missing(arg_del, error, capfd, tmp_path):
    cfg = pt_cfg.replace(arg_del, 'dummy')
    path = tmp_path / 'pt_missing.cfg'
    path.write_text(cfg)

    pyrat = pb.run(str(path))
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_pressure'" \
           in captured.out
    assert error in captured.out

