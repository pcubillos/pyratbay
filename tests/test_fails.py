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


@pytest.mark.parametrize('runmode',
    ['None', 'invalid'])
def test_init_runmode(capfd, runmode):
    pyrat = pb.init('fail_runmode_{:s}.cfg'.format(runmode.lower()))
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'argum.py', function: 'checkinputs', line:",
            "Invalid runmode ({:s}). Select from: ['tli', 'pt', 'atmosphere',".
             format(runmode)]
    for cap in caps:
        assert cap in captured.out


@pytest.mark.parametrize('call', ['command_line', 'interpreter'])
def test_run_tli(capfd, call):
    if call == 'interpreter':
        pyrat = pb.run('fail_tli_tlifile.cfg')
        assert pyrat is None
    elif call == 'command_line':
        subprocess.call('../pbay.py -c fail_tli_tlifile.cfg'.split())
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'run', line:",
            "No output TLI file specified."]
    for cap in caps:
        assert cap in captured.out


@pytest.mark.parametrize('call,runmode',
    [('command_line', 'opacity'),
     ('command_line', 'mcmc'),
     ('interpreter', 'opacity'),
     ('interpreter', 'mcmc')])
def test_run_opacity_extfile(capfd, call, runmode):
    config = 'fail_{:s}_extfile.cfg'.format(runmode)
    if call == 'interpreter':
        pyrat = pb.run(config)
        assert pyrat is None
    elif call == 'command_line':
        subprocess.call('../pbay.py -c {:s}'.format(config).split())
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'run', line:",
            "Unspecified extinction-coefficient file (extfile)."]
    for cap in caps:
        assert cap in captured.out


@pytest.mark.parametrize('call', ['command_line', 'interpreter'])
def test_run_mcmc_mcmcfile(capfd, call):
    config = 'fail_mcmc_mcmcfile.cfg'
    if call == 'interpreter':
        pyrat = pb.run(config)
        assert pyrat is None
    elif call == 'command_line':
        subprocess.call('../pbay.py -c {:s}'.format(config).split())
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'run', line:",
            "No MCMC file specified."]
    for cap in caps:
        assert cap in captured.out

