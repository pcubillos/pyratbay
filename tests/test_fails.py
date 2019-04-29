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


@pytest.mark.parametrize('param',
    ['verb', 'wnosamp', 'nlayers', 'ncpu', 'nDop', 'nLor', 'quadrature',
     'nsamples', 'nchains', 'burnin', 'thinning', 'resume'])
def test_invalid_integer_all_params(tmp_path, capfd, param):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={param:'abc'})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert ("Invalid data type for {:s}, could not convert string to "
            "integer: 'abc'".format(param)) in captured.out


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


@pytest.mark.parametrize('param',
    ['wnlow', 'wnhigh', 'wnstep', 'resolution', 'xsolar', 'tmin', 'tmax',
     'tstep', 'ethresh', 'vextent', 'Dmin', 'Dmax', 'Lmin', 'Lmax', 'DLratio',
     'fpatchy', 'maxdepth', 'qcap', 'tlow', 'thigh', 'grbreak', 'grnmin',
     'gstar', 'tstar', 'gplanet', 'tint'])
def test_invalid_float_all_params(tmp_path, capfd, param):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={param:'abc'})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert ("Invalid data type for {:s}, could not convert string to "
            "float: 'abc'".format(param)) in captured.out


@pytest.mark.parametrize('param',
    ['runmode', 'rayleigh', 'hazes', 'alkali', 'path', 'tmodel', 'retflag'])
def test_invalid_choice(tmp_path, capfd, param, invalid):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={param:'invalid'})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert invalid[param] in captured.out


@pytest.mark.parametrize('param',
    ['starspec', 'kurucz', 'marcs', 'phoenix', 'filter',
     'dblist', 'molfile', 'csfile'])
def test_file_not_found(tmp_path, capfd, param, invalid_file):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={param:'nope.dat'})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert invalid_file[param] in captured.out


@pytest.mark.parametrize('param',
    ['atmfile', 'tlifile', 'extfile', 'mcmcfile', 'outspec', 'ptfile',
     'logfile'])
def test_invalid_file_path(tmp_path, capfd, param, invalid_path):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={param:'nope/file.dat'})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert invalid_path[param] in captured.out
    # params from 'test_file_not_found' do not raise folder error since
    # they catch file not found first.

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Parameter boundaries:

@pytest.mark.parametrize('param, value',
    [('wllow',   ' -1.0 um'),
     ('wlhigh',  ' -1.0 um'),
     ('wnhigh',  ' -1.0'),
     ('wnlow',   ' -1.0'),
     ('wnstep',  ' -1.0'),
     ('resolution', ' -100.0'),
     ('nlayers',  ' 1'),
     ('ptop',    ' -1.0 bar'),
     ('pbottom', ' -1.0 bar'),
     ('refpressure', ' -1.0 bar'),
     ('radhigh', ' -1.0 rjup'),
     ('radstep', ' -100.0 km'),
     ('pbottom', ' -1.0 bar'),
     ('mplanet', ' -1.0 mjup'),
     ('rplanet', ' -1.0 rjup'),
     ('gplanet', ' -1000.0'),
     ('tint',    ' -100.0'),
     ('smaxis',  ' -0.01 au'),
     ('rstar',   ' -1.0 rsun'),
     ('mstar',   ' -1.0 msun'),
     ('gstar',   ' -1000.0'),
     ('tstar',   ' -5000.0'),
     ('Dmin',    ' -1e-6'),
     ('Dmax',    ' -1e-1'),
     ('Lmin',    ' -1e-6'),
     ('Lmax',    ' -1e-1'),
     ('DLratio', ' -0.1'),
     ('tmin',    ' -100'),
     ('tmax',    ' -100'),
     ('tstep',   ' -100'),
     ('ethresh', ' -1e15'),
     ('qcap',    ' -0.5'),
     ('grnmin',  ' -0.5'),
     ('nsamples', ' 0'),
     ('burnin',   ' 0'),
    ])
def test_greater_than(tmp_path, capfd, param, value):
    cfg = make_config(tmp_path, ROOT+'tests/pt_tcea.cfg',
        reset={param:value})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert "({:s}) must be > ".format(param) in captured.out


@pytest.mark.parametrize('param',
    ['verb', 'wnosamp', 'nDop', 'nLor', 'thinning', 'nchains', 'ncpu',
     'quadrature', 'grbreak', 'radlow', 'fpatchy', 'maxdepth', 'vextent'])
def test_greater_equal(tmp_path, capfd, param):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={param:'-10'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert "({:s}) must be >= ".format(param) in captured.out


@pytest.mark.parametrize('param', ['verb'])
def test_lower_than(tmp_path, capfd, param):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={param:'10'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert "({:s}) must be < ".format(param) in captured.out


@pytest.mark.parametrize('param', ['fpatchy', 'qcap'])
def test_lower_equal(tmp_path, capfd, param):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        reset={param:'1.1'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert "({:s}) must be <= ".format(param) in captured.out


def test_tcea_missing_mass_units(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/pt_tcea.cfg',
        reset={'mplanet':'1.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Invalid units 'None' for parameter mplanet." in captured.out

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# pt runmode fails:
@pytest.mark.parametrize('param', ['nlayers', 'ptop', 'pbottom'])
def test_pt_pressure_missing(tmp_path, capfd, undefined, param):
    cfg = make_config(tmp_path, ROOT+'tests/pt_isothermal.cfg',
        remove=[param])
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'driver.py', function: 'check_pressure'" \
           in captured.out
    assert undefined[param] in captured.out


# This is valid for any get_param() input:
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
def test_pt_tpars_mismatch(tmp_path, capfd, cfile, error):
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

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# atmosphere runmode fails:

@pytest.mark.parametrize('param',
    ['atmfile', 'species'])
def test_uniform_missing(tmp_path, capfd, param, undefined):
    cfg = make_config(tmp_path, ROOT+'tests/atmosphere_uniform_test.cfg',
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_atm'" \
           in captured.out
    assert undefined[param] in captured.out


def test_uniform_uniform_mismatch(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/atmosphere_uniform_test.cfg',
        reset={'uniform':'0.85 0.15'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_atm'" in captured.out
    assert "Number of uniform abundances (2) does not match the number " \
           "of species (7)." in captured.out


@pytest.mark.parametrize('param',
    ['atmfile', 'species', 'elements'])
def test_tea_missing(tmp_path, capfd, param, undefined):
    cfg = make_config(tmp_path, ROOT+'tests/atmosphere_tea_test.cfg',
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_atm'" in captured.out
    assert undefined[param] in captured.out

