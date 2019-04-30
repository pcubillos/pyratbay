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
    ['runmode', 'rayleigh', 'hazes', 'alkali', 'path',
     'tmodel', 'molmodel', 'retflag'])
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

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# spectrum runmode fails (setup):

@pytest.mark.parametrize('param',
    ['wllow', 'wlhigh', 'wnstep', 'wnosamp'])
def test_spectrum_missing(tmp_path, capfd, param, undefined_spec):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'makesample.py', function: 'make_wavenumber'" \
           in captured.out
    assert undefined_spec[param] in captured.out


def test_spectrum_inconsistent_wl_bounds(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'wllow':'2.0 um', 'wlhigh':'1.0 um'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'makesample.py', function: 'make_wavenumber'" \
           in captured.out
    assert 'Wavenumber low boundary (10000.0 cm-1) must be larger than the ' \
           'high boundary\n(5000.0 cm-1).' in captured.out


@pytest.mark.parametrize('param',
    ['rstar', 'path', 'outspec'])
def test_spectrum_transmission_missing(tmp_path, capfd, param, undefined_spec):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert undefined_spec[param] in captured.out


@pytest.mark.parametrize('param',
    ['mplanet', 'rplanet', 'gplanet'])
def test_spectrum_hydrostatic_equilibrium(tmp_path, capfd, param):
    keep = ['mplanet', 'rplanet', 'gplanet']
    keep.remove(param)
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=keep)
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'makesample.py', function: 'make_atmprofiles'" \
           in captured.out
    assert 'Cannot compute hydrostatic equilibrium.  Must define ' \
           'at least two of\nmplanet, rplanet, or gplanet.' in captured.out


def test_spectrum_refpressure(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=['refpressure'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'makesample.py', function: 'make_atmprofiles'" \
           in captured.out
    assert 'Cannot compute hydrostatic equilibrium.  Undefined reference ' \
           'pressure level\n(refpressure).' in captured.out


def test_spectrum_inconsistent_mass_radius_gravity(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'gplanet':'1400.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'run'" in captured.out
    assert 'All mplanet, rplanet, and gplanet were provided, but values ' \
           'are inconsistent\n(>5%): g(M,R) =  1487.2 cm s-2 and ' \
           'gplanet =  1400.0 cm s-2.'  in captured.out


@pytest.mark.parametrize('param',
    ['tlifile',])
def test_spectrum_invalid_file(tmp_path, capfd, param, invalid_file):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={param:'nope.dat'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert invalid_file[param] in captured.out


@pytest.mark.parametrize('vmin,vmax',
    [('Dmin', 'Dmax'),
     ('Lmin', 'Lmax')])
def test_spectrum_inconsistent_voigt_bounds(tmp_path, capfd, vmin, vmax):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={vmin:'1e5', vmax:'1e4'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert '{:s} (10000 cm-1) must be > {:s} (100000 cm-1).'. \
           format(vmax,vmin) in captured.out


def test_spectrum_rpars_mismatch(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'rpars':'1.0 1.0 1.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert 'Number of input Rayleigh parameters (3) does not match the ' \
           'number of\nrequired model parameters (2).' in captured.out


def test_spectrum_hpars_mismatch(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'hpars':'1.0 1.0 1.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert 'Number of input haze parameters (3) does not match the number ' \
           'of required\nmodel parameters (1).' in captured.out


@pytest.mark.parametrize('value',
    ['10 60 90', '0 30 60 100', '0 30 90 60'])
def test_spectrum_raygrid(tmp_path, capfd, invalid_raygrid, value):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'raygrid':value})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert invalid_raygrid[value] in captured.out


def test_spectrum_uncert_mismatch(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'data':'1.0 2.0', 'uncert':'0.1'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert 'Number of data uncertainty values (1) does not match the ' \
           'number of data\npoints (2).' in captured.out


def test_spectrum_filter_mismatch(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'data':'1.0 2.0',
               'filter':ROOT+'tests/filters/filter_test_WFC3_G141_1.133um.dat'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert 'Number of filter bands (1) does not match the number of ' \
           'data points (2).' in captured.out


@pytest.mark.parametrize('param', ['rstar', 'tstar', 'smaxis'])
def test_spectrum_tcea_parameters(tmp_path, capfd, param, undefined):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=[param],
        reset={'path':'eclipse', 'tmodel':'TCEA'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert undefined[param] in captured.out


@pytest.mark.parametrize('param', ['rplanet', 'mplanet'])
def test_spectrum_tcea_gplanet(tmp_path, capfd, param, undefined):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        remove=[param, 'gplanet'],
        reset={'path':'eclipse', 'tmodel':'TCEA'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert undefined['gplanet'] in captured.out


def test_bulk_not_in_atm(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'bulk':'N2'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert "These bulk species are not present in the atmosphere: ['N2']." in captured.out


def test_molfree_mismatch(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'molmodel':'vert'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert "molmodel is set, but there are no molfree." in captured.out


def test_molfree_mismatch2(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'molmodel':'vert vert', 'molfree':'H2O'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert "There should be one molfree for each molmodel:" in captured.out
    assert "molmodel: ['vert', 'vert']" in captured.out
    assert "molfree: ['H2O']" in captured.out


def test_molfree_mismatch3(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'molmodel':'vert', 'molfree':'N2'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert "These species are not present in the atmosphere: ['N2']." in captured.out


def test_bulk_molfree_overlap(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'molmodel':'vert', 'molfree':'H2', 'bulk':'H2'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert "These species were marked as both bulk and variable-abundance: ['H2']." in captured.out


@pytest.mark.parametrize('param', ['tstar', 'gstar'])
def test_kurucz_missing_pars(tmp_path, capfd, param, undefined):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'kurucz':'fp00k0odfnew.pck'},
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert undefined[param] in captured.out


def test_spectrum_opacity_tmin(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'extfile':str(tmp_path/'new_opacity.dat'),
               'tmin':'10.0', 'tmax':'1000.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'extinction.py', function: 'calc_extinction'" \
           in captured.out
    assert ("Requested extinction-coefficient table temperature "
            "(tmin=10.0 K) below the\nlowest available TLI temperature "
            "(70.0 K).") in captured.out


def test_spectrum_opacity_tmax(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/spectrum_transmission_test.cfg',
        reset={'extfile':str(tmp_path/'new_opacity.dat'),
               'tmin':'1000.0', 'tmax':'5000.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'extinction.py', function: 'calc_extinction'" \
           in captured.out
    assert ("Requested extinction-coefficient table temperature "
            "(tmax=5000.0 K) above the\nhighest available TLI temperature "
            "(3000.0 K).") in captured.out


@pytest.mark.skip
def test_molecule_not_in_molfile():
    pass

@pytest.mark.skip
def test_incompatible_tli():
    pass

@pytest.mark.skip
def test_unbounded_ptop_pbottom():
    pass

@pytest.mark.skip
def test_crosssec_mol_not_in_atm():
    pass

