# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

import os
import subprocess
import pytest

from conftest import make_config

import pyratbay as pb
import pyratbay.io as io
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check runmode and command_line vs interpreter runs:
@pytest.mark.parametrize('runmode', ['None', 'invalid'])
@pytest.mark.parametrize('call',    ['command_line', 'interpreter'])
def test_run_runmode(tmp_path, capfd, runmode, call):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'runmode':runmode})
    if call == 'interpreter':
        pyrat = pb.run(cfg)
        assert pyrat is None
    elif call == 'command_line':
        subprocess.call(f'pbay -c {cfg}'.split())
    captured = capfd.readouterr()
    caps = ["Error in module: 'parser.py', function: 'parse'",
           f"Invalid running mode (runmode): {runmode}. Select from: "
            "['tli', 'atmosphere',"]
    for cap in caps:
        assert cap in captured.out


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check output for each case is defined:
def test_run_tli(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/tli_hitran_1.1-1.7um_test.cfg',
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
    cfg = make_config(tmp_path, f'{ROOT}tests/configs/{cfile}',
        remove=['extfile'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'run'",
            "Undefined extinction-coefficient file (extfile)."]
    for cap in caps:
        assert cap in captured.out


def test_run_mcmc_mcmcfile(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/mcmc_transmission_test.cfg',
        remove=['mcmcfile'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    caps = ["Error in module: 'driver.py', function: 'run'",
            "Undefined MCMC file (mcmcfile)."]
    for cap in caps:
        assert cap in captured.out


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check input units:
@pytest.mark.parametrize('param, var',
    [('wlunits', 'wavelength'),
     ('runits', 'radius'),
     ('punits', 'pressure'),
     ('dunits', 'data')])
def test_invalid_units(tmp_path, capfd, param, var):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'invalid'})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert "Invalid {:s} units ({:s}): invalid".format(var, param) \
           in captured.out


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check integer and float data types:
@pytest.mark.parametrize('param, value',
    [('nlayers', '10.5'),
     ('nlayers', '10 20'),
     ('nlayers', 'a')])
def test_invalid_integer_type(tmp_path, capfd, param, value):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:value})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert ("Invalid data type for {:s}, could not convert string to "
            "integer: '{:s}'".format(param, value)) in captured.out


@pytest.mark.parametrize('param',
    ['verb', 'wnosamp', 'nlayers', 'ncpu', 'ndop', 'nlor', 'quadrature',
     'nsamples', 'nchains', 'burnin', 'thinning', 'resume'])
def test_invalid_integer_all_params(tmp_path, capfd, param):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
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
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_tcea.cfg',
        reset={param:value})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert ("Invalid data type for {:s}, could not convert string to "
            "float: '{:s}'".format(param, value)) in captured.out


@pytest.mark.parametrize('param',
    ['wnlow', 'wnhigh', 'wnstep', 'resolution', 'xsolar', 'tmin', 'tmax',
     'tstep', 'ethresh', 'vextent', 'dmin', 'dmax', 'lmin', 'lmax', 'dlratio',
     'fpatchy', 'maxdepth', 'qcap', 'tlow', 'thigh', 'grbreak', 'grnmin',
     'gstar', 'tstar', 'gplanet', 'tint'])
def test_invalid_float_all_params(tmp_path, capfd, param):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'abc'})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert ("Invalid data type for {:s}, could not convert string to "
            "float: 'abc'".format(param)) in captured.out


@pytest.mark.parametrize('param',
    ['runmode', 'rayleigh', 'clouds', 'alkali', 'rt_path',
     'tmodel', 'molmodel', 'retflag'])
def test_invalid_choice(tmp_path, capfd, param, invalid):
    cfg = make_config(
        tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'invalid'})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert invalid[param] in captured.out


@pytest.mark.parametrize('param',
    ['starspec', 'kurucz', 'marcs', 'phoenix', 'filters',
     'dblist', 'molfile', 'csfile'])
def test_file_not_found(tmp_path, capfd, param, invalid_file):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'nope.dat'})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert invalid_file[param] in captured.out


@pytest.mark.parametrize('param',
    ['atmfile', 'tlifile', 'extfile', 'mcmcfile', 'specfile', 'ptfile',
     'logfile'])
def test_invalid_file_path(tmp_path, capfd, param, invalid_path):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'nope/file.dat'})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert invalid_path[param] in captured.out
    # params from 'test_file_not_found' do not raise folder error since
    # they catch file not found first.


def test_missing_mass_units(tmp_path, capfd):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_tcea.cfg',
        reset={'mplanet':'1.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Invalid units 'None' for parameter mplanet." in captured.out


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
     ('smaxis',  ' -0.01 au'),
     ('rstar',   ' -1.0 rsun'),
     ('mstar',   ' -1.0 msun'),
     ('gstar',   ' -1000.0'),
     ('tstar',   ' -5000.0'),
     ('dmin',    ' -1e-6'),
     ('dmax',    ' -1e-1'),
     ('lmin',    ' -1e-6'),
     ('lmax',    ' -1e-1'),
     ('dlratio', ' -0.1'),
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
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_tcea.cfg',
        reset={param:value})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert "({:s}) must be > ".format(param) in captured.out


@pytest.mark.parametrize('param',
    ['wnosamp', 'ndop', 'nlor', 'thinning', 'nchains', 'ncpu', 'tint',
     'quadrature', 'grbreak', 'fpatchy', 'maxdepth', 'vextent'])
def test_greater_equal(tmp_path, capfd, param):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'-10'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert "({:s}) must be >= ".format(param) in captured.out


@pytest.mark.parametrize('param', ['verb'])
def test_lower_than(tmp_path, capfd, param):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'10'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert "({:s}) must be < ".format(param) in captured.out


@pytest.mark.parametrize('param', ['fpatchy', 'qcap'])
def test_lower_equal(tmp_path, capfd, param):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'1.1'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert "({:s}) must be <= ".format(param) in captured.out


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# atmosphere (temperature profiles) runmode fails:
@pytest.mark.parametrize('param', ['nlayers', 'ptop', 'pbottom'])
def test_pt_pressure_missing(tmp_path, capfd, undefined, param):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        remove=[param])
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'driver.py', function: 'check_pressure'" \
           in captured.out
    assert undefined[param] in captured.out


# This is valid for any get_param() input:
@pytest.mark.parametrize('value', ['a', '10.0 bar 30.0'])
@pytest.mark.parametrize('param', ['ptop', 'pbottom'])
def test_pressure_invalid_type(tmp_path, capfd, param, value):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:value})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert "Invalid value '{:s}' for parameter {:s}.". \
           format(value, param) in captured.out


@pytest.mark.parametrize('param', ['ptop', 'pbottom'])
def test_pressure_invalid_units(tmp_path, capfd, param):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'10.0 20.0'})
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()
    assert pyrat is None
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert f"Invalid units for value '10.0 20.0' for parameter {param}." \
           in captured.out


@pytest.mark.parametrize('param', ['tmodel', 'tpars'])
def test_pt_temperature_missing(tmp_path, capfd, param, undefined):
    cfg = make_config(tmp_path, ROOT+'tests/configs/pt_isothermal.cfg',
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_temp'" \
           in captured.out
    assert undefined[param] in captured.out


@pytest.mark.parametrize('tmodel, npars',
    [('isothermal', 1),
     ('tcea', 6),
     ('madhu', 6)])
def test_pt_tpars_mismatch(tmp_path, capfd, tmodel, npars):
    cfg = make_config(tmp_path, f'{ROOT}tests/configs/pt_{tmodel}.cfg',
        reset={'tpars':'100.0 200.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'atmosphere.py', function: 'temperature'" \
           in captured.out
    assert f"Wrong number of parameters (2) for the {tmodel} temperature " \
           f"model ({npars})" in captured.out


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# atmosphere (chemistry) runmode fails:

def test_missing_atmfile(tmp_path, capfd, undefined):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_uniform_test.cfg',
        remove=['atmfile'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_atm'" in captured.out
    assert undefined['atmfile'] in captured.out


@pytest.mark.parametrize('chem', ['uniform', 'tea'])
def test_atmosphere_missing_species(tmp_path, capfd, undefined, chem):
    cfg = make_config(tmp_path,
        f'{ROOT}tests/configs/atmosphere_{chem}_test.cfg',
        remove=['species'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_atm'" in captured.out
    assert undefined['species'] in captured.out


def test_atmosphere_uniform_missing_uniform(tmp_path, capfd, undefined):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_uniform_test.cfg',
        remove=['uniform'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_atm'" in captured.out
    assert undefined['uniform'] in captured.out


def test_atmosphere_tea_missing_elements(tmp_path, capfd, undefined):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_tea_test.cfg',
        remove=['elements'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_atm'" in captured.out
    assert undefined['elements'] in captured.out


def test_atmosphere_uniform_mismatch_uniform(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_uniform_test.cfg',
        reset={'uniform':'0.85 0.15'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_atm'" in captured.out
    assert "Number of uniform abundances (2) does not match the number " \
           "of species (7)." in captured.out


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# atmosphere (altitude) runmode fails:

def test_atmosphere_hydro_missing_refpressure(tmp_path, capfd, undefined):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_hydro_test.cfg',
        remove=['refpressure'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_altitude'" \
        in captured.out
    assert undefined['refpressure'] in captured.out


def test_atmosphere_hydro_missing_all_planet_props(tmp_path, capfd, undefined):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_hydro_test.cfg',
        remove=['mplanet', 'rplanet', 'gplanet'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_altitude'" \
        in captured.out
    assert ("Cannot compute hydrostatic-equilibrium radius profile.  Must\n"
            "    define at least two of mplanet, rplanet, or gplanet.") \
        in captured.out

@pytest.mark.parametrize('param', ['mplanet', 'rplanet', 'gplanet'])
def test_atmosphere_hydro_missing_two_props(tmp_path, capfd, param):
    params = {
        'mplanet': '1.0 mjup',
        'rplanet': '1.0 rjup',
        'gplanet': '2479.0'
    }
    missing = list(params)
    missing.remove(param)
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_hydro_test.cfg',
        remove=missing,
        reset={param:params[param]})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_altitude'" \
        in captured.out
    assert ("Cannot compute hydrostatic-equilibrium radius profile.  Must\n"
            f"    define either {missing[0]} or {missing[1]}.") in captured.out


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# spectrum runmode fails (setup):

@pytest.mark.parametrize('param',
    ['wllow', 'wlhigh',
     'pbottom', 'ptop', 'refpressure',
     'mstar', 'rstar', 'smaxis',
    ])
def test_spectrum_missing_units(tmp_path, capfd, param):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={param:'1.1'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'parser.py', function: 'parse'" in captured.out
    assert f"Invalid units 'None' for parameter {param}." in captured.out


def test_spectrum_inconsistent_wl_bounds(tmp_path, capfd):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'wllow':'2.0 um', 'wlhigh':'1.0 um'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'makesample.py', function: 'make_wavenumber'" \
           in captured.out
    assert 'Wavenumber low boundary (10000.0 cm-1) must be larger than the ' \
           'high boundary\n(5000.0 cm-1).' in captured.out


@pytest.mark.parametrize('param',
    ['rstar', 'rt_path', 'specfile'])
def test_spectrum_transmission_missing(tmp_path, capfd, param, undefined_spec):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert undefined_spec[param] in captured.out


def test_spectrum_missing_chemistry_new_atmfile(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'atmfile':'{ROOT}tests/inputs/atmosphere_new_test.atm',
               'ptop':'1e-6 bar',
               'pbottom':'100.0 bar',
               'nlayers':'81',
               'tmodel':'tcea',
               'tpars':'-4.84 -0.8 -0.8 0.5 1200.0 100.0',
            })
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_atm'" \
           in captured.out
    assert 'Undefined chemistry model (chemistry).' in captured.out


def test_spectrum_missing_chemistry_no_atmfile(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['atmfile'],
        reset={'ptop':'1e-6 bar',
               'pbottom':'100.0 bar',
               'nlayers':'81',
               'tmodel':'tcea',
               'tpars':'-4.84 -0.8 -0.8 0.5 1200.0 100.0',
            })
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'check_atm'" \
           in captured.out
    assert 'Undefined chemistry model (chemistry).' in captured.out


def test_spectrum_no_radius(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        remove=['radmodel'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'makesample.py', function: 'make_atmprofiles'" \
           in captured.out
    assert 'Cannot compute radius profile.  Need to set a radius model ' \
           '(radmodel) or\nprovide an input radius array in the ' \
           'atmospheric file.' in captured.out


@pytest.mark.parametrize('param',
    ['mplanet', 'rplanet', 'gplanet'])
@pytest.mark.parametrize('atm',
    [f'{ROOT}/tests/inputs/atmosphere_uniform_test.atm',
     f'{ROOT}/tests/inputs/atmosphere_uniform_radius.atm'])
def test_spectrum_hydro_MRGplanet(tmp_path, capfd, param, atm):
    keep = ['mplanet', 'rplanet', 'gplanet']
    keep.remove(param)
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'atmfile':atm},
        remove=keep)
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'makesample.py', function: 'make_atmprofiles'" \
           in captured.out
    assert 'Cannot compute hydrostatic-equilibrium radius profile.  Must ' \
           'define at least\ntwo of mplanet, rplanet, or gplanet.' \
           in captured.out


@pytest.mark.parametrize('atm',
    [f'{ROOT}/tests/inputs/atmosphere_uniform_test.atm',
     f'{ROOT}/tests/inputs/atmosphere_uniform_radius.atm'])
def test_spectrum_hydro_refpressure(tmp_path, capfd, atm):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'atmfile':atm},
        remove=['refpressure'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'makesample.py', function: 'make_atmprofiles'" \
           in captured.out
    assert 'Cannot compute hydrostatic-equilibrium radius profile.  ' \
           'Undefined reference\npressure level (refpressure).' in captured.out


@pytest.mark.parametrize('value', ['1.00e-09 bar', '1.00e+03 bar'])
@pytest.mark.parametrize('param', ['pbottom', 'ptop'])
def test_spectrum_unbounded_pressures(tmp_path, capfd, param, value):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={param:value})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'makesample.py', function: 'make_atmprofiles'" \
           in captured.out
    assert ('{}-pressure boundary ({}={}) lies outside of the\n'
            'atmospheric-file range 1.00e-06--1.00e+02 bar.'.
            format(param[1:].capitalize(), param, value)) in captured.out


def test_spectrum_invalid_pressure_ranges(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'ptop':'1.0e-02 bar', 'pbottom':'1.0e-03 bar'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'makesample.py', function: 'make_atmprofiles'" \
           in captured.out
    assert ('Bottom-layer pressure (1.00e-03 bar) must be higher than the '
            'top-layer\npressure (1.00e-02 bar).') in captured.out


def test_spectrum_inconsistent_mass_radius_gravity(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'gplanet':'1400.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'driver.py', function: 'run'" in captured.out
    assert 'All mplanet, rplanet, and gplanet were provided, but values ' \
           'are inconsistent\n(>5%): g(M,R) =  1487.3 cm s-2 and ' \
           'gplanet =  1400.0 cm s-2.'  in captured.out


@pytest.mark.parametrize('param',
    ['tlifile',])
def test_spectrum_invalid_file(tmp_path, capfd, param, invalid_file):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={param:'nope.dat'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert invalid_file[param] in captured.out


@pytest.mark.parametrize('vmin,vmax',
    [('dmin', 'dmax'),
     ('lmin', 'lmax')])
def test_spectrum_inconsistent_voigt_bounds(tmp_path, capfd, vmin, vmax):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={vmin:'1e5', vmax:'1e4'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert '{:s} (10000 cm-1) must be > {:s} (100000 cm-1).'. \
           format(vmax,vmin) in captured.out


def test_spectrum_rpars_mismatch(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'rpars':'1.0 1.0 1.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert 'Number of input Rayleigh parameters (3) does not match the ' \
           'number of\nrequired model parameters (2).' in captured.out


def test_spectrum_cpars_mismatch(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'cpars':'1.0 1.0 1.0'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert 'Number of input cloud parameters (3) does not match the number ' \
           'of required\nmodel parameters (1).' in captured.out


@pytest.mark.parametrize('value',
    ['10 60 90', '0 30 60 100', '0 30 90 60'])
def test_spectrum_raygrid(tmp_path, capfd, invalid_raygrid, value):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'raygrid':value})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert invalid_raygrid[value] in captured.out


def test_spectrum_uncert_mismatch(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'data':'1.0 2.0', 'uncert':'0.1'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert 'Number of data uncertainty values (1) does not match the ' \
           'number of data\npoints (2).' in captured.out


def test_spectrum_filters_mismatch(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'data':'1.0 2.0',
            'filters':ROOT+'tests/filters/filter_test_WFC3_G141_1.133um.dat'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert 'Number of filter bands (1) does not match the number of ' \
           'data points (2).' in captured.out


def test_spectrum_params_misfit(tmp_path, capfd):
    # Without evaulating params:
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'tmodel':'tcea',
               'retflag':'temp',
               'params':'-4.67 -0.8'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert "The number of input fitting parameters (params, 2) does not " \
           "match\n    the number of required parameters (6)." in captured.out


def test_eval_params_misfit(tmp_path, capfd):
    # Without evaulating params:
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'tmodel':'tcea',
               'retflag':'temp'})
    pyrat = pb.run(cfg)
    pyrat.eval([-4.67, -0.8])
    captured = capfd.readouterr()
    assert "The number of input fitting parameters (2) does not " \
           "match\n    the number of required parameters (6)." in captured.out


def test_bulk_not_in_atm(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'bulk':'N2'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert "These bulk species are not present in the atmosphere: ['N2']." \
           in captured.out


def test_molfree_mismatch(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'molmodel':'vert'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert "molmodel is set, but there are no molfree." in captured.out


def test_molfree_mismatch2(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'molmodel':'vert vert', 'molfree':'H2O'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert "There should be one molfree for each molmodel:" in captured.out
    assert "molmodel: ['vert', 'vert']" in captured.out
    assert "molfree: ['H2O']" in captured.out


def test_molfree_mismatch3(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'molmodel':'vert', 'molfree':'N2'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert "These species are not present in the atmosphere: ['N2']." in captured.out


def test_bulk_molfree_overlap(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'molmodel':'vert', 'molfree':'H2', 'bulk':'H2'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert "These species were marked as both bulk and variable-abundance: ['H2']." in captured.out


@pytest.mark.parametrize('param', ['tstar', 'gstar'])
def test_kurucz_missing_pars(tmp_path, capfd, param, undefined):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'kurucz':f'{ROOT}/tests/inputs/fp00k0odfnew.pck'},
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert undefined[param] in captured.out


def test_spectrum_opacity_invalid_tmin(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'extfile':str(tmp_path/'new_opacity.dat'),
               'tmin':'0.1', 'tmax':'1000.0', 'tstep':'900'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'extinction.py', function: 'calc_extinction'" \
           in captured.out
    assert ("Requested extinction-coefficient table temperature "
            "(tmin=0.1 K) below the\nlowest available TLI temperature "
            "(1.0 K).") in captured.out


@pytest.mark.parametrize('param',
    ['tmodel', 'clouds', 'rayleigh', 'molmodel', 'bulk'])
def test_spectrum_missing_retflag_models(tmp_path, capfd, param,undefined_mcmc):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'retflag':'temp mol ray cloud',
               'tmodel':'isothermal',
               'clouds':'deck',
               'rayleigh':'lecavelier',
               'molmodel':'vert', 'molfree':'H2O', 'bulk':'H2'},
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert undefined_mcmc[param] in captured.out


def test_spectrum_opacity_invalid_tmax(tmp_path, capfd):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'extfile':str(tmp_path/'new_opacity.dat'),
               'tmin':'1000.0', 'tmax':'6000.0', 'tstep':'100'})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'extinction.py', function: 'calc_extinction'" \
           in captured.out
    assert ("Requested extinction-coefficient table temperature "
            "(tmax=6000.0 K) above the\nhighest available TLI temperature "
            "(5000.0 K).") in captured.out


@pytest.mark.parametrize('param', ['tmin', 'tmax', 'tstep', 'tlifile'])
def test_spectrum_opacity_missing(tmp_path, capfd, param, undefined_opacity):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'extfile':str(tmp_path/'new_opacity.dat'),
               'tmin':'300.0', 'tmax':'3000.0', 'tstep':'900'},
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert undefined_opacity[param] in captured.out


def test_molecule_not_in_molfile(tmp_path, capfd):
    # Modify atm:
    units, species, press, temp, q, rad = \
        io.read_atm(ROOT+'tests/inputs/atmosphere_uniform_test.atm')
    press = press * pb.tools.u(units[0])
    species[-1] = 'X'
    new_atm = str(tmp_path/'new_atmosphere_uniform_test.atm')
    io.write_atm(new_atm, press, temp, species, q, punits=units[0])

    cfg = make_config(tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={'atmfile':new_atm})
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'read_atm.py', function: 'get_constants'" \
           in captured.out
    assert "These species: ['X'] are not listed in the molecules info file" \
           in captured.out


@pytest.mark.skip
def test_incompatible_tli():
    pass

@pytest.mark.skip
def test_crosssec_mol_not_in_atm():
    pass


@pytest.mark.parametrize('param', ['tmin', 'tmax', 'tstep', 'tlifile'])
def test_opacity_missing(tmp_path, capfd, param, undefined_opacity):
    cfg = make_config(tmp_path, ROOT+'tests/configs/opacity_test.cfg',
        reset={'extfile':str(tmp_path/'new_opacity.dat')},
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'check_spectrum'" \
           in captured.out
    assert undefined_opacity[param] in captured.out


@pytest.mark.parametrize('param',
    ['retflag', 'params', 'data', 'uncert', 'filters', 'rstar',
     'sampler', 'nsamples', 'burnin', 'nchains'])
def test_mcmc_missing(tmp_path, capfd, param, undefined_mcmc):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/mcmc_transmission_test.cfg',
        reset={
            'rt_path': 'emission',
            'kurucz': f'{ROOT}tests/inputs/fp00k0odfnew.pck'},
        remove=[param])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert \
        "Error in module: 'argum.py', function: 'check_spectrum'" \
        in captured.out
    assert undefined_mcmc[param] in captured.out


def test_mcmc_missing_starspec(tmp_path, capfd):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/mcmc_transmission_test.cfg',
        reset={'rt_path':'emission', 'retflag':'mol'},
        remove=['tstar', 'tmodel'])
    pyrat = pb.run(cfg)
    assert pyrat is None
    captured = capfd.readouterr()
    assert "Error in module: 'argum.py', function: 'setup'" in captured.out
    assert ('Undefined stellar flux model.  Set starspec, kurucz, or tstar '
            '(for a\nblackbody spectrum).') in captured.out


@pytest.mark.parametrize('param', ['tlifile', 'csfile', 'extfile'])
def test_spectrum_temperature_bounds(tmp_path, capfd, param, invalid_temp):
    remove = [par for par in ['tlifile', 'csfile', 'extfile'] if par != param]
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={
            'tmodel': 'isothermal',
            'tpars': '6000.0',
            'extfile': f'{ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz'},
        remove=remove)
    pyrat = pb.run(cfg)
    assert pyrat is not None
    captured = capfd.readouterr()
    assert invalid_temp[param] in captured.out
    assert pyrat.spec.spectrum is None

