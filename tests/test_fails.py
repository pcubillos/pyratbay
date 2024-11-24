# Copyright (c) 2021-2023 Patricio Cubillos
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import subprocess
import pytest
import re

from conftest import make_config

import pyratbay as pb
import pyratbay.io as io
from pyratbay.constants import ROOT

os.chdir(ROOT+'tests')

check_spectrum_err = "Error in module: 'argum.py', function: 'check_spectrum'"
setup_err = "Error in module: 'argum.py', function: 'setup'"


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check runmode and command_line vs interpreter runs:
@pytest.mark.parametrize('runmode', ['None', 'invalid'])
@pytest.mark.parametrize('call', ['command_line', 'interpreter'])
def test_run_runmode(tmp_path, capfd, runmode, call):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'runmode':runmode},
    )
    error = (
        f"Invalid running mode (runmode): '{runmode}'. Select from: "
    )
    if call == 'interpreter':
        with pytest.raises(ValueError, match=re.escape(error)):
            pyrat = pb.run(cfg)
    elif call == 'command_line':
        subprocess.call(f'pbay -c {cfg}'.split())
        captured = capfd.readouterr()
        assert error in captured.err


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check output for each case is defined:
def test_run_mcmc_mcmcfile(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/mcmc_transmission_test.cfg',
        remove=['mcmcfile'],
        reset={'runmode': 'spectrum'},
    )
    error = "Undefined retrieval file (mcmcfile)"
    with pytest.raises(ValueError, match=re.escape(error)):
        pyrat = pb.Pyrat(cfg)
        pyrat.retrieval()


@pytest.mark.skip(reason='TBD: fix input atmosphere, add radius')
def test_run_opacity_extfile_missing(tmp_path):
    reset = {
        'runmode': 'spectrum',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/opacity_test.cfg',
        reset=reset,
        remove=['extfile'],
    )
    error = re.escape(
        'Undefined output opacity file (extfile) needed to compute '
        'opacity table'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.Pyrat(cfg)
        pyrat.compute_opacity()


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check input units:
@pytest.mark.parametrize('param, var',
    [('wlunits', 'wavelength'),
     ('runits', 'radius'),
     ('punits', 'pressure'),
     ('dunits', 'data')])
def test_invalid_units(tmp_path, param, var):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'invalid'},
    )
    error = f"Invalid {var} units ({param}): invalid"
    with pytest.raises(ValueError, match=re.escape(error)):
        pyrat = pb.run(cfg)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check integer and float data types:
@pytest.mark.parametrize('param, value',
    [('nlayers', '10.5'),
     ('nlayers', '10 20'),
     ('nlayers', 'a')])
def test_invalid_integer_type(tmp_path, param, value):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:value},
    )
    error = (
        f"Invalid data type for {param}, could not convert string to "
        f"integer: '{value}'")
    with pytest.raises(ValueError, match=re.escape(error)):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('param',
    ['verb', 'wnosamp', 'nlayers', 'ndop', 'nlor', 'quadrature',
     'nsamples', 'nchains', 'burnin', 'thinning'])
def test_invalid_integer_all_params(tmp_path, param):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'abc'},
    )
    error = (
        f"Invalid data type for {param}, could not convert string to "
        "integer: 'abc'")
    with pytest.raises(ValueError, match=re.escape(error)):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('param, value',
    [('tstar', '100 200'),
     ('tstar', 'a')])
def test_invalid_float_type(tmp_path, param, value):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_guillot.cfg',
        reset={param:value},
    )
    error = (
        f"Invalid data type for {param}, could not convert string to "
        f"float: '{value}'")
    with pytest.raises(ValueError, match=re.escape(error)):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('param',
    ['wnlow', 'wnhigh', 'wnstep', 'resolution', 'xsolar', 'tmin', 'tmax',
     'tstep', 'ethresh', 'vextent', 'dmin', 'dmax', 'lmin', 'lmax', 'dlratio',
     'fpatchy', 'maxdepth', 'qcap', 'tlow', 'thigh', 'grbreak', 'grnmin',
     'log_gstar', 'tstar', 'gplanet', 'tint'])
def test_invalid_float_all_params(tmp_path, param):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'abc'},
    )
    error = (
        f"Invalid data type for {param}, could not convert string to "
        "float: 'abc'")
    with pytest.raises(ValueError, match=re.escape(error)):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize(
    'param, reason',
    [
        ('runmode',  "running mode"),
        ('rayleigh', "Rayleigh model"),
        ('clouds',   "cloud model"),
        ('alkali',   "alkali model"),
        ('rt_path',  "radiative-transfer observing geometry"),
        ('tmodel',   "temperature model"),
        ('retflag',  "retrieval flag"),
        ('sampler',  "posterior sampler"),
    ]
)
def test_invalid_choice(tmp_path, param, reason):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'invalid'},
    )
    error = (
        re.escape(f"Invalid {reason} ({param}): 'invalid'.")
        + " Select from: *"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('param',
    ['starspec', 'kurucz', 'marcs', 'phoenix', 'filters',
     'dblist', 'molfile', 'csfile'])
def test_file_not_found(tmp_path, param, invalid_file):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'nope.dat'},
    )
    error = re.escape(invalid_file[param]) + '*.'
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_invalid_logfile_path(tmp_path, invalid_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={'logfile':'ok/nope/file.dat'},
    )
    error = re.escape('No such file or directory')
    with pytest.raises(FileNotFoundError, match=error):
        pyrat = pb.run(cfg)



@pytest.mark.parametrize(
    'param',
    [
        'output_atmfile',
        'tlifile',
        'extfile',
        'mcmcfile',
        'specfile',
        'ptfile',
    ],
)
def test_invalid_file_path(tmp_path, param, invalid_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'nope/file.dat'},
    )
    error = re.escape(invalid_path[param])
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)
    # params from 'test_file_not_found' do not raise folder error since
    # they catch file not found first.


def test_missing_mass_units(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_guillot.cfg',
        reset={'mplanet':'1.0'},
    )
    error = "Invalid units 'None' for parameter mplanet"
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


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
     ('pbottom', ' -1.0 bar'),
     ('mplanet', ' -1.0 mjup'),
     ('rplanet', ' -1.0 rjup'),
     ('gplanet', ' -1000.0'),
     ('smaxis',  ' -0.01 au'),
     ('rstar',   ' -1.0 rsun'),
     ('mstar',   ' -1.0 msun'),
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
def test_greater_than(tmp_path, param, value):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_guillot.cfg',
        reset={param:value},
    )
    error = re.escape(f"({param}) must be > ") + '*.'
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('param',
    ['wnosamp', 'ndop', 'nlor', 'thinning', 'nchains', 'ncpu', 'tint',
     'quadrature', 'grbreak', 'fpatchy', 'maxdepth', 'vextent'])
def test_greater_equal(tmp_path, param):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'-10'},
    )
    error = re.escape(f"({param}) must be >= ") + '*.'
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('param', ['verb'])
def test_lower_than(tmp_path, param):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'10'},
    )
    error = re.escape(f"({param}) must be < ") + '*.'
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('param', ['fpatchy', 'qcap'])
def test_lower_equal(tmp_path, param):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'1.1'},
    )
    error = re.escape(f"({param}) must be <= ") + '*.'
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# atmosphere (temperature profiles) runmode fails:
@pytest.mark.parametrize('param', ['nlayers', 'ptop', 'pbottom'])
def test_pt_pressure_missing(tmp_path, param):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        remove=['atmfile', param],
    )
    error = re.escape(
        'Cannot compute pressure profile, either set {ptop, pbottom, '
        'nlayers} parameters, or provide an input PT profile (ptfile) '
        'or atmospheric file (atmfile)'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


# This is valid for any get_param() input:
@pytest.mark.parametrize('value', ['a', '10.0 bar 30.0'])
@pytest.mark.parametrize('param', ['ptop', 'pbottom'])
def test_pressure_invalid_type(tmp_path, param, value):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:value},
    )
    error = f"Invalid value '{value}' for parameter {param}"
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('param', ['ptop', 'pbottom'])
def test_pressure_invalid_units(tmp_path, param):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        reset={param:'10.0 20.0'},
    )
    error = f"Invalid units for value '10.0 20.0' for parameter {param}"
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('param', ['tmodel', 'tpars'])
def test_pt_temperature_missing(tmp_path, param):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/pt_isothermal.cfg',
        remove=['atmfile', 'tmodel'],
    )
    error = re.escape(
        'Cannot compute temperature profile, either set a temperature model '
        '(tmodelname) and parameters (tpars), or provide an input PT '
        'profile (ptfile) or atmospheric file (atmfile)'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('tmodel, npars',
    [('isothermal', 1),
     ('guillot', 6),
     ('madhu', 6)])
def test_pt_tpars_mismatch(tmp_path, tmodel, npars):
    cfg = make_config(
        tmp_path,
        f'{ROOT}tests/configs/pt_{tmodel}.cfg',
        reset={'tpars':'100.0 200.0'},
    )
    error = re.escape(
        'Number of temperature parameters (2) does not match the '
        f'required number of parameters ({npars}) of the {tmodel} model'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# atmosphere (chemistry) runmode fails:

@pytest.mark.parametrize('chem', ['uniform', 'tea'])
def test_atmosphere_missing_species(tmp_path, chem):
    cfg = make_config(
        tmp_path,
        f'{ROOT}tests/configs/atmosphere_{chem}_test.cfg',
        remove=['species'],
    )
    error = re.escape('Undefined atmospheric species list (species)')
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_atmosphere_uniform_missing_uniform(tmp_path):
    cfg = make_config(tmp_path,
        ROOT+'tests/configs/atmosphere_uniform_test.cfg',
        remove=['uniform'])
    error = re.escape(
        'Undefined list of uniform volume mixing ratios (uniform) '
        'for uniform chemistry model'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_atmosphere_uniform_mismatch_uniform(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/atmosphere_uniform_test.cfg',
        reset={'uniform':'0.85 0.15'},
    )
    error = re.escape(
        'Number of uniform abundances (2) does not match the number '
        'of species (7)'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# atmosphere (altitude) runmode fails:

def test_atmosphere_hydro_missing_all_planet_props(tmp_path):
    cfg = make_config(
        tmp_path,
        f'{ROOT}tests/configs/atmosphere_hydro_test.cfg',
        remove=['rplanet', 'mplanet', 'gplanet', 'refpressure'],
    )
    error = re.escape(
        'Cannot compute hydrostatic-equilibrium radius profile.\n'
        'Undefined planet radius (rplanet).\n'
        'Undefined planet mass (mplanet) or surface gravity (gplanet).\n'
        'Undefined reference pressure level (refpressure).'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_atmosphere_hydro_missing_refpressure(tmp_path):
    cfg = make_config(
        tmp_path,
        f'{ROOT}tests/configs/atmosphere_hydro_test.cfg',
        remove=['refpressure'],
    )
    error = re.escape(
        'Cannot compute hydrostatic-equilibrium radius profile.\n'
        'Undefined reference pressure level (refpressure).'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_atmosphere_hydro_missing_mass_gravity(tmp_path):
    cfg = make_config(
        tmp_path,
        f'{ROOT}tests/configs/atmosphere_hydro_test.cfg',
        remove=['mplanet', 'gplanet'],
    )
    error = re.escape(
        'Cannot compute hydrostatic-equilibrium radius profile.\n'
        'Undefined planet mass (mplanet) or surface gravity (gplanet).'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# spectrum runmode fails (setup):

@pytest.mark.parametrize(
    'param',
    [
        'pbottom',
        'ptop',
        'refpressure',
        'mstar',
        'rstar',
        'smaxis',
    ])
def test_spectrum_missing_units(tmp_path, param):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={param:'1.1'},
    )
    error = f"Invalid units 'None' for parameter {param}"
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_missing_wl_units(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'wllow': '1.1', 'wlhigh': '2.0'},
    )
    error = "Invalid units 'None' for parameter wllow"
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_inconsistent_wl_bounds(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'wllow':'2.0 um', 'wlhigh':'1.0 um'},
    )
    error = re.escape(
        'Wavenumber low boundary (10000.0 cm-1) must be larger than the '
        'high boundary (5000.0 cm-1)')
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize(
    'param',
    ['rstar', 'rt_path'],
)
def test_spectrum_transmission_missing(tmp_path, param, undefined_spec):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        remove=[param],
    )
    error = re.escape(undefined_spec[param])
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_missing_chemistry(tmp_path):
    reset = {
        'ptop': '1e-6 bar',
        'pbottom': '100.0 bar',
        'nlayers': '81',
        'tmodel': 'guillot',
        'tpars': '-4.84 -0.8 -0.8 0.5 1200.0 100.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
        remove=['atmfile', 'radmodel'],
    )
    error = re.escape(
        'Missing atmospheric volume mixing ratios. Need to either read '
        'an input profile or compute one via a chemistry model (chemistry)'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_no_radius(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        remove=['radmodel'],
    )
    error = re.escape(
        'Missing atmospheric radius profile.  Need to either read an '
        'input profile or compute one via the radmodel argument'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('atm',
    [f'{ROOT}/tests/inputs/atmosphere_uniform_test.atm',
     f'{ROOT}/tests/inputs/atmosphere_uniform_radius.atm'])
def test_spectrum_hydro_missing_MGplanet(tmp_path, atm):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'atmfile': atm},
        remove=['mplanet', 'gplanet'],
    )
    error = re.escape(
        'Cannot compute hydrostatic-equilibrium radius profile.\n'
        'Undefined planet mass (mplanet) or surface gravity (gplanet).'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('atm',
    [f'{ROOT}/tests/inputs/atmosphere_uniform_test.atm',
     f'{ROOT}/tests/inputs/atmosphere_uniform_radius.atm'])
def test_spectrum_hydro_missing_rplanet(tmp_path, atm):
    #keep = ['mplanet', 'gplanet']
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'atmfile': atm},
        remove=['rplanet'],
    )
    error = re.escape(
        'Cannot compute hydrostatic-equilibrium radius profile.\n'
        'Undefined planet radius (rplanet).'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize(
    'atm',
    [
        f'{ROOT}/tests/inputs/atmosphere_uniform_test.atm',
        f'{ROOT}/tests/inputs/atmosphere_uniform_radius.atm',
    ]
)
def test_spectrum_hydro_refpressure(tmp_path, atm):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'atmfile':atm},
        remove=['refpressure'],
    )
    error = re.escape(
        'Cannot compute hydrostatic-equilibrium radius profile.\n'
        'Undefined reference pressure level (refpressure)'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.skip(reason='Do I want to reset boundaries?')
@pytest.mark.parametrize('value', ['1.00e-09 bar', '1.00e+03 bar'])
@pytest.mark.parametrize('param', ['pbottom', 'ptop'])
def test_spectrum_unbounded_pressures(tmp_path, param, value):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={param:value},
    )
    error = re.escape(
        f'{param[1:].capitalize()}-pressure boundary ({param}={value}) '
        'lies outside of the atmospheric file range 1.00e-06--1.00e+02 bar'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_invalid_pressure_ranges(tmp_path):
    reset = {
        'ptop': '1.0e-02 bar',
        'pbottom': '1.0e-03 bar',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        'Bottom-layer pressure (1.00e-03 bar) must be higher than the '
        'top-layer pressure (1.00e-02 bar)'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('param',
    ['tlifile',])
def test_spectrum_invalid_file(tmp_path, param, invalid_file):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={param:'nope.dat'},
    )
    error = re.escape(invalid_file[param]) + '*.'
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('vmin,vmax',
    [('dmin', 'dmax'),
     ('lmin', 'lmax')])
def test_spectrum_inconsistent_voigt_bounds(tmp_path, vmin, vmax):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset={vmin:'1e5', vmax:'1e4'},
    )
    error = re.escape(
        f'Voigt {vmax} (1.0000e+04 cm-1) must be > {vmin} (1.0000e+05 cm-1)'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_rpars_mismatch(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'rpars':'1.0 1.0 1.0'},
    )
    error = re.escape(
        'Number of input Rayleigh parameters (3) does not match the '
        'number of required model parameters (2)'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_cpars_mismatch(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'cpars':'1.0 1.0 1.0'},
    )
    error = re.escape(
        'Number of input cloud parameters (3) does not match the number '
        'of required model parameters (1)'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('value',
    ['10 60 90',
     '0 30 60 100',
     '0 30 90 60',
    ]
)
def test_spectrum_raygrid(tmp_path, invalid_raygrid, value):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'raygrid':value},
    )
    error = re.escape(invalid_raygrid[value])
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_line_sample_missing_species(tmp_path):
    reset = {
        'species': 'H2  He  H   Na  CH4  CO  CO2',
        'chemistry': 'tea',
        'extfile': f'{ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        remove=['tlifile'],
        reset=reset,
    )
    error = re.escape(
        "Species ['H2O'], required for opacity model line sampling, "
        "are not present in the atmosphere"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_line_by_line_missing_species(tmp_path):
    reset = {
        'species': 'H2  He  H   Na  CH4  CO  CO2',
        'chemistry': 'tea',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset=reset,
    )
    error = re.escape(
        "The species 'H2O' for isotopes ['161' '181' '171' '162' '182' "
        "'172' '262' '282' '272'] is not present in the atmosphere"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_alkali_missing_species(tmp_path):
    reset = {
        'alkali': 'potassium_vdw',
        'atmfile': '{ROOT}tests/inputs/atmosphere_uniform_no_potassium.atm',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        "Species ['K'], required for opacity model potassium_vdw, "
        "are not present in the atmosphere"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_cia_missing_species(tmp_path):
    reset = {
        'species': 'H2  H   Na  H2O CH4  CO  CO2',
        'chemistry': 'tea',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        "Species ['He'], required for opacity model CIA H2-He, "
        "are not present in the atmosphere"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_rayleigh_missing_species(tmp_path):
    reset = {
        'species': 'H2  H   Na  H2O CH4  CO  CO2',
        'chemistry': 'tea',
        'rayleigh': 'lecavelier dalgarno_He',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        remove=['csfile'],
        reset=reset,
    )
    error = re.escape(
        "Species ['He'], required for opacity model dalgarno_He, "
        "are not present in the atmosphere"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_h_ion_missing_species(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        remove=['tlifile', 'csfile', 'alkali', 'clouds'],
        reset={'h_ion': 'h_ion_john1988'},
    )
    error = re.escape(
        "Species ['e-'], required for opacity model H- "
        "bound-free/free-free, are not present in the atmosphere"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_uncert_mismatch(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'data':'1.0 2.0', 'uncert':'0.1'},
    )
    error = re.escape(
        'Number of data uncertainty values (1) does not match the '
        'number of data points (2)'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_filters_mismatch(tmp_path):
    reset = {
        'data': '1.0 2.0',
        'filters': ROOT+'tests/filters/filter_test_WFC3_G141_1.133um.dat',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        'Number of filter bands (1) does not match the number of '
        'data points (2)'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)

@pytest.mark.parametrize(
    'ncolumns',
    [1, 3, 4, 6, 9],
)
def test_spectrum_invalid_retrieval_params_entry(tmp_path, ncolumns):
    entry = 'T_iso 1500.0 300.0 3500.0 10.0 900.0 100.0 100.0 1.0'.split()
    ret_pars = " ".join(entry[0:ncolumns])
    reset = {
        'tmodel':'isothermal',
        'retrieval_params': ret_pars,
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        f"Invalid number of fields for retrieval_params entry\n'{ret_pars}'"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_invalid_retrieval_params_pname(tmp_path):
    reset = {
        'tmodel':'isothermal',
        'retrieval_params': 'log_H2O -3.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        "Invalid retrieval parameter 'log_H2O'. Possible values are:\n"
        "['log_p_ref', 'R_planet', 'M_planet', 'f_patchy', 'T_eff', "
        "'f_dilution', 'T_iso', 'log_k_ray', 'alpha_ray', 'log_p_cl']"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_repeated_retrieval_params(tmp_path):
    reset = {
        'retrieval_params': 'R_planet 1.0\nR_planet 1.5',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape("Repeated parameter names: ['R_planet']")
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_insuficient_retrieval_params_temp(tmp_path):
    reset = {
        'tmodel': 'isothermal',
        'retrieval_params': 'R_planet -3.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape('Not all temperature parameters were defined (tpars)')
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_insuficient_retrieval_params_mol(tmp_path):
    reset = {
        'chemistry': 'tea',
        'vmr_vars': '[M/H]',
        'retrieval_params': 'R_planet -3.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape('Not all vmr parameter values were defined (vmr_vars)')
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_insuficient_retrieval_params_cloud(tmp_path):
    reset = {
        'cloud': 'deck',
        'retrieval_params': 'R_planet -3.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
        remove=['cpars'],
    )
    error = re.escape("Undefined parameter values for cloud model 'deck'")
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_insuficient_retrieval_params_rayleigh(tmp_path):
    reset = {
        'rayleigh': 'lecavelier',
        'retrieval_params': 'R_planet -3.0',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
        remove=['rpars'],
    )
    error = re.escape("Undefined parameter values for cloud model 'lecavelier'")
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_spectrum_params_misfit(tmp_path):
    # Without evaulating params:
    reset = {
        'tmodel':'guillot',
        'retrieval_params': 'T_iso 1500.0 ',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        """Invalid retrieval parameter 'T_iso'. Possible values are:\n"""
        """['log_p_ref', 'R_planet', 'M_planet', 'f_patchy', 'T_eff', """
        """'f_dilution', "log_kappa'", 'log_gamma1', 'log_gamma2', 'alpha', """
        """'T_irr', 'T_int', 'log_k_ray', 'alpha_ray', 'log_p_cl']"""
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_eval_params_misfit(tmp_path, capfd):
    # Without evaulating params:
    reset = {
        'tmodel': 'isothermal',
        'retrieval_params': 'T_iso 1500.0 ',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    pyrat = pb.run(cfg)
    captured = capfd.readouterr()

    # Now, only a warning
    warning = (
        "The number of input fitting parameters (2) does not match\n"
        "    the number of required parameters (1)"
    )
    pyrat.eval([1200.0, 0.5])
    captured = capfd.readouterr()
    assert warning in captured.out


def test_bulk_not_in_atm(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'bulk':'N2'},
    )
    error = re.escape(
        "These bulk species are not present in the atmosphere: ['N2']"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_vmr_vars_missing_vmr(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'vmr_vars': 'log_N2'},
    )
    error = re.escape(
        "Invalid vmr_vars variable 'log_N2', species N2 is not in the atmosphere"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_bulk_molfree_overlap(tmp_path):
    reset = {
        'vmr_vars': 'log_H2',
        'bulk': 'H2',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        "These species were marked as both bulk and variable-abundance: ['H2']"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_invalid_equil_without_tea(tmp_path):
    reset = {
        'vmr_vars': '[M/H]',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape("vmr_vars variable '[M/H]' requires chemistry=tea")
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_invalid_equil_vmr_vars(tmp_path):
    reset = {
        'chemistry': 'tea',
        'vmr_vars': 'zen',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape("Unrecognized VMR model (vmr_vars): 'zen'")
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_invalid_metal_vmr_vars(tmp_path):
    reset = {
        'chemistry': 'tea',
        'vmr_vars': '[X/H]',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        "Invalid vmr_vars variable '[X/H]', "
        "element 'X' is not in the atmosphere"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_invalid_ratio_vmr_vars(tmp_path):
    reset = {
        'chemistry': 'tea',
        'vmr_vars': 'C/O/X',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        "Invalid vmr_vars variable 'C/O/X', "
        "element 'O/X' is not in the atmosphere"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize('param', ['tstar', 'log_gstar'])
def test_kurucz_missing_pars(tmp_path, param):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'kurucz':f'{ROOT}/tests/inputs/mock_fp00k0odfnew.pck'},
        remove=[param],
    )
    if param == 'tstar':
        error = re.escape('Undefined stellar temperature (tstar)')
    elif param == 'log_gstar':
        error = re.escape('Undefined stellar gravity (log_gstar)')
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.filterwarnings("ignore: The 'retflag' argument")
@pytest.mark.parametrize('param',
    ['tmodel', 'clouds', 'rayleigh', 'vmr_vars', 'bulk'])
def test_spectrum_missing_retflag_models(tmp_path, param, undefined_mcmc):
    reset = {
        'retflag': 'temp mol ray cloud',
        'tmodel': 'isothermal',
        'clouds': 'deck',
        'rayleigh': 'lecavelier',
        'vmr_vars': 'log_H2O',
        'bulk': 'H2',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        remove=[param],
        reset=reset,
    )
    error = undefined_mcmc[param]
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_compute_opacity_invalid_tmin(tmp_path):
    reset = {
        'runmode': 'opacity',
        'extfile': str(tmp_path/'new_opacity.npz'),
        'tmin': '0.1',
        'tmax': '1000.0',
        'tstep': '900',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset=reset,
    )
    error = re.escape(
        "Requested cross-section table temperature (tmin=0.1 K) "
        "below the lowest available TLI temperature (1.0 K)"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_compute_opacity_invalid_tmax(tmp_path):
    reset = {
        'runmode': 'opacity',
        'extfile': str(tmp_path/'new_opacity.npz'),
        'tmin':'1000.0',
        'tmax':'9000.0',
        'tstep':'100',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset=reset,
    )
    error = re.escape(
        "Requested cross-section table temperature (tmax=9000.0 K) "
        "above the highest available TLI temperature (6000.0 K)"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_read_opacity_missing_extfile(tmp_path):
    efile = str(tmp_path/'non_existent_exttable_test_300-3000K_1.1-1.7um.npz')
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'extfile':efile, 'wnstep': '1.0'},
    )
    error = re.escape(
        f"These input cross-section files are missing: ['{efile}']"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_compute_opacity_make_multiple_extfiles(tmp_path):
    extfiles = [
        str(tmp_path/'new_opacity1.npz'),
        str(tmp_path/'new_opacity2.npz'),
    ]
    reset = {
        'runmode': 'opacity',
        'extfile': '\n  '.join(extfiles),
        'tmin':'1000.0',
        'tmax':'6000.0',
        'tstep':'100',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_test.cfg',
        reset=reset,
    )
    error = re.escape(
        'Computing opacity table, but there was more than one '
        'output opacity file (extfile)'
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.skip(reason='TBI')
def test_read_opacity_mismatched_sizes(tmp_path):
    extfiles = [
        str(tmp_path/'new_opacity.npz'),
        str(tmp_path/'new_opacity.npz'),
    ]
    reset = {
        'extfile': '\n    '.join(extfiles),
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        f"Shape of the cross-section file '{extfiles[1]}' "
        "does not match with previous file shapes."
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.skip(reason='TBI')
def test_read_opacity_mismatched_values(tmp_path):
    efiles = [
        f'{ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz',
        str(tmp_path/'new_opacity.npz'),
    ]
    reset = {
        'extfile': '\n    '.join(efiles),
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
    )
    error = re.escape(
        f"Tabulated temperature values in file '{efiles[1]}' "
        "do not match with previous arrays"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.parametrize(
    'param',
    ['tmin', 'tmax', 'tstep', 'tlifile'],
)
def test_opacity_missing(tmp_path, param):
    reset = {'extfile': str(tmp_path/'new_opacity.npz')}
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/opacity_test.cfg',
        reset=reset,
        remove=[param],
    )
    if param == 'tlifile':
        err = 'input TLI files (tlifile)'
    if param == 'tmin':
        err = 'lower temperature boundary (tmin)'
    if param == 'tmax':
        err = 'upper temperature boundary (tmax)'
    if param == 'tstep':
        err = 'temperature sampling step (tstep)'
    error = re.escape(f'Undefined {err} needed to compute opacity table')

    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_molecule_not_in_molfile(tmp_path):
    # Modify atm:
    units, species, press, temp, q, rad = \
        io.read_atm(ROOT+'tests/inputs/atmosphere_uniform_test.atm')
    press = press * pb.tools.u(units[0])
    species[-1] = 'X'
    new_atm = str(tmp_path/'new_atmosphere_uniform_test.atm')
    io.write_atm(new_atm, press, temp, species, q, punits=units[0])

    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset={'atmfile':new_atm},
    )
    error = re.escape(
        "These species: ['X'] are not listed in the molecules info file:"
    )
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


@pytest.mark.skip
def test_incompatible_tli():
    pass

@pytest.mark.skip
def test_crosssec_mol_not_in_atm():
    pass


@pytest.mark.parametrize(
    'param',
    [
        'data',
        'uncert',
        'filters',
        'rstar',
        'sampler',
        'nsamples',
        'burnin',
        'nchains',
    ]
)
def test_mcmc_missing(tmp_path, param, undefined_mcmc):
    reset = {
        'rt_path': 'emission',
        'kurucz': f'{ROOT}tests/inputs/mock_fp00k0odfnew.pck',
        'log_gstar': '4.5',
        'sampler': 'snooker',
    }
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/mcmc_transmission_test.cfg',
        reset=reset,
        remove=[param],
    )
    error = re.escape(undefined_mcmc[param])
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)


def test_mcmc_missing_starspec(tmp_path):
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/mcmc_transmission_test.cfg',
        reset={'rt_path':'emission'},
        remove=['tstar', 'tmodel'],
    )
    error = re.escape(
        'Undefined stellar flux model.  Set starspec, kurucz, or tstar '
        '(for a blackbody spectrum)')
    with pytest.raises(ValueError, match=error):
        pyrat = pb.run(cfg)



@pytest.mark.skip(reason='TBD: warning or error...')
@pytest.mark.parametrize('param', ['tlifile', 'csfile', 'extfile'])
def test_spectrum_temperature_bounds(tmp_path, capfd, param, invalid_temp):
    reset = {
        'tmodel': 'isothermal',
        'tpars': '6000.0',
        'extfile': f'{ROOT}tests/outputs/exttable_test_300-3000K_1.1-1.7um.npz',
    }
    remove = [
        par
        for par in ['tlifile', 'csfile', 'extfile']
        if par != param
    ]
    cfg = make_config(
        tmp_path,
        ROOT+'tests/configs/spectrum_transmission_extfile.cfg',
        reset=reset,
        remove=remove,
    )
    pyrat = pb.run(cfg)
    assert pyrat is not None
    captured = capfd.readouterr()
    assert invalid_temp[param] in captured.out
    assert pyrat.spec.spectrum is None

