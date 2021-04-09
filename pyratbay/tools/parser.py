# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'Namespace',
    'parse',
    'parse_str',
    'parse_int',
    'parse_float',
    'parse_array',
    ]

import os
import argparse
from datetime import date
import configparser

import numpy as np
import mc3.utils as mu

from . import tools as pt
from .. import constants as pc
from ..VERSION import __version__


class Namespace(argparse.Namespace):
    """A container object to hold variables."""
    def __init__(self, args=None, log=None):
        if args is not None:
            super(Namespace, self).__init__(**args)
        if log is None:
            self._log = mu.Log(logname=None, verb=2, width=80)
        else:
            self._log = log

    def get_path(self, pname, desc='', exists=False):
        """
        Extract pname file path (or list of paths) from Namespace,
        return the canonical path.

        Examples
        --------
        >>> import pyratbay.tools as pt
        >>> ns = pt.Namespace({'f1':'file1', 'f23':['file2', 'file3']})
        >>> # Get path of a single file:
        >>> ns.get_path('f1')
        >>> # Get path of a list of files:
        >>> ns.get_path('f23')
        >>> # Attempt to get non-existing file:
        >>> ns.get_path('f1', desc='Configuration', exists=True)
        """
        value = getattr(self, pname)
        if value is None:
            return None

        if isinstance(value, str):
            values = [value]
            is_list = False
        else:
            values = value
            is_list = True

        values = [os.path.realpath(val.replace('{ROOT}', pc.ROOT))
                  for val in values]

        for val in values:
            if exists and not os.path.isfile(val):
                self._log.error(
                    f"{desc} file ({pname}) does not exist: '{val}'",
                    tracklev=-3)
            if not os.path.exists(os.path.dirname(val)):
                self._log.error(
                    f"Folder for {desc} file ({pname}) does not exist: '{val}'",
                    tracklev=-3)

        if not is_list:
            return values[0]
        return values

    def get_choice(self, pname, desc, choices, take_none=True):
        value = getattr(self, pname)
        if value is None and take_none:
            return None

        if isinstance(value, list):
            values = value
            is_list = True
        else:
            values = [value]
            is_list = False

        for value in values:
            if value not in choices:
                self._log.error(f"Invalid {desc} ({pname}): {value}. "
                    f"Select from: {choices}.", tracklev=-3)
        if not is_list:
            return values[0]
        return values

    def get_default(self, pname, desc, default=None, wflag=False,
                    gt=None, ge=None, lt=None, le=None):
        """
        Extract pname variable from Namespace; if None, return
        default.  If any of gt, ge, lt, or le is not None, run
        greater/lower/equal checks.

        Parameters
        ----------
        pname: String
            Parameter name.
        desc: String
            Parameter description.
        default: Any
            Parameter default value.
        gt: Float
            If not None, check output is greater than gt.
        ge: Float
            If not None, check output is greater-equal than gt.
        lt: Float
            If not None, check output is lower than gt.
        le: Float
            If not None, check output is lower-equal than gt.
        """
        value = getattr(self, pname)
        if value is None and default is not None:
            if wflag:
                self._log.warning(f'{desc} ({pname}) defaulted to: {default}')
            value = default

        if value is None:
            return None

        if gt is not None and value <= gt:
            self._log.error(f'{desc} ({pname}) must be > {gt}', tracklev=-3)
        if ge is not None and value < ge:
            self._log.error(f'{desc} ({pname}) must be >= {ge}', tracklev=-3)
        if lt is not None and lt <= value:
            self._log.error(f'{desc} ({pname}) must be < {lt}', tracklev=-3)
        if le is not None and le < value:
            self._log.error(f'{desc} ({pname}) must be <= {le}', tracklev=-3)
        return value

    def get_units(self, pname):
        """
        Extract units from a value input.
        Return None if value does not have units or has an invalid format.

        Parameters
        ----------
        pname: String
            Parameter name.

        Returns
        -------
        units: String
        """
        value = getattr(self, pname)
        if not isinstance(value, str):
            return None

        par = value.split()
        if len(par) != 2:
            return None

        units = par[1]
        return units

    def get_param(self, pname, units, desc, gt=None, ge=None):
        try:
            value = pt.get_param(getattr(self, pname), units)
        except ValueError as error:
            self._log.error(f'{error} for parameter {pname}.', tracklev=-3)

        if value is None:
            return None

        if gt is not None and value <= gt:
            self._log.error(f'{desc} ({pname}) must be > {gt}', tracklev=-3)
        if ge is not None and value < ge:
            self._log.error('{desc} ({pname}) must be >= {ge}', tracklev=-3)
        return value


def parse_str(args, param):
    """Parse a string parameter into args."""
    if param not in args:
        args[param] = None
    else:
        args[param] = str(args[param])


def parse_int(args, param):
    """
    Convert a dictionary's parameter from string to integer.
    Raise ValueError if the operation is not possible.
    Set parameter to None if it was not in the dictinary.

    Parameters
    ----------
    args: dict
        Dictionary where to operate.
    param: String
        Parameter to cast to int.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> inputs = ['10', '-10', '+10', '10.0', '1e1',
    >>>           '10.5', 'None', 'True', 'inf', '10 20']
    >>> args = {f'par{i}':val for i,val in enumerate(inputs)}
    >>> for i,var in enumerate(inputs):
    >>>     try:
    >>>         par = f'par{i}'
    >>>         pt.parse_int(args, par)
    >>>         print(f"{par}: '{var}' -> {args[par]}")
    >>>     except ValueError as e:
    >>>         print(e)
    par0: '10' -> 10
    par1: '-10' -> -10
    par2: '+10' -> 10
    par3: '10.0' -> 10
    par4: '1e1' -> 10
    Invalid data type for par5, could not convert string to integer: '10.5'
    Invalid data type for par6, could not convert string to integer: 'None'
    Invalid data type for par7, could not convert string to integer: 'True'
    Invalid data type for par8, could not convert string to integer: 'inf'
    Invalid data type for par9, could not convert string to integer: '10 20'
    """
    if param not in args:
        args[param] = None
        return

    try:
        val = np.double(args[param])
    except:
        raise ValueError(
            f"Invalid data type for {param}, could not convert string "
            f"to integer: '{args[param]}'")
    if not np.isfinite(val) or int(val) != val:
        raise ValueError(
            f"Invalid data type for {param}, could not convert string "
            f"to integer: '{args[param]}'")
    args[param] = int(val)



def parse_float(args, param):
    """
    Convert a dictionary's parameter from string to float.
    Raise ValueError if the operation is not possible.
    Set parameter to None if it was not in the dictinary.

    Parameters
    ----------
    args: dict
        Dictionary where to operate.
    param: String
        Parameter to cast to float.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> inputs = ['10', '-10', '+10', '10.5', '1e1', 'inf', 'nan',
    >>>           'None', 'True', '10 20']
    >>> args = {f'par{i}':val for i,val in enumerate(inputs)}
    >>> for i,var in enumerate(inputs):
    >>>     try:
    >>>         par = f'par{i}'
    >>>         pt.parse_float(args, par)
    >>>         print(f"{par}: '{var}' -> {args[par]}")
    >>>     except ValueError as e:
    >>>         print(e)
    par0: '10' -> 10.0
    par1: '-10' -> -10.0
    par2: '+10' -> 10.0
    par3: '10.5' -> 10.5
    par4: '1e5' -> 10.0
    par5: 'inf' -> inf
    par6: 'nan' -> nan
    Invalid data type for par7, could not convert string to float: 'None'
    Invalid data type for par8, could not convert string to float: 'True'
    Invalid data type for par9, could not convert string to float: '10 20'
    """
    if param not in args:
        args[param] = None
        return

    try:
        val = np.double(args[param])
    except:
        raise ValueError(
            f"Invalid data type for {param}, could not convert string "
            f"to float: '{args[param]}'")
    args[param] = val


def parse_array(args, param):
    r"""
    Convert a dictionary's parameter from string to iterable.
    If possible cast into a float numpy array; otherwise,
    set as a list of strings.
    Assume any blank character delimits the elements in the string.
    Set parameter to None if it was not in the dictinary.

    Parameters
    ----------
    args: dict
        Dictionary where to operate.
    param: String
        Parameter to cast to array.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> inputs = ['10 20', '10.0 20.0', 'a b', 'a\n b']
    >>> args = {f'par{i}':val for i,val in enumerate(inputs)}
    >>> for i,var in enumerate(inputs):
    >>>     par = f'par{i}'
    >>>     pt.parse_array(args, par)
    >>>     print(f"{par}: {repr(var)} -> {repr(args[par])}")
    par0: '10 20' -> array([10., 20.])
    par1: '10.0 20.0' -> array([10., 20.])
    par2: 'a b' -> ['a', 'b']
    par3: 'a\n b' -> ['a', 'b']
    """
    if param not in args:
        args[param] = None
        return

    val = args[param].split()
    try:
        val = np.asarray(val, np.double)
    except:
        pass
    args[param] = val


def parse(pyrat, cfile, no_logfile=False, mute=False):
    """
    Read the command line arguments.

    Parameters
    ----------
    cfile: String
        A Pyrat Bay configuration file.
    no_logfile: Bool
        If True, enforce not to write outputs to a log file
        (e.g., to prevent overwritting log of a previous run).
    mute: Bool
        If True, enforce verb to take a value of -1.
    """
    with pt.log_error():
        if not os.path.isfile(cfile):
            raise ValueError(f"Configuration file '{cfile}' does not exist.")
        config = configparser.ConfigParser()
        config.optionxform = str  # Enable case-sensitive variable names
        config.read([cfile])
        if "pyrat" not in config.sections():
            raise ValueError(
                f"\nInvalid configuration file: '{cfile}', no [pyrat] section.")
    args = dict(config.items("pyrat"))

    # Parse data type:
    with pt.log_error():
        parse_str(args,   'runmode')
        parse_int(args,   'ncpu')
        parse_int(args,   'verb')
        parse_array(args, 'dblist')
        parse_array(args, 'pflist')
        parse_array(args, 'dbtype')
        parse_array(args, 'tlifile')
        parse_array(args, 'csfile')
        parse_str(args,   'molfile')
        parse_array(args, 'extfile')
        # Spectrum sampling options:
        parse_str(args,   'wlunits')
        parse_str(args,   'wllow')
        parse_str(args,   'wlhigh')
        parse_float(args, 'wnlow')
        parse_float(args, 'wnhigh')
        parse_float(args, 'wnstep')
        parse_int(args,   'wnosamp')
        parse_float(args, 'resolution')
        # Atmospheric sampling options:
        parse_str(args,   'tmodel')
        parse_array(args, 'tpars')
        parse_str(args,   'radlow')
        parse_str(args,   'radhigh')
        parse_str(args,   'radstep')
        parse_str(args,   'runits')
        parse_str(args,   'punits')
        parse_int(args,   'nlayers')
        parse_str(args,   'ptop')
        parse_str(args,   'pbottom')
        parse_str(args,   'atmfile')
        parse_str(args,   'radmodel')
        # Variables for chemistry calculations
        parse_str(args,   'chemistry')
        parse_array(args, 'species')
        parse_array(args, 'uniform')
        parse_str(args,   'ptfile')
        parse_str(args,   'solar')
        parse_float(args, 'xsolar')
        parse_array(args, 'escale')
        parse_array(args, 'elements')
        # Extinction options:
        parse_float(args, 'tmin')
        parse_float(args, 'tmax')
        parse_float(args, 'tstep')
        parse_float(args, 'ethresh')
        # Voigt-profile options:
        parse_float(args, 'vextent')
        parse_float(args, 'vcutoff')
        parse_float(args, 'dmin')
        parse_float(args, 'dmax')
        parse_int(args,   'ndop')
        parse_float(args, 'lmin')
        parse_float(args, 'lmax')
        parse_int(args,   'nlor')
        parse_float(args, 'dlratio')
        # Hazes and clouds options:
        parse_array(args, 'clouds')
        parse_array(args, 'cpars')
        parse_array(args, 'rayleigh')
        parse_array(args, 'rpars')
        parse_float(args, 'fpatchy')
        parse_array(args, 'alkali')
        parse_float(args, 'alkali_cutoff')
        # Optical depth options:
        parse_str(args,   'rt_path')
        parse_float(args, 'maxdepth')
        parse_array(args, 'raygrid')
        parse_int(args,   'quadrature')
        # Data options:
        parse_str(args,   'dunits')
        parse_array(args, 'data')
        parse_array(args, 'uncert')
        parse_array(args, 'filters')
        # Abundances:
        parse_array(args, 'molmodel')
        parse_array(args, 'molfree')
        parse_array(args, 'molpars')
        parse_array(args, 'bulk')
        # Retrieval options:
        parse_str(args,   'mcmcfile')
        parse_array(args, 'retflag')
        parse_float(args, 'qcap')
        parse_array(args, 'params')
        parse_array(args, 'pstep')
        parse_float(args, 'tlow')
        parse_float(args, 'thigh')
        parse_array(args, 'pmin')
        parse_array(args, 'pmax')
        parse_array(args, 'prior')
        parse_array(args, 'priorlow')
        parse_array(args, 'priorup')
        parse_str(args,   'sampler')
        parse_int(args,   'nsamples')
        parse_int(args,   'nchains')
        parse_int(args,   'burnin')
        parse_int(args,   'thinning')
        parse_float(args, 'grbreak')
        parse_float(args, 'grnmin')
        parse_int(args,   'resume')      # False, action='store_true')
        # Stellar models:
        parse_str(args,   'starspec')
        parse_str(args,   'kurucz')
        parse_str(args,   'marcs')
        parse_str(args,   'phoenix')
        # System parameters:
        parse_str(args,   'rstar')
        parse_float(args, 'gstar')
        parse_float(args, 'tstar')
        parse_str(args,   'mstar')
        parse_str(args,   'rplanet')
        parse_str(args,   'refpressure')
        parse_str(args,   'mplanet')
        parse_str(args,   'mpunits')
        parse_float(args, 'gplanet')
        parse_str(args,   'smaxis')
        parse_float(args, 'tint')
        # Outputs:
        parse_str(args,   'specfile')
        parse_str(args,   'logfile')
        parse_array(args, 'logxticks')
        parse_array(args, 'yran')

    # Cast into a Namespace to make my life easier:
    args = Namespace(args)
    args.configfile = cfile

    if mute:
        args.verb = -1
    pyrat.verb = args.get_default('verb', 'Verbosity', 2, lt=5)
    runmode = pyrat.runmode = args.get_choice(
        'runmode', 'running mode', pc.rmodes, take_none=False)

    # Define logfile name and initialize log object:
    pyrat.lt.tlifile = args.get_path('tlifile', 'TLI')
    pyrat.atm.atmfile = args.get_path('atmfile', 'Atmospheric')
    pyrat.spec.specfile = args.get_path('specfile', 'Spectrum')
    pyrat.ex.extfile = args.get_path('extfile', 'Extinction-coefficient')
    pyrat.ret.mcmcfile = args.get_path('mcmcfile', 'MCMC')

    outfile_dict = {
        'tli': args.tlifile,
        'atmosphere': args.atmfile,
        'spectrum': args.specfile,
        'opacity': args.extfile,
        'mcmc': args.mcmcfile,
        }
    outfile = outfile_dict[args.runmode]
    if args.logfile is None and outfile is not None:
        if args.runmode in ['tli', 'opacity']:
            outfile = outfile[0]
        args.logfile = os.path.splitext(outfile)[0] + '.log'

    args.logfile = args.get_path('logfile', 'Log')

    # Override logfile if requested:
    if no_logfile:
        args.logfile = None

    args.logfile = pt.path(args.logfile)
    log = pyrat.log = mu.Log(
        logname=args.logfile, verb=pyrat.verb, width=80, append=args.resume)
    args._log = log

    # Welcome message:
    log.head(
       f"{log.sep}\n"
        "  Python Radiative Transfer in a Bayesian framework (Pyrat Bay).\n"
       f"  Version {__version__}.\n"
       f"  Copyright (c) 2021-{date.today().year} Patricio Cubillos.\n"
        "  Pyrat Bay is open-source software under the GNU GPLv2 lincense "
          "(see LICENSE).\n"
       f"{log.sep}\n\n")

    log.head(f"Read command-line arguments from configuration file: '{cfile}'")

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Parse valid inputs and defaults:
    pyrat.inputs = args

    phy  = pyrat.phy
    spec = pyrat.spec
    atm  = pyrat.atm

    pyrat.lt.dblist = args.get_path('dblist', 'Opacity database', exists=True)
    pyrat.mol.molfile = args.get_path('molfile', 'Molecular data', exists=True)
    pyrat.cs.files = args.get_path('csfile', 'Cross-section', exists=True)
    pyrat.atm.ptfile = args.get_path('ptfile', 'Pressure-temperature')

    spec.wlunits = args.get_default('wlunits', 'Wavelength units')
    if spec.wlunits is not None and not hasattr(pc, spec.wlunits):
        log.error(f'Invalid wavelength units (wlunits): {spec.wlunits}')
    spec.wllow = args.get_param(
        'wllow', spec.wlunits, 'Wavelength lower boundary', gt=0.0)
    spec.wlhigh = args.get_param(
        'wlhigh', spec.wlunits, 'Wavelength higher boundary', gt=0.0)
    if spec.wlunits is None:
        spec.wlunits = args.get_units('wllow')

    spec.wnlow  = args.get_default(
        'wnlow', 'Wavenumber lower boundary',  gt=0.0)
    spec.wnhigh = args.get_default(
        'wnhigh', 'Wavenumber higher boundary', gt=0.0)

    spec.wnstep = args.get_default(
        'wnstep', 'Wavenumber sampling step', gt=0.0)
    spec.wnosamp = args.get_default(
        'wnosamp', 'Wavenumber oversampling factor', ge=1)
    spec.resolution = args.get_default(
        'resolution', 'Spectral resolution', gt=0.0)

    atm.runits = args.get_default('runits', 'Planetary-radius units')
    if atm.runits is not None and not hasattr(pc, atm.runits):
        log.error(f'Invalid radius units (runits): {atm.runits}')
    phy.rplanet = args.get_param(
        'rplanet', atm.runits, 'Planetary radius', gt=0.0)
    if atm.runits is None:
        atm.runits = args.get_units('rplanet')

    atm.nlayers = args.get_default(
        'nlayers', 'Number of atmospheric layers', gt=1)

    atm.rmodelname = args.get_choice(
        'radmodel', 'Radius-profile model', pc.radmodels)

    # Pressure inputs:
    atm.punits = args.get_default('punits', 'Pressure units')
    if atm.punits is not None and not hasattr(pc, atm.punits):
        log.error(f'Invalid pressure units (punits): {atm.punits}')

    atm.pbottom = args.get_param(
        'pbottom', atm.punits, 'Pressure at bottom of atmosphere', gt=0.0)
    atm.ptop = args.get_param(
        'ptop', atm.punits, 'Pressure at top of atmosphere', gt=0.0)
    atm.refpressure = args.get_param(
        'refpressure', atm.punits, 'Planetary reference pressure level', gt=0.0)

    if atm.punits is None and atm.pbottom is not None:
        atm.punits = args.get_units('pbottom')
    elif atm.punits is None and atm.ptop is not None:
        atm.punits = args.get_units('ptop')
    elif atm.punits is None and atm.refpressure is not None:
        atm.punits = args.get_units('refpressure')
    # else, set atm.punits from atmospheric file in read_atm().

    # Radius boundaries:
    atm.radlow = args.get_param(
        'radlow', atm.runits, 'Radius at bottom of atmosphere', ge=0.0)
    atm.radhigh = args.get_param(
        'radhigh', atm.runits, 'Radius at top of atmosphere', gt=0.0)
    atm.radstep = args.get_param(
        'radstep', atm.runits, 'Radius sampling step', gt=0.0)

    # Chemistry:
    atm.chemistry = args.get_choice(
       'chemistry', 'Chemical model', pc.chemmodels)

    escale = args.get_default(
       'escale', 'Elemental abundance scaling factors', [])
    atm.escale = {
        atom: float(fscale)
        for atom,fscale in zip(escale[::2], escale[1::2])
        }

    # System physical parameters:
    phy.mpunits = args.get_default('mpunits', 'Planetary-mass units')
    if phy.mpunits is not None and not hasattr(pc, phy.mpunits):
        log.error(f'Invalid planet mass units (mpunits): {phy.mpunits}')
    phy.mplanet = args.get_param(
        'mplanet', phy.mpunits, 'Planetary mass', gt=0.0)
    if phy.mpunits is None:
        phy.mpunits = args.get_units('mplanet')

    phy.gplanet = args.get_default(
        'gplanet', 'Planetary surface gravity (cm s-2)', gt=0.0)
    phy.tint = args.get_default(
        'tint', 'Planetary internal temperature', 100.0, ge=0.0)

    phy.smaxis = args.get_param(
        'smaxis', None, 'Orbital semi-major axis', gt=0.0)
    phy.rstar = args.get_param(
        'rstar', None, 'Stellar radius', gt=0.0)
    phy.mstar = args.get_param(
        'mstar', None, 'Stellar mass', gt=0.0)
    phy.gstar = args.get_default(
        'gstar', 'Stellar surface gravity', gt=0.0)
    phy.tstar = args.get_default(
        'tstar', 'Stellar effective temperature (K)', gt=0.0)

    pyrat.voigt.extent = args.get_default(
        'vextent', 'Voigt profile extent in HWHM', 100.0, ge=1.0)
    pyrat.voigt.cutoff = args.get_default(
        'vcutoff', 'Voigt profile cutoff in cm-1', 25.0, ge=0.0)
    pyrat.voigt.ndop = args.get_default(
        'ndop', 'Number of Doppler-width samples', 40, ge=1)
    pyrat.voigt.dmin = args.get_default(
        'dmin', 'Minimum Doppler HWHM (cm-1)', gt=0.0)
    pyrat.voigt.dmax = args.get_default(
        'dmax', 'Maximum Doppler HWHM (cm-1)', gt=0.0)
    pyrat.voigt.nlor = args.get_default(
        'nlor', 'Number of Lorentz-width samples', 40, ge=1)
    pyrat.voigt.lmin = args.get_default(
        'lmin', 'Minimum Lorentz HWHM (cm-1)', gt=0.0)
    pyrat.voigt.lmax = args.get_default(
        'lmax', 'Maximum Lorentz HWHM (cm-1)', gt=0.0)
    pyrat.voigt.dlratio = args.get_default(
        'dlratio', 'Doppler/Lorentz-width ratio threshold', 0.1, gt=0.0)

    pyrat.ex.tmin = args.get_param(
        'tmin', 'kelvin', 'Minimum temperature of opacity grid', gt=0.0)
    pyrat.ex.tmax = args.get_param(
        'tmax', 'kelvin', 'Maximum temperature of opacity grid', gt=0.0)
    pyrat.ex.tstep = args.get_default(
        'tstep', "Opacity grid's temperature sampling step in K", gt=0.0)

    pyrat.rayleigh.pars = args.rpars
    pyrat.cloud.pars     = args.cpars
    pyrat.rayleigh.model_names = args.get_choice(
        'rayleigh', 'Rayleigh model', pc.rmodels)
    pyrat.cloud.model_names = args.get_choice(
        'clouds', 'cloud model', pc.cmodels)
    pyrat.alkali.model_names = args.get_choice(
        'alkali', 'alkali model', pc.amodels)
    pyrat.alkali.cutoff = args.get_default(
        'alkali_cutoff',
        'Alkali profiles hard cutoff from line center (cm-1)', 4500.0, gt=0.0)
    pyrat.cloud.fpatchy = args.get_default(
        'fpatchy', 'Patchy-cloud fraction', ge=0.0, le=1.0)

    pyrat.od.rt_path = args.get_choice(
        'rt_path', 'radiative-transfer observing geometry', pc.rt_paths)
    pyrat.spec._rt_path = pyrat.od.rt_path

    pyrat.ex.ethresh = args.get_default(
        'ethresh', 'Extinction-cofficient threshold', 1e-15, gt=0.0)
    pyrat.od.maxdepth = args.get_default(
        'maxdepth', 'Maximum optical-depth', 10.0, ge=0.0)

    phy.starspec = args.get_path('starspec', 'Stellar spectrum', exists=True)
    phy.kurucz   = args.get_path('kurucz',   'Kurucz model',     exists=True)
    phy.marcs    = args.get_path('marcs',    'MARCS model',      exists=True)
    phy.phoenix  = args.get_path('phoenix',  'PHOENIX model',    exists=True)

    spec.raygrid = args.get_default(
        'raygrid', 'Emission raygrid (deg)', np.array([0, 20, 40, 60, 80.]))
    spec.quadrature = args.get_default(
        'quadrature', 'Number of Gaussian-quadrature points', ge=1)

    pyrat.obs.units = args.get_default(
        'dunits', 'Data units', 'none', wflag=args.data is not None)
    if not hasattr(pc, pyrat.obs.units):
        log.error(f'Invalid data units (dunits): {pyrat.obs.units}')
    pyrat.obs.data = args.get_param('data', pyrat.obs.units, 'Data')
    pyrat.obs.uncert = args.get_param(
        'uncert', pyrat.obs.units, 'Uncertainties')
    pyrat.obs.filters = args.get_path(
        'filters', 'Filter pass-bands', exists=True)

    pyrat.ret.retflag = args.get_choice(
        'retflag', 'retrieval flag', pc.retflags)
    if pyrat.ret.retflag is None:
        pyrat.ret.retflag = []

    pyrat.ret.params   = args.params
    pyrat.ret.pstep    = args.pstep
    pyrat.ret.pmin     = args.pmin
    pyrat.ret.pmax     = args.pmax
    pyrat.ret.prior    = args.prior
    pyrat.ret.priorlow = args.priorlow
    pyrat.ret.priorup = args.priorup

    pyrat.ret.qcap = args.get_default(
        'qcap', 'Metals volume-mixing-ratio cap', 1.0, gt=0.0, le=1.0)
    pyrat.ret.params = args.params
    pyrat.ret.pstep = args.pstep
    pyrat.ret.tlow = args.get_default(
        'tlow', 'Retrieval low-temperature (K) bound', 0,
        wflag=(runmode=='mcmc'))
    pyrat.ret.thigh = args.get_default(
        'thigh', 'Retrieval high-temperature (K) bound', np.inf,
        wflag=(runmode=='mcmc'))
    pyrat.ret.sampler = args.sampler
    pyrat.ret.nsamples = args.get_default(
        'nsamples', 'Number of MCMC samples', gt=0)
    pyrat.ret.burnin = args.get_default(
        'burnin', 'Number of burn-in samples per chain', gt=0)
    pyrat.ret.thinning = args.get_default(
        'thinning', 'MCMC posterior thinning', 1, ge=1)
    pyrat.ret.nchains = args.get_default(
        'nchains', 'Number of MCMC parallel chains', ge=1)
    pyrat.ret.grbreak = args.get_default(
        'grbreak', 'Gelman-Rubin convergence criteria', 0.0, ge=0)
    pyrat.ret.grnmin = args.get_default(
        'grnmin', 'Gelman-Rubin convergence fraction', 0.5, gt=0.0)

    atm.molmodel = args.get_choice(
        'molmodel', 'molecular-abundance model', pc.molmodels)
    atm.molfree = args.molfree
    atm.molpars = args.molpars
    atm.bulk = args.bulk
    atm.tmodelname = args.get_choice('tmodel', 'temperature model', pc.tmodels)
    atm.tpars = args.tpars
    pyrat.ncpu = args.get_default('ncpu', 'Number of processors', 1, ge=1)

    return
