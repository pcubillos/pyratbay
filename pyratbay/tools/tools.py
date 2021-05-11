# Copyright (c) 2021 Patricio Cubillos
# Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

__all__ = [
    'log_error',
    'cd',
    'tmp_reset',
    'binsearch',
    'divisors',
    'unpack',
    'u',
    'get_param',
    'ifirst', 'ilast',
    'isfile',
    'file_exists',
    'path',
    'Formatted_Write',
    'make_tea',
    'Timer',
    'get_exomol_mol',
    'cia_hitran', 'cia_borysow',
    'radius_to_depth',
    'depth_to_radius',
    'ignore_system_exit',
    ]


import os
import re
import struct
import time
import numbers
import string
import textwrap
import itertools
import functools
from collections import Iterable
from contextlib import contextmanager
import configparser

import numpy as np
import mc3.utils as mu

from .. import constants as pc
from .. import io as io
from ..lib import _indices


@contextmanager
def log_error(log=None, error=None):
    """Capture exceptions into a log.error() call."""
    try:
        yield
    except Exception as e:
        if log is None:
            log = mu.Log(logname=None, verb=1, width=80)
        if error is None:
            error = str(e)
        log.error(error, tracklev=-4)


@contextmanager
def cd(newdir):
    """
    Context manager for changing the current working directory.
    Taken from here: https://stackoverflow.com/questions/431684/
    """
    olddir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(olddir)


def recursive_setattr(obj, attr, val):
    """Recursive setattr, see https://stackoverflow.com/questions/31174295"""
    pre, _, post = attr.rpartition('.')
    return setattr(recursive_getattr(obj, pre) if pre else obj, post, val)


def recursive_getattr(obj, attr):
    """Recursive getattr, see https://stackoverflow.com/questions/31174295"""
    def _getattr(obj, attr):
        return getattr(obj, attr)
    return functools.reduce(_getattr, [obj] + attr.split('.'))


@contextmanager
def tmp_reset(obj, *attrs, **tmp_attrs):
    """
    Temporarily remove attributes from an object.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> o   = type('obj', (object,), {'x':1.0, 'y':2.0})
    >>> obj = type('obj', (object,), {'z':3.0, 'w':4.0, 'o':o})
    >>> # All listed arguments are set to None:
    >>> with pt.tmp_reset(obj, 'o.x', 'z'):
    >>>     print(obj.o.x, obj.o.y, obj.z, obj.w)
    (None, 2.0, None, 4.0)
    >>> # Keyword arguments can be set to a value, but cannot be recursive:
    >>> with pt.tmp_reset(obj, 'o.x', z=10):
    >>>     print(obj.o.x, obj.o.y, obj.z, obj.w)
    (None, 2.0, 10, 4.0)
    """
    orig_attrs = {}
    for attr in attrs:
        orig_attrs[attr] = recursive_getattr(obj, attr)
        recursive_setattr(obj, attr, None)
    for attr, tmp_val in tmp_attrs.items():
        orig_attrs[attr] = recursive_getattr(obj, attr)
        recursive_setattr(obj, attr, tmp_val)
    yield

    for attr, orig_val in orig_attrs.items():
        recursive_setattr(obj, attr, orig_val)


def binsearch(tli, wnumber, rec0, nrec, upper=True):
    r"""
    Do a binary+linear search in TLI dbfile for record with wavenumber
    immediately less equal to wnumber (if upper is True), or greater
    equal to wnumber (if upper) is False (considering duplicate values
    in tli file).

    Parameters
    ----------
    tli: File object
        TLI file where to search.
    wnumber: Scalar
        Target wavenumber in cm-1.
    rec0: Integer
        File position of first wavenumber record.
    nrec: Integer
        Number of wavenumber records.
    upper: Boolean
        If True, consider wnumber as an upper boundary. If False,
        consider wnumber as a lower boundary.

    Returns
    -------
    irec: Integer
        Index of record nearest to target. Return -1 if out of bounds.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> import struct
    >>> # Mock a TLI file:
    >>> wn = [0.0, 1.0, 1.0, 1.0, 2.0, 2.0]
    >>> with open('tli_demo.dat', 'wb') as tli:
    >>>     tli.write(struct.pack(str(len(wn))+"d", *wn))
    >>> # Now do bin searches for upper and lower boundaries:
    >>> with open('tli_demo.dat', 'rb') as tli:
    >>>     bs_lower = [pt.binsearch(tli, target, 0, len(wn), upper=False)
    >>>                 for target in [-1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5]]
    >>>     bs_upper = [pt.binsearch(tli, target, 0, len(wn), upper=True)
    >>>                 for target in [-1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5]]
    >>> print(bs_lower, bs_upper, sep='\n')
    [0, 0, 1, 1, 4, 4, -1]
    [-1, 0, 0, 3, 3, 5, 5]
    """
    if nrec <= 0:
        raise ValueError('Requested binsearch over a zero a zero-sized array.')

    # Initialize indices and current record:
    irec = ilo = 0
    ihi = nrec - 1
    tli.seek(rec0, 0)
    current = first = struct.unpack('d', tli.read(8))[0]
    tli.seek(rec0 + ihi*pc.dreclen, 0)
    last = struct.unpack('d', tli.read(8))[0]

    # Out of bounds:
    if wnumber < first and upper:
        return -1
    if last < wnumber and not upper:
        return -1

    # Binary search:
    while ihi - ilo > 1:
        irec = (ihi + ilo) // 2
        tli.seek(rec0 + irec*pc.dreclen, 0)
        current = struct.unpack('d', tli.read(8))[0]

        if current > wnumber:
            ihi = irec
        else:
            ilo = irec

    # Linear search:
    if upper and current > wnumber:
        return irec - 1
    elif not upper and current < wnumber:
        return irec + 1
    elif upper:
        while current <= wnumber:
            irec += 1
            if irec > nrec-1:
                return nrec-1
            tli.seek(rec0 + irec*pc.dreclen, 0)
            current = struct.unpack('d', tli.read(8))[0]
        return irec - 1
    else:
        while current >= wnumber:
            irec -= 1
            if irec < 0:
                return 0
            tli.seek(rec0 + irec*pc.dreclen, 0)
            current = struct.unpack('d', tli.read(8))[0]
        return irec + 1


def divisors(number):
    """
    Find all the integer divisors of number.
    """
    divs = []
    for i in np.arange(1, number/2+1):
      if number % i == 0:
        divs.append(i)
    divs.append(number)
    return np.asarray(divs, np.int)


def unpack(file, n, dtype):
    r"""
    Wrapper for struct unpack.

    Parameters
    ----------
    file: File object
        File object to read from.
    n: Integer
        Number of elements to read from file.
    dtype: String
        Data type of the bytes read.

    Returns
    -------
    output: Scalar, tuple, or string
        If dtype is 's' return the string (decoded as UTF-8).
        If there is a single element to read, return the scalar value.
        Else, return a tuple with the elements read.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> import struct
    >>> import numpy as np
    >>> # Store a string and numbers in a binary file:
    >>> with open('delete_me.dat', 'wb') as bfile:
    >>>     bfile.write(struct.pack('3s', 'H2O'.encode('utf-8')))
    >>>     bfile.write(struct.pack('h', 3))
    >>>     bfile.write(struct.pack('3f', np.pi, np.e, np.inf))

    >>> # Unpack them:
    >>> with open('delete_me.dat', 'rb') as bfile:
    >>>     string = pt.unpack(bfile, 3, 's')
    >>>     number = pt.unpack(bfile, 1, 'h')
    >>>     values = pt.unpack(bfile, 3, 'f')

    >>> # See outputs:
    >>> print(string, number, values, sep='\n')
    H2O
    3
    (3.1415927410125732, 2.7182817459106445, inf)
    """
    # Calculate number of bytes and read:
    size = struct.calcsize(f'{n}{dtype}')
    output = struct.unpack(f'{n}{dtype}', file.read(size))

    if dtype == 's':
        return output[0].decode('utf-8')
    elif n == 1:
        return output[0]
    return output


def u(units):
    """
    Get the conversion factor (to the CGS system) for units.

    Parameters
    ----------
    units: String
        Name of units.

    Returns
    -------
    value: Float
        Value of input units in CGS units.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> for units in ['cm', 'm', 'rearth', 'rjup', 'au']:
    >>>     print(f'{units} = {pt.u(units)} cm')
    cm = 1.0 cm
    m = 100.0 cm
    rearth = 637810000.0 cm
    rjup = 7149200000.0 cm
    au = 14959787069100.0 cm
    """
    # Accept only valid units:
    if not hasattr(pc, units):
        raise ValueError(
            f"Units '{units}' does not exist in pyratbay.constants.")
    return getattr(pc, units)


def get_param(param, units='none', gt=None, ge=None):
    """
    Read a parameter that may or may not have units.
    If it doesn't, default to the 'units' input argument.

    Parameters
    ----------
    param: String, Float, integer, or ndarray
        The parameter value (which may contain the units).
    units: String
        The default units for the parameter.
    gt: Float
        If not None, check output is greater than gt.
    ge: Float
        If not None, check output is greater-equal than gt.

    Returns
    -------
    value: Float or integer

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> # One meter in cm:
    >>> pt.get_param('1.0 m')
    100.0

    >>> # Alternatively, specify units in second argument:
    >>> pt.get_param(1.0, 'm')
    100.0

    >>> # Units in 'param' take precedence over 'unit':
    >>> pt.get_param('1.0 m', 'km')
    100.0

    >>> # Request returned value to be positive:
    >>> pt.get_param('-30.0 kelvin', gt=0.0)
    ValueError: Value -30.0 must be > 0.0.
    """
    if param is None:
        return None

    # Split the parameter if it has a white-space:
    if isinstance(param, str):
        par = param.split()
        if len(par) > 2:
            raise ValueError(f"Invalid value '{param}'")
        if len(par) == 2:
            units = par[1]
            if not hasattr(pc, units):
                raise ValueError(f"Invalid units for value '{param}'")
        try:
            value = np.float(par[0])
        except:
            raise ValueError(f"Invalid value '{param}'")
    else:
        value = param

    # Use given units:
    if isinstance(param, (numbers.Number, np.ndarray)) \
            or (isinstance(param, str) and len(par) == 1):
        if units is None or not hasattr(pc, units):
            raise ValueError(f"Invalid units '{units}'")

    # Apply the units:
    value *= u(units)

    if gt is not None and value <= gt:
        raise ValueError(f'Value {value} must be > {gt}')
    if ge is not None and value < ge:
        raise ValueError(f'Value {value} must be >= {ge}')

    return value


def ifirst(data, default_ret=-1):
    """
    Get the first index where data is True or 1.

    Parameters
    ----------
    data: 1D bool/integer iterable
        An array of bools or integers.
    default_ret: Integer
        Default returned value when no value in data is True or 1.

    Returns
    -------
    first: integer
       First index where data == True or 1.  Return default_ret otherwise.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> import numpy as np
    >>> print(pt.ifirst([1,0,0]))
    0
    >>> print(pt.ifirst(np.arange(5)>2.5))
    3
    >>> print(pt.ifirst([False, True, True]))
    1
    >>> print(pt.ifirst([False, False, False]))
    -1
    >>> print(pt.ifirst([False, False, False], default_ret=0))
    0
    """
    return _indices.ifirst(np.asarray(data, np.int), default_ret)


def ilast(data, default_ret=-1):
    """
    Get the last index where data is 1 or True.

    Parameters
    ----------
    data: 1D bool/integer iterable
        An array of bools or integers.
    default_ret: Integer
        Default returned value when no value in data is True or 1.

    Returns
    -------
    last: integer
       Last index where data == 1 or True.  Return default_ret otherwise.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> import numpy as np
    >>> print(pt.ilast([1,0,0]))
    0
    >>> print(pt.ilast(np.arange(5)<2.5))
    2
    >>> print(pt.ilast([False, True, True]))
    2
    >>> print(pt.ilast([False, False, False]))
    -1
    >>> print(pt.ilast([False, False, False], default_ret=0))
    0
    """
    return _indices.ilast(np.asarray(data, np.int), default_ret)


def isfile(path):
    """
    Check whether a path (or list of paths) is a regular file.

    Parameters
    ----------
    path:  String or list
        Path(s) to check.

    Returns
    -------
    status: Integer
        If path is None, return -1.
        If any path is not a regular file, return 0.
        If all paths are a regular file, return 1.

    Examples (for Python 2.7, import from pathlib2)
    --------
    >>> import pyratbay.tools as pt
    >>> from pathlib import Path
    >>> # Mock couple files:
    >>> file1, file2 = './tmp_file1.deleteme', './tmp_file2.deleteme'
    >>> Path(file1).touch()
    >>> Path(file2).touch()
    >>> # Input is None:
    >>> print(pt.isfile(None))
    -1
    >>> # All input files exist:
    >>> print(pt.isfile(file1))
    1
    >>> print(pt.isfile([file1]))
    1
    >>> print(pt.isfile([file1, file2]))
    1
    >>> # At least one input does not exist:
    >>> print(pt.isfile('nofile'))
    0
    >>> print(pt.isfile(['nofile']))
    0
    >>> print(pt.isfile([file1, 'nofile']))
    0
    """
    # None exception:
    if path is None:
        return -1

    if isinstance(path, str):
        paths = [path]
    else:
        paths = path

    # Regular file or not:
    return int(all(os.path.isfile(path) for path in paths))


def file_exists(pname, desc, value):
    """
    Check that a file or list of files (value) exist.  If not None
    and file(s) do not exist, raise a ValueError.

    Parameters
    ----------
    pname: String
        Parameter name.
    desc: String
        Parameter description.
    value: String or list of strings
        File path(s) to check.

    Examples (for Python 2.7, import from pathlib2)
    --------
    >>> import pyratbay.tools as pt
    >>> from pathlib import Path
    >>> # None is OK:
    >>> pt.file_exists('none', 'None input', None)
    >>> # Create a file, check it exists:
    >>> Path('./new_tmp_file.dat').touch()
    >>> pt.file_exists('testfile', 'Test', 'new_tmp_file.dat')
    >>> # Non-existing file throws error:
    >>> pt.file_exists('testfile', 'Test', 'no_file.dat')
    ValueError: Test file (testfile) does not exist: 'no_file.dat'
    """
    if value is None:
        return

    if isinstance(value, str):
        values = [value]
    else:
        values = value

    for value in values:
        if not os.path.isfile(value):
            raise ValueError(f"{desc} file ({pname}) does not exist: '{value}'")


def path(filename):
    """
    Ensure file names have non-null path

    Parameters
    ----------
    filename: String
        A file name.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> print(pt.path('file.txt'))
    ./file.txt
    >>> print(pt.path('./file.txt'))
    ./file.txt
    >>> print(pt.path('/home/user/file.txt'))
    /home/user/file.txt
    """
    if filename is None:
        return None
    path, filename = os.path.split(filename)
    if path == '':
        path = '.'
    return f'{path}/{filename}'


class Formatted_Write(string.Formatter):
    """
    Write (and keep) formatted, wrapped text to string.

    Following PEP3101, this class subclasses Formatter to handle
    None when a specific format is set.

    Examples
    --------
    >>> import numpy as np
    >>> import pyratbay.tools as pt
    >>> fmt = pt.Formatted_Write()
    >>> rstar = np.pi/3.14
    >>> fmt.write('Stellar radius (rstar, rsun):  {:.2f}', rstar)
    >>> fmt.write('Stellar radius (rstar, rsun):  {:.2f}', None)
    >>> fmt.write('Stellar radius (rstar, rsun):  {}',     rstar)
    >>> fmt.write('Stellar radius (rstar, rsun):  {}',     None)
    >>> print(fmt.text)
    Stellar radius (rstar, rsun):  1.00
    Stellar radius (rstar, rsun):  None
    Stellar radius (rstar, rsun):  1.0005072145190423
    Stellar radius (rstar, rsun):  None
    """
    def __init__(self, indent=0, si=4, fmt=None, edge=None, lw=80, prec=None):
        """
        Parameters
        ----------
        indent: Integer
            Number of blanks for indentation in first line.
        si: Integer
            Number of blanks for indentation in subsequent lines.
        fmt: dict of callables.
            Default formatting for numpy arrays (as in formatting in
            np.printoptions).
        edge: Integer
            Default number of array items in summary at beginning/end
            (as in edgeitems in np.printoptions).
        lw: Integer
            Default number of characters per line (as in linewidth in
            np.printoptions).
        prec: Integer
            Default precision for floating point values (as in precision
            in np.printoptions).
        """
        self.text = ''
        self.indent = indent
        self.si = si
        # Numpy formatting:
        self.fmt  = fmt
        self.edge = edge
        self.lw   = lw
        self.prec = prec

    def format_field(self, value, spec):
        if value is None:
            return 'None'
        return super(Formatted_Write, self).format_field(value, spec)

    def write(self, text, *format, **numpy_fmt):
        """
        Write formatted text.
        See __init__ arguments for avaiable numpy_fmt items.
        """
        numpy_fmt.setdefault('fmt',  self.fmt)
        numpy_fmt.setdefault('edge', self.edge)
        numpy_fmt.setdefault('lw',   self.lw)
        numpy_fmt.setdefault('prec', self.prec)

        if numpy_fmt['edge'] is None:
            threshold = None
        else:
            threshold = 2*numpy_fmt['edge']

        fmt = {
            'formatter': numpy_fmt['fmt'],
            'edgeitems': numpy_fmt['edge'],
            'threshold': threshold,
            'linewidth': numpy_fmt['lw'],
            'precision': numpy_fmt['prec'],
            }
        with np.printoptions(**fmt):
            text = super(Formatted_Write, self).format(text, *format)

        indspace = ' '*self.indent
        if self.si is None:
            sindspace = indspace
        else:
            sindspace = ' '*self.si

        for line in text.splitlines():
            self.text += textwrap.fill(
                line, break_long_words=False, initial_indent=indspace,
                subsequent_indent=sindspace, width=80)
            self.text += '\n'


def make_tea(maxiter=200, savefiles=False, times=False, location_TEA=None,
    abun_file=None, location_out='./TEA', verb=1, ncpu=1):
    """
    Make a TEA configuration file.

    Parameters
    ----------
    TBD
    """
    if location_TEA is None:
        location_TEA = os.path.realpath(pc.ROOT + 'pyratbay/TEA/')

    # Open new configparser:
    config = configparser.ConfigParser()
    config.add_section('TEA')
    config.set('TEA', 'maxiter',      str(maxiter))
    config.set('TEA', 'savefiles',    str(savefiles))
    config.set('TEA', 'times',        str(times))
    config.set('TEA', 'location_TEA', str(location_TEA))
    config.set('TEA', 'location_out', str(location_out))
    config.set('TEA', 'abun_file',    str(abun_file))
    config.set('TEA', 'verb',         str(verb))
    config.set('TEA', 'ncpu',         str(ncpu))

    # For completeness:
    config.add_section('PRE-ATM')
    config.set('PRE-ATM', 'PT_file',        'None')
    config.set('PRE-ATM', 'pre_atm_name',   'None')
    config.set('PRE-ATM', 'input_elem',     'None')
    config.set('PRE-ATM', 'output_species', 'None')

    # Write TEA configuration file:
    with open('TEA.cfg', 'w') as config_file:
        config.write(config_file)


class Timer(object):
    """
    Timer to get the time (in seconds) since the last call.
    """
    def __init__(self):
        self.t0 = time.time()

    def clock(self):
        tnew = time.time()
        delta = tnew - self.t0
        self.t0 = tnew
        return delta


def get_exomol_mol(dbfile):
    """
    Parse an exomol file to extract the molecule and isotope name.

    Parameters
    ----------
    dbfile: String
        An exomol line-list file (must follow ExoMol naming convention).

    Returns
    -------
    molecule: String
        Name of the molecule.
    isotope: String
        Name of the isotope (See Tennyson et al. 2016, jmosp, 327).

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> filenames = [
    >>>     '1H2-16O__POKAZATEL__00400-00500.trans.bz2',
    >>>     '1H-2H-16O__VTT__00250-00500.trans.bz2',
    >>>     '12C-16O2__HITEMP.pf',
    >>>     '12C-16O-18O__Zak.par',
    >>>     '12C-1H4__YT10to10__01100-01200.trans.bz2',
    >>>     '12C-1H3-2H__MockName__01100-01200.trans.bz2'
    >>>    ]
    >>> for db in filenames:
    >>>     print(pt.get_exomol_mol(db))
    ('H2O', '116')
    ('H2O', '126')
    ('CO2', '266')
    ('CO2', '268')
    ('CH4', '21111')
    ('CH4', '21112')
    """
    atoms = os.path.split(dbfile)[1].split('_')[0].split('-')
    elements = []
    isotope  = ''
    for atom in atoms:
        match = re.match(r"([0-9]+)([a-z]+)([0-9]*)", atom, re.I)
        N = 1 if match.group(3) == '' else int(match.group(3))
        elements += N * [match.group(2)]
        isotope  += match.group(1)[-1:] * N

    composition = [list(g[1]) for g in itertools.groupby(elements)]
    molecule = ''.join([
        c[0] + str(len(c))*(len(c)>1)
        for c in composition])

    return molecule, isotope


def cia_hitran(ciafile, tstep=1, wstep=1):
    """
    Re-write a HITRAN CIA file into Pyrat Bay format.
    See Richard et al. (2012) and https://www.cfa.harvard.edu/HITRAN/

    Parameters
    ----------
    ciafile: String
        A HITRAN CIA file.
    tstep: Integer
        Slicing step size along temperature dimension.
    wstep: Integer
        Slicing step size along wavenumber dimension.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> # Before moving on, download a HITRAN CIA files from the link above.
    >>> ciafile = 'H2-H2_2011.cia'
    >>> pt.cia_hitran(ciafile, tstep=2, wstep=10)
    """
    # Extract CS data:
    with open(ciafile, 'r') as f:
        info = f.readline().strip().split()
        species = info[0].split('-')
        temps, data, wave = [], [], []
        wnmin, wnmax = -1, -1
        f.seek(0)
        for line in f:
            if line.strip().startswith('-'.join(species)):
                info = line.strip().split()
                # if wn ranges differ, trigger new set
                if float(info[1]) != wnmin or float(info[2]) != wnmax:
                    wnmin = float(info[1])
                    wnmax = float(info[2])
                temp = float(info[4])
                nwave = int  (info[3])
                wn = np.zeros(nwave, np.double)
                cs = np.zeros(nwave, np.double)
                i = 0
                continue
            # else, read in opacities
            wn[i], cs[i] = line.split()[0:2]
            i += 1
            if i == nwave:
                temps.append(temp)
                # Thin the arrays in wavenumber if requested:
                data.append(cs[::wstep])
                wave.append(wn[::wstep])

    # Identify sets:
    temps = np.array(temps)
    ntemps = len(temps)
    i = 0
    while i < ntemps:
        wn = wave[i]
        j = i
        while j < ntemps and len(wave[j])==len(wn) and np.all(wave[j]-wn==0):
            j += 1
        temp = temps[i:j:tstep]
        # Set cm-1 amagat-2 units:
        cs = np.array(data[i:j])[::tstep] * pc.amagat**2

        pair = '-'.join(species)
        wl_min = 1.0/(wn[-1]*pc.um)
        wl_max = 1.0/(wn[0]*pc.um)

        csfile = (
            f'CIA_HITRAN_{pair}_{wl_min:.1f}-{wl_max:.1f}um_'
            f'{temp[0]:04.0f}-{temp[-1]:04.0f}K.dat')
        header = (
            f'# This file contains the reformated {pair} CIA data from\n'
            f'# HITRAN file: {ciafile}\n\n')
        io.write_cs(csfile, cs, species, temp, wn, header)
        i = j


def cia_borysow(ciafile, species1, species2):
    """
    Re-write a Borysow CIA file into Pyrat Bay format.
    See http://www.astro.ku.dk/~aborysow/programs/

    Parameters
    ----------
    ciafile: String
        A HITRAN CIA file.
    species1: String
        First CIA species.
    species2: String
        Second CIA species.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> # Before moving on, download a HITRAN CIA files from the link above.
    >>> ciafile = 'ciah2he_dh_quantmech'
    >>> pt.cia_borysow(ciafile, 'H2', 'He')
    """
    data = np.loadtxt(ciafile, skiprows=3)
    wn = data[:,0]
    cs = data[:,1:].T

    with open(ciafile) as f:
        line = f.readline()
        temp = f.readline().split()[1:]
    temp = [float(t.replace('K','')) for t in temp]

    species = [species1, species2]

    pair = '-'.join(species)
    wl_min = 1.0/(wn[-1]*pc.um)
    wl_max = 1.0/(wn[0]*pc.um)
    file_name = os.path.basename(ciafile)

    csfile = (
        f'CIA_Borysow_{pair}_{wl_min:.1f}-{wl_max:.1f}um_'
        f'{temp[0]:04.0f}-{temp[-1]:04.0f}K.dat')
    header = (
        f'# This file contains the reformated {pair} CIA data from:\n'
        f'# http://www.astro.ku.dk/~aborysow/programs/{file_name}\n\n')
    io.write_cs(csfile, cs, species, temp, wn, header)


def radius_to_depth(rprs, rprs_err):
    r"""
    Compute transit depth (and uncertainties) from input
    planet=to-star radius-ratio, with error propagation.

    Parameters
    ----------
    rprs: Float or float iterable
        Planet-to-star radius ratio.
    rprs_err: Float or float iterable
        Uncertainties of the radius ratios.

    Returns
    -------
    depth: Float or float ndarray
        Transit depth for given radius ratio.
    depth_err: Float or float ndarray
        Uncertainties of the transit depth.

    Examples
    --------
    >>> import numpy as np
    >>> import pyratbay.tools as pt
    >>> rprs = 1.2
    >>> rprs_err = 0.25
    >>> depth, depth_err = pt.radius_to_depth(rprs, rprs_err)
    >>> print(f'Depth = {depth} +/- {depth_err}')
    Depth = 1.44 +/- 0.6

    >>> rprs = [1.2, 1.5]
    >>> rprs_err = [0.25, 0.3]
    >>> depth, depth_err = pt.radius_to_depth(rprs, rprs_err)
    >>> print('Depth    Uncert\n' +
    >>>     '\n'.join([f'{d} +/- {de:.1f}' for d,de in zip(depth, depth_err)]))
    Depth    Uncert
    1.44 +/- 0.6
    2.25 +/- 0.9
    """
    if not isinstance(rprs, Iterable):
        pass
    elif not isinstance(rprs, np.ndarray):
        rprs = np.array(rprs)
        rprs_err = np.array(rprs_err)

    depth = rprs**2.0
    depth_err = 2.0 * rprs * rprs_err

    return depth, depth_err


def depth_to_radius(depth, depth_err):
    r"""
    Compute planet-to-star radius ratio (and uncertainties) from
    input transit depth, with error propagation.

    Parameters
    ----------
    depth: Float or float iterable
        Transit depth.
    depth_err: Float or float iterable
        Uncertainties of the transit depth.

    Returns
    -------
    rprs: Float or float ndarray
        Planet-to-star radius ratio.
    rprs_err: Float or float ndarray
        Uncertainties of the radius ratio rprs.

    Examples
    --------
    >>> import numpy as np
    >>> import pyratbay.tools as pt
    >>> depth = 1.44
    >>> depth_err = 0.6
    >>> rprs, rprs_err = pt.depth_to_radius(depth, depth_err)
    >>> print(f'Rp/Rs = {rprs} +/- {rprs_err}')
    Rp/Rs = 1.2 +/- 0.25

    >>> depth = [1.44, 2.25]
    >>> depth_err = [0.6, 0.9]
    >>> rprs, rprs_err = pt.depth_to_radius(depth, depth_err)
    >>> print('Rp/Rs   Uncert\n'
    >>>     + '\n'.join([f'{r} +/- {re}' for r,re in zip(rprs, rprs_err)]))
    Rp/Rs   Uncert
    1.2 +/- 0.25
    1.5 +/- 0.3
    """
    if not isinstance(depth, Iterable):
        pass
    elif not isinstance(depth, np.ndarray):
        depth = np.array(depth)
        depth_err = np.array(depth_err)

    rprs = np.sqrt(depth)
    rprs_err = 0.5 * depth_err / rprs
    return rprs, rprs_err


def ignore_system_exit(func):
    """Decorator to ignore SystemExit exceptions."""
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except SystemExit:
            return None
    return new_func

