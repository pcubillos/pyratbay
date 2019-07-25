.. _api:

API
===


pyratbay
________


.. py:module:: pyratbay

.. py:class:: Pyrat(cfile)

.. code-block:: pycon

    Main Pyrat object.

.. code-block:: pycon

    Parse the command-line arguments into the pyrat object.

    Parameters
    ----------
    args: Namespace
        Object storing user-input attributes to initialize Pyrat.
    log: Log object
        An MCcubed.utils.Log instance to log screen outputs to file.

    Examples
    --------
    >>> import pyratbay as pb
    >>> # There are two ways to initialize a Pyrat object:
    >>> # Initialize only:
    >>> pyrat = pb.init('spectrum_transmission.cfg')
    >>> # This is equivalent to:
    >>> args, log = pb.tools.parse('spectrum_transmission.cfg')
    >>> pyrat = pb.Pyrat(args, log)

    >>> # Initialize and compute a spectrum:
    >>> pyrat = pb.run('spectrum_transmission.cfg')

.. py:function:: run(cfile, init=False)
.. code-block:: pycon

    Pyrat Bay initialization driver.

    Parameters
    ----------
    cfile: String
        A Pyrat Bay configuration file.
    init: Bool
        If True, only initialize a Pyrat object (no spectra calculation).
        This is useful when computing spectra interactively.


pyratbay.constants
__________________


.. py:module:: pyratbay.constants

.. py:data:: h
.. code-block:: pycon

  6.62607004e-27

.. py:data:: k
.. code-block:: pycon

  1.38064852e-16

.. py:data:: c
.. code-block:: pycon

  29979245800.0

.. py:data:: G
.. code-block:: pycon

  6.67408e-08

.. py:data:: eV
.. code-block:: pycon

  8065.49179

.. py:data:: A
.. code-block:: pycon

  1e-08

.. py:data:: nm
.. code-block:: pycon

  1e-07

.. py:data:: um
.. code-block:: pycon

  0.0001

.. py:data:: mm
.. code-block:: pycon

  0.1

.. py:data:: cm
.. code-block:: pycon

  1.0

.. py:data:: m
.. code-block:: pycon

  100.0

.. py:data:: km
.. code-block:: pycon

  100000.0

.. py:data:: au
.. code-block:: pycon

  14959787069100.0

.. py:data:: pc
.. code-block:: pycon

  3.085677581305729e+18

.. py:data:: rearth
.. code-block:: pycon

  637810000.0

.. py:data:: rjup
.. code-block:: pycon

  7149200000.0

.. py:data:: rsun
.. code-block:: pycon

  69570000000.0

.. py:data:: barye
.. code-block:: pycon

  1.0

.. py:data:: mbar
.. code-block:: pycon

  1000.0

.. py:data:: pascal
.. code-block:: pycon

  10.0

.. py:data:: bar
.. code-block:: pycon

  1000000.0

.. py:data:: atm
.. code-block:: pycon

  1010000.0

.. py:data:: gram
.. code-block:: pycon

  1.0

.. py:data:: kg
.. code-block:: pycon

  1000.0

.. py:data:: mearth
.. code-block:: pycon

  5.9724e+27

.. py:data:: mjup
.. code-block:: pycon

  1.8982e+30

.. py:data:: msun
.. code-block:: pycon

  1.9885e+33

.. py:data:: amu
.. code-block:: pycon

  1.66053904e-24

.. py:data:: me
.. code-block:: pycon

  9.10938356e-28

.. py:data:: kelvin
.. code-block:: pycon

  1.0

.. py:data:: amagat
.. code-block:: pycon

  2.6867811e+19

.. py:data:: e
.. code-block:: pycon

  4.803205e-10

.. py:data:: percent
.. code-block:: pycon

  0.01

.. py:data:: ppt
.. code-block:: pycon

  0.001

.. py:data:: ppm
.. code-block:: pycon

  1e-06

.. py:data:: none
.. code-block:: pycon

  1

.. py:data:: C1
.. code-block:: pycon

  1129583489393.1277

.. py:data:: C2
.. code-block:: pycon

  1.4387773538277204

.. py:data:: C3
.. code-block:: pycon

  8.852821819282496e-13

.. py:data:: tlireclen
.. code-block:: pycon

  26

.. py:data:: dreclen
.. code-block:: pycon

  8

.. py:data:: ireclen
.. code-block:: pycon

  4

.. py:data:: sreclen
.. code-block:: pycon

  2

.. py:data:: ROOT
.. code-block:: pycon

  '/Users/pato/Dropbox/IWF/projects/2014_pyratbay/2016-01-08_develop/pyratbay/'

.. py:data:: dbases
.. code-block:: pycon

  ['hitran', 'exomol', 'repack', 'pands', 'tioschwenke', 'voplez', 'vald']

.. py:data:: rmodes
.. code-block:: pycon

  ['tli', 'pt', 'atmosphere', 'opacity', 'spectrum', 'mcmc']

.. py:data:: retflags
.. code-block:: pycon

  ['temp', 'rad', 'mol', 'ray', 'cloud', 'patchy']

.. py:data:: tmodels
.. code-block:: pycon

  ['isothermal', 'tcea', 'madhu_inv', 'madhu_noinv']

.. py:data:: molmodels
.. code-block:: pycon

  ['vert', 'scale']

.. py:data:: amodels
.. code-block:: pycon

  ['sodium_vdw', 'potassium_vdw']

.. py:data:: rmodels
.. code-block:: pycon

  ['dalgarno_H', 'dalgarno_H2', 'dalgarno_He', 'lecavelier']

.. py:data:: cmodels
.. code-block:: pycon

  ['deck', 'ccsgray']


pyratbay.io
___________


.. py:module:: pyratbay.io

.. py:function:: save_pyrat(pyrat, pfile=None)
.. code-block:: pycon

    Save a pyrat instance into a pickle file.

    Parameters
    ----------
    pyrat: A Pyrat instance
        Object to save.
    pfile: String
        Name of output file.  Default to the pyrat logname (changing
        the extension to '.pickle').

.. py:function:: load_pyrat(pfile)
.. code-block:: pycon

    Load a pyrat instance from a pickle file.

    Parameters
    ----------
    pfile: String
        Name of input pickle file.

    Returns
    -------
    pyrat: A Pyrat instance
        Loaded object.

.. py:function:: write_spectrum(wl, spectrum, filename, type, wlunits='um')
.. code-block:: pycon

    Write a Pyrat spectrum to file.

    Parameters
    ----------
    wl: 1D float iterable
        Wavelength array in cm units.
    spectrum: 1D float iterable
        Spectrum array. (rp/rs)**2 for transmission (unitless),
        planetary flux for emission (erg s-1 cm-2 cm units).
    filename: String
        Output file name.
    type: String
        Data type:
        'transit' for transmission,
        'eclipse' for emission,
        'filter' for a instrumental filter transmission.
    wlunits: String
        Output units for wavelength.

    Examples
    --------
    >>> # See read_spectrum() examples.

.. py:function:: read_spectrum(filename, wn=True)
.. code-block:: pycon

    Read a Pyrat spectrum file, a plain text file with two-columns: the
    wavelength and signal.  If wn is true, this function converts
    wavelength to wavenumber in cm-1.  The very last comment line sets
    the wavelength units (the first string following a blank, e.g., the
    string '# um' sets the wavelength units as microns).
    If the units are not defined, assume wavelength units are microns.

    Parameters
    ----------
    filename: String
       Path to output Transit spectrum file to read.
    wn: Boolean
       If True convert wavelength to wavenumber.

    Return
    ------
    wave: 1D float ndarray
       The spectrum's wavenumber (in cm units) or wavelength array (in
       the input file's units).
    spectrum: 1D float ndarray
       The spectrum in the input file.

    Examples
    --------
    >>> import pyratbay.io as io
    >>> # Write a spectrum to file:
    >>> nwave = 7
    >>> wl = np.linspace(1.1, 1.7, nwave) * 1e-4
    >>> spectrum = np.ones(nwave)
    >>> io.write_spectrum(wl, spectrum,
    >>>     filename='sample_spectrum.dat', type='transit', wlunits='um')
    >>> # Take a look at the output file:
    >>> with open('sample_spectrum.dat', 'r') as f:
    >>>     print("".join(f.readlines()))
    # Wavelength        (Rp/Rs)**2
    #         um          unitless
         1.10000   1.000000000e+00
         1.20000   1.000000000e+00
         1.30000   1.000000000e+00
         1.40000   1.000000000e+00
         1.50000   1.000000000e+00
         1.60000   1.000000000e+00
         1.70000   1.000000000e+00
    >>> # Now, read from file (getting wavenumber array):
    >>> wn, flux = io.read_spectrum('sample_spectrum.dat')
    >>> print(wn)
    [9090.90909091 8333.33333333 7692.30769231 7142.85714286 6666.66666667
     6250.         5882.35294118]
    >>> print(flux)
    [1. 1. 1. 1. 1. 1. 1.]
    >>> # Read from file (getting wavelength array):
    >>> wl, flux = io.read_spectrum('sample_spectrum.dat', wn=False)
    >>> print(wl)
    [1.1 1.2 1.3 1.4 1.5 1.6 1.7]
    >>> print(flux)
    [1. 1. 1. 1. 1. 1. 1.]

.. py:function:: write_opacity(ofile, molID, temp, press, wn, etable)
.. code-block:: pycon

    Write an opacity table as a binary file.

    Parameters
    ----------
    ofile: String
        Path to a Pyrat Bay opacity file.
    molID: 1D integer ndarray
        molecule ID.
    temp: 1D float ndarray
        Temperature (Kelvin degree).
    press: 1D float ndarray
        Pressure (barye).
    wn: 1D float ndarray
        Wavenumber (cm-1).
    etable: 4D float ndarray
        Tabulated opacities (cm-1).

.. py:function:: read_opacity(ofile)
.. code-block:: pycon

    Read an opacity table from file.

    Parameters
    ----------
    ofile: String
        Path to a Pyrat Bay opacity file.

    Returns
    -------
    sizes: 4-element integer tuple
        Sizes of the dimensions of the opacity table:
        (nmol, ntemp, nlayers, nwave)
    arrays: 4-element 1D ndarray tuple
        The dimensions of the opacity table:
        - molecule ID (integer, unitless, see inputs/molecules.dat)
        - temperature (float, Kelvin)
        - pressure    (float, barye)
        - wavenumber  (float, cm-1)
    etable: 4D float ndarray tuple
        The tabulated opacities (cm-1), of shape [nmol, ntemp, nlayers, nwave].

.. py:function:: write_pf(pffile, pf, isotopes, temp, header=None)
.. code-block:: pycon

    Write a partition-function file in Pyrat Bay format.

    Parameters
    ----------
    pffile: String
        Output partition-function file.
    pf: 2D float iterable
        Partition-function data (of shape [niso, ntemp]).
    isotopes: 1D string iterable
        Isotope names.
    temp: 1D float iterable
        Temperature array.
    header: String
        A header for the partition-function file (must be as comments).

    Examples
    --------
    >>> # See read_pf() examples.

.. py:function:: read_pf(pffile)
.. code-block:: pycon

    Read a partition-function file.

    Parameters
    ----------
    pffile: String
        Partition function file to read.

    Returns
    -------
    pf: 2D float ndarray
        The partition function data (of shape [niso, ntemp]).
    isotopes: List of strings
         The names of the tabulated isotopes.
    temp: 1D float ndarray
        Array with temperature sample.

    Examples
    --------
    >>> import pyratbay.io as io
    >>> # Generate some mock PF data and write to file:
    >>> pffile = 'PF_Exomol_NH3.dat'
    >>> isotopes = ['4111', '5111']
    >>> temp   = np.linspace(10,100,4)
    >>> pf     = np.array([np.logspace(0,3,4),
    >>>                    np.logspace(1,4,4)])
    >>> header = '# Mock partition function for NH3.\n'
    >>> io.write_pf(pffile, pf, isotopes, temp, header)

    >>> # Now, read it back:
    >>> pf, iso, temp = io.read_pf(pffile)
    >>> for item in [iso, temp, pf]:
    >>>     print(item)
    ['4111' '5111']
    [ 10.  40.  70. 100.]
    [[1.e+00 1.e+01 1.e+02 1.e+03]
     [1.e+01 1.e+02 1.e+03 1.e+04]]

.. py:function:: write_cs(csfile, cs, species, temp, wn, header=None)
.. code-block:: pycon

    Write a cross-section file in Pyrat Bay format.

    Parameters
    ----------
    csfile: String
        Output cross-section file.
    cs: 2D float iterable
        Cross-section opacity in units of cm-1 amagat^-N, with N the
        number of species, of shape [ntemp, nwave].
    species: 1D string iterable
        Species names.
    temp: 1D float iterable
        Temperature array in Kelvin degree.
    wn: 1D float iterable
        Wavenumber array in cm-1.
    header: String
        A header for the cross-section file (must be as comments).

    Examples
    --------
    >>> # See read_cs() examples.

.. py:function:: read_cs(csfile)
.. code-block:: pycon

    Read a cross-section file.

    Parameters
    ----------
    csfile: String
        Partition function file to read.

    Returns
    -------
    cs: 2D float ndarray
        Cross-section opacity in units of cm-1 amagat^-N, with N the
        number of species, of shape [ntemp, nwave].
    species: 1D string list
        Species names.
    temp: 1D float ndarray
        Temperature array in Kelvin degree.
    wn: 1D float ndarray
        Wavenumber array in cm-1.

    Examples
    --------
    >>> import pyratbay.io as io
    >>> # Generate some mock PF data and write to file:
    >>> csfile = 'CS_Mock-HITRAN_H2-H2.dat'
    >>> species = ['H2', 'H2']
    >>> temp = np.linspace(100, 1000, 3)
    >>> wn   = np.arange(10, 15, 1.0)
    >>> cs   = np.array([np.logspace( 0,-4,5),
    >>>                  np.logspace(-1,-5,5),
    >>>                  np.logspace(-2,-6,5)])
    >>> header = '# Mock cross-section for H2-H2.\n'
    >>> io.write_cs(csfile, cs, species, temp, wn, header)
    >>> # Now, read it back:
    >>> cs, species, temp, wn = io.read_cs(csfile)
    >>> for item in [species, temp, wn, cs]:
    >>>     print(item)
    ['H2', 'H2']
    [ 100.  550. 1000.]
    [10. 11. 12. 13. 14.]
    [[1.e+00 1.e-01 1.e-02 1.e-03 1.e-04]
     [1.e-01 1.e-02 1.e-03 1.e-04 1.e-05]
     [1.e-02 1.e-03 1.e-04 1.e-05 1.e-06]]

.. py:function:: read_pt(ptfile)
.. code-block:: pycon

    Read a pressure and temperature profile from a file.

    Parameters
    ----------
    ptfile: String
        Input file with pressure (in bars, first column) and temperature
        profiles (in Kelvin degree, second column).

    Returns
    -------
    pressure: 1D float ndarray
        Pressure profile in barye.
    temperature: 1D float ndarray
        Temperature profile in Kelvin.

    Examples
    --------
    >>> import pyratbay.io as io
    >>> ptfile = 'pt_profile.dat'
    >>> temp  = np.array([100.0, 150.0, 200.0, 175.0, 150.0])
    >>> press = np.array([1e-6,  1e-4,  1e-2,  1e0,   1e2])
    >>> with open(ptfile, 'w') as f:
    >>>     for p,t in zip(press, temp):
    >>>         f.write('{:.3e}  {:5.1f}\n'.format(p, t))
    >>> pressure, temperature = io.read_pt(ptfile)
    >>> for p,t in zip(pressure, temperature):
    >>>     print('{:.1e} barye  {:5.1f} K'.format(p, t))
    1.0e+00 barye  100.0 K
    1.0e+02 barye  150.0 K
    1.0e+04 barye  200.0 K
    1.0e+06 barye  175.0 K
    1.0e+08 barye  150.0 K


pyratbay.tools
______________


.. py:module:: pyratbay.tools

.. py:function:: log_error(log=None, error=None)
.. code-block:: pycon

    Capture exceptions into a log.error() call.

.. py:function:: tmp_reset(obj, *attrs, **tmp_attrs)
.. code-block:: pycon

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

.. py:function:: binsearch(tli, wnumber, rec0, nrec, upper=True)
.. code-block:: pycon

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

.. py:function:: divisors(number)
.. code-block:: pycon

    Find all the integer divisors of number.

.. py:function:: unpack(file, n, dtype)
.. code-block:: pycon

    Wrapper for struct unpack.

    Parameters
    ----------
    file: File object
       File object to read from.
    n: Integer
       Number of elements to read from file.
    dtype: String
       Data type of the bytes read.

    Returns:
    --------
    output: Scalar, tuple, or string
       If dtype is 's' return the string.
       If there is a single element to read, return the scalar value.
       Else, return a tuple with the elements read.

.. py:function:: u(units)
.. code-block:: pycon

    Get the conversion factor (to the CGS system) for units.

    Parameters
    ----------
    units: String
       Name of units

.. py:function:: get_param(pname, value, units, log=None, gt=None, ge=None, tracklev=-3)
.. code-block:: pycon

    Read a parameter that may have units.
    If it doesn't, default to the 'units' input argument.

    Parameters
    ----------
    pname: String
        The parameter name.
    value: String, Float, integer, or ndarray
        The parameter value (which may contain the units).
    units: String
        The default units for the parameter.
    log: Log object
        Screen-output log handler.
    gt: Float
        If not None, check output is greater than gt.
    ge: Float
        If not None, check output is greater-equal than gt.
    tracklev: Integer
        Error track level.

    Returns
    -------
    value: Float or integer

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> for line in ['One km in cm:',
    >>>              pt.get_param('size', 1.0, 'km'),
    >>>              "units in 'param' take precedence over 'unit':",
    >>>              pt.get_param('size', '10 cm', 'km')]
    >>>     print(line)
    One km in cm:
    100000.0
    units in 'param' take precedence over 'unit':
    10
    # Cast to integer:
    10

.. py:function:: ifirst(data, default_ret=-1)
.. code-block:: pycon

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

.. py:function:: ilast(data, default_ret=-1)
.. code-block:: pycon

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

.. py:function:: isfile(path)
.. code-block:: pycon

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

.. py:function:: file_exists(pname, desc, value)
.. code-block:: pycon

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

.. py:function:: path(filename)
.. code-block:: pycon

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

.. py:class:: Formatted_Write(indent=0, si=4, fmt=None, edge=None, lw=80, prec=None)

.. code-block:: pycon

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

.. code-block:: pycon

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

.. py:function:: make_tea(cfile=None, maxiter=100, save_headers=False, save_outputs=False, doprint=False, times=False, location_TEA=None, abun_file=None, location_out='./TEA')
.. code-block:: pycon

    Make a TEA configuration file.

    Parameters
    ----------
    cfile: String
        Input configuration file to get arguments for TEA config file.

.. py:function:: clock(t0=1560844128.89099)
.. code-block:: pycon

    Timer generator yields the time (in seconds) since the last call.

.. py:function:: get_exomol_mol(dbfile)
.. code-block:: pycon

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

.. py:function:: pf_exomol(pf_files)
.. code-block:: pycon

    Re-format Exomol partition-function files for use with Pyrat Bay.

    Parameters
    ----------
    pf_files: String or List of strings
       Input Exomol partition-function filenames.  If there are
       multiple isotopes, all of them must correspond to the same
       molecule.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> # A single file:
    >>> pt.pf_exomol('14N-1H3__MockBYTe.pf')
    Written partition-function file:
      'PF_Exomol_NH3.dat'
    for molecule NH3, with isotopes ['4111'],
    and temperature range 1--1600 K.

    >>> # Multiple files (isotopes) for a molecule:
    >>> pt.pf_exomol(['14N-1H3__MockBYTe.pf', '15N-1H3__MockBYTe-15.pf'])
    Warning: Length of PF files do not match.  Zero-padding the shorter array(s).
    Written partition-function file:
      'PF_Exomol_NH3.dat'
    for molecule NH3, with isotopes ['4111', '5111'],
    and temperature range 1--2000 K.

.. py:function:: pf_kurucz(pf_file)
.. code-block:: pycon

    Re-format Kurucz's partition-function files for use with Pyrat Bay.

    Parameters
    ----------
    pf_file: String
        Input partition-function from Kurucz webpage.  Currently only H2O
        and TiO are available (probably there's no need for any other support).
        Files can be downloaded from these links:
          http://kurucz.harvard.edu/molecules/h2o/h2opartfn.dat
          http://kurucz.harvard.edu/molecules/tio/tiopart.dat

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> # Before moving on, download Kurucz's PF files from links above.
    >>> pf_files = ['h2opartfn.dat', 'tiopart.dat']
    >>> for pf_file in pf_files:
    >>>     pt.pf_kurucz(pf_file)
    Written partition-function file:
      'PF_kurucz_H2O.dat'
    for molecule H2O, with isotopes ['1H1H16O', '1H1H17O', '1H1H18O', '1H2H16O'],
    and temperature range 10--6000 K.

    Written partition-function file:
      'PF_kurucz_TiO.dat'
    for molecule TiO, with isotopes ['66', '76', '86', '96', '06'],
    and temperature range 10--6000 K.

.. py:function:: cia_hitran(ciafile, tstep=1, wstep=1)
.. code-block:: pycon

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

.. py:function:: cia_borysow(ciafile, species1, species2)
.. code-block:: pycon

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

.. py:function:: tophat(wl0, width, margin=None, dlambda=None, resolution=None, ffile=None)
.. code-block:: pycon

      Generate a top-hat filter function, with transmission = 1.0 from
      wl0-width/2 to wl0+width/2, and an extra margin with transmission
      = 0.0 at each end.

      Parameters
      ----------
      ffile: String
          Name of the output file.
      wl0:  Float
          Filter central wavelength in microns.
      width: Float
          Filter width in microns.
      margin: Float
          Margin (in microns) with zero-valued transmission, to append
          at each end of the filter.
      dlambda: Float
          Spectral sampling rate in microns.
      resolution: Float
          Spectral sampling resolution (used if dlambda is None).
      ffile: String
          If not None, save filter to file.

      Examples
      --------
      >>> import pyratbay.tools as pt
      >>> wl0     = 1.50
      >>> width   = 0.50
      >>> margin  = 0.10
      >>> dlambda = 0.05
      >>> wl, trans = pt.tophat(wl0, width, margin, dlambda)
      >>> print(wl, trans, sep='
    ')
      [1.15 1.2  1.25 1.3  1.35 1.4  1.45 1.5  1.55 1.6  1.65 1.7  1.75 1.8
       1.85]
      [0. 0. 0. 1. 1. 1. 1. 1. 1. 1. 1. 1. 0. 0. 0.]
  

.. py:function:: resample(signal, wn, specwn, normalize=False)
.. code-block:: pycon

    Resample signal from wn to specwn wavenumber sampling using a linear
    interpolation.

    Parameters
    ----------
    signal: 1D ndarray
        A spectral signal sampled at wn.
    wn: 1D ndarray
        Signal's wavenumber sampling, in cm-1 units.
    specwn: 1D ndarray
        Wavenumber sampling to resample into, in cm-1 units.
    normalize: Bool
        If True, normalized the output resampled signal to integrate to
        1.0 (note that a normalized curve when specwn is a decreasing
        function results in negative values for resampled).

    Returns
    -------
    resampled: 1D ndarray
        The interpolated signal.
    wnidx: 1D ndarray
        The indices of specwn covered by the input signal.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> import numpy as np
    >>> wn     = np.linspace(1.3, 1.7, 11)
    >>> signal = np.array(np.abs(wn-1.5)<0.1, np.double)
    >>> specwn = np.linspace(1, 2, 51)
    >>> resampled, wnidx = pt.resample(signal, wn, specwn)
    >>> print(wnidx, specwn[wnidx], resampled, sep='\n')
    [16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34]
    [1.32 1.34 1.36 1.38 1.4  1.42 1.44 1.46 1.48 1.5  1.52 1.54 1.56 1.58
     1.6  1.62 1.64 1.66 1.68]
    [0.  0.  0.  0.  0.5 1.  1.  1.  1.  1.  1.  1.  1.  1.  0.5 0.  0.  0.
     0. ]

.. py:function:: band_integrate(spectrum, specwn, bandtrans, bandwn)
.. code-block:: pycon

    Integrate a spectrum over the band transmission.

    Parameters
    ----------
    spectrum: 1D float iterable
        Spectral signal to be integrated.
    specwn: 1D float iterable
        Wavenumber of spectrum in cm-1.
    bandtrans: 1D float iterable
        List of normalized interpolated band transmission values in each filter.
    bandwn: 1D float iterable

    Returns
    -------
    bflux: 1D float list
        Band-integrated values.

    Example
    -------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pyratbay.tools     as pt
    >>> import pyratbay.io        as io
    >>> import pyratbay.starspec  as ps
    >>> import pyratbay.constants as pc
    >>> # Load Spitzer IRAC filters:
    >>> wn1, irac1 = io.read_spectrum(pc.ROOT+'inputs/filters/spitzer_irac1_sa.dat')
    >>> wn2, irac2 = io.read_spectrum(pc.ROOT+'inputs/filters/spitzer_irac2_sa.dat')
    >>> # Spectrum to integrate:
    >>> wn = np.arange(1500, 5000.1, 1.0)
    >>> sflux = ps.bbflux(wn, 1800.0)
    >>> # Integrate over single filter:
    >>> bandflux = pt.band_integrate(sflux, wn, irac1, wn1)
    >>> # Integrate over multiple:
    >>> bandfluxes = pt.band_integrate(sflux, wn, [irac1,irac2], [wn1, wn2])
    >>> # Plot the results:
    >>> meanwn = [np.mean(wn1), np.mean(wn2)]
    >>> width = 0.5*(np.amax(wn1)-np.amin(wn1)), 0.5*(np.amax(wn2)-np.amin(wn2))
    >>> plt.figure(1)
    >>> plt.clf()
    >>> plt.semilogy(wn, sflux, 'k')
    >>> plt.plot(wn1, (irac1+1)*4e4, 'red')
    >>> plt.plot(wn2, (irac2+1)*4e4, 'blue')
    >>> plt.errorbar(meanwn[0], bandflux, xerr=width[0], fmt='o', color='red')
    >>> plt.errorbar(meanwn, bandfluxes, xerr=width, fmt='o', color='none',
    >>>              mec='k', ecolor='k')
    >>> plt.xlim(np.amin(wn), np.amax(wn))
    >>> plt.ylim(4e4, 1.2e5)
    >>> plt.xlabel('Wavenumber  (cm$^{-1}$)')
    >>> plt.ylabel(r'Flux  (erg s$^{-1}$ cm$^{-2}$ cm)')

.. py:function:: ignore_system_exit(func)
.. code-block:: pycon

    Decorator to ignore SystemExit exceptions.

.. py:function:: cf(optdepth, pressure, B)
.. code-block:: pycon

    Evaluate the contribution function equation as in Knutson et al. (2009)
    ApJ, 690, 822; Equation (2).

    Parameters
    ----------
    optdepth: 2D float ndarray
       Optical depth at each layer and wavenumber [nlayers, nwave].
    pressure: 1D float ndarray
       Atmospheric pressure array [nlayers].
    B: 2D float ndarray
       Plank emission at each layer and wavenumber [nlayers, nwave].

    Returns
    -------
    cf: 2D float ndarray
       The contribution function at each layer and wavenumber [nlayers, nwave].

.. py:function:: transmittance(optdepth, ideep)
.. code-block:: pycon

    Compute the transmittance spectrum for the impact-parameter
    raypaths of a transmission model.

    Parameters
    ----------
    optdepth: 2D float ndarray
       Optical depth at each layer and wavenumber [nlayers, nwave].
    ideep: 1D integer ndarray
       Impact-parameter indices of deepest-calculated optical depth
       at each wavenumber.

.. py:function:: bandcf(cf, bandtrans, wn, bandidx)
.. code-block:: pycon

    Compute band-averaged contribution functions or transmittances.

    Parameters
    ----------
    cf: 2D float ndarray
       The contribution function or transmittance [nlayers, nwave]
    bandtrans: List of 1D ndarrays
       List of band transmission curves.
    wn: 1D float ndarray
       The wavenumber sampling (in cm-1).
    bandidx: List of 1D ndarrays
       List of wavenumber-index arrays for each band transmission curve.

    Returns
    -------
    bandcf: 2D float ndarray
       The band-integrated contribution functions.

.. py:class:: Namespace(args=None, log=None)

.. code-block:: pycon

    A container object to hold variables.

.. code-block:: pycon

    Initialize self.  See help(type(self)) for accurate signature.

.. py:function:: parse(pyrat, cfile)
.. code-block:: pycon

    Read the command line arguments.

    Parameters
    ----------
    cfile: String
        A Pyrat Bay configuration file.

    Returns
    -------
    args: Namespace
        Object storing the attributes defined in this function, with
        the values given in cfile.
    log: Log object
        An MCcubed.utils.Log instance to log screen outputs to file.

.. py:function:: parse_str(args, param)
.. code-block:: pycon

    Parse a string parameter into args.

.. py:function:: parse_int(args, param)
.. code-block:: pycon

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
    >>> args = {'par{}'.format(i):val for i,val in enumerate(inputs)}
    >>> for i,var in enumerate(inputs):
    >>>     try:
    >>>         par = 'par{}'.format(i)
    >>>         pt.parse_int(args, par)
    >>>         print("{:s}: '{:s}' -> {}".format(par, var, args[par]))
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

.. py:function:: parse_float(args, param)
.. code-block:: pycon

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
    >>> args = {'par{}'.format(i):val for i,val in enumerate(inputs)}
    >>> for i,var in enumerate(inputs):
    >>>     try:
    >>>         par = 'par{}'.format(i)
    >>>         pt.parse_float(args, par)
    >>>         print("{:s}: '{:s}' -> {}".format(par, var, args[par]))
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

.. py:function:: parse_array(args, param)
.. code-block:: pycon

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
    >>> args = {'par{}'.format(i):val for i,val in enumerate(inputs)}
    >>> for i,var in enumerate(inputs):
    >>>     par = 'par{}'.format(i)
    >>>     pt.parse_array(args, par)
    >>>     print("{:s}: {:s} -> {}".format(par, repr(var), repr(args[par])))
    par0: '10 20' -> array([10., 20.])
    par1: '10.0 20.0' -> array([10., 20.])
    par2: 'a b' -> ['a', 'b']
    par3: 'a\n b' -> ['a', 'b']


pyratbay.blackbody
__________________


.. py:module:: pyratbay.blackbody

.. py:data:: Bwn
.. code-block:: pycon

  <built-in function Bwn>

.. py:data:: Bwn2D
.. code-block:: pycon

  <built-in function Bwn2D>


pyratbay.broadening
___________________


.. py:module:: pyratbay.broadening

.. py:class:: Lorentz(x0=0.0, hwhm=1.0, scale=1.0)

.. code-block:: pycon

    1D Lorentz profile model.

    Parameters
    ----------
    x0: Float
       Profile center location.
    hwhm: Float
       Profile's half-width at half maximum.
    scale: Float
       Scale of the profile (scale=1 returns a profile with integral=1.0).

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pyratbay.broadening as b
    >>> lor = b.Lorentz(x0=0.0, hwhm=2.5, scale=1.0)
    >>> # Half-width at half maximum is ~2.5:
    >>> x = np.linspace(-10.0, 10.0, 100001)
    >>> print(0.5 * np.ptp(x[lor(x)>0.5*np.amax(lor(x))]))
    2.4998
    >>> # Integral is ~ 1.0:
    >>> x = np.linspace(-5000.0, 5000.0, 100001)
    >>> print(np.trapz(lor(x), x))
    0.999681690140321
    >>> # Take a look at a Lorenzt profile:
    >>> x = linspace(-10, 10, 101)
    >>> plt.plot(x, lor(x))

.. code-block:: pycon

    Initialize self.  See help(type(self)) for accurate signature.

.. py:class:: Gauss(x0=0.0, hwhm=1.0, scale=1.0)

.. code-block:: pycon

    1D Gaussian profile model.

    Parameters
    ----------
    x0: Float
       Profile center location.
    hwhm: Float
       Profile's half-width at half maximum.
    scale: Float
       Scale of the profile (scale=1 returns a profile with integral=1.0).

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pyratbay.broadening as b
    >>> gauss = b.Gauss(x0=0.0, hwhm=2.5, scale=1.0)
    >>> # Half-width at half maximum is ~2.5:
    >>> x = np.linspace(-10.0, 10.0, 100001)
    >>> print(0.5 * np.ptp(x[gauss(x)>0.5*np.amax(gauss(x))]))
    2.4998
    >>> # Integral is ~ 1.0:
    >>> x = np.linspace(-5000.0, 5000.0, 100001)
    >>> print(np.trapz(gauss(x), x))
    1.0
    >>> # Take a look at a Lorenzt profile:
    >>> x = linspace(-10, 10, 101)
    >>> plt.plot(x, gauss(x))

.. code-block:: pycon

    Initialize self.  See help(type(self)) for accurate signature.

.. py:class:: Voigt(x0=0.0, hwhmL=1.0, hwhmG=1.0, scale=1.0)

.. code-block:: pycon

    1D Voigt profile model.

    Parameters
    ----------
    x0: Float
       Line center location.
    hwhmL: Float
       Half-width at half maximum of the Lorentz distribution.
    hwhmG: Float
       Half-width at half maximum of the Gaussian distribution.
    scale: Float
       Scale of the profile (scale=1 returns a profile with integral=1.0).

    Example
    -------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pyratbay.broadening as b
    >>> Nl = 5
    >>> Nw = 10.0
    >>> hG = 1.0
    >>> HL = np.logspace(-2, 2, Nl)
    >>> l = b.Lorentz(x0=0.0)
    >>> d = b.Gauss  (x0=0.0, hwhm=hG)
    >>> v = b.Voigt  (x0=0.0, hwhmG=hG)

    >>> plt.figure(11, (6,6))
    >>> plt.clf()
    >>> plt.subplots_adjust(0.15, 0.1, 0.95, 0.95, wspace=0, hspace=0)
    >>> for i in np.arange(Nl):
    >>>   hL = HL[i]
    >>>   ax = plt.subplot(Nl, 1, 1+i)
    >>>   v.hwhmL = hL
    >>>   l.hwhm  = hL
    >>>   width = 0.5346*hL + np.sqrt(0.2166*hL**2+hG**2)
    >>>   x = np.arange(-Nw*width, Nw*width, width/1000.0)
    >>>   plt.plot(x/width, l(x), lw=2.0, color="b",         label="Lorentz")
    >>>   plt.plot(x/width, d(x), lw=2.0, color="limegreen", label="Doppler")
    >>>   plt.plot(x/width, v(x), lw=2.0, color="orange",    label="Voigt",
    >>>            dashes=(8,2))
    >>>   plt.ylim(np.amin([l(x), v(x)]), 3*np.amax([l(x), v(x), d(x)]))
    >>>   ax.set_yscale("log")
    >>>   plt.text(0.025, 0.75, r"$\rm HW_L/HW_G={:4g}$".format(hL/hG),
    >>>            transform=ax.transAxes)
    >>>   plt.xlim(-Nw, Nw)
    >>>   plt.xlabel(r"$\rm x/HW_V$", fontsize=12)
    >>>   plt.ylabel(r"$\rm Profile$")
    >>>   if i != Nl-1:
    >>>       ax.set_xticklabels([""])
    >>>   if i == 0:
    >>>       plt.legend(loc="upper right", fontsize=11)

.. code-block:: pycon

    Initialize self.  See help(type(self)) for accurate signature.

.. py:function:: min_widths(min_temp, min_wn, max_mass, dlratio=0.1)
.. code-block:: pycon

    Estimate the minimum Doppler and Lorentz half-widths at half maximum
    (cm-1) for a given atmosphere.

    Parameters
    ----------
    min_temp: Float
        Minimum atmospheric tmperature (Kelvin degrees).
    min_wn: Float
        Minimum spectral wavenumber (cm-1).
    max_mass: Float
        Maximum mass of molecule/isotope (amu).
    dlratio: Float
        Doppler--Lorentz width ratio.

    Returns
    -------
    dmin: Float
        Minimum Doppler HWHM (cm-1).
    lmin: Float
        Minimum Lorentz HWHM (cm-1).

    Examples
    --------
    >>> wn = np.amin(pyrat.spec.wn)
    >>> temperature = np.amin(pyrat.atm.temp)
    >>> mols = np.unique(pyrat.iso.imol)
    >>> mols = mols[np.where(mols>=0)]
    >>> molmass = np.amax(pyrat.mol.mass[mols])
    >>> # TBD

.. py:function:: max_widths(min_temp, max_temp, max_wn, min_mass, max_rad, max_press)
.. code-block:: pycon

    Estimate the maximum Doppler and Lorentz half-widths at half maximum
    (cm-1) for a given atmosphere.

    Parameters
    ----------
    min_temp: Float
        Minimum atmospheric tmperature (Kelvin degrees).
    max_temp: Float
        Maximum atmospheric tmperature (Kelvin degrees).
    max_wn: Float
        Maximum spectral wavenumber (cm-1).
    min_mass: Float
        Minimum mass of molecule/isotope (amu).
    max_press: Float
        Maximum collision radius of molecule/isotope (Angstrom).
    max_press: Float
        Maximum atmospheric pressure (barye).

    Returns
    -------
    dmax: Float
        Maximum Doppler HWHM (cm-1).
    lmax: Float
        Maximum Lorentz HWHM (cm-1).

    Examples
    --------
    >>> pmax = np.amax(pyrat.atm.press)
    >>> wn = np.amax(pyrat.spec.wn)
    >>> temperature = np.amax(pyrat.atm.temp)
    >>> mols = np.unique(pyrat.iso.imol) # Molecules with transitions
    >>> mols = mols[np.where(mols>=0)]   # Remove -1's
    >>> mmin = np.amin(pyrat.mol.mass[mols])
    >>> rmax = np.amax(pyrat.mol.radius[mols])
    >>> # TBD


pyratbay.lineread
_________________


.. py:module:: pyratbay.lineread

.. py:function:: makeTLI(dblist, pflist, dbtype, tlifile, wllow, wlhigh, wlunits, log)
.. code-block:: pycon

    Driver function to create a TLI file.

    Parameters
    ----------
    dblist: List of strings
        Opacity databases to read.
    pflist: List of strings
        Partition function for each of the databases.
    dbtype: List of strings
        Type of each database.
    tlifile: String
        Output TLI file name.
    wllow: String or float
        Lower wavelength boundary to consider. If float, assume units
        from wlunits input.  Otherwise, wllow sets the value and units
        (for example: '1.0 um').
    wlhigh: String or float
        High wavelength boundary to consider. If float, assume units
        from wlunits input.  Otherwise, wlhigh sets the value and units.
    wlunits: String
        Wavelength units (when not specified in wllow nor wlhigh).
    log: Log object
        An MCcubed.utils.Log instance to log screen outputs to file.


pyratbay.lineread.database
__________________________


.. py:module:: pyratbay.lineread.database

.. py:class:: hitran(dbfile, pffile, log)

.. code-block:: pycon

    HITRAN/HITEMP database reader.

.. code-block:: pycon

    Initialize HITRAN database object.

    Parameters
    ----------
    dbfile: String
        File with the Database line-transition info.
    pffile: String
        File with the partition function.
    log: Log object
        An MCcubed.utils.Log instance to log screen outputs to file.

.. py:class:: pands(dbfile, pffile, log)

.. code-block:: pycon

    Partridge & Schwenke (1997) H2O database reader.

.. code-block:: pycon

    Initialize P&S database object.

    Parameters
    ----------
    dbfile: String
        File with the Database line-transition info.
    pffile: String
        File with the partition function.
    log: Log object
        An MCcubed.utils.Log instance to log screen outputs to file.

.. py:class:: tioschwenke(dbfile, pffile, log)

.. code-block:: pycon

    Notes:
    ------
    Linelist and partition function downloaded from:
      http://kurucz.harvard.edu/molecules/tio/tioschwenke.bin
      http://kurucz.harvard.edu/molecules/tio/tiopart.dat

    There might be a problem with the linebreak character of the partition
    function.  One way to fix is, on vim do: :%s/    /    /g

.. code-block:: pycon

    Initialize self.  See help(type(self)) for accurate signature.

.. py:class:: voplez(dbfile, pffile, log)

.. code-block:: pycon

    Download the linelist from:

.. code-block:: pycon

    Initializer.

.. py:class:: vald(dbfile, pffile, log)

.. code-block:: pycon

    Notes
    -----
      Download linelist from: http://vald.astro.uu.se/~vald/php/vald.php
         Selecting 'Extract Element' and 'Short format'.
      Download partition functions from:

.. code-block:: pycon

    Initialize Basic data for the Database.

    Parameters
    ----------
    dbfile: String
       File with the Database line-transition info.
    pffile: String
       File with the partition function.
    log: File
       File object to store the log.

.. py:class:: exomol(dbfile, pffile, log)

.. code-block:: pycon

    Exomol database reader.

.. code-block:: pycon

    Initialize Exomol database object.

    Parameters
    ----------
    dbfile: String
        File with the Database line-transition info.
    pffile: String
        File with the partition function.
    log: Log object
        An MCcubed.utils.Log instance to log screen outputs to file.

.. py:class:: repack(dbfile, pffile, log)

.. code-block:: pycon

    Repack database reader.

.. code-block:: pycon

    Initialize Exomol database object.

    Parameters
    ----------
    dbfile: String
        File with the Database line-transition info.
    pffile: String
        File with the partition function.
    log: Log object
        An MCcubed.utils.Log instance to log screen outputs to file.


pyratbay.plots
______________


.. py:module:: pyratbay.plots

.. py:function:: spectrum(spectrum, wavelength, path, data=None, uncert=None, bandwl=None, bandflux=None, bandtrans=None, bandidx=None, starflux=None, rprs=None, label='model', bounds=None, logxticks=None, gaussbin=2.0, yran=None, filename=None, fignum=501)
.. code-block:: pycon

    Plot a transmission or emission model spectrum with (optional) data
    points with error bars and band-integrated model.

    Parameters
    ----------
    spectrum: 1D float ndarray
        Planetary spectrum evaluated at wavelength.
    wavelength: 1D float ndarray
        The wavelength of the model in microns.
    path: String
        Observing-geometry path: transit or eclipse.
    data: 1D float ndarray
        Observing data points at each bandwl.
    uncert: 1D float ndarray
        Uncertainties of the data points.
    bandwl: 1D float ndarray
        The mean wavelength for each band/data point.
    bandflux: 1D float ndarray
        Band-integrated model spectrum at each bandwl.
    bandtrans: List of 1D float ndarrays
        Transmission curve for each band.
    bandidx: List of 1D float ndarrays.
        The indices in wavelength for each bandtrans.
    starflux: 1D float ndarray
        Stellar spectrum evaluated at wavelength.
    rprs: Float
        Planet-to-star radius ratio.
    logxticks: 1D float ndarray
        If not None, switch the X-axis scale from linear to log, and set
        the X-axis ticks at the locations given by logxticks.
    gaussbin: Integer
        Standard deviation for Gaussian-kernel smoothing (in number of samples).
    yran: 1D float ndarray
        Figure's Y-axis boundaries.
    filename: String
        If not None, save figure to filename.
    fignum: Integer
        Figure number.

    Returns
    -------
    ax: AxesSubplot instance
        The matplotlib Axes of the figure.

.. py:function:: cf(bandcf, bandwl, path, pressure, radius, rtop=0, filename=None, filters=None, fignum=-21)
.. code-block:: pycon

    Plot the band-integrated normalized contribution functions
    (emission) or transmittance (transmission).

    Parameters
    ----------
    bandcf: 2D float ndarray
        Band-integrated contribution functions [nfilters, nlayers].
    bandwl: 1D float ndarray
        Mean wavelength of the bands in microns.
    path: String
        Observing geometry (transit or eclipse).
    pressure: 1D float ndarray
        Layer's pressure array (barye units).
    radius: 1D float ndarray
        Layer's impact parameter array (cm units).
    rtop: Integer
        Index of topmost valid layer.
    filename: String
        Filename of the output figure.
    filters: 1D string ndarray
        Name of the filter bands (optional).
    fignum: Integer
        Figure number.

    Notes
    -----
    - The dashed lines denote the 0.16 and 0.84 percentiles of the
      cumulative contribution function or the transmittance (i.e.,
      the boundaries of the central 68% of the respective curves).
    - If there are more than 80 filters, this code will thin the
      displayed filter names.

.. py:function:: posterior_pt(posterior, tmodel, targs, tpars, ifree, pressure, bestpars=None, filename=None)
.. code-block:: pycon

    Plot the posterior PT profile.

    Parameters
    ----------
    posterior: 2D float ndarray
        MCMC posterior distribution for tmodel (of shape [nparams, nfree]).
    tmodel: Callable
        Temperature-profile model.
    tpars: 1D float ndarray
        Temperature-profile parameters (including fixed parameters).
    ifree: 1D bool ndarray
        Mask of free (True) and fixed (False) parameters in tpars.
        The number of free parameters must match nfree in posterior.
    targs: List
        List of additional arguments for tmodel.
    pressure: 1D float ndarray
        The atmospheric pressure profile in barye.
    bestpars: 1D float ndarray
        Best-fitting temperature-profile parameters.
    filename: String
        If not None, save figure to filename.


pyratbay.starspec
_________________


.. py:module:: pyratbay.starspec

.. py:function:: read_kurucz(filename, temp=None, logg=None)
.. code-block:: pycon

    Extract stellar flux models from a Kurucz file.
    Kurucz model files can be found at http://kurucz.harvard.edu/grids.html

    Parameters
    ----------
    filename: String
       Name of a Kurucz model file.
    temp: Float
       Requested surface temperature for the Kurucz model.
       If temp and logg are not None, return the model with the closest
       surface temperature and gravity.
    logg: Float
       Requested log10 of the surface gravity for the Kurucz model
       (where g is in cgs units).

    Returns
    -------
    flux: 1D or 2D float ndarray
       If temp and logg are not None, a 1D array with the kurucz surface
       flux per unit wavenumber (erg s-1 cm-2 cm) of the closest model to
       the input temperature and gravity.
       Else, a 2D array with all kurucz models in file, of shape
       [nmodels, nwave].
    wavenumber: 1D ndarray
       Wavenumber sampling of the flux models (in cm-1 units).
    ktemp: Scalar or 1D float ndarray
       Surface temperature of the output models (in Kelvin degrees).
    klogg: Scalar or 1D float ndarray
       log10 of the stellar surface gravity of the output models (in cm s-2).
    continuum: 2D ndarray
       The models' fluxes with no line absorption.  Same units and
       shape of flux. Returned only if temp and logg are None.

    Examples
    --------
    >>> import pyratbay.starspec  as ps
    >>> import pyratbay.constants as pc
    >>> import numpy as np
    >>> # Download a Kurucz stellar model file from:
    >>> # http://kurucz.harvard.edu/grids/gridp00odfnew/fp00k0odfnew.pck
    >>> # Read a single model from the kurucz file:
    >>> kfile = 'fp00k0odfnew.pck'
    >>> tsun = 5770.0  # Sun's surface temperature
    >>> gsun = 4.44    # Sun's surface gravity (log)
    >>> flux, wn, ktemp, klogg = ps.read_kurucz(kfile, tsun, gsun)
    >>> # Compute brightness at 1 AU from a 1 Rsun radius star:
    >>> s = np.trapz(flux, wn) * (pc.rsun/pc.au)**2
    >>> print("Solar constant [T={:.0f} K, logg={:.1f}]:  S = {:.1f} W m-2".
    >>>       format(ktemp, klogg, s * 1e-3))
    Solar constant [T=5750 K, logg=4.5]:  S = 1340.0 W m-2
    >>> # Pretty close to the solar constant: ~1361 W m-2

    >>> # Read the whole set of models in file:
    >>> # (in this case, ktemp and klogg are 1D arrays)
    >>> fluxes, wn, ktemp, klogg, continua = ps.read_kurucz(kfile)

.. py:function:: bbflux(wn, teff)
.. code-block:: pycon

    Compute the emission flux of a blackbody at temperature Teff
    in wavenumber space.

    Parameters
    ----------
    wn: 1D float iterable
       Wavenumber array where to evaluate the flux (cm-1).
    teff: Float
       The effective temperature (Kelvin).

    Return
    ------
    flux: 1D float ndarray
       blackbody flux (erg s-1 cm-2 cm) at wn.

    Examples
    --------
    >>> import pyratbay.starspec  as ps
    >>> import pyratbay.constants as pc
    >>> import numpy as np
    >>> tsun = 5772.0
    >>> wn = np.logspace(-1, 5, 30000)
    >>> flux = ps.bbflux(wn, tsun)
    >>> # Solar constant:
    >>> s = np.trapz(flux, wn) * (pc.rsun/pc.au)**2
    >>> print("Solar constant (Teff={:.0f}K): S = {:.1f} W m-2\n"
    >>>       "Wien's displacement law: wn(flux_max) = {:.1f} cm-1\n"
    >>>       "             5.879E10 Hz/K * Teff / c = {:.1f} cm-1".
    >>>       format(tsun, s*1e-3, wn[np.argmax(flux)], 5.879e10*tsun/pc.c))
    Solar constant (Teff=5772K): S = 1361.2 W m-2
    Wien's displacement law: wn(flux_max) = 11318.0 cm-1
                 5.879E10 Hz/K * Teff / c = 11319.0 cm-1


pyratbay.atmosphere
___________________


.. py:module:: pyratbay.atmosphere

.. py:function:: writeatm(atmfile, pressure, temperature, species, abundances, punits, header, radius=None, runits=None)
.. code-block:: pycon

    Write an atmospheric file following the Pyrat format.

    Parameters
    ----------
    atmfile: String
       Name of output atmospheric file.
    pressure: 1D float ndarray
       Monotonously decreasing pressure profile (in barye).
    temperature: 1D float ndarray
       Temperature profile for pressure layers (in Kelvin).
    species: 1D string ndarray
       List of atmospheric species.
    abundances: 2D float ndarray
       The species mole mixing ratio (of shape [nlayers,nspecies]).
    punits:  String
       Pressure units of output.
    header:  String
       Header message (comment) to include at the top of the file.
    radius: 1D float ndarray
       Monotonously increasing radius profile (in cm).
    runits:  String
       Radius units of output.

.. py:function:: readatm(atmfile)
.. code-block:: pycon

    Read a Pyrat atmospheric file.

    Parameters
    ----------
    atmfile: String
       File path to a Pyrat Bay's atmospheric file.

    Returns
    -------
    units: 4-element string tuple
        Units for pressure, temperature, abundance, and radius as given
        in the atmospheric file.
    species: 1D string ndarray
        The list of species names read from the atmospheric file (of
        size nspec).
    press: 1D float ndarray
        The atmospheric pressure profile (of size nlayers). The
        file's @PRESSURE keyword indicates the ouptput units.
    temp: 1D float ndarray
        The atmospheric temperature profile (of size nlayers). The
        file's @TEMPERATURE keyword indicates the ouptput units.
    q: 2D float ndarray
        The mixing ratio profiles of the atmospheric species (of size
        [nlayers,nspec]).  The file's @ABUNDANCE indicates the output
        units.
    radius: 1D float ndarray
        The atmospheric altiture profile (of size nlayers).  None if the
        atmospheric file does not contain a radius profile.
        The file's @RADIUS keyword indicates the output units.

.. py:function:: makeatomic(solar, afile, xsolar=1.0, swap=None)
.. code-block:: pycon

    Generate an (atomic) elemental-abundances file by scaling the
    solar abundances from Asplund et al. (2009).
    http://adsabs.harvard.edu/abs/2009ARA%26A..47..481A

    Parameters
    ----------
    solar: String
       Input solar abundances filename.
    afile: String
       Output filename of modified atomic abundances.
    xsolar: Integer
       Multiplication factor for metal abundances (everything
       except H and He).
    swap: 2-element string tuple
        Swap the abundances of the given elements.

.. py:function:: readatomic(afile)
.. code-block:: pycon

    Read an elemental (atomic) composition file.

    Parameters
    ----------
    afile: String
      File with atomic composition.

    Returns
    -------
    anum: 1D integer ndarray
       Atomic number (except for Deuterium, which has anum=0).
    symbol: 1D string ndarray
       Elemental chemical symbol.
    dex: 1D float ndarray
       Logarithmic number-abundance, scaled to log(H) = 12.
    name: 1D string ndarray
       Element names.
    mass: 1D float ndarray
       Elemental mass in amu.

    Uncredited developers
    ---------------------
    Jasmina Blecic

.. py:function:: makepreatm(pressure, temp, afile, elements, species, patm)
.. code-block:: pycon

    Generate a pre-atm file for TEA containing the elemental abundances
    at each atmospheric layer.

    Parameters
    ----------
    pressure: String
       Pressure atmospheric profile (bar).
    temp: 1D float array
       Temperature atmospheric profile (in K).
    afile: String
       Name of the elemental abundances file.
    elements: List of strings
       List of input elemental species.
    species: List of strings
       List of output molecular species.
    patm: String
       Output pre-atmospheric filename.

    Uncredited developers
    ---------------------
    Jasmina Blecic

.. py:function:: TEA2pyrat(teafile, atmfile, req_species)
.. code-block:: pycon

    Format a TEA atmospheric file into a Pyrat atmospheric file.

    Paramters
    ---------
    teafile:  String
       Input TEA atmospheric file.
    atmfile:  String
       Output Pyrat atmospheric file.
    req_species: List of strings
       The requested species for output.

.. py:function:: readmol(file)
.. code-block:: pycon

    Read a molecules file to extract their ID, symbol, mass, and diameter.

    Parameters
    ----------
    file: String
       The molecule file path.

    Returns
    -------
    molID: 1D integer ndarray
       The molecules' ID.
    symbol: 1D string ndarray
       The molecule's name.
    mass: 1D float ndarray
       The mass of the molecules (in g mol-1).
    diam: 1D float ndarray
       The collisional diameter of the molecules (in Angstrom).

    Notes
    -----
    In all truthfulness, these are species, not only molecules, as the
    file also contain elemental particles.

.. py:function:: pressure(ptop, pbottom, nlayers, units='bar', log=None, verb=0)
.. code-block:: pycon

    Compute a log-scale pressure profile.

    Parameters
    ----------
    ptop: String or Float
       Pressure at the top of the atmosphere. If string, may contain the units.
    pbottom: String or Float
       Pressure at the bottom of the atmosphere. If string, may contain the units.
    nlayers: Integer
       Number of pressure layers.
    units: String
       Pressure input units (if not defined in ptop, pbottom).
       Available units are: barye, mbar, pascal, bar (default), and atm.
    log: Log object
       Screen-output log handler.
    verb: Integer
       Verbosity level (when log is None). Print out when verb > 0.

    Returns
    -------
    press: 1D float ndarray
       The pressure profile (in barye units).

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants  as pc
    >>> nlayers = 9
    >>> # These are all equivalent:
    >>> p1 = pa.pressure(ptop=1e-6,   pbottom=1e2, nlayers=nlayers)
    >>> p2 = pa.pressure(1e-6,        1e2,         nlayers, 'bar')
    >>> p3 = pa.pressure('1e-6 bar', '1e2 bar',    nlayers)
    >>> p4 = pa.pressure(1e-6*pc.bar, 1e2*pc.bar,  nlayers, 'barye')
    >>> print(p1/pc.bar)
    [1.e-06 1.e-05 1.e-04 1.e-03 1.e-02 1.e-01 1.e+00 1.e+01 1.e+02]

.. py:function:: temperature(tmodel, pressure=None, rstar=None, tstar=None, tint=100.0, gplanet=None, smaxis=None, runits='cm', nlayers=None, log=None, tparams=None)
.. code-block:: pycon

    Temperature profile wrapper.

    Parameters
    ----------
    tmodel: String
        Name of the temperature model.
    pressure: 1D float ndarray
        Atmospheric pressure profile in barye units.
    rstar: String or float
        Stellar radius. If string, may contain the units.
    tstar: String or float
        Stellar temperature in kelvin degrees.
    tint: String or float
        Planetary internal temperature in kelvin degrees.
    gplanet: String or float
        Planetary atmospheric temperature in cm s-2.
    smaxis: String or float
        Orbital semi-major axis. If string, may contain the units.
    runits: String
        Default units for rstar and smaxis.  Available units are: A, nm, um,
        mm, cm (default), m, km, au, pc, rearth, rjup, rsun.
    nlayers: Integer
        Number of pressure layers.
    log: Log object
        Screen-output log handler.
    tparams: 1D float ndarray
        Temperature model parameters. If None, return a tuple with the
        temperature model, its arguments, and the number or required parameters.

    Returns
    -------
    If tparams is not None:
        temperature: 1D float ndarray
            The evaluated atmospheric temperature profile.
    If tparams is None:
        Tmodel: Callable
            The atmospheric temperature model.
        targs: List
            The list of additional arguments (besides the model parameters).
        ntpars: Integer
            The expected number of model parameters.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> nlayers = 11
    >>> # Isothermal profile:
    >>> temp_iso = pa.temperature("isothermal", tparams=1500.0, nlayers=nlayers)
    >>> print(temp_iso)
    [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.]
    >>> # Three-channel Eddington-approximation profile:
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, "bar")
    >>> tparams = np.array([-1.5, -0.8, -0.8, 0.5, 1.0])
    >>> temp = pa.temperature('tcea', pressure, rstar="0.756 rsun", tstar=5040,
            tint=100.0, gplanet=2200.0, smaxis="0.031 au", tparams=tparams)
    >>> print(temp)
    [1047.04157312 1047.04189805 1047.04531644 1047.08118784 1047.45648563
     1051.34469989 1088.69956369 1311.86379107 1640.12857767 1660.02396061
     1665.30121021]

.. py:function:: uniform(pressure, temperature, species, abundances, punits='bar', log=None, atmfile=None)
.. code-block:: pycon

    Generate an atmospheric file with uniform abundances.
    Save it into atmfile.

    Parameters
    ----------
    pressure: 1D float ndarray
        Monotonously decreasing pressure profile (in punits).
    temperature: 1D float ndarray
        Temperature profile for pressure layers (in Kelvin).
    species: 1D string ndarray
        List of atmospheric species.
    abundances: 1D float ndarray
        The species mole mixing ratio.
    punits:  String
       Pressure units.
    log: Log object
        Screen-output log handler.
    atmfile: String
        If not None, filename to save atmospheric model.

    Returns
    -------
    qprofiles: 2D Float ndarray
        Abundance profiles of shape [nlayers,nspecies]

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> atmfile = "atm_test.dat"
    >>> nlayers = 11
    >>> punits  = 'bar'
    >>> pressure    = pa.pressure(1e-8, 1e2, nlayers, punits)
    >>> temperature = pa.tmodels.isothermal(1500.0, nlayers)
    >>> species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    >>> abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    >>> qprofiles = pa.abundances(atmfile, pressure, temperature, species,
    >>>                           quniform=abundances, punits=punits)
    >>> print(qprofiles)
    [[8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]
     [8.496e-01 1.500e-01 1.000e-04 1.000e-04 1.000e-08 1.000e-04]]

.. py:function:: abundances(atmfile, pressure, temperature, species, elements=None, quniform=None, punits='bar', xsolar=1.0, solar='/Users/pato/Dropbox/IWF/projects/2014_pyratbay/2016-01-08_develop/pyratbay/inputs/AsplundEtal2009.txt', log=None, verb=0)
.. code-block:: pycon

    Wrapper to compute atmospheric abundaces for given pressure and
    temperature profiles with either uniform abundances or TEA.

    Parameters
    ----------
    atmfile: String
       Output file where to save the atmospheric model.
    pressure: 1D float ndarray
       Atmospheric pressure profile (barye).
    temperature: 1D float ndarray
       Atmospheric temperature profile (Kelvin).
    species: 1D string list
       Atmospheric composition.
    elements: 1D strings list
       Elemental composition.
    quniform: 1D float ndarray
       If not None, the species abundances at all layers.
    punits: String
       Output pressure units.
    xsolar: Float
       Metallicity enhancement factor.
    solar: String
       Solar elemental abundances file.
    log: Log object
       Screen-output log handler.
    verb: Integer
       Verbosity level.

    Returns
    -------
    q: 2D float ndarray
       Atmospheric abundances (volume mixing fraction).

    Example
    -------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants  as pc
    >>> atmfile = "pbtea.atm"
    >>> nlayers = 100
    >>> press = np.logspace(-8, 3, nlayers) * pc.bar
    >>> temp  = 900+500/(1+np.exp(-(np.log10(press)+1.5)*1.5))
    >>> species = ['H2O', 'CH4', 'CO', 'CO2', 'NH3', 'C2H2', 'C2H4', 'HCN',
    >>>            'N2', 'H2', 'H', 'He']
    >>> # Automatically get 'elements' necessary from the list of species:
    >>> Q = pa.abundances("pbtea.atm", press, temp, species)

.. py:function:: hydro_g(pressure, temperature, mu, g, p0=None, r0=None)
.. code-block:: pycon

    Calculate radii using the hydrostatic-equilibrium equation considering
    a constant gravity.

    Parameters
    ----------
    pressure: 1D float ndarray
       Atmospheric pressure for each layer (in barye).
    temperature: 1D float ndarray
       Atmospheric temperature for each layer (in K).
    mu: 1D float ndarray
       Mean molecular mass for each layer (in g mol-1).
    g: Float
       Atmospheric gravity (in cm s-2).
    p0: Float
       Reference pressure level (in barye) where radius(p0) = r0.
    r0: Float
       Reference radius level (in cm) corresponding to p0.

    Returns
    -------
    radius: 1D float ndarray
       Radius for each layer (in cm).

    Notes
    -----
    If the reference values (p0 and r0) are not given, set radius = 0.0
    at the bottom of the atmosphere.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants as pc
    >>> nlayers = 11
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    >>> temperature = pa.tmodels.isothermal(1500.0, nlayers)
    >>> mu = np.tile(2.3, nlayers)
    >>> g = pc.G * pc.mjup / pc.rjup**2
    >>> r0 = 1.0 * pc.rjup
    >>> p0 = 1.0 * pc.bar
    >>> # Radius profile in Jupiter radii:
    >>> radius = pa.hydro_g(pressure, temperature, mu, g, p0, r0) / pc.rjup
    >>> print(radius)
    [1.0563673  1.04932138 1.04227547 1.03522956 1.02818365 1.02113774
     1.01409182 1.00704591 1.         0.99295409 0.98590818]

.. py:function:: hydro_m(pressure, temperature, mu, M, p0, r0)
.. code-block:: pycon

    Calculate radii using the hydrostatic-equilibrium equation considering
    a variable gravity: g=G*M/r**2

    Parameters
    ----------
    pressure: 1D float ndarray
       Atmospheric pressure for each layer (in barye).
    temperature: 1D float ndarray
       Atmospheric temperature for each layer (in K).
    mu: 1D float ndarray
       Mean molecular mass for each layer (in g mol-1).
    M: Float
       Planetary mass (in g).
    p0: Float
       Reference pressure level (in barye) where radius(p0) = r0.
    r0: Float
       Reference radius level (in cm) corresponding to p0.

    Returns
    -------
    radius: 1D float ndarray
       Radius for each layer (in cm).

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants  as pc
    >>> nlayers = 11
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    >>> temperature = pa.tmodels.isothermal(1500.0, nlayers)
    >>> mu = np.tile(2.3, nlayers)
    >>> Mp = 1.0 * pc.mjup
    >>> r0 = 1.0 * pc.rjup
    >>> p0 = 1.0 * pc.bar
    >>> # Radius profile in Jupiter radii:
    >>> radius = pa.hydro_m(pressure, temperature, mu, Mp, p0, r0) / pc.rjup
    >>> print(radius)
    [1.05973436 1.05188019 1.04414158 1.036516   1.029001   1.02159419
     1.01429324 1.00709591 1.         0.99300339 0.986104  ]

.. py:function:: stoich(species)
.. code-block:: pycon

      Extract the elemental composition from a list of species.

      Parameters
      ----------
      species: 1D string list
         List of species.

      Returns
      -------
      elements: 1D string list
         List of elements contained in species list.
      stoich: 2D integer ndarray
         Stoichiometric elemental values for each species (number of elements).

      Examples
      --------
      >>> import pyratbay.atmosphere as pa
      >>> species = ['H2', 'He', 'H2O', 'CO', 'CO2', 'CH4']
      >>> elements, stoichs = pa.stoich(species)
      >>> print('{}
    {}'.format(elements, stoichs))
      ['C', 'H', 'He', 'O']
      [[0 2 0 0]
       [0 0 1 0]
       [0 2 0 1]
       [1 0 0 1]
       [1 0 0 2]
       [1 4 0 0]]
  

.. py:function:: meanweight(abundances, species, molfile='/Users/pato/Dropbox/IWF/projects/2014_pyratbay/2016-01-08_develop/pyratbay/inputs/molecules.dat')
.. code-block:: pycon

    Calculate the mean molecular weight (a.k.a. mean molecular mass)
    for the given abundances composition.

    Parameters
    ----------
    abundances: 2D float iterable
        Species mol-mixing-fraction array of shape [nlayers,nmol].
    species: 1D string iterable
        Species names.
    molfile: String
        A molecules file with the species info.

    Returns
    -------
    mu: 1D float ndarray
        Mean molecular weight at each layer for the input abundances.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> species     = ['H2', 'He', 'H2O', 'CO', 'CO2', 'CH4']
    >>> abundances  = [[0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]]
    >>> mu = pa.meanweight(abundances, species)
    >>> print(mu)
    [2.31928918]

.. py:function:: IGLdensity(abundances, pressure, temperature)
.. code-block:: pycon

    Use the Ideal gas law to calculate the density in molecules cm-3.

    Parameters
    ----------
    abundances: 2D float ndarray
       Species volume mixing ratio profiles, of shape [nlayers,nmol].
    pressure: 1D ndarray
       Atmospheric pressure profile (in barye units).
    temperature: 1D ndarray
       Atmospheric temperature (in kelvin).

    Returns
    -------
    density: 2D float ndarray
       Atmospheric density in molecules per centimeter^3.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> atmfile = "uniform_test.atm"
    >>> nlayers = 11
    >>> pressure    = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    >>> temperature = np.tile(1500.0, nlayers)
    >>> species     = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    >>> abundances  = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    >>> qprofiles = pa.uniform(pressure, temperature, species, abundances)
    >>> dens = pa.IGLdensity(qprofiles, pressure, temperature)
    >>> print(dens[0])
    [4.10241993e+10 7.24297303e+09 4.82864869e+06 4.82864869e+06
     4.82864869e+02 4.82864869e+06]

.. py:function:: qcapcheck(Q, qcap, ibulk)
.. code-block:: pycon

    Check if the cummulative abundance of traces exceed qcap.

    Parameters
    ----------
    Q: 2D float ndarray
       Mole mixing ratio of the species in the atmosphere [Nlayers, Nspecies].
    qcap: Float
       Cap threshold for cummulative trace abundances.
    ibulk: 1D integer ndarray
       Indices of the bulk species to calculate the mixing ratio.

    Returns
    -------
    qcapcheck: Bool
       Flag indicating whether trace abundances sum more than qcap.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> # Make an atmosphere:
    >>> pressure    = pa.pressure(ptop=1e-8, pbottom=1e2, nlayers=11, units='bar')
    >>> temperature = np.tile(1500.0, 11)
    >>> species     = ["H2", "He", "H2O"]
    >>> abundances  = [0.8495, 0.15, 5e-4]
    >>> qprofiles = pa.uniform(pressure, temperature, species, abundances)
    >>> ibulk = [0,1]
    >>> # Sum of all metals (H2O) is not above qcap:
    >>> qcap = 1e-3
    >>> print(pa.qcapcheck(qprofiles, qcap, ibulk))
    False
    >>> # Sum of all metals (H2O) is exceedes qcap:
    >>> qcap = 1e-4
    >>> print(pa.qcapcheck(qprofiles, qcap, ibulk))
    True

.. py:function:: balance(Q, ibulk, ratio, invsrat)
.. code-block:: pycon

    Balance the mole mixing ratios of the bulk species, Q[ibulk],
    such that Sum(Q) = 1.0 at each level.

    Parameters
    ----------
    Q: 2D float ndarray
       Mole mixing ratio of the species in the atmosphere [Nlayers, Nspecies].
    ibulk: 1D integer ndarray
       Indices of the bulk species to calculate the mixing ratio.
    ratio: 2D float ndarray
       Abundance ratio between species indexed by ibulk.
    invsrat: 1D float ndarray
       Inverse of the sum of the ratios (at each layer).

    Notes
    -----
    Let the bulk abundance species be the remainder of the sum of the trace
    species:
       Q_{\rm bulk} = \sum Q_j = 1.0 - \sum Q_{\rm trace}.
    This code assumes that the abundance ratio among bulk species
    remains constant in each layer:
       {\rm ratio}_j = Q_j/Q_0.
    The balanced abundance of the bulk species is then:
       Q_j = \frac{{\rm ratio}_j * Q_{\rm bulk}} {\sum {\rm ratio}}.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> q = np.tile([0.8, 0.2, 0.5], (5,1))
    >>> q[4] = 0.5, 0.5, 0.5
    >>> ibulk = [0, 1]
    >>> bratio, invsrat = pa.ratio(q, ibulk)
    >>> pa.balance(q, ibulk, bratio, invsrat)
    >>> print(np.sum(q,axis=1))
    [ 1.  1.  1.  1.  1.]
    >>> print(q[:,1]/q[:,0])
    [ 0.25  0.25  0.25  0.25  1.  ]

.. py:function:: ratio(Q, ibulk)
.. code-block:: pycon

    Calculate the abundance ratios of the species indexed by ibulk, relative
    to the first species in the list.

    Parameters
    ----------
    Q: 2D float ndarray
       Mole mixing ratio of the species in the atmosphere [Nlayers, Nspecies].
    ibulk: 1D integer ndarray
       Indices of the species to calculate the ratio.

    Returns
    -------
    bratio: 2D float ndarray
       Abundance ratio between species indexed by ibulk.
    invsrat: 1D float ndarray
       Inverse of the sum of the ratios (at each layer).

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> q = np.tile([0.8, 0.2], (5,1))
    >>> q[4] = 0.5, 0.5
    >>> ibulk = [0, 1]
    >>> bratio, invsrat = pa.ratio(q, ibulk)
    >>> print(bratio)
    [[ 1.    0.25]
     [ 1.    0.25]
     [ 1.    0.25]
     [ 1.    0.25]
     [ 1.    1.  ]]
    >>> print(invsrat)
    [ 0.8  0.8  0.8  0.8  0.5]

.. py:function:: qscale(Q, spec, molmodel, molfree, molpars, bulk, qsat=None, iscale=None, ibulk=None, bratio=None, invsrat=None)
.. code-block:: pycon

    Scale specified species abundances and balance bulk abundances to
    conserve sum(Q)=1 in each layer.

    Parameters
    ----------
    Q:  2D float ndarray
       Mole mixing ratio of the species in the atmosphere [Nlayers, Nspecies].
    spec:  1D string ndarray
       Names of the species in the atmosphere.
    molmodel: 1D string ndarray
       Model to vary the species abundances.
    molfree:  1D string ndarray
       Names of the species to vary their abundances.
    molpars:  1D float ndarray
       Scaling factor (dex) for each species in molfree.
    bulk:  1D string ndarray
       Names of the bulk (dominant) species.
    qsat:  Float
       Maximum allowed combined abundance for trace species.
    iscale:  1D integer ndarray
       Indices of molfree species in Q.
    ibulk:  1D integer ndarray
       Indices of the bulk species in Q.
    bratio:  2D float ndarray
       Abundance ratios between the bulk species (relative to bulk[0]).
    invsrat:  1D float ndarray
       Inverse of the sum of the ratios (at each layer).

    Returns
    -------
    q:  2D float ndarray
       The modified atmospheric abundance profiles.

    Notes
    -----
    iscale, ibulk, bratio, and invsrat are optional parameters to speed up
       the routine.
    I'm not completely happy with qsat yet.


pyratbay.atmosphere.tmodels
___________________________


.. py:module:: pyratbay.atmosphere.tmodels

.. py:function:: tcea(tparams, pressure, rstar, tstar, tint, gplanet, smaxis, runits='cm')
.. code-block:: pycon

    Compute Three-channel Eddington Approximation (TCEA) temperature
    profile model.

    tparams: 1D iterable
        TCEA model parameters:
        log10(kappa):  Planck thermal IR opacity in units cm^2/gr
        log10(gamma1): Visible-to-thermal stream Planck mean opacity ratio
        log10(gamma2): Visible-to-thermal stream Planck mean opacity ratio
        alpha: Visible-stream partition (0.0--1.0)
        beta:  'catch-all' for albedo, emissivity, and day--night
               redistribution (on the order of unity)
    pressure: 1D float ndarray
        Atmospheric pressure profile in barye units.
    rstar: String or float
        Stellar radius (default in cm). If string, may specify units.
    tstar: String or float
        Stellar temperature in Kelvin degrees.
    tint: String or float
        Planetary internal temperature in Kelvin degrees.
    gplanet: String or float
        Planetary atmospheric temperature in cm s-2.
    smaxis: String or float
        Orbital semi-major axis (default in cm). If string, may specify units.
    runits: String
        Default units for rstar and smaxis.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants  as pc
    >>> nlayers = 11
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, 'bar')
    >>> rstar   = 0.756 * pc.rsun
    >>> tstar   = 5040.0  # K
    >>> tint    = 100.0   # K
    >>> gplanet = 2200.0   # cm s-2
    >>> smaxis  = 0.031 * pc.au
    >>> tparams = [-1.5, -0.8, -0.8, 0.5, 1.0]
    >>> temp = pa.tmodels.tcea(tparams, pressure, rstar, tstar, tint,
                               gplanet, smaxis)
    >>> print(temp)
    [1047.04157312 1047.04189805 1047.04531644 1047.08118784 1047.45648563
     1051.34469989 1088.69956369 1311.86379107 1640.12857767 1660.02396061
     1665.30121021]

.. py:function:: isothermal(tparams, nlayers)
.. code-block:: pycon

    Compute isothermal temperature profile.

    Parameters
    ----------
    tparams: scalar or iterable
        Temperature of the isothermal profile (Kelvin degree).
    nlayers: Integer
        Number of layers in temperature profile.

    Returns
    -------
    temp: 1D float ndarray
        Temperature profile.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> nlayers = 8
    >>> temp = pa.tmodels.isothermal(1500.0, nlayers)
    >>> print(temp)
    [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.]
    >>> temp = pa.tmodels.isothermal(np.array([1500.0]), nlayers)
    >>> print(temp)
    [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.]

.. py:function:: madhu_inv(params, pressure)
.. code-block:: pycon

    Calculates PT profile for inversion case based on Equation (2) from
    Madhusudhan & Seager 2009.

    Parameters
    ----------
    pressure: 1D float ndarray
        Pressure array, needs to be equally spaced in log space from bottom
        to top of the atmosphere. Must be given in bars.
    params: 1D float ndarray
        Temperature model parameters:
            a1 - float, exponential factor in Layer 1,
            a2 - float, exponential factor in Layer 2,
            p1 - floa, pressure boundary between Layer 1 and 2 (in bars).
            p2 - float, pressure in the middle of tLayer 2
            p3 - float, pressure boundary between Layers 2 and 3 (in bars).
            T3 - float, temperature in the Layer 3.

    Returns
    -------
    temp: 1D float ndarray
        Gaussian smoothed temperatures.

    Example
    -------
    >>> import pyratbay.atmosphere as pa
    >>> # array of pressures, equally spaced in log space
    >>> press = pa.pressure(1e-5, 1e2, 100, 'bar')
    >>> # Params [a1,  a2,   p1,     p2,   p3,  T3]
    >>> params = 0.51, 0.25, 4.5e-3, 0.01, 1.0, 1500.0
    >>> temp = pa.tmodels.madhu_inv(params, press/pc.bar)

.. py:function:: madhu_noinv(params, pressure)
.. code-block:: pycon

    Calculates PT profile for inversion case based on Equation (2) from
    Madhusudhan & Seager 2009.

    Parameters
    ----------
    pressure: 1D float ndarray
        Pressure array, needs to be equally spaced in log space from bottom
        to top of the atmosphere. Must be given in bars.
    params: 1D float ndarray
        Temperature model parameters:
            a1 - float, exponential factor in Layer 1,
            a2 - float, exponential factor in Layer 2,
            p1 - floa, pressure boundary between Layer 1 and 2 (in bars).
            p3 - float, pressure boundary between Layers 2 and 3 (in bars).
            T3 - float, temperature in the Layer 3.

    Returns
    -------
    T_smooth:  1D array of floats, Gaussian smoothed temperatures,
               no kinks on Layer boundaries

    Example
    -------
    >>> import pyratbay.atmosphere as pa
    >>> # array of pressures, equally spaced in log space
    >>> press = pa.pressure(1e-5, 1e2, 100, 'bar')
    >>> # Params [a1,  a2,   p1,     p3,  T3]
    >>> params = 0.51, 0.25, 4.5e-3, 1.0, 1500.0
    >>> temp = pa.tmodels.madhu_noinv(params, press/pc.bar)


pyratbay.atmosphere.clouds
__________________________


.. py:module:: pyratbay.atmosphere.clouds

.. py:class:: CCSgray()

.. code-block:: pycon

    Constant cross-section gray cloud model.

.. code-block:: pycon

    Initialize self.  See help(type(self)) for accurate signature.

.. py:class:: Deck()

.. code-block:: pycon

    Instantly opaque gray cloud deck at given pressure.

.. code-block:: pycon

    Initialize self.  See help(type(self)) for accurate signature.

.. py:function:: get_model(name)
.. code-block:: pycon

    Get a cloud model by its name.


pyratbay.atmosphere.rayleigh
____________________________


.. py:module:: pyratbay.atmosphere.rayleigh

.. py:class:: Dalgarno(mol)

.. code-block:: pycon

    Rayleigh-scattering model from Dalgarno (1962), Kurucz (1970), and
    Dalgarno & Williams (1962).

.. code-block:: pycon

    Parameters
    ----------
    mol: String
       The species, which can be H, He, or H2.

.. py:class:: Lecavelier()

.. code-block:: pycon

    Rayleigh-scattering model from Lecavelier des Etangs et al. (2008).
    AA, 485, 865.

.. code-block:: pycon

    Initialize self.  See help(type(self)) for accurate signature.

.. py:function:: get_model(name)
.. code-block:: pycon

    Get a Rayleigh model by its name.


pyratbay.atmosphere.alkali
__________________________


.. py:module:: pyratbay.atmosphere.alkali

.. py:class:: SodiumVdW()

.. code-block:: pycon

    Sodium Van der Waals model (Burrows et al. 2000, ApJ, 531).

.. code-block:: pycon

    Initialize self.  See help(type(self)) for accurate signature.

.. py:class:: PotassiumVdW()

.. code-block:: pycon

    Potassium Van der Waals model (Burrows et al. 2000, ApJ, 531).

.. code-block:: pycon

    Initialize self.  See help(type(self)) for accurate signature.

.. py:function:: get_model(name)
.. code-block:: pycon

    Get an alkali model by its name.

