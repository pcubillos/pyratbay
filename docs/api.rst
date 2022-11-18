.. _api:

API
===


pyratbay
________


.. py:module:: pyratbay

.. py:class:: Pyrat(cfile, no_logfile=False, mute=False)

    .. code-block:: pycon

        Main Pyrat object.


        Parse the command-line arguments into the pyrat object.

        Parameters
        ----------
        cfile: String
            A Pyrat Bay configuration file.
        no_logfile: Bool
            If True, enforce not to write outputs to a log file
            (e.g., to prevent overwritting log of a previous run).
        mute: Bool
            If True, enforce verb to take a value of -1.

        Examples
        --------
        >>> import pyratbay as pb
        >>> # Initialize and execute task:
        >>> pyrat = pb.run('spectrum_transmission.cfg')

        >>> # Initialize only:
        >>> pyrat = pb.Pyrat('spectrum_transmission.cfg')
        >>> # Then, setup internal varible for spectra evaluation:
        >>> pyrat.setup_spectrum()

    .. py:method:: band_integrate(spectrum=None)
    .. code-block:: pycon

        Band-integrate transmission spectrum (transit) or planet-to-star
        flux ratio (eclipse) over transmission band passes.

    .. py:method:: eval(params, retmodel=True, verbose=False)
    .. code-block:: pycon

        Fitting routine for MCMC.

        Parameters
        ----------
        params: 1D float iterable
           Array of fitting parameters that define the atmosphere.
        retmodel: Bool
           Flag to include the model spectra in the return.
        verbose: Bool
           Flag to print out if a run failed.

        Returns
        -------
        spectrum: 1D float ndarray
           The output model spectra.  Returned only if retmodel=True.
        bandflux: 1D float ndarray
           The waveband-integrated spectrum values.

    .. py:method:: get_ec(layer)
    .. code-block:: pycon

        Extract extinction-coefficient contribution (in cm-1) from each
        component of the atmosphere at the requested layer.

        Parameters
        ----------
        layer: Integer
           The index of the atmospheric layer where to extract the EC.

        Returns
        -------
        ec: 2D float ndarray
           An array of shape [ncomponents, nwave] with the EC spectra
           (in cm-1) from each component of the atmosphere.
        label: List of strings
           The names of each atmospheric component that contributed to EC.

    .. py:method:: hydro(pressure, temperature, mu, g, mass, p0, r0)
    .. code-block:: pycon

        Hydrostatic-equilibrium driver.
        Depending on self.atm.rmodelname, select between the g=GM/r**2
        (hydro_m) or constant-g (hydro_g) formula to compute
        the hydrostatic-equilibrium radii of the planet layers.

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
        mass: Float
           Planetary mass (in g).
        p0: Float
           Reference pressure level (in barye) where radius(p0) = r0.
        r0: Float
           Reference radius level (in cm) corresponding to p0.

    .. py:method:: percentile_spectrum(nmax=None)
    .. code-block:: pycon

        Compute spectrum posterior percentiles.

    .. py:method:: plot_spectrum(spec='model', **kwargs)
    .. code-block:: pycon

        Plot spectrum.

        Parameters
        ----------
        spec: String
            Flag indicating which model to plot.  By default plot the
            latest evaulated model (spec='model').  Other options are
            'best' or 'median' to plot the posterior best-fit or median
            model, in which case, the code will plot the 1- and 2-sigma
            boundaries if they have been computed (see
            self.percentile_spectrum).
        kwargs: dict
            Dictionary of arguments to pass into plots.spectrum().
            See help(pyratbay.plots.spectrum).

        Returns
        -------
        ax: AxesSubplot instance
            The matplotlib Axes of the figure.

    .. py:method:: plot_temperature(**kwargs)
    .. code-block:: pycon

        Plot temperature profile.
        If self.ret.posterior exitst, plot the best fit, median, and
        the '1sigma/2sigma' boundaries of the temperature posterior
        distribution.

        Parameters
        ----------
        kwargs: dict
            Dictionary of arguments to pass into plots.temperature().
            See help(pyratbay.plots.temperature).

        Returns
        -------
        ax: AxesSubplot instance
            The matplotlib Axes of the figure.

    .. py:method:: radiative_equilibrium(nsamples=None, continue_run=False, convection=False)
    .. code-block:: pycon

        Compute radiative-thermochemical equilibrium atmosphere.
        Currently there is no convergence criteria implemented,
        some 100--300 iterations are typically sufficient to converge
        to a stable temperature-profile solution.

        Parameters
        ----------
        nsamples: Integer
            Number of radiative-equilibrium iterations to run.
        continue_run: Bool
            If True, continue from a previous radiative-equilibrimu run.
        convection: Bool
            If True, skip convective flux calculation in the radiative
            equilibrium calculation.

        Returns
        -------
        There are no returned values, but this method updates the
        temperature profile (self.atm.temp) and abundances (self.atm.q)
        with the values from the last radiative-equilibrium iteration.

        This method also defines pyrat.atm.radeq_temps, a 2D array
        containing all temperature-profile iterations.

    .. py:method:: run(temp=None, abund=None, radius=None)
    .. code-block:: pycon

        Evaluate a Pyrat spectroscopic model

        Parameters
        ----------
        temp: 1D float ndarray
            Updated atmospheric temperature profile in Kelvin, of size nlayers.
        abund: 2D float ndarray
            Updated atmospheric abundances profile by number density, of
            shape [nlayers, nmol].
        radius: 1D float ndarray
            Updated atmospheric altitude profile in cm, of size nlayers.

    .. py:method:: set_filters()
    .. code-block:: pycon

        Set observational variables (pyrat.obs) based on given parameters.

.. py:function:: run(cfile, run_step='run', no_logfile=False)
.. code-block:: pycon

    Pyrat Bay initialization driver.

    Parameters
    ----------
    cfile: String
        A Pyrat Bay configuration file.
    run_step: String
        If 'dry': only read the configuration file into a Pyrat object
        If 'init': initialize a Pyrat object (but no spectra calculation).
        If 'run': run all the way (default).
    no_logfile: Bool
        If True, enforce not to write outputs to a log file
        (e.g., to prevent overwritting log of a previous run).


pyratbay.constants
__________________


.. py:module:: pyratbay.constants

.. py:data:: h
.. code-block:: pycon

  6.62607015e-27

.. py:data:: k
.. code-block:: pycon

  1.380649e-16

.. py:data:: c
.. code-block:: pycon

  29979245800.0

.. py:data:: G
.. code-block:: pycon

  6.674299999999999e-08

.. py:data:: sigma
.. code-block:: pycon

  5.670374419e-05

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

  14959787070000.0

.. py:data:: pc
.. code-block:: pycon

  3.0856775814913674e+18

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

  1.6605390666e-24

.. py:data:: me
.. code-block:: pycon

  9.1093837015e-28

.. py:data:: kelvin
.. code-block:: pycon

  1.0

.. py:data:: sec
.. code-block:: pycon

  1.0

.. py:data:: min
.. code-block:: pycon

  60.0

.. py:data:: hour
.. code-block:: pycon

  3600.0

.. py:data:: day
.. code-block:: pycon

  86400.0

.. py:data:: amagat
.. code-block:: pycon

  2.686780111e+19

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

  1129583353135.2844

.. py:data:: C2
.. code-block:: pycon

  1.4387768775039338

.. py:data:: C3
.. code-block:: pycon

  8.852821681767784e-13

.. py:data:: ROOT
.. code-block:: pycon

  os.path.realpath(os.path.dirname(__file__) + '/../..') + '/'

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

.. py:data:: dbases
.. code-block:: pycon

  ['Hitran', 'Exomol', 'Repack', 'Pands', 'Tioschwenke', 'Voplez', 'Vald']

.. py:data:: rmodes
.. code-block:: pycon

  ['tli', 'atmosphere', 'opacity', 'spectrum', 'radeq', 'mcmc']

.. py:data:: transmission_rt
.. code-block:: pycon

  ['transit']

.. py:data:: emission_rt
.. code-block:: pycon

  ['emission', 'emission_two_stream']

.. py:data:: rt_paths
.. code-block:: pycon

  ['transit', 'emission', 'emission_two_stream']

.. py:data:: retflags
.. code-block:: pycon

  ['temp', 'rad', 'press', 'mol', 'ray', 'cloud', 'patchy', 'mass']

.. py:data:: tmodels
.. code-block:: pycon

  ['isothermal', 'guillot', 'madhu']

.. py:data:: chemmodels
.. code-block:: pycon

  ['uniform', 'tea']

.. py:data:: radmodels
.. code-block:: pycon

  ['hydro_m', 'hydro_g']

.. py:data:: molmodels
.. code-block:: pycon

  ['vert', 'scale', 'equil']

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

.. py:function:: write_atm(atmfile, pressure, temperature, species=None, abundances=None, radius=None, punits='bar', runits=None, header=None)
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
    radius: 1D float ndarray
        Monotonously increasing radius profile (in cm).
    punits:  String
        Pressure units of output.
    runits:  String
        Radius units of output.
    header:  String
        Header message (comment) to include at the top of the file.

    Examples
    --------
    >>> import numpy as np
    >>> import pyratbay.io as io
    >>> import pyratbay.atmosphere as pa

    >>> atmfile = 'WASP-00b.atm'
    >>> nlayers = 5
    >>> pressure    = pa.pressure('1e-8 bar', '1e2 bar', nlayers)
    >>> temperature = pa.tmodels.Isothermal(nlayers)(1500.0)
    >>> species     = "H2 He H2O".split()
    >>> abundances  = [0.8499, 0.15, 1e-4]
    >>> qprofiles = pa.uniform(pressure, temperature, species, abundances)
    >>> io.write_atm(atmfile, pressure, temperature, species, qprofiles,
    >>>     punits='bar', header='# Example atmospheric file:\n')
    >>> # Print output file:
    >>> with open(atmfile, 'r') as f:
    >>>     print(f.read())
    # Example atmospheric file:
    # Pressure units:
    @PRESSURE
    bar
    # Temperatures units:
    @TEMPERATURE
    kelvin
    # Abundance units (mixing ratio):
    @ABUNDANCE
    volume

    # Atmospheric composition:
    @SPECIES
    H2  He  H2O

    # Pressure  Temperature  H2            He            H2O
    @DATA
    1.0000e-08     1500.000  8.499000e-01  1.500000e-01  1.000000e-04
    3.1623e-06     1500.000  8.499000e-01  1.500000e-01  1.000000e-04
    1.0000e-03     1500.000  8.499000e-01  1.500000e-01  1.000000e-04
    3.1623e-01     1500.000  8.499000e-01  1.500000e-01  1.000000e-04
    1.0000e+02     1500.000  8.499000e-01  1.500000e-01  1.000000e-04

.. py:function:: read_atm(atmfile)
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

    Examples
    --------
    >>> # Continuing example from io.write_atm():
    >>> import pyratbay.io as io

    >>> atmfile = 'WASP-00b.atm'
    >>> units, specs, pressure, temp, q, rad = io.read_atm(atmfile)
    >>> print(units, specs, pressure, temp, q, rad, sep='\n')
    ('bar', 'kelvin', 'number', None)
    ['H2' 'He' 'H2O']
    [1.0000e-08 3.1623e-06 1.0000e-03 3.1623e-01 1.0000e+02]
    [1500. 1500. 1500. 1500. 1500.]
    [[8.499e-01 1.500e-01 1.000e-04]
     [8.499e-01 1.500e-01 1.000e-04]
     [8.499e-01 1.500e-01 1.000e-04]
     [8.499e-01 1.500e-01 1.000e-04]
     [8.499e-01 1.500e-01 1.000e-04]]
    None

.. py:function:: write_spectrum(wl, spectrum, filename, type, wlunits='um')
.. code-block:: pycon

    Write a spectrum to file.

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
        - 'transit' for transmission
        - 'emission' for emission
        - 'filter' for a instrumental filter transmission
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

.. py:function:: write_opacity(ofile, species, temp, press, wn, opacity)
.. code-block:: pycon

    Write an opacity table as a binary file.

    Parameters
    ----------
    ofile: String
        Output filename where to save the opacity data.
        File extension must be .npz
    species: 1D string iterable
        Species names.
    temp: 1D float ndarray
        Temperature array (Kelvin degree).
    press: 1D float ndarray
        Pressure array (barye).
    wn: 1D float ndarray
        Wavenumber array (cm-1).
    opacity: 4D float ndarray
        Tabulated opacities (cm2 molecule-1) of shape
        [nspec, ntemp, nlayers, nwave].

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
        (nspec, ntemp, nlayers, nwave)
    arrays: 4-element 1D ndarray tuple
        The dimensions of the opacity table:
        - species     (string, the species names)
        - temperature (float, Kelvin)
        - pressure    (float, barye)
        - wavenumber  (float, cm-1)
    opacity: 4D float ndarray tuple
        The tabulated opacities (cm2 molecule-1), of shape
        [nspec, ntemp, nlayers, nwave].

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

.. py:function:: read_observations(obs_file)
.. code-block:: pycon

    Read an observations file.

    Parameters
    ----------
    obs_file: String
        Path to file containing observations info, see Notes below.

    Returns
    -------
    filters: List
        Filter passband objects.
    data: 1D string list
        The transit or eclipse depths for each filter.
    uncert: 1D float ndarray
        The depth uncertainties.

    Notes
    -----
    An obs_file contains passband info (indicated by a '@DATA' flag),
    one line per passband. The passband info could be:
    (1) a path to a file containing the spectral response, or
    (2) a tophat filter defined by a central wavelength, half-width,
    and optionally a name.

    A @DEPTH_UNITS flag sets the depth and uncert units (which can be
    set to: none, percent, ppt, ppm).
    This flag also indicates that there's data and uncerts to read
    as two columns before the passband info.

    Comment and blank lines are ignored,
    central-wavelength and half-width units are always microns.

    Examples
    --------
    >>> import pyratbay.io as io
    >>> # File including depths and uncertainties:
    >>> obs_file = 'observations.dat'
    >>> bands, data, uncert = io.read_observations(obs_file)

    >>> # File including only the passband info:
    >>> obs_file = 'filters.dat'
    >>> bands = io.read_observations(obs_file)

.. py:function:: read_atomic(afile)
.. code-block:: pycon

    Read an elemental (atomic) composition file.

    Parameters
    ----------
    afile: String
        File with atomic composition.

    Returns
    -------
    atomic_num: 1D integer ndarray
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

.. py:function:: read_molecs(file)
.. code-block:: pycon

    Read a molecules file to extract their names, masses, and radii.
    The output also includes the ions denoted by a '-' and '+'
    character appended at the end of the species names.

    Parameters
    ----------
    file: String
        The molecule file path.

    Returns
    -------
    names: 1D string ndarray
        The molecules' names.
    masses: 1D float ndarray
        The mass of the molecules (in g mol-1).
    radii: 1D float ndarray
        The collisional radius of the molecules (in angstrom).

    Notes
    -----
    In all truthfulness, these are species, not only molecules, as the
    file also contain elemental particles.

    Examples
    --------
    >>> import pyratbay.io as io
    >>> import pyratbay.constants as pc
    >>> names, masses, radii = io.read_molecs(
    >>>     pc.ROOT+'pyratbay/data/molecules.dat')
    >>> names = list(names)
    >>> print(f"H2O: mass = {masses[names.index('H2O')]} g mol-1, "
    >>>       f"radius = {radii[names.index('H2O')]} angstrom.")
    H2O: mass = 18.015 g mol-1, radius = 1.6 Angstrom.

.. py:function:: read_isotopes(file)
.. code-block:: pycon

    Read an isotopes file to extract their molecule, hitran name,
    exomol name, isotopic ratio, and mass.

    Parameters
    ----------
    file: String
        The isotope file path.

    Returns
    -------
    mol_ID: 1D integer ndarray
        HITRAN molecule ID.
    mol: 1D string ndarray
        Molecule names.
    hitran_iso: 1D string ndarray
        Isotope name as in HITRAN database.
    exomol_iso: 1D string ndarray
        Isotope name based on exomol database.
    iso_ratio: 1D float ndarray
        Isotopic ratios.
    iso_mass: 1D float ndarray
        The mass of the molecules (in g mol-1).

    Examples
    --------
    >>> import pyratbay.io as io
    >>> import pyratbay.constants as pc
    >>> ID, mol, hit_iso, exo_iso, ratio, mass = \
    >>>     io.read_isotopes(pc.ROOT+'pyratbay/data/isotopes.dat')
    >>> print("H2O isotopes:\n iso    iso    isotopic  mass"
    >>>                    "\n hitran exomol ratio     g/mol")
    >>> for i in range(len(mol)):
    >>>     if mol[i] == 'H2O':
    >>>         print(f" {hit_iso[i]:6} {exo_iso[i]:6} "
    >>>               f"{ratio[i]:.3e} {mass[i]:.4f}")
    H2O isotopes:
    iso    iso    isotopic  mass
    hitran exomol ratio     g/mol
    161    116    9.973e-01 18.0106
    181    118    1.999e-03 20.0148
    171    117    3.719e-04 19.0148
    162    126    3.107e-04 19.0168
    182    000    6.230e-07 21.0211
    172    000    1.158e-07 20.0211
    262    226    2.420e-08 20.0210
    282    000    0.000e+00 22.0000
    272    000    0.000e+00 21.0000

.. py:function:: import_xs(filename, source, read_all=True, ofile=None)
.. code-block:: pycon

    Read a cross-section opacity file from an external source.

    Parameters
    ----------
    filename: String
        The opacity pickle file to read.
    source: String
        The cross-section source: exomol or taurex (see note below).
    read_all: Bool
        If True, extract all contents in the file: cross-section,
        pressure, temperature, and wavenumber.
        If False, extract only the cross-section data.
    ofile: String
        If not None, store Exomol XS data into a Pyratbay opacity
        format.

    Returns
    -------
    xs: 3D float ndarray
        Opacity cross-section in cm2 molecule-1.
        with shape [npress, ntemp, nwave].
    pressure: 1D float ndarray
        Pressure sample of the opacity file (in barye units)
    temperature: 1D float ndarray
        Temperature sample of the opacity file (in Kelvin degrees units).
    wavenumber: 1D float ndarray
        Wavenumber sample of the opacity file (in cm-1 units).
    species: String
        The species name.

    Notes
    -----
    - exomol cross sections (Chubb et al. 2020, AA) can be found here:
    http://www.exomol.com/data/data-types/opacity/
    - taurex cross sections (Al-Refaie et al. 2019) can be found here:
    https://taurex3-public.readthedocs.io/en/latest/user/taurex/quickstart.html

    Examples
    --------
    >>> # For this example, you'll need to have/download the following
    >>> # file into the current folder:
    >>> # http://www.exomol.com/db/H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL__R15000_0.3-50mu.xsec.TauREx.h5
    >>> import pyratbay.io as io
    >>> filename = '1H2-16O__POKAZATEL__R15000_0.3-50mu.xsec.TauREx.h5'
    >>> xs_H2O, press, temp, wn, species = io.import_xs(filename, 'exomol')

.. py:function:: import_tea(teafile, atmfile, req_species=None)
.. code-block:: pycon

    Format a TEA atmospheric file into a Pyrat atmospheric file.

    Paramters
    ---------
    teafile:  String
        Input TEA atmospheric file.
    atmfile:  String
        Output Pyrat atmospheric file.
    req_species: List of strings
        The requested species for output.  If None, request all species
        in teafile.

.. py:function:: export_pandexo(pyrat, baseline, transit_duration, Vmag=None, Jmag=None, Hmag=None, Kmag=None, metal=0.0, instrument=None, n_transits=1, resolution=None, noise_floor=0.0, sat_level=80.0, save_file=True)
.. code-block:: pycon

    Parameters
    ----------
    pyrat: A Pyrat instance
        Pyrat object from which to extract the system physical properties.
    baseline: Float or string
        Total observing time in sec (float) or with given units (string).
    transit_duration: Float or string
        Transit/eclipse duration in sec (float) or with given units (string).
    metal: Float
        Stellar metallicity as log10(Fe/H).
    Vmag: Float
        Stellar magnitude in the Johnson V band.
        Only one of Vmag, Jmag, Hmag, or Kmag should be defined.
    Jmag: Float
        Stellar magnitude in the Johnson J band.
        Only one of Vmag, Jmag, Hmag, or Kmag should be defined.
    Hmag: Float
        Stellar magnitude in the Johnson H band.
        Only one of Vmag, Jmag, Hmag, or Kmag should be defined.
    Kmag: Float
        Stellar magnitude in the Johnson Kband.
        Only one of Vmag, Jmag, Hmag, or Kmag should be defined.
    instrument: String or list of strings or dict
        Observing instrument to simulate.
        If None, this function returns the input dictionary.
    n_transits: Integer
        Number of transits/eclipses.
    resolution: Float
        Approximate output spectral sampling R = 0.5*lambda/delta-lambda.
    sat_level: Float
        Saturation level in percent of full well.
    noise_floor: Float or string
        Noise-floor level in ppm at all wavelengths (if float) or
        wavelength dependent (if string, filepath).
    save_file: Bool or string
        If string, store pandexo output pickle file with this filename.
        If True, store pandexo output with default name based on
        the pyrat object's output filename.

    Returns
    -------
    pandexo_sim: dict
        Output from pandexo.engine.justdoit.run_pandexo().
        Note this dict has R=None, noccultations=1 (as suggested in pandexo).
    wavelengths: List of 1D float arrays
        Wavelengths of simulated observed spectra for each instrument.
        Returned only if instrument is not None.
    spectra: List of 1D float arrays
        Simulated observed spectra for each instrument.
        Returned only if instrument is not None.
    uncertainties: List of 1D float arrays
        Uncertainties of simulated observed spectra for each instrument.
        Returned only if instrument is not None.

    Examples
    --------
    >>> import pyratbay as pb
    >>> import pyratbay.io as io

    >>> pyrat = pb.run('demo_spectrum-transmission.cfg')
    >>> instrument = 'NIRCam F322W2'
    >>> #instrument = jdi.load_mode_dict(instrument)
    >>> baseline = '4.0 hour'
    >>> transit_duration = '2.0 hour'
    >>> resolution = 100.0
    >>> n_transits = 2
    >>> Jmag = 8.0
    >>> metal = 0.0

    >>> pandexo_sim, wls, spectra, uncerts = io.export_pandexo(
    >>>     pyrat, baseline, transit_duration,
    >>>     n_transits=n_transits,
    >>>     resolution=resolution,
    >>>     instrument=instrument,
    >>>     Jmag=Jmag,
    >>>     metal=metal)


pyratbay.tools
______________


.. py:module:: pyratbay.tools

.. py:function:: log_error(log=None, error=None)
.. code-block:: pycon

    Capture exceptions into a log.error() call.

.. py:function:: cd(newdir)
.. code-block:: pycon

    Context manager for changing the current working directory.
    Taken from here: https://stackoverflow.com/questions/431684/

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

.. py:function:: u(units)
.. code-block:: pycon

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

.. py:function:: get_param(param, units='none', gt=None, ge=None)
.. code-block:: pycon

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

    .. py:method:: write(text, *format, **numpy_fmt)
    .. code-block:: pycon

        Write formatted text.
        See __init__ arguments for avaiable numpy_fmt items.

.. py:class:: Timer()

    .. code-block:: pycon

        Timer to get the time (in seconds) since the last call.


        Initialize self.  See help(type(self)) for accurate signature.

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

.. py:function:: radius_to_depth(rprs, rprs_err)
.. code-block:: pycon

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

.. py:function:: depth_to_radius(depth, depth_err)
.. code-block:: pycon

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

.. py:function:: ignore_system_exit(func)
.. code-block:: pycon

    Decorator to ignore SystemExit exceptions.

.. py:class:: Namespace(args=None, log=None)

    .. code-block:: pycon

        A container object to hold variables.


        Initialize self.  See help(type(self)) for accurate signature.

    .. py:method:: get_default(pname, desc, default=None, wflag=False, gt=None, ge=None, lt=None, le=None)
    .. code-block:: pycon

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

    .. py:method:: get_path(pname, desc='', exists=False)
    .. code-block:: pycon

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

    .. py:method:: get_units(pname)
    .. code-block:: pycon

        Extract units from a value input.
        Return None if value does not have units or has an invalid format.

        Parameters
        ----------
        pname: String
            Parameter name.

        Returns
        -------
        units: String

.. py:function:: parse(pyrat, cfile, no_logfile=False, mute=False)
.. code-block:: pycon

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
    >>> args = {f'par{i}':val for i,val in enumerate(inputs)}
    >>> for i,var in enumerate(inputs):
    >>>     par = f'par{i}'
    >>>     pt.parse_array(args, par)
    >>>     print(f"{par}: {repr(var)} -> {repr(args[par])}")
    par0: '10 20' -> array([10., 20.])
    par1: '10.0 20.0' -> array([10., 20.])
    par2: 'a b' -> ['a', 'b']
    par3: 'a\n b' -> ['a', 'b']


pyratbay.opacity
________________


.. py:module:: pyratbay.opacity

.. py:function:: make_tli(dblist, pflist, dbtype, tlifile, wllow, wlhigh, wlunits, log)
.. code-block:: pycon

    Create a transition-line-information (TLI) file.

    Parameters
    ----------
    dblist: List of strings
        Opacity databases to read.
    pflist: List of strings
        Partition function for each of the databases.
    dbtype: List of strings
        Database type of each database.
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
        An mc3.utils.Log instance to log screen outputs to file.


pyratbay.opacity.linelist
_________________________


.. py:module:: pyratbay.opacity.linelist

.. py:class:: Hitran(dbfile, pffile, log)

    .. code-block:: pycon

        HITRAN/HITEMP database reader.


        Initialize HITRAN database object.

        Parameters
        ----------
        dbfile: String
            File with the Database line-transition info.
        pffile: String
            File with the partition function.
        log: Log object
            An mc3.utils.Log instance to log screen outputs to file.

    .. py:method:: binsearch(dbfile, wave, ilo, ihi, searchup=True)
    .. code-block:: pycon

        Do a binary (and then linear) search for wavelength/wavenumber in
        file 'dbfile' between record positions ilo and ihi.

        Parameters
        ----------
        dbfile: File object
           File where to search.
        wave: Scalar
           Target wavelength/wavenumber (as given in each specific database).
        ilo: Integer
           Lowest index record to search.
        ihi: Integer
           highest index record to search.
        searchup: Boolean
           Search up (True) or down (False) the records for duplicate results
           after the binary search.

        Returns:
        --------
        irec:  Integer
           Record index for wave.

    .. py:method:: dbread(iwn, fwn, verb)
    .. code-block:: pycon

        Read line-transition info between wavenumbers iwn and fwn.

        Parameters
        ----------
        iwn: Float
            Lower wavenumber boundary in cm-1.
        fwn: Float
            Upper wavenumber boundary in cm-1.
        verb: Integer
            Verbosity threshold.

        Returns
        -------
        wnumber: 1D float ndarray
            Line-transition central wavenumber (cm-1).
        gf: 1D float ndarray
            gf value (unitless).
        elow: 1D float ndarray
            Lower-state energy (cm-1).
        isoID: 1D integer ndarray
            Isotope index.

    .. py:method:: get_iso(molname, dbtype)
    .. code-block:: pycon

        Get isotopic info from isotopes.dat file.

        Parameters
        ----------
        mol: String
            If not None, extract data based on this molecule name.
        dbtype: String
            Database type (for isotope names).

        Returns
        -------
        molID: Integer
            HITRAN molecule ID.
        isotopes: List of strings
            Isotopes names.
        mass: List of floats
            Masses for each isotope.
        isoratio: List of integers
            Isotopic terrestrial abundance ratio.

    .. py:method:: getpf(verbose=0)
    .. code-block:: pycon

        Compute partition function for specified source.

        Returns
        -------
        temp: 1D float ndarray
            Array with temperature sample.
        PF: 2D float ndarray
            The partition function data for each isotope at each temperature.
        isotopes: List of strings
            The names of the tabulated isotopes

    .. py:method:: readwave(dbfile, irec)
    .. code-block:: pycon

        Read irec-th wavenumber record from FILE dbfile.

        Parameters
        ----------
        dbfile: File object
            File where to extract the wavenumber.
        irec: Integer
            Index of record.

        Returns
        -------
        wavenumber: Float
            Wavenumber value in cm-1.

.. py:class:: Exomol(dbfile, pffile, log)

    .. code-block:: pycon

        Exomol database reader.


        Initialize Exomol database object.

        Parameters
        ----------
        dbfile: String
            File with the Database line-transition info.
        pffile: String
            File with the partition function.
        log: Log object
            An mc3.utils.Log instance to log screen outputs to file.

    .. py:method:: binsearch(dbfile, wave, ilo, ihi, searchup=True)
    .. code-block:: pycon

        Do a binary (and then linear) search for wavelength/wavenumber in
        file 'dbfile' between record positions ilo and ihi.

        Parameters
        ----------
        dbfile: File object
           File where to search.
        wave: Scalar
           Target wavelength/wavenumber (as given in each specific database).
        ilo: Integer
           Lowest index record to search.
        ihi: Integer
           highest index record to search.
        searchup: Boolean
           Search up (True) or down (False) the records for duplicate results
           after the binary search.

        Returns:
        --------
        irec:  Integer
           Record index for wave.

    .. py:method:: dbread(iwn, fwn, verb)
    .. code-block:: pycon

        Read line-transition info between wavenumbers iwn and fwn.

        Parameters
        ----------
        iwn: Float
            Lower wavenumber boundary in cm-1.
        fwn: Float
            Upper wavenumber boundary in cm-1.
        verb: Integer
            Verbosity threshold.

        Returns
        -------
        wnumber: 1D float ndarray
            Line-transition central wavenumber (cm-1).
        gf: 1D float ndarray
            gf value (unitless).
        elow: 1D float ndarray
            Lower-state energy (cm-1).
        isoID: 1D integer ndarray
            Isotope index.

    .. py:method:: get_iso(molname, dbtype)
    .. code-block:: pycon

        Get isotopic info from isotopes.dat file.

        Parameters
        ----------
        mol: String
            If not None, extract data based on this molecule name.
        dbtype: String
            Database type (for isotope names).

        Returns
        -------
        molID: Integer
            HITRAN molecule ID.
        isotopes: List of strings
            Isotopes names.
        mass: List of floats
            Masses for each isotope.
        isoratio: List of integers
            Isotopic terrestrial abundance ratio.

    .. py:method:: getpf(verbose=0)
    .. code-block:: pycon

        Compute partition function for specified source.

        Returns
        -------
        temp: 1D float ndarray
            Array with temperature sample.
        PF: 2D float ndarray
            The partition function data for each isotope at each temperature.
        isotopes: List of strings
            The names of the tabulated isotopes

    .. py:method:: readwave(dbfile, irec)
    .. code-block:: pycon

        Read irec-th wavenumber record from FILE dbfile.

        Parameters
        ----------
        dbfile: File object
            File where to extract the wavenumber.
        irec: Integer
            Index of record.

        Returns
        -------
        wavenumber: Float
            Wavenumber value in cm-1.

.. py:class:: Repack(dbfile, pffile, log)

    .. code-block:: pycon

        Repack database reader.


        Initialize Exomol database object.

        Parameters
        ----------
        dbfile: String
            File with the Database line-transition info.
        pffile: String
            File with the partition function.
        log: Log object
            An mc3.utils.Log instance to log screen outputs to file.

    .. py:method:: binsearch(dbfile, wave, ilo, ihi, searchup=True)
    .. code-block:: pycon

        Do a binary (and then linear) search for wavelength/wavenumber in
        file 'dbfile' between record positions ilo and ihi.

        Parameters
        ----------
        dbfile: File object
           File where to search.
        wave: Scalar
           Target wavelength/wavenumber (as given in each specific database).
        ilo: Integer
           Lowest index record to search.
        ihi: Integer
           highest index record to search.
        searchup: Boolean
           Search up (True) or down (False) the records for duplicate results
           after the binary search.

        Returns:
        --------
        irec:  Integer
           Record index for wave.

    .. py:method:: dbread(iwn, fwn, verb)
    .. code-block:: pycon

        Read line-transition info between wavenumbers iwn and fwn.

        Parameters
        ----------
        iwn: Float
            Lower wavenumber boundary in cm-1.
        fwn: Float
            Upper wavenumber boundary in cm-1.
        verb: Integer
            Verbosity threshold.

        Returns
        -------
        wnumber: 1D float ndarray
            Line-transition central wavenumber (cm-1).
        gf: 1D float ndarray
            gf value (unitless).
        elow: 1D float ndarray
            Lower-state energy (cm-1).
        isoID: 1D integer ndarray
            Isotope index.

    .. py:method:: get_iso(molname, dbtype)
    .. code-block:: pycon

        Get isotopic info from isotopes.dat file.

        Parameters
        ----------
        mol: String
            If not None, extract data based on this molecule name.
        dbtype: String
            Database type (for isotope names).

        Returns
        -------
        molID: Integer
            HITRAN molecule ID.
        isotopes: List of strings
            Isotopes names.
        mass: List of floats
            Masses for each isotope.
        isoratio: List of integers
            Isotopic terrestrial abundance ratio.

    .. py:method:: getpf(verbose=0)
    .. code-block:: pycon

        Compute partition function for specified source.

        Returns
        -------
        temp: 1D float ndarray
            Array with temperature sample.
        PF: 2D float ndarray
            The partition function data for each isotope at each temperature.
        isotopes: List of strings
            The names of the tabulated isotopes

    .. py:method:: readwave(dbfile, irec)
    .. code-block:: pycon

        Read irec-th wavenumber record from FILE dbfile.

        Parameters
        ----------
        dbfile: File object
            File where to extract the wavenumber.
        irec: Integer
            Index of record.

        Returns
        -------
        wavenumber: Float
            Wavenumber value in cm-1.

.. py:class:: Pands(dbfile, pffile, log)

    .. code-block:: pycon

        Partridge & Schwenke (1997) H2O database reader.


        Initialize P&S database object.

        Parameters
        ----------
        dbfile: String
            File with the Database line-transition info.
        pffile: String
            File with the partition function.
        log: Log object
            An mc3.utils.Log instance to log screen outputs to file.

    .. py:method:: binsearch(dbfile, wave, ilo, ihi, searchup=True)
    .. code-block:: pycon

        Do a binary (and then linear) search for wavelength/wavenumber in
        file 'dbfile' between record positions ilo and ihi.

        Parameters
        ----------
        dbfile: File object
           File where to search.
        wave: Scalar
           Target wavelength/wavenumber (as given in each specific database).
        ilo: Integer
           Lowest index record to search.
        ihi: Integer
           highest index record to search.
        searchup: Boolean
           Search up (True) or down (False) the records for duplicate results
           after the binary search.

        Returns:
        --------
        irec:  Integer
           Record index for wave.

    .. py:method:: dbread(iwn, fwn, verb)
    .. code-block:: pycon

        Read line-transition info between wavenumbers iwn and fwn.

        Parameters
        ----------
        iwn: Float
            Lower wavenumber boundary in cm-1.
        fwn: Float
            Upper wavenumber boundary in cm-1.
        verb: Integer
            Verbosity threshold.

        Returns
        -------
        wnumber: 1D float ndarray
            Line-transition central wavenumber (cm-1).
        gf: 1D float ndarray
            gf value (unitless).
        elow: 1D float ndarray
            Lower-state energy (cm-1).
        isoID: 1D integer ndarray
            Isotope index.

    .. py:method:: get_iso(molname, dbtype)
    .. code-block:: pycon

        Get isotopic info from isotopes.dat file.

        Parameters
        ----------
        mol: String
            If not None, extract data based on this molecule name.
        dbtype: String
            Database type (for isotope names).

        Returns
        -------
        molID: Integer
            HITRAN molecule ID.
        isotopes: List of strings
            Isotopes names.
        mass: List of floats
            Masses for each isotope.
        isoratio: List of integers
            Isotopic terrestrial abundance ratio.

    .. py:method:: getpf(verbose=0)
    .. code-block:: pycon

        Compute partition function for specified source.

        Returns
        -------
        temp: 1D float ndarray
            Array with temperature sample.
        PF: 2D float ndarray
            The partition function data for each isotope at each temperature.
        isotopes: List of strings
            The names of the tabulated isotopes

    .. py:method:: readwave(dbfile, irec)
    .. code-block:: pycon

        Read irec-th wavelength record from FILE dbfile.

        Parameters
        ----------
        dbfile: File object
            File where to extract the wavelength.
        irec: Integer
            Index of record.

        Returns
        -------
        recwl: Unsigned integer
            Wavelength value as given in the P&S binary file.

.. py:class:: Tioschwenke(dbfile, pffile, log)

    .. code-block:: pycon

        Notes:
        ------
        Linelist and partition function downloaded from:
          http://kurucz.harvard.edu/molecules/tio/tioschwenke.bin
          http://kurucz.harvard.edu/molecules/tio/tiopart.dat

        There might be a problem with the linebreak character of the partition
        function.  One way to fix is, on vim do: :%s/        /        /g


        Initialize self.  See help(type(self)) for accurate signature.

    .. py:method:: binsearch(dbfile, wave, ilo, ihi, searchup=True)
    .. code-block:: pycon

        Do a binary (and then linear) search for wavelength/wavenumber in
        file 'dbfile' between record positions ilo and ihi.

        Parameters
        ----------
        dbfile: File object
           File where to search.
        wave: Scalar
           Target wavelength/wavenumber (as given in each specific database).
        ilo: Integer
           Lowest index record to search.
        ihi: Integer
           highest index record to search.
        searchup: Boolean
           Search up (True) or down (False) the records for duplicate results
           after the binary search.

        Returns:
        --------
        irec:  Integer
           Record index for wave.

    .. py:method:: dbread(iwn, fwn, verb)
    .. code-block:: pycon

        Read the Schwenke TiO database.

        Parameters
        ----------
        iwn: Scalar
            Initial wavenumber limit (in cm-1).
        fwn: Scalar
            Final wavenumber limit (in cm-1).
        verb: Integer
            Verbosity threshold.

        Returns
        -------
        wnumber: 1D float ndarray
            Line-transition central wavenumber (centimeter-1).
        gf: 1D float ndarray
            gf value (unitless).
        elow: 1D float ndarray
            Lower-state energy (centimeter-1).
        isoID: 2D integer ndarray
            Isotope index (0, 1, 2, 3, ...).

    .. py:method:: get_iso(molname, dbtype)
    .. code-block:: pycon

        Get isotopic info from isotopes.dat file.

        Parameters
        ----------
        mol: String
            If not None, extract data based on this molecule name.
        dbtype: String
            Database type (for isotope names).

        Returns
        -------
        molID: Integer
            HITRAN molecule ID.
        isotopes: List of strings
            Isotopes names.
        mass: List of floats
            Masses for each isotope.
        isoratio: List of integers
            Isotopic terrestrial abundance ratio.

    .. py:method:: getpf(verbose=0)
    .. code-block:: pycon

        Compute partition function for specified source.

        Returns
        -------
        temp: 1D float ndarray
            Array with temperature sample.
        PF: 2D float ndarray
            The partition function data for each isotope at each temperature.
        isotopes: List of strings
            The names of the tabulated isotopes

    .. py:method:: readwave(dbfile, irec)
    .. code-block:: pycon

        Read wavelength parameter from irec record in dbfile database.

        Parameters
        ----------
        dbfile: File object
           File where to extract the wavelength.
        irec: Integer
           Index of record.

        Returns
        -------
        rec_wl: integer
           Wavelength value at record irec, as given in dbfile database.

.. py:class:: Voplez(dbfile, pffile, log)

    .. code-block:: pycon

        Download the linelist from:


        Initializer.

    .. py:method:: binsearch(dbfile, wave, ilo, ihi, searchup=True)
    .. code-block:: pycon

        Do a binary (and then linear) search for wavelength/wavenumber in
        file 'dbfile' between record positions ilo and ihi.

        Parameters
        ----------
        dbfile: File object
           File where to search.
        wave: Scalar
           Target wavelength/wavenumber (as given in each specific database).
        ilo: Integer
           Lowest index record to search.
        ihi: Integer
           highest index record to search.
        searchup: Boolean
           Search up (True) or down (False) the records for duplicate results
           after the binary search.

        Returns:
        --------
        irec:  Integer
           Record index for wave.

    .. py:method:: dbread(iwn, fwn, verb)
    .. code-block:: pycon

        Read the B. Plez VO database between the wavelengths iwl and fwl.

        Parameters:
        -----------
        iwn: Scalar
           Initial wavenumber limit (in cm-1).
        fwn: Scalar
           Final wavenumber limit (in cm-1).
        verb: Integer
           Verbosity threshold.

        Returns:
        --------
        wnumber: 1D float ndarray
          Line-transition central wavenumber (centimeter-1).
        gf: 1D float ndarray
          gf value (unitless).
        elow: 1D float ndarray
          Lower-state energy (centimeter-1).
        isoID: 2D integer ndarray
          Isotope index (0, 1, 2, 3, ...).

        Developers:
        -----------
        Patricio Cubillos (UCF).
        Sarah Blumenthal (UCF).

        Notes:
        ------
        The Plez VO database is an ASCII format.
        The line transitions are sorted in increasing wavelength (micron) order.

    .. py:method:: get_iso(molname, dbtype)
    .. code-block:: pycon

        Get isotopic info from isotopes.dat file.

        Parameters
        ----------
        mol: String
            If not None, extract data based on this molecule name.
        dbtype: String
            Database type (for isotope names).

        Returns
        -------
        molID: Integer
            HITRAN molecule ID.
        isotopes: List of strings
            Isotopes names.
        mass: List of floats
            Masses for each isotope.
        isoratio: List of integers
            Isotopic terrestrial abundance ratio.

    .. py:method:: getpf(verbose=0)
    .. code-block:: pycon

        Compute partition function for specified source.

        Returns
        -------
        temp: 1D float ndarray
            Array with temperature sample.
        PF: 2D float ndarray
            The partition function data for each isotope at each temperature.
        isotopes: List of strings
            The names of the tabulated isotopes

    .. py:method:: readwave(dbfile, irec)
    .. code-block:: pycon

        Extract the wavelength from record irec.

        Parameters:
        -----------
        dbfile: File pointer
           Pointer to file being read.
        irec: Integer
           Index of record to read.

        Returns:
        --------
        wl: Float
           The wavelength in microns for record irec.

.. py:class:: Vald(dbfile, pffile, log)

    .. code-block:: pycon

        Notes
        -----
          Download linelist from: http://vald.astro.uu.se/~vald/php/vald.php
             Selecting 'Extract Element' and 'Short format'.
          Download partition functions from:


        Initialize Basic data for the Database.

        Parameters
        ----------
        dbfile: String
            File with the Database line-transition info.
        pffile: String
            File with the partition function.
        log: File
            File object to store the log.

    .. py:method:: binsearch(dbfile, wave, ilo, ihi, searchup=True)
    .. code-block:: pycon

        Do a binary (and then linear) search for wavelength/wavenumber in
        file 'dbfile' between record positions ilo and ihi.

        Parameters
        ----------
        dbfile: File object
           File where to search.
        wave: Scalar
           Target wavelength/wavenumber (as given in each specific database).
        ilo: Integer
           Lowest index record to search.
        ihi: Integer
           highest index record to search.
        searchup: Boolean
           Search up (True) or down (False) the records for duplicate results
           after the binary search.

        Returns:
        --------
        irec:  Integer
           Record index for wave.

    .. py:method:: dbread(iwn, fwn, verb)
    .. code-block:: pycon

        Read a VALD database.

        Parameters
        ----------
        iwn: Scalar
            Initial wavenumber limit (in cm-1).
        fwn: Scalar
            Final wavenumber limit (in cm-1).
        verb: Integer
            Verbosity threshold.

        Returns
        -------
        wnumber: 1D float ndarray
            Line-transition central wavenumber (cm-1).
        gf: 1D float ndarray
            gf value (unitless).
        elow: 1D float ndarray
            Lower-state energy (cm-1).
        isoID: 2D integer ndarray
          Isotope index.

    .. py:method:: get_iso(molname, dbtype)
    .. code-block:: pycon

        Get isotopic info from isotopes.dat file.

        Parameters
        ----------
        mol: String
            If not None, extract data based on this molecule name.
        dbtype: String
            Database type (for isotope names).

        Returns
        -------
        molID: Integer
            HITRAN molecule ID.
        isotopes: List of strings
            Isotopes names.
        mass: List of floats
            Masses for each isotope.
        isoratio: List of integers
            Isotopic terrestrial abundance ratio.

    .. py:method:: getinfo()
    .. code-block:: pycon

        Doc me.

    .. py:method:: getpf(verbose=0)
    .. code-block:: pycon

        Compute partition function for specified source.

        Returns
        -------
        temp: 1D float ndarray
            Array with temperature sample.
        PF: 2D float ndarray
            The partition function data for each isotope at each temperature.
        isotopes: List of strings
            The names of the tabulated isotopes

    .. py:method:: readwave(dbfile, irec)
    .. code-block:: pycon

        Read irec-th wavenumber record from FILE dbfile.

        Parameters
        ----------
        dbfile: File object
            File where to extract the wavelength.
        irec: Integer
            Index of record.

        Returns
        -------
        wavenumber: Unsigned integer
            Wavenumber value in cm-1.


pyratbay.opacity.partitions
___________________________


.. py:module:: pyratbay.opacity.partitions

.. py:function:: get_tips_molname(molID)
.. code-block:: pycon

    Get the TIPS molecule name for given molecule ID.

    Parameters
    ----------
    molID: Integer
        HITRAN molecule ID. See for example: https://hitran.org/lbl/

    Returns
    -------
    molname: String
        Name of molecule.

    Examples
    --------
    >>> import pyratbay.opacity.partitions as pf
    >>> print(pf.get_tips_molname(1), pf.get_tips_molname(6))
    H2O CH4

.. py:function:: tips(molecule, isotopes=None, outfile=None, db_type='as_tips')
.. code-block:: pycon

    Extract TIPS 2021 partition-function values for given molecule.
    If requested, write the partition-function into a file for use
    with Pyrat Bay.
    References:
        Gamache et al. (2017), JQSRT, 203, 70.
        Gamache et al. (2021), JQSRT, 271, 107713.

    Parameters
    ----------
    molecule: String
        Name of the molecule.
    isotopes: String or list of strings
        If not None, only extract the requested isotopes.
    outfile: String
        If not None, save output to file.
        If outfile == 'default', save output to file named as
        PF_tips_molecule.dat
    db_type: String
        If db_type == 'as_exomol', return isotopic names following
        the exomol notation.

    Returns
    -------
    pf: 2D float ndarray
        TIPS partition function for input molecule.
    isotopes: 1D string list
        List of isotopes.
    temp: 1D float ndarray
        Partition-function temperature samples (K).

    Examples
    --------
    >>> import pyratbay.opacity.partitions as pf
    >>> pf_data, isotopes, temp = pf.tips('H2O', outfile='default')

    Written partition-function file:
      'PF_tips_H2O.dat'
    for molecule H2O, with isotopes ['161', '181', '171', '162', '182', '172', '262', '282', '272'],
    and temperature range 1--5000 K.

.. py:function:: exomol(pf_files, outfile=None)
.. code-block:: pycon

    Extract ExoMol partition-function values from input files.
    If requested, write the partition-function into a file for use
    with Pyrat Bay.

    Parameters
    ----------
    pf_files: String or List of strings
        Input Exomol partition-function filenames.  If there are
        multiple isotopes, all of them must correspond to the same
        molecule.
    outfile: String
        If not None, save output to file.
        If outfile == 'default', save output to file named as
        PF_exomol_molecule.dat

    Returns
    -------
    pf: 2D float ndarray
        TIPS partition function for input molecule.
    isotopes: 1D string list
        List of isotopes.
    temp: 1D float ndarray
        Partition-function temperature samples (K).

    Examples
    --------
    >>> # First, download ExoMol data to current dictory, e.g.:
    >>> # wget http://www.exomol.com/db/NH3/14N-1H3/BYTe/14N-1H3__BYTe.pf
    >>> # wget http://www.exomol.com/db/NH3/15N-1H3/BYTe-15/15N-1H3__BYTe-15.pf
    >>> import pyratbay.opacity.partitions as pf
    >>> # A single file:
    >>> pf_data, isotopes, temp = pf.exomol('14N-1H3__BYTe.pf',
    >>>     outfile='default')
    Written partition-function file:
      'PF_exomol_NH3.dat'
    for molecule NH3, with isotopes ['4111'],
    and temperature range 1--1600 K.

    >>> # Multiple files (isotopes) for a molecule:
    >>> pf_data, isotopes, temp = pf.exomol(
    >>>     ['14N-1H3__BYTe.pf', '15N-1H3__BYTe-15.pf'], outfile='default')

    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Warning:
        Length of PF files do not match.  Trimming to shorter size.
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    Written partition-function file:
      'PF_exomol_NH3.dat'
    for molecule NH3, with isotopes ['4111', '5111'],
    and temperature range 1--1600 K.

.. py:function:: kurucz(pf_file, outfile=None, type_flag='as_exomol')
.. code-block:: pycon

    Extract Kurucz partition-function values from input file.
    If requested, write the partition-function into a file for use
    with Pyrat Bay.

    Parameters
    ----------
    pf_file: String
        Input partition-function from Kurucz webpage.  Currently only H2O
        and TiO are available (probably there's no need for any other support).
        Files can be downloaded from these links:
          http://kurucz.harvard.edu/molecules/h2o/h2opartfn.dat
          http://kurucz.harvard.edu/molecules/tio/tiopart.dat
    outfile: String
        If not None, save output to file.
        If outfile == 'default', save output to file named as
        PF_kurucz_molecule.dat

    Returns
    -------
    pf: 2D float ndarray
        TIPS partition function for input molecule.
    isotopes: 1D string list
        List of isotopes.
    temp: 1D float ndarray
        Partition-function temperature samples (K).

    Examples
    --------
    >>> # First, download kurucz data to current dictory, e.g.:
    >>> # wget http://kurucz.harvard.edu/molecules/h2o/h2opartfn.dat
    >>> # wget http://kurucz.harvard.edu/molecules/tio/tiopart.dat

    >>> import pyratbay.opacity.partitions as pf
    >>> pf_data, isotopes, temp = pf.kurucz('h2opartfn.dat', outfile='default')

    Written partition-function file:
      'PF_kurucz_H2O.dat'
    for molecule H2O, with isotopes ['1H1H16O', '1H1H17O', '1H1H18O', '1H2H16O'],
    and temperature range 10--6000 K.

    >>> pf_data, isotopes, temp = pf.kurucz('tiopart.dat', outfile='default')

    Written partition-function file:
      'PF_kurucz_TiO.dat'
    for molecule TiO, with isotopes ['66', '76', '86', '96', '06'],
    and temperature range 10--6000 K.


pyratbay.opacity.broadening
___________________________


.. py:module:: pyratbay.opacity.broadening

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
        >>> import pyratbay.opacity.broadening as b
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


        Initialize self.  See help(type(self)) for accurate signature.

    .. py:method:: eval(x)
    .. code-block:: pycon

        Compute Lorentz profile over the specified coordinates range.

        Parameters
        ----------
        x: 1D float ndarray
           Input coordinates where to evaluate the profile.

        Returns
        -------
        l: 1D float ndarray
           The line profile at the x locations.

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
        >>> import pyratbay.opacity.broadening as b
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


        Initialize self.  See help(type(self)) for accurate signature.

    .. py:method:: eval(x)
    .. code-block:: pycon

        Compute Gaussian profile over the specified coordinates range.

        Parameters
        ----------
        x: 1D float ndarray
           Input coordinates where to evaluate the profile.

        Returns
        -------
        g: 1D float ndarray
           The line profile at the x locations.

.. py:class:: Voigt(x0=0.0, hwhm_L=1.0, hwhm_G=1.0, scale=1.0)

    .. code-block:: pycon

        1D Voigt profile model.

        Parameters
        ----------
        x0: Float
            Line center location.
        hwhm_L: Float
            Half-width at half maximum of the Lorentz distribution.
        hwhm_G: Float
            Half-width at half maximum of the Gaussian distribution.
        scale: Float
            Scale of the profile (scale=1 returns a profile with integral=1.0).

        Examples
        --------
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> import pyratbay.opacity.broadening as b

        >>> hwhm_G = 1.0
        >>> hwhm_L = 1.0
        >>> voigt = b.Voigt(x0=0.0, hwhm_L=hwhm_L, hwhm_G=hwhm_G)

        >>> plt.figure('A Voigt profile', (6,4))
        >>> plt.clf()
        >>> ax = plt.subplot(1, 1, 1)
        >>> x = np.linspace(-15.0, 15.0, 3001)
        >>> plt.plot(x, voigt(x), lw=2.0, color="orange")
        >>> plt.xlim(np.amin(x), np.amax(x))
        >>> plt.xlabel(r"x", fontsize=12)
        >>> plt.ylabel(r"Voigt profile", fontsize=12)

        >>> # Compare a range of Voigt, Lorentz, and Doppler profiles:
        >>> lorentz = b.Lorentz(x0=0.0, hwhm=hwhm_L)
        >>> doppler = b.Gauss(x0=0.0, hwhm=hwhm_G)
        >>> nplots = 5
        >>> HWHM_L = np.logspace(-2, 2, nplots)
        >>> nwidths = 10.0
        >>> plt.figure(11, (6,6))
        >>> plt.clf()
        >>> plt.subplots_adjust(0.15, 0.1, 0.95, 0.95, wspace=0, hspace=0)
        >>> for i,hwhm_L in enumerate(HWHM_L):
        >>>     ax = plt.subplot(nplots, 1, 1+i)
        >>>     voigt.hwhm_L = lorentz.hwhm = hwhm_L
        >>>     width = 0.5346*hwhm_L + np.sqrt(0.2166*hwhm_L**2+hwhm_G**2)
        >>>     x = np.arange(-nwidths*width, nwidths*width, width/1000.0)
        >>>     plt.plot(x/width, lorentz(x), lw=2.0, color="blue", label="Lorentz")
        >>>     plt.plot(
        >>>         x/width, doppler(x), lw=2.0, color="limegreen", label="Doppler")
        >>>     plt.plot(
        >>>         x/width, voigt(x), lw=2.0, color="orange", label="Voigt",
        >>>         dashes=(4,1))
        >>>     ymin = np.amin([lorentz(x), voigt(x)])
        >>>     ymax = np.amax([lorentz(x), voigt(x), doppler(x)])
        >>>     plt.ylim(ymin, 3*ymax)
        >>>     ax.set_yscale("log")
        >>>     plt.text(
        >>>         0.025, 0.75, rf"$\rm HWHM_L/HWHM_G={hwhm_L/hwhm_G:4g}$",
        >>>         transform=ax.transAxes)
        >>>     plt.xlim(-nwidths, nwidths)
        >>>     plt.xlabel(r"$\rm x/HWHM_V$", fontsize=12)
        >>>     plt.ylabel("Profile")
        >>>     if i != nplots-1:
        >>>         ax.set_xticklabels([])
        >>>     if i == 0:
        >>>         plt.legend(loc="upper right", fontsize=11)


        Initialize self.  See help(type(self)) for accurate signature.

    .. py:method:: eval(x)
    .. code-block:: pycon

        Evaluate the Voigt profile at the specified coordinates range.

        Parameters
        ----------
        x: 1D float ndarray
            Input coordinates where to evaluate the profile.

        Returns
        -------
        v: 1D float ndarray
            The line profile at the x locations.

.. py:function:: doppler_hwhm(temperature, mass, wn)
.. code-block:: pycon

    Get Doppler half-width at half maximum broadening.

    Parameters
    ----------
    temperature: Float scalar or ndarray
        Atmospheric temperature (Kelvin degree).
    mass: Float scalar or ndarray
        Mass of the species (AMU).
    wn: Float scalar or ndarray
        Wavenumber (cm-1).

    Returns
    -------
    dop_hwhm: Float scalar or ndarray
        The Doppler half-width at half maximum broadening (cm-1).

    Note
    ----
    All inputs must have compatible data shapes to be broadcastable.

    Examples
    --------
    >>> import pyratbay.opacity.broadening as b
    >>> # Doppler HWHM at 1000K and 1 micron, for H2O and CO2:
    >>> temperature = 1000.0
    >>> wn = 10000.0
    >>> mass = np.array([18.0, 44.0])
    >>> dop_hw = b.doppler_hwhm(temperature, mass, wn)
    >>> print(f'Doppler broadening:\n H2O        CO2\n{dop_hw}')
    Doppler broadening:
     H2O        CO2
    [0.02669241 0.01707253]

.. py:function:: lorentz_hwhm(temperature, pressure, masses, radii, vol_mix_ratio, imol)
.. code-block:: pycon

    Get Lorentz half-width at half maximum broadening.

    Parameters
    ----------
    temperature: Float scalar or ndarray
        Atmospheric temperature (Kelvin degree).
    pressure: Float scalar or ndarray
        Atmospheric pressure (barye).
    masses: 1D float ndarray
        Masses of atmospheric species (AMU).
    radii: 1D float ndarray
        Collision radius of atmospheric species (cm).
    vol_mix_ratio: 1D float ndarray
        Volume mixing ratio of atmospheric species.
    imol: Integer
        Index of species to calculate the HWHM (in masses/radii arrays).

    Returns
    -------
    lor_hwhm: Float scalar or ndarray
        The Lorentz half-width at half maximum broadening (cm-1).

    Note
    ----
    The temperature, pressure, and imol inputs must have compatible
    shapes to be broadcastable.

    Examples
    --------
    >>> import pyratbay.opacity.broadening as b
    >>> import pyratbay.constants as pc
    >>> # Lorenz HWHM at 1000K and 1 bar, for H2O and CO2:
    >>> temperature = 1000.0
    >>> pressure = 1.0 * pc.bar
    >>> #                  H2O   CO2   H2    He
    >>> masses = np.array([18.0, 44.0, 2.0,  4.0])
    >>> radii  = np.array([1.6,  1.9,  1.45, 1.4]) * pc.A
    >>> vmr    = np.array([1e-4, 1e-4, 0.85, 0.15])
    >>> imol = np.array([0, 1])
    >>> lor_hw = b.lorentz_hwhm(temperature, pressure, masses, radii, vmr, imol)
    >>> print(f'Lorentz broadening:\n H2O        CO2\n{lor_hw}')
    Lorentz broadening:
     H2O        CO2
    [0.03691111 0.04308068]

.. py:function:: min_widths(min_temp, max_temp, min_wn, max_mass, min_rad, min_press)
.. code-block:: pycon

      Estimate the minimum Doppler and Lorentz half-widths at half maximum
      (cm-1) for an H2-dominated atmosphere.

      Parameters
      ----------
      min_temp: Float
          Minimum atmospheric tmperature (Kelvin degrees).
      max_temp: Float
          Maximum atmospheric tmperature (Kelvin degrees).
      min_wn: Float
          Minimum spectral wavenumber (cm-1).
      max_mass: Float
          Maximum mass of molecule/isotope (amu).
      min_rad: Float
          Minimum collisional radius (cm).
      min_press: Float
          Minimum atmospheric pressure (barye).

      Returns
      -------
      dmin: Float
          Minimum Doppler HWHM (cm-1).
      lmin: Float
          Minimum Lorentz HWHM (cm-1).

      Examples
      --------
      >>> import pyratbay.opacity.broadening as b
      >>> import pyratbay.constants as pc
      >>> min_temp =  100.0
      >>> max_temp = 3000.0
      >>> min_wn   = 1.0/(10.0*pc.um)
      >>> max_mass = 18.015    # H2O molecule
      >>> min_rad  = 1.6*pc.A  # H2O molecule
      >>> min_press = 1e-5 * pc.bar
      >>> dmin, lmin = b.min_widths(min_temp, max_temp, min_wn, max_mass,
      >>>     min_rad, min_press)
      >>> print('Minimum Doppler half width: {:.2e} cm-1
    '
      >>>       'Minimum Lorentz half width: {:.2e} cm-1'.format(dmin,lmin))
  

.. py:function:: max_widths(min_temp, max_temp, max_wn, min_mass, max_rad, max_press)
.. code-block:: pycon

      Estimate the maximum Doppler and Lorentz half-widths at half maximum
      (cm-1) for an H2-dominated atmosphere.

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
      max_rad: Float
          Maximum collisional radius (cm).
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
      >>> import pyratbay.opacity.broadening as b
      >>> import pyratbay.constants as pc
      >>> min_temp =  100.0
      >>> max_temp = 3000.0
      >>> max_wn   = 1.0/(1.0*pc.um)
      >>> min_mass = 18.015    # H2O molecule
      >>> max_rad  = 1.6*pc.A  # H2O molecule
      >>> max_press = 100.0*pc.bar
      >>> dmax, lmax = b.max_widths(min_temp, max_temp, max_wn, min_mass,
      >>>     max_rad, max_press)
      >>> print('Maximum Doppler half width: {:.2e} cm-1
    '
      >>>       'Maximum Lorentz half width: {:.2e} cm-1'.format(dmax,lmax))
  


pyratbay.plots
______________


.. py:module:: pyratbay.plots

.. py:function:: alphatize(colors, alpha, bg='w')
.. code-block:: pycon

    Get rgb representation of a color as if it had the specified alpha.

    Parameters
    ----------
    colors: color or iterable of colors
        The color to alphatize.
    alpha: Float
        Alpha value to apply.
    bg: color
        Background color.

    Returns
    -------
    rgb: RGB or list of RGB color arrays
        The RGB representation of the alphatized color (or list of colors).

    Examples
    --------
    >>> import pyrabay.plots as pp
    >>> pp.alphatize('r', 0.5)
    array([1. , 0.5, 0.5])
    >>> pp.alphatize(['r', 'b'], 0.8)
    [array([1. , 0.2, 0.2]), array([0.2, 0.2, 1. ])]

.. py:function:: spectrum(spectrum, wavelength, rt_path, data=None, uncert=None, bandwl=None, bandflux=None, bandtrans=None, bandidx=None, starflux=None, rprs=None, label='model', bounds=None, logxticks=None, gaussbin=2.0, yran=None, filename=None, fignum=501, axis=None)
.. code-block:: pycon

    Plot a transmission or emission model spectrum with (optional) data
    points with error bars and band-integrated model.

    Parameters
    ----------
    spectrum: 1D float ndarray
        Planetary spectrum evaluated at wavelength.
    wavelength: 1D float ndarray
        The wavelength of the model in microns.
    rt_path: String
        Radiative-transfer observing geometry (transit, eclipse, or emission).
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
    label: String
        Label for spectrum curve.
    bounds: Tuple
        Tuple with -2, -1, +1, and, +2 sigma boundaries of spectrum.
        If not None, plot shaded area between +/-1sigma and +/-2sigma
        boundaries.
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
    axis: AxesSubplot instance
        The matplotlib Axes of the figure.

    Returns
    -------
    ax: AxesSubplot instance
        The matplotlib Axes of the figure.

.. py:function:: contribution(contrib_func, wl, rt_path, pressure, radius, rtop=0, filename=None, filters=None, fignum=-21)
.. code-block:: pycon

    Plot the band-integrated normalized contribution functions
    (emission) or transmittance (transmission).

    Parameters
    ----------
    contrib_func: 2D float ndarray
        Band-integrated contribution functions [nfilters, nlayers].
    wl: 1D float ndarray
        Mean wavelength of the bands in microns.
    rt_path: String
        Radiative-transfer observing geometry (emission or transit).
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

    Returns
    -------
    ax: AxesSubplot instance
        The matplotlib Axes of the figure.

    Notes
    -----
    - The dashed lines denote the 0.16 and 0.84 percentiles of the
      cumulative contribution function or the transmittance (i.e.,
      the boundaries of the central 68% of the respective curves).
    - If there are more than 80 filters, this code will thin the
      displayed filter names.

.. py:function:: temperature(pressure, profiles=None, labels=None, colors=None, bounds=None, punits='bar', ax=None, filename=None, theme='blue', alpha=[0.8, 0.6], fs=13, lw=2.0, fignum=504)
.. code-block:: pycon

    Plot temperature profiles.

    Parameters
    ----------
    pressure: 1D float ndarray
        The atmospheric pressure profile in barye.
    profiles: iterable of 1D float ndarrays
        Temperature profiles to plot.
    labels: 1D string iterable
        Labels for temperature profiles.
    colors: 1D string iterable.
        Colors for temperature profiles.
    bounds: Tuple
        Tuple with -1sigma, +1sigma, -2sigma, and +2sigma temperature
        boundaries.
        If not None, plot shaded area between +/-1sigma and +/-2sigma
        boundaries.
    punits: String
        Pressure units for output plot (input units are always barye).
    ax: AxesSubplot instance
        If not None, plot into the given axis.
    filename: String
        If not None, save plot to given file name.
    theme: String
        The histograms' color theme for bounds regions.
        Only 'blue' and 'orange' themes are valid at the moment.
        Alternatively, provide a two-element iterable to provide the colors.
    alpha: 2-element float iterable
        Alpha transparency for bounds regions.
    fs: Float
        Labels font sizes.
    lw: Float
        Lines width.
    fignum: Integer
        Figure's number (ignored if axis is not None).

    Returns
    -------
    ax: AxesSubplot instance
        The matplotlib Axes of the figure.

.. py:function:: abundance(vol_mix_ratios, pressure, species, highlight=None, xlim=None, punits='bar', colors=None, dashes=None, filename=None, lw=2.0, fignum=505, fs=13, legend_fs=None, ax=None)
.. code-block:: pycon

    Plot atmospheric volume-mixing-ratio abundances.

    Parameters
    ----------
    vol_mix_ratios: 2D float ndarray
        Atmospheric volume mixing ratios to plot [nlayers,nspecies].
    pressure: 1D float ndarray
        Atmospheric pressure [nlayers], units are given by punits argument.
    species: 1D string iterable
        Atmospheric species names [nspecies].
    highlight: 1D string iterable
        List of species names to highlight.  Non-highlighed species are
        plotted with alpha=0.4, below the highligted species, and are
        not considered to set the default xlim (e.g., might not be shown
        if their abundances are too low).
        If None, all input species are highlighted.
    xlim: 2-element float iterable
        Volume mixing ratio plotting boundaries.
    punits: String
        Pressure units.
    colors: 1D string iterable
        List of colors to use.
        - If len(colors) >= len(species), colors are assigned to each
          species irrespective of highlight.
        - If len(colors) < len(species), the display will cycle the
          color list using solid, long-dashed, short-dashed, and dotted
          line styles (all highlight species being displayed before the rest).
        - If colors == 'default', use pyratbay.plots.default_colors
          dict to assign colors.
        - If colors is None, use matplotlib's default color cycler.
    dashes: 1D dash-sequence iterable
        List of line-styles for each species, irrespective of highlight.
        len(dashes) has to be equal to len(species).
        Alternatively, dashes can by a dash-sequence Cycler.
    filename: String
        If not None, save plot to given file name.
    lw: Float
        Lines width.
    fignum: Integer
        Figure's number (ignored if axis is not None).
    fs: Float
        Labels font sizes.
    legend_fs: Float
        Legend font size.  If legend_fs is None, default to fs-2.
        If legend_fs <= 0, do not plot a legend.
    ax: AxesSubplot instance
        If not None, plot into the given axis.

    Returns
    -------
    ax: AxesSubplot instance
        The matplotlib Axes of the figure.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.plots as pp

    >>> nlayers = 51
    >>> pressure = pa.pressure('1e-6 bar', '1e2 bar', nlayers)
    >>> temperature = pa.temperature('isothermal', pressure,  params=1000.0)
    >>> species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 TiO VO H2 H He Na K'.split()
    >>> vmr = pa.chemistry('tea', pressure, temperature, species).vmr
    >>> ax = pp.abundance(
    >>>     vmr, pressure, species, colors='default',
    >>>     highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())

.. py:data:: default_colors
.. code-block:: pycon

  {'H2O': 'navy', 'CO2': 'red', 'CO': 'limegreen', 'CH4': 'orange', 'H2': 'deepskyblue', 'He': 'seagreen', 'HCN': '0.7', 'NH3': 'magenta', 'C2H2': 'brown', 'C2H4': 'pink', 'N2': 'gold', 'H': 'olive', 'TiO': 'black', 'VO': 'peru', 'Na': 'darkviolet', 'K': 'cornflowerblue'}


pyratbay.spectrum
_________________


.. py:module:: pyratbay.spectrum

.. py:function:: blackbody_wn(...)
.. code-block:: pycon

    Calculate the Planck emission function in wavenumber space:
       Bnu(T) = 2 h c**2 nu**3 / (exp(hc*nu/kT)-1),             
    with units of erg s-1 sr-1 cm-2 cm.                         
                                                            
    Parameters                                                  
    ----------                                                  
    wn: 1D float ndarray                                        
        Wavenumber spectrum (cm-1).                             
    temp: Float                                                 
        Temperature (Kelvin).                                   
    B: 1D float ndarray [optional]                              
        Array to store the Planck emission.                     
                                                            
    Returns                                                     
    -------                                                     
    (If B was not provided as input:)                           
    B: 1D float ndarray                                         
        Planck emission function at wn (erg s-1 sr-1 cm-2 cm).

.. py:function:: blackbody_wn_2D(...)
.. code-block:: pycon

    Compute the Planck emission function in wavenumber space:       
       Bnu(T) = 2 h c**2 nu**3 / (exp(hc*nu/kT)-1),                  
    with units of erg s-1 sr-1 cm-2 cm.                              
                                                                 
    Parameters                                                       
    ----------                                                       
    wn: 1D float ndarray                                             
        Wavenumber spectrum (cm-1).                                  
    temp: 1D float ndarray                                           
        Temperature at each layer (K).                               
    B: 2D float ndarray [optional]                                   
        Array to store the Planck emission of shape [nlayers, nwave].
    last: 1D integer ndarray [optional]                              
        Indices of last layer to evaluate at each wavenumber.        
                                                                 
    Returns                                                          
    -------                                                          
    (If B was not provided as input:)                                
    B: 2D float ndarray                                              
        Planck emission function at wn (erg s-1 sr-1 cm-2 cm).

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
    >>> import pyratbay.spectrum as ps
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
    >>> import pyratbay.spectrum as ps
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

.. py:class:: PassBand(filter_file, wn=None)

    .. code-block:: pycon

        A Filter passband object.


        Parameters
        ----------
        filter_file: String
            Path to filter file containing wavelength (um) and passband
            response in two columns.
            Comment and blank lines are ignored.

        Examples
        --------
        >>> import pyratbay.spectrum as ps
        >>> import pyratbay.constants as pc
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np

        >>> filter_file = f'{pc.ROOT}pyratbay/data/filters/spitzer_irac2_sa.dat'
        >>> band = ps.PassBand(filter_file)

        >>> # Evaluate over a wavelength array (um):
        >>> wl = np.arange(3.5, 5.5, 0.001)
        >>> out_wl, out_response = band(wl)

        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.plot(out_wl, out_response)
        >>> plt.plot(band.wl, band.response)  # Same variables
        >>> # Note wl differs from band.wl, but original array can be used as:
        >>> plt.plot(wl[band.idx], band.response)

        >>> # Evaluate over a wavenumber array:
        >>> wn = 1e4 / wl
        >>> band(wn=wn)
        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.plot(band.wn, band.response, dashes=(5,3))
        >>> plt.plot(wn[band.idx], out_response)

    .. py:method:: save_filter(save_file)
    .. code-block:: pycon

        Write filter response function data to file, into two columns
        the wavelength (um) and the response function.

        Parameters
        ----------
        save_file: String
            File where to save the filter data.

.. py:class:: Tophat(wl0, half_width, name='tophat')

    .. code-block:: pycon

        A Filter passband object with a tophat-shaped passband.


        Parameters
        ----------
        wl0: Float
            The passband's central wavelength (um units).
        half_width: Float
            The passband's half-width (um units).
        name: Str
            A user-defined name for the filter when calling str(self),
            e.g., to identify the instrument provenance of this filter.

        Examples
        --------
        >>> import pyratbay.spectrum as ps
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np

        >>> hat = ps.Tophat(4.5, 0.5)

        >>> # Evaluate over a wavelength array (um units):
        >>> wl = np.arange(3.5, 5.5, 0.001)
        >>> out_wl, out_response = hat(wl)

        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.plot(out_wl, out_response)
        >>> plt.plot(hat.wl, hat.response)  # Same variables
        >>> # Note wl differs from hat.wl, but original array can be used as:
        >>> plt.plot(wl[hat.idx], hat.response)

        >>> # Evaluate over a wavenumber array:
        >>> wn = 1e4 / wl
        >>> hat(wn=wn)
        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.plot(hat.wn, hat.response, dashes=(5,3))
        >>> plt.plot(wn[hat.idx], out_response)

    .. py:method:: save_filter(save_file)
    .. code-block:: pycon

        Write filter response function data to file, into two columns
        the wavelength (um) and the response function.

        Parameters
        ----------
        save_file: String
            File where to save the filter data.

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
    >>> import pyratbay.spectrum as ps
    >>> wl0     = 1.50
    >>> width   = 0.50
    >>> margin  = 0.10
    >>> dlambda = 0.05
    >>> wl, trans = ps.tophat(wl0, width, margin, dlambda)
    >>> print(wl, trans, sep='\n')
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
    >>> import pyratbay.spectrum as ps
    >>> import numpy as np
    >>> wn     = np.linspace(1.3, 1.7, 11)
    >>> signal = np.array(np.abs(wn-1.5)<0.1, np.double)
    >>> specwn = np.linspace(1, 2, 51)
    >>> resampled, wnidx = ps.resample(signal, wn, specwn)
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

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pyratbay.io as io
    >>> import pyratbay.spectrum as ps
    >>> import pyratbay.constants as pc
    >>> # Load Spitzer IRAC filters:
    >>> wn1, irac1 = io.read_spectrum(
    >>>     pc.ROOT+'pyratbay/data/filters/spitzer_irac1_sa.dat')
    >>> wn2, irac2 = io.read_spectrum(
    >>>     pc.ROOT+'pyratbay/data/filters/spitzer_irac2_sa.dat')
    >>> # Spectrum to integrate:
    >>> wn = np.arange(1500, 5000.1, 1.0)
    >>> sflux = ps.bbflux(wn, 1800.0)
    >>> # Integrate over single filter:
    >>> bandflux = ps.band_integrate(sflux, wn, irac1, wn1)
    >>> # Integrate over multiple:
    >>> bandfluxes = ps.band_integrate(sflux, wn, [irac1,irac2], [wn1, wn2])
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

.. py:function:: contribution_function(optdepth, pressure, B)
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
        The contribution function at each layer and wavenumber
        of shape [nlayers, nwave].

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

.. py:function:: band_cf(cf, bandtrans, wn, bandidx)
.. code-block:: pycon

    Compute band-averaged contribution functions or transmittances.

    Parameters
    ----------
    cf: 2D float ndarray
        The contribution function or transmittance of
        shape [nlayers, nwave].
    bandtrans: List of 1D ndarrays
        List of band transmission curves.
    wn: 1D float ndarray
        The wavenumber sampling (in cm-1).
    bandidx: List of 1D ndarrays
        List of wavenumber-index arrays for each band transmission curve.

    Returns
    -------
    bandcf: 2D float ndarray
        The band-integrated contribution functions of
        shape [nlayers, nbands].

.. py:function:: convective_flux(pressure, temperature, cp, gravity, mu, rho, alpha=1.0, beta=1.0)
.. code-block:: pycon

    Estimate the convective flux for an atmosphere following mixing-
    length theory as described in 'Modern Astrophysics (Carrol & Ostlie).

    Parameters
    ----------
    pressure: 1D float array
        Atmospheric pressure profile (barye).
    temperature: 1D float array
        Atmospheric temperature profile (kelvin degree).
    cp: 1D float array
        Atmospheric specific heat capacity at constant pressure
        (erg K-1 mol-1).
    gravity: 1D float array
        Atmospheric gravity profile (cm s-2).
    mu: 1D float array
        Atmospheric mean molecular mass profile (g mol-1).
    rho: 1D float array
        Atmospheric mass-density profile (g cm-3).
    alpha: Float
        Mixing-length scaling parameter: alpha = l/H, with l the
        mixing length and H the pressure scale height.
    beta: Float
        Free parameter from the average kinetic energy velocity
        estimation. Should take a value between 0 and 1.

    Returns
    -------
    F_conv: 1D float array
        Estimated convective flux (erg s-1 cm-2).
        This is not zero only where the actual temperature gradient
        is larger than the adiabatic gradient.

.. py:function:: radiative_equilibrium(pressure, temperature, nsamples, chem_model, two_stream_rt, wavenumber, spec, atm, radeq_temps=None, convection=False, tmin=0.0, tmax=6000.0, f_scale=None, mplanet=None, mol_mass=None)
.. code-block:: pycon

    Compute radiative-thermochemical equilibrium atmosphere.
    Currently there is no convergence criteria implemented,
    some 100--300 iterations are typically sufficient to converge
    to a stable temperature-profile solution.

    Parameters
    ----------
    nsamples: Integer
        Number of radiative-equilibrium iterations to run.
    continue_run: Bool
        If True, continue from a previous radiative-equil. run.
    no_convection: Bool
        If True, skip convective flux calculation in the radiative
        equilibrium calculation.

    Returns
    -------
    There are no returned values, but this method updates the
    temperature profile (self.atm.temp) and abundances (self.atm.q)
    with the values from the last radiative-equilibrium iteration.

    This method also defines pyrat.atm.radeq_temps, a 2D array
    containing all temperature-profile iterations.


pyratbay.atmosphere
___________________


.. py:module:: pyratbay.atmosphere

.. py:function:: pressure(ptop, pbottom, nlayers, units='bar', log=None, verb=0)
.. code-block:: pycon

    Compute a log-scale pressure profile.

    Parameters
    ----------
    ptop: String or Float
        Pressure at the top of the atmosphere. If string, may contain units.
    pbottom: String or Float
       Pressure at the bottom of the atmosphere. If string, may contain units.
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
    >>> import pyratbay.constants as pc

    >>> nlayers = 9
    >>> # These are all equivalent:
    >>> p1 = pa.pressure(ptop=1e-6, pbottom=1e2, nlayers=nlayers)
    >>> p2 = pa.pressure(1e-6, 1e2, nlayers, 'bar')
    >>> p3 = pa.pressure('1e-6 bar', '1e2 bar', nlayers)
    >>> p4 = pa.pressure(1e-6*pc.bar, 1e2*pc.bar, nlayers, 'barye')
    >>> print(p1/pc.bar)
    [1.e-06 1.e-05 1.e-04 1.e-03 1.e-02 1.e-01 1.e+00 1.e+01 1.e+02]

.. py:function:: temperature(tmodel, pressure=None, nlayers=None, log=None, params=None)
.. code-block:: pycon

    Temperature profile wrapper.

    Parameters
    ----------
    tmodel: String
        Name of the temperature model.
    pressure: 1D float ndarray
        Atmospheric pressure profile in barye units.
    nlayers: Integer
        Number of pressure layers.
    log: Log object
        Screen-output log handler.
    params: 1D float ndarray
        Temperature model parameters. If None, return a tuple with the
        temperature model, its arguments, and the number or required parameters.

    Returns
    -------
    If params is not None:
        temperature: 1D float ndarray
            The evaluated atmospheric temperature profile.
    If params is None:
        temp_model: Callable
            The atmospheric temperature model.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa

    >>> nlayers = 11
    >>> # Isothermal profile:
    >>> temp_iso = pa.temperature("isothermal", params=1500.0, nlayers=nlayers)
    >>> print(temp_iso)
    [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.]

    >>> # Three-channel Eddington-approximation profile:
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, "bar")
    >>> params = np.array([-4.84, -0.8, -0.8, 0.5, 1200.0, 100.0])
    >>> temp = pa.temperature('guillot', pressure, params=params)
    >>> print(temp)
    [1046.89057381 1046.89090056 1046.89433798 1046.93040895 1047.30779086
     1051.21739055 1088.76131307 1312.57904127 1640.18896334 1659.78818839
     1665.09706555]

.. py:function:: chemistry(chem_model, pressure, temperature, species, metallicity=0.0, e_abundances={}, e_scale={}, e_ratio={}, q_uniform=None, solar_file=None, log=None, verb=1, atmfile=None, punits='bar')
.. code-block:: pycon

    Compute atmospheric abundaces for given pressure and
    temperature profiles with either uniform abundances or TEA.

    Parameters
    ----------
    chem_model: String
        Name of chemistry model.
    pressure: 1D float ndarray
        Atmospheric pressure profile (barye).
    temperature: 1D float ndarray
        Atmospheric temperature profile (Kelvin).
    species: 1D string list
        Output atmospheric composition.
    metallicity: Float
        Metallicity enhancement factor in dex units relative to solar.
    e_abundances: Dictionary
        Custom elemental abundances.
        The dict contains the name of the element and their custom
        abundance in dex units relative to H=12.0.
        These values override metallicity.
    e_scale: Dictionary
        Scaling abundance factor for specified atoms by the respective
        values (in dex units, in addition to metallicity scaling).
        E.g. (3x solar): e_scale = {'C': np.log10(3.0)}
    e_ratio: Dictionary
        Custom elemental abundances scaled relative to another element.
        The dict contains the pair of elements joined by an underscore
        and their ratio in dex units, e.g., for a C/O ratio of 0.8 set
        e_ratio = {'C_O': 0.8}.
        These values modify the abundances after metallicity and
        e_scale have been applied.
    solar_file: String
        Input solar elemental abundances file (default Asplund et al. 2021).
    log: Log object
        Screen-output log handler.
    verb: Integer
        Verbosity level.
    atmfile: String
        If not None, output file where to save the atmospheric model.
    punits: String
        Output pressure units.

    Returns
    -------
    model: Callable
        The atmospheric chemistry-network model.

    Example
    -------
    >>> import pyratbay.atmosphere as pa
    >>> import matplotlib.pyplot as plt
    >>> import pyratbay.plots as pp

    >>> nlayers = 100
    >>> T0 = 1500.0
    >>> pressure = pa.pressure('1e-8 bar', '1e3 bar', nlayers)
    >>> temperature = pa.temperature('isothermal', pressure, nlayers, params=T0)
    >>> species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 H2 H He H+ e-'.split()
    >>> # Equilibrium abundances model:
    >>> chem_model = 'tea'
    >>> chem_network = pa.chemistry(
    >>>     chem_model, pressure, temperature, species,)

    >>> q_uniform = np.array([
    >>>     5e-4, 3e-5, 2e-4, 1e-8,  1e-6, 1e-14, 1e-13, 5e-10,
    >>>     1e-4, 0.85, 5e-3, 0.14,  3e-23, 1e-23])
    >>> chem_model = 'uniform'
    >>> chem_network_unif = pa.chemistry(
    >>>     chem_model, pressure, temperature, species, q_uniform=q_uniform)

    >>> # Plot the results:
    >>> ax1 = pp.abundance(
    >>>     chem_network.vmr, pressure, chem_network.species,
    >>>     colors='default', xlim=[1e-30, 3.0])

    >>> ax2 = pp.abundance(
    >>>     chem_network_unif.vmr, pressure, chem_network_unif.species,
    >>>     colors='default', xlim=[1e-30, 3.0], fignum=506)

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
    >>> nlayers = 11
    >>> punits = 'bar'
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, punits)
    >>> tmodel = pa.tmodels.Isothermal(nlayers)
    >>> species    = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    >>> abundances = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    >>> qprofiles = pa.uniform(pressure, tmodel(1500.0), species,
    >>>     abundances=abundances, punits=punits)
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

.. py:function:: abundance(pressure, temperature, species, elements=None, quniform=None, atmfile=None, punits='bar', metallicity=0.0, e_abundances={}, e_scale={}, e_ratio={}, solar_file=None, log=None, verb=1, ncpu=None, xsolar=None, escale=None)
.. code-block:: pycon

    Compute atmospheric abundaces for given pressure and
    temperature profiles with either uniform abundances or TEA.

    Parameters
    ----------
    pressure: 1D float ndarray
        Atmospheric pressure profile (barye).
    temperature: 1D float ndarray
        Atmospheric temperature profile (Kelvin).
    species: 1D string list
        Output atmospheric composition.
    elements: 1D strings list
        Input elemental composition (default to minimum list of elements
        required to form species).
    quniform: 1D float ndarray
        If not None, the output species abundances (isobaric).
    atmfile: String
        If not None, output file where to save the atmospheric model.
    punits: String
        Output pressure units.
    metallicity: Float
        Metallicity enhancement factor in dex units relative to solar.
    e_abundances: Dictionary
        Custom elemental abundances.
        The dict contains the name of the element and their custom
        abundance in dex units relative to H=12.0.
        These values override metallicity.
    e_scale: Dict
        Scaling abundance factor for specified atoms by the respective
        values (in dex units, in addition to metallicity scaling).
        E.g. (3x solar): e_scale = {'C': np.log10(3.0)}
    solar_file: String
        Input solar elemental abundances file (default Asplund et al. 2021).
    log: Log object
        Screen-output log handler.
    verb: Integer
        Verbosity level.

    ncpu: Integer [DEPRECATED]
        Number of parallel CPUs to use in TEA calculation.
    xsolar: Float [DEPRECATED]
        Metallicity enhancement factor.  Deprecated, use metallicity instead.
    escale: Dict [DEPRECATED]
        Multiplication factor for specified atoms (dict's keys)
        by the respective values (on top of the xsolar scaling).
        Deprecated, use e_scale instead.

    Returns
    -------
    vmr: 2D float ndarray
       Atmospheric volume mixing fraction abundances of shape
       [nlayers, nspecies].

    Example
    -------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants as pc

    >>> nlayers = 100
    >>> press = np.logspace(-8, 3, nlayers) * pc.bar
    >>> temp  = np.tile(900.0, nlayers)
    >>> species = 'H2O CH4 CO CO2 NH3 C2H2 C2H4 HCN N2 H2 H He'.split()
    >>> # Thermochemical equilibrium abundances for requested species:
    >>> vmr = pa.abundance(press, temp, species)

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
    >>> temperature = pa.tmodels.Isothermal(nlayers)(1500.0)
    >>> mu = np.tile(2.3, nlayers)
    >>> g = pc.G * pc.mjup / pc.rjup**2
    >>> r0 = 1.0 * pc.rjup
    >>> p0 = 1.0 * pc.bar
    >>> # Radius profile in Jupiter radii:
    >>> radius = pa.hydro_g(pressure, temperature, mu, g, p0, r0) / pc.rjup
    >>> print(radius)
    [1.0563673  1.04932138 1.04227547 1.03522956 1.02818365 1.02113774
     1.01409182 1.00704591 1.         0.99295409 0.98590818]

.. py:function:: hydro_m(pressure, temperature, mu, mass, p0, r0)
.. code-block:: pycon

    Calculate radii using the hydrostatic-equilibrium equation considering
    a variable gravity: g(r) = G*mass/r**2

    Parameters
    ----------
    pressure: 1D float ndarray
        Atmospheric pressure for each layer (in barye).
    temperature: 1D float ndarray
        Atmospheric temperature for each layer (in K).
    mu: 1D float ndarray
        Mean molecular mass for each layer (in g mol-1).
    mass: Float
        Object's mass (in g).
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
    It is possible that this hydrostatic solution diverges when an
    atmosphere is too puffy.  In such cases, some returned radii
    will have np.inf values (at the top layers).

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants  as pc

    >>> nlayers = 11
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    >>> temperature = pa.tmodels.Isothermal(nlayers)(1500.0)
    >>> mu = np.tile(2.3, nlayers)
    >>> mplanet = 1.0 * pc.mjup
    >>> r0 = 1.0 * pc.rjup
    >>> p0 = 1.0 * pc.bar
    >>> # Radius profile in Jupiter radii:
    >>> radius = pa.hydro_m(pressure, temperature, mu, mplanet, p0, r0)/pc.rjup
    >>> print(radius)
    [1.05973436 1.05188019 1.04414158 1.036516   1.029001   1.02159419
     1.01429324 1.00709591 1.         0.99300339 0.986104  ]

.. py:function:: rhill(smaxis, mplanet, mstar)
.. code-block:: pycon

    Compute the Hill radius.  If any argument is None, return inf.

    Parameters
    ----------
    smaxis: Float
        Orbital semi-major axis (in cm).
    mplanet: Float
        Planetary mass (in g).
    mstar: Float
        Stellar mass (in g).

    Returns
    -------
    rhill: Float
        Hill radius of planet.

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
        >>> print(f'{elements}
    {stoichs}')
        ['C', 'H', 'He', 'O']
        [[0 2 0 0]
         [0 0 1 0]
         [0 2 0 1]
         [1 0 0 1]
         [1 0 0 2]
         [1 4 0 0]]
    

.. py:function:: mean_weight(abundances, species=None, molfile=None, mass=None)
.. code-block:: pycon

    Calculate the mean molecular weight (a.k.a. mean molecular mass)
    for the given abundances composition.

    Parameters
    ----------
    abundances: 2D float iterable
        Species volume mixing fraction, of shape [nlayers,nmol].
    species: 1D string iterable
        Species names.
    molfile: String
        A molecules file with the species info.  If None, use
        pyratbay's default molecules.dat file.
    mass: 1D float ndarray
        The mass for each one of the species (g mol-1).
        If not None, this variable takes precedence over molfile.

    Returns
    -------
    mu: 1D float ndarray
        Mean molecular weight at each layer for the input abundances.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> species = 'H2 He H2O CO CO2 CH4'.split()
    >>> abundances = [[0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]]
    >>> mu = pa.mean_weight(abundances, species)
    >>> print(mu)
    [2.31939114]

.. py:function:: ideal_gas_density(abundances, pressure, temperature)
.. code-block:: pycon

    Use the Ideal gas law to calculate number density in molecules cm-3.

    Parameters
    ----------
    abundances: 2D float ndarray
        Species volume mixing fraction, of shape [nlayers,nmol].
    pressure: 1D ndarray
        Atmospheric pressure profile (in barye units).
    temperature: 1D ndarray
        Atmospheric temperature profile (in kelvin).

    Returns
    -------
    density: 2D float ndarray
        Atmospheric density in molecules cm-3.

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
    >>> dens = pa.ideal_gas_density(qprofiles, pressure, temperature)
    >>> print(dens[0])
    [4.10241993e+10 7.24297303e+09 4.82864869e+06 4.82864869e+06
     4.82864869e+02 4.82864869e+06]

.. py:function:: equilibrium_temp(tstar, rstar, smaxis, A=0.0, f=1.0, tstar_unc=0.0, rstar_unc=0.0, smaxis_unc=0.0)
.. code-block:: pycon

    Calculate equilibrium temperature and uncertainty.

    Parameters
    ----------
    tstar: Scalar
        Effective temperature of host star (in kelvin degrees).
    rstar: Scalar
        Radius of host star (in cm).
    smaxis: Scalar
        Orbital semi-major axis (in cm).
    A: Scalar
        Planetary bond albedo.
    f: Scalar
        Planetary energy redistribution factor:
        f=0.5  no redistribution (total dayside reemission)
        f=1.0  good redistribution (4pi reemission)
    tstar_unc: Scalar
        Effective temperature uncertainty (in kelvin degrees).
    rstar_unc: Scalar
        Stellar radius uncertainty (in cm).
    smaxis_unc: Scalar
        Semi-major axis uncertainty (in cm).

    Returns
    -------
    teq: 1D ndarray
        Planet equilibrium temperature in kelvin degrees.
    teq_unc: 1D ndarray
        Equilibrium temperature uncertainty.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants as pc
    >>> import numpy as np
    >>> # HD 209458b (Stassun et al. 2017):
    >>> tstar  = 6091.0
    >>> rstar  = 1.19 * pc.rsun
    >>> smaxis = 0.04747 * pc.au
    >>> tstar_unc  = 10.0
    >>> rstar_unc  = 0.02 * pc.rsun
    >>> smaxis_unc = 0.00046 * pc.au
    >>> A = 0.3
    >>> f = np.array([1.0, 0.5])
    >>> teq, teq_unc = pa.equilibrium_temp(tstar, rstar, smaxis, A, f,
    >>>     tstar_unc, rstar_unc, smaxis_unc)
    >>> print(f'HD 209458b T_eq =\n    '
    >>>    f'{teq[0]:.1f} +/- {teq_unc[0]:.1f} K (day--night redistributed)\n'
    >>>    f'    {teq[1]:.1f} +/- {teq_unc[1]:.1f} K (instant re-emission)')
    HD 209458b T_eq =
        1345.1 +/- 13.2 K (day--night redistribution)
        1599.6 +/- 15.7 K (instant re-emission)

.. py:function:: transit_path(radius, nskip=0)
.. code-block:: pycon

    Calculate the distances between layers for a set of rays grazing
    an atmosphere defined by concentric spheres at radii radius.
    Assume that the grazing rays have an impact parameter = radius.
    See for example, Fig 1 of Molliere et al. (2019), AA, 627, 67.

    Parameters
    ----------
    radius: 1D float ndarray
        Atmospheric radius profile array (from top to bottom).
    nskip: Integer
        Number of layers to skip from the top of the atmosphere.

    Returns
    -------
    path: List of 1D float ndarrays
        List where each element is a 1D array of the paths at an
        impact parameter defined.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> nlayers = 5
    >>> radius = np.linspace(5.0, 1.0, nlayers)
    >>> path = pa.transit_path(radius)
    >>> # First path grazes the outer layer (thus, empty path)
    >>> # Second path is sqrt(5.0**2 - 4.0**2), and so on
    >>> for p in path:
    >>>     print(p)
    []
    [3.]
    [1.35424869 2.64575131]
    [1.11847408 1.22803364 2.23606798]
    [1.02599614 1.04455622 1.09637632 1.73205081]
    >>> # Now, ignore the top layer:
    >>> path = pa.transit_path(radius, nskip=1)
    >>> for p in path:
    >>>     print(p)
    []
    []
    [2.64575131]
    [1.22803364 2.23606798]
    [1.04455622 1.09637632 1.73205081]

.. py:function:: temperature_posterior(posterior, tmodel, tpars, ifree, pressure)
.. code-block:: pycon

    Compute the median and inter-quantiles regions (68% and 95%)
    of a temperature profile posterior.

    Parameters
    ----------
    posterior: 2D float ndarray [nsamples, nfree]
        A posterior distribution for tmodel.
    tmodel: Callable
        Temperature-profile model.
    tpars: 1D float iterable [npars]
        Temperature-profile parameters (including fixed parameters).
    ifree: 1D bool iterable [npars]
        Mask of free (True) and fixed (False) parameters in tpars.
        The number of free parameters must match nfree in posterior.
    pressure: 1D float ndarray [nlayers]
        The atmospheric pressure profile in barye.

    Returns
    -------
    median: 1D float ndarray [nlayers]
        The matplotlib Axes of the figure.
    low1: 1D float ndarray [nlayers]
        Lower temperature boundary of the 68%-interquantile range.
    high1: 1D float ndarray [nlayers]
        Upper temperature boundary of the 68%-interquantile range.
    low2: 1D float ndarray [nlayers]
        Lower temperature boundary of the 95%-interquantile range.
    high2: 1D float ndarray [nlayers]
        Upper temperature boundary of the 95%-interquantile range.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants as pc
    >>> import pyratbay.plots as pp
    >>> import numpy as np

    >>> # Non-inverted temperature profile:
    >>> pressure = pa.pressure('1e-6 bar', '1e2 bar', nlayers=100)
    >>> tmodel = pa.tmodels.Guillot(pressure)
    >>> tpars = np.array([-4.0, -1.0, 0.0, 0.0, 1000.0, 0.0])
    >>> # Simulate posterior where kappa' and gamma are variable:
    >>> nsamples = 5000
    >>> ifree = [True, True, False, False, False, False]
    >>> posterior = np.array([
    >>>     np.random.normal(tpars[0], 0.5, nsamples),
    >>>     np.random.normal(tpars[1], 0.1, nsamples)]).T
    >>> tpost = pa.temperature_posterior(
    >>>     posterior, tmodel, tpars, ifree, pressure)
    >>> ax = pp.temperature(pressure, profiles=tpost[0], bounds=tpost[1:])

.. py:function:: qcapcheck(vmr, qcap, ibulk)
.. code-block:: pycon

    Check if the cummulative abundance of traces exceed qcap.

    Parameters
    ----------
    vmr: 2D float ndarray
        Volume mixing ratio of atmospheric species [nlayers, nspecies].
    qcap: Float
        Cap threshold for cummulative trace abundances.
    ibulk: 1D integer ndarray
        Indices of the bulk species to calculate the mixing ratio.

    Returns
    -------
    vmr_cap_flag: Bool
        Flag indicating whether trace abundances sum more than qcap.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> # Make an atmosphere:
    >>> pressure = pa.pressure(ptop=1e-8, pbottom=1e2, nlayers=11, units='bar')
    >>> temperature = np.tile(1500.0, 11)
    >>> species = ["H2", "He", "H2O"]
    >>> abundances = [0.8495, 0.15, 5e-4]
    >>> qprofiles = pa.uniform(pressure, temperature, species, abundances)
    >>> ibulk = [0,1]
    >>> # Sum of all metals (H2O) does not exceed qcap:
    >>> qcap = 1e-3
    >>> print(pa.qcapcheck(qprofiles, qcap, ibulk))
    False
    >>> # Sum of all metals (H2O) exceedes qcap:
    >>> qcap = 1e-4
    >>> print(pa.qcapcheck(qprofiles, qcap, ibulk))
    True

.. py:function:: balance(vmr, ibulk, ratio, invsrat)
.. code-block:: pycon

    Balance the volume mixing ratios of bulk species, vmr[ibulk],
    such that sum(vmr) = 1.0 at each level.

    Parameters
    ----------
    vmr: 2D float ndarray
        Volume mixing ratio of atmospheric  species [nlayers, nspecies].
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
        vmr_bulk = sum vmr_j = 1.0 - sum vmr_trace.
    This code assumes that the abundance ratio among bulk species
    remains constant in each layer:
        {\rm ratio}_j = vmr_j/vmr_0.
    The balanced abundance of the bulk species is then:
        vmr_j = \frac{{\rm ratio}_j * vmr_{\rm bulk}} {\sum {\rm ratio}}.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>>
    >>> vmr = np.tile([0.8, 0.2, 0.1], (5,1))
    >>> ibulk = [0, 1]
    >>> bratio, invsrat = pa.ratio(vmr, ibulk)
    >>> pa.balance(vmr, ibulk, bratio, invsrat)
    >>> # Balanced VMRs:
    >>> print(vmr[0])
    [0.72 0.18 0.1 ]
    >>> # Sum of VMRs equals one at each layer:
    >>> print(np.sum(vmr, axis=1))
    [1. 1. 1. 1. 1.]
    >>> # Ratio of 'bulk' species remains constant:
    >>> print(vmr[:,1]/vmr[:,0])
    [0.25 0.25 0.25 0.25 0.25]

.. py:function:: ratio(vmr, ibulk)
.. code-block:: pycon

    Calculate the abundance ratios of the species indexed by ibulk, relative
    to the first species in the list.

    Parameters
    ----------
    vmr: 2D float ndarray
        Volume mixing ratio of atmospheric species [nlayers, nspecies].
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
    >>> vmr = np.tile([0.8, 0.2], (5,1))
    >>> ibulk = [0, 1]
    >>> bratio, invsrat = pa.ratio(vmr, ibulk)
    >>> print(bratio)
    [[ 1.    0.25]
     [ 1.    0.25]
     [ 1.    0.25]
     [ 1.    0.25]
     [ 1.    0.25]]
    >>> print(invsrat)
    [ 0.8  0.8  0.8  0.8  0.8]

.. py:function:: qscale(vmr, species, molmodel, molfree, molpars, bulk, qsat=None, iscale=None, ibulk=None, bratio=None, invsrat=None)
.. code-block:: pycon

    Scale specified species abundances and balance bulk abundances to
    conserve sum(vmr)=1 in each layer.

    Parameters
    ----------
    vmr: 2D float ndarray
        Volume mixing ratio of atmospheric species [nlayers, nspecies].
    species: 1D string ndarray
        Names of the species in the atmosphere.
    molmodel: 1D string ndarray
        Model to vary the species abundances.
    molfree: 1D string ndarray
        Names of the species to vary their abundances.
    molpars: 1D float ndarray
        Scaling factor (dex) for each species in molfree.
    bulk: 1D string ndarray
        Names of the bulk (dominant) species.
    qsat: Float
        Maximum allowed combined abundance for trace species.
    iscale: 1D integer ndarray
        Indices of molfree species in vmr.
    ibulk: 1D integer ndarray
        Indices of bulk species in vmr.
    bratio: 2D float ndarray
        Abundance ratios between the bulk species (relative to bulk[0]).
    invsrat: 1D float ndarray
        Inverse of the sum of the ratios (at each layer).

    Returns
    -------
    scaled_vmr: 2D float ndarray
       The modified atmospheric VMR profiles.

    Notes
    -----
    iscale, ibulk, bratio, and invsrat are optional parameters to
    speed up the routine.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa

    >>> nlayers = 2
    >>> vmr = np.tile([0.85, 0.15, 1e-4, 1e-4], (nlayers,1))
    >>> species = ['H2', 'He' ,'H2O', 'CH4']
    >>> # Set the H2O abundance to 1.0e-3:
    >>> molmodel = ['vert']
    >>> molfree = ['H2O']
    >>> molpars = [-3.0]
    >>> bulk = ['H2', 'He']
    >>> scaled_vmr = pa.qscale(
    >>>     vmr, species, molmodel, molfree, molpars, bulk,
    >>> )
    >>> print(scaled_vmr)
    [[8.49065e-01 1.49835e-01 1.00000e-03 1.00000e-04]
     [8.49065e-01 1.49835e-01 1.00000e-03 1.00000e-04]]


pyratbay.atmosphere.tmodels
___________________________


.. py:module:: pyratbay.atmosphere.tmodels

.. py:class:: Isothermal(pressure)

    .. code-block:: pycon

        Isothermal temperature profile model.


        Parameters
        ----------
        pressure: 1D float iterable
            Pressure array where to evaluate the temperature profile.

.. py:class:: Guillot(pressure, gravity=None)

    .. code-block:: pycon

        Guillot (2010) temperature profile based on the three-channel
        Eddington approximation, as described Line et al. (2013)


        Parameters
        ----------
        pressure: 1D float ndarray
            Atmospheric pressure profile (barye).
        gravity: 1D float ndarray or scalar
            Atmospheric gravity profile (cm s-2).
            If None, assume a constant gravity of 1 cm s-2, in which
            case, one should regard the kappa parameter as
            kappa' = kappa/gravity.

        Note that the input gravity can be a scalar value used at all
        atmospheric pressures (as it has been used so far in the
        literature.  However, from a parametric point of view, this is
        redundant, as it only acts as a scaling factor for kappa.
        Ideally, one would wish to input a pressure-dependent gravity,
        but such profile would need to be derived from a hydrostatic
        equilibrium calculation, for example.  Unfortunately, HE cannot
        be solved without knowing the temperature, thus making this a
        circular problem (shrug emoji).

.. py:class:: Guillot(pressure, gravity=None)

    .. code-block:: pycon

        Guillot (2010) temperature profile based on the three-channel
        Eddington approximation, as described Line et al. (2013)


        Parameters
        ----------
        pressure: 1D float ndarray
            Atmospheric pressure profile (barye).
        gravity: 1D float ndarray or scalar
            Atmospheric gravity profile (cm s-2).
            If None, assume a constant gravity of 1 cm s-2, in which
            case, one should regard the kappa parameter as
            kappa' = kappa/gravity.

        Note that the input gravity can be a scalar value used at all
        atmospheric pressures (as it has been used so far in the
        literature.  However, from a parametric point of view, this is
        redundant, as it only acts as a scaling factor for kappa.
        Ideally, one would wish to input a pressure-dependent gravity,
        but such profile would need to be derived from a hydrostatic
        equilibrium calculation, for example.  Unfortunately, HE cannot
        be solved without knowing the temperature, thus making this a
        circular problem (shrug emoji).

.. py:class:: Madhu(pressure)

    .. code-block:: pycon

        Temperature profile model by Madhusudhan & Seager (2009)


        Parameters
        ----------
        pressure: 1D float ndarray
            Pressure array in barye.

.. py:function:: get_model(name, *args, **kwargs)
.. code-block:: pycon

    Get a temperature-profile model by its name.


pyratbay.atmosphere.clouds
__________________________


.. py:module:: pyratbay.atmosphere.clouds

.. py:class:: CCSgray()

    .. code-block:: pycon

        Constant cross-section gray cloud model.


        Initialize self.  See help(type(self)) for accurate signature.

    .. py:method:: extinction(wn, pressure)
    .. code-block:: pycon

        Calculate a uniform gray-cloud cross section in cm2 molec-1:
           cross section = s0 * 10**pars[0],
        between layers with pressure 10**pars[1] -- 10**pars[2] bar
        (top and bottom layers, respectively).
        s0 is the H2 Rayleigh cross section at 0.35 um.

        Parameters
        ----------
        wn:  1D float ndarray
           Wavenumber array in cm-1.

.. py:class:: Deck()

    .. code-block:: pycon

        Instantly opaque gray cloud deck at given pressure.


        Initialize self.  See help(type(self)) for accurate signature.

    .. py:method:: extinction(pressure, radius, temp)
    .. code-block:: pycon

        Calculate gray-cloud deck model that's optically thin
        above ptop, and becomes instantly opaque at ptop, with
        ptop (bar) = 10**pars[0].

        Parameters
        ----------
        pressure: 1D float ndarray
            Atmospheric pressure profile (in barye).
        radius: 1D float ndarray
            Atmospheric radius profile (in cm).
        temp: 1D float ndarray
            Atmospheric temperature profile (in Kelvin degree).

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


        Parameters
        ----------
        mol: String
           The species, which can be H, He, or H2.

.. py:class:: Lecavelier()

    .. code-block:: pycon

        Rayleigh-scattering model from Lecavelier des Etangs et al. (2008).
        AA, 485, 865.


        Initialize self.  See help(type(self)) for accurate signature.

    .. py:method:: extinction(wn)
    .. code-block:: pycon

        Calculate the Rayleigh cross section in cm2 molec-1:
            cross section = f_ray * s0 * (lambda/l0)**alpha_ray,
        parameterized as params = [log10(f_ray), alpha_ray).

        Parameters
        ----------
        wn:  1D float ndarray
            Wavenumber array in cm-1.

.. py:function:: get_model(name)
.. code-block:: pycon

    Get a Rayleigh model by its name.


pyratbay.atmosphere.alkali
__________________________


.. py:module:: pyratbay.atmosphere.alkali

.. py:class:: SodiumVdW(cutoff=4500.0)

    .. code-block:: pycon

        Sodium Van der Waals model (Burrows et al. 2000, ApJ, 531).


        Initialize self.  See help(type(self)) for accurate signature.

    .. py:method:: absorption(press, temp, wn)
    .. code-block:: pycon

        Evaluate alkali model's opacity cross section (cm2 molecule-1).

.. py:class:: PotassiumVdW(cutoff=4500.0)

    .. code-block:: pycon

        Potassium Van der Waals model (Burrows et al. 2000, ApJ, 531).


        Initialize self.  See help(type(self)) for accurate signature.

    .. py:method:: absorption(press, temp, wn)
    .. code-block:: pycon

        Evaluate alkali model's opacity cross section (cm2 molecule-1).

.. py:function:: get_model(name, *args)
.. code-block:: pycon

    Get an alkali model by its name.


pyratbay.pyrat.read_atm
_______________________


.. py:module:: pyratbay.pyrat.read_atm


pyratbay.pyrat.crosssec
_______________________


.. py:module:: pyratbay.pyrat.crosssec


pyratbay.pyrat.clouds
_____________________


.. py:module:: pyratbay.pyrat.clouds


pyratbay.pyrat.rayleigh
_______________________


.. py:module:: pyratbay.pyrat.rayleigh


pyratbay.pyrat.alkali
_____________________


.. py:module:: pyratbay.pyrat.alkali


pyratbay.pyrat.optical_depth
____________________________


.. py:module:: pyratbay.pyrat.optical_depth


pyratbay.pyrat.spectrum
_______________________


.. py:module:: pyratbay.pyrat.spectrum

.. py:function:: spectrum(pyrat)
.. code-block:: pycon

    Spectrum calculation driver.

.. py:function:: modulation(pyrat)
.. code-block:: pycon

    Calculate modulation spectrum for transit geometry.

.. py:function:: intensity(pyrat)
.. code-block:: pycon

    Calculate the intensity spectrum [units] for eclipse geometry.

.. py:function:: flux(pyrat)
.. code-block:: pycon

    Calculate the hemisphere-integrated flux spectrum [units] for eclipse
    geometry.

.. py:function:: two_stream(pyrat)
.. code-block:: pycon

    Two-stream approximation radiative transfer
    following Heng et al. (2014)

    This function defines downward (flux_down) and uppward fluxes
    (flux_up) into pyrat.spec, and sets the emission spectrum as the
    uppward flux at the top of the atmosphere (flux_up[0]):

    flux_up: 2D float ndarray
        Upward flux spectrum through each layer under the two-stream
        approximation (erg s-1 cm-2 cm).
    flux_down: 2D float ndarray
        Downward flux spectrum through each layer under the two-stream
        approximation (erg s-1 cm-2 cm).

