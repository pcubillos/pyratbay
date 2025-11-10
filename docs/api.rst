.. _api:

API
===


pyratbay
________


.. py:module:: pyratbay

.. py:class:: Pyrat(cfg_file, log=True, mute=False)

    .. code-block:: pycon

        Parse the command-line arguments into the pyrat object.

        Parameters
        ----------
        cfg_file: String
            A Pyrat Bay configuration file.
        log: Bool
            Flag to save screen outputs to file (True) or not (False)
            (e.g., to prevent overwritting log of a previous run).
        mute: Bool
            If True, enforce verb to take a value of -1.

        Examples
        --------
        >>> import pyratbay as pb
        >>> pyrat = pb.run('spectrum_transmission.cfg')

    .. py:method:: band_contribution()
    .. code-block:: pycon

        Compute contribution functions or transmittance at each band.

    .. py:method:: band_integrate()
    .. code-block:: pycon

        Band-integrate transmission spectrum (transit) or planet-to-star
        flux ratio (eclipse) over transmission band passes.

    .. py:method:: compute_opacity()
    .. code-block:: pycon

        Calculate opacity (cm2 molecule-1) tabulated over
        temperature, pressure, and wavenumber arrays

    .. py:method:: eval(params, retmodel=True, skip=[])
    .. code-block:: pycon

        Fitting routine for atmospheric retrieval

        Parameters
        ----------
        params: 1D float iterable
            Array of fitting parameters that define the atmosphere.
        retmodel: Bool
            Flag to include the model spectra in the return.
        skip: List of strings
            If set, the opacity from the model names or line-sample species
            listed here will be neglected.

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

    .. py:method:: optical_depth()
    .. code-block:: pycon

        Calculate the optical depth.

    .. py:method:: plot_spectrum(spec='model', **kwargs)
    .. code-block:: pycon

        Plot spectrum.

        Parameters
        ----------
        spec: String
            Flag indicating which model to plot.  By default plot the
            latest evaulated model (spec='model').  Another option is
            'best', to plot the posterior best-fit (after a retrieval
            posterior run).
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
        temperature profile (self.atm.temp) and abundances (self.atm.vmr)
        with the values from the last radiative-equilibrium iteration.

        This method also defines self.atm.radeq_temps, a 2D array
        containing all temperature-profile iterations.

    .. py:method:: retrieval()
    .. code-block:: pycon

        Run an MCMC or nested-sampling atmospheric retrieval.

    .. py:method:: run(temp=None, vmr=None, radius=None, skip=[])
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
        skip: List of strings
            If set, the opacity from the model names or line-sample species
            listed here will be neglected.

.. py:class:: Atmosphere(inputs, wn=None, log=None)

    .. code-block:: pycon

        Initialize an atmospheric model object
        There are four main properties to compute (in this order):
        The pressure, the temperature, the volume mixing ratios (VMRs),
        and the radius profiles.

        Properties can be read from input profiles, computed by models,
        or skipped, depending on the configuration file.

        The rules are simple:
        - if there is a model in the config file, calculate the property
        - else if there is an input atmfile or ptfile, read properties from file
        - else, skip the calculation
        - if calculate p, any further reads (T,VMR,r) will interpolate

    .. py:method:: calc_profiles(temp=None, vmr=None, radius=None, abund=None)
    .. code-block:: pycon

        Update temperature, abundances, and radius profiles

        Parameters
        ----------
        temp: 1D float ndarray
            Layer's temperature profile (Kelvin) sorted from top to bottom.
        vmr: 2D float ndarray
            Species mole mixing ratio profiles [nlayers, nmol].
        radius: 1D float ndarray
            Layer's altitude profile (in cm), same order as temp.

    .. py:method:: parse_abundance_parameters()
    .. code-block:: pycon

        Sort out variables related to the VMR modeling

    .. py:method:: rad_model(pressure, temperature, mu, mplanet, gplanet, p0, r0)
    .. code-block:: pycon

        Calculate radius profile in hydrostatic equilibrium.

        Depending on self.rmodelname, select between a g(r)=GM/r**2
        (hydro_m) or constant-g (hydro_g) formula to compute
        the hydrostatic-equilibrium radii at each layer.

        Parameters
        ----------
        pressure: 1D float ndarray
            Atmospheric pressure for each layer (in bar).
        temperature: 1D float ndarray
            Atmospheric temperature for each layer (in K).
        mu: 1D float ndarray
            Mean molecular mass for each layer (in g mol-1).
        mplanet: Float
            Planetary mass in g (ignored for hydro_g).
        gplanet: Float
            Atmospheric gravity in cm s-2 (ignored for hydro_m).
        p0: Float
            Reference pressure level (in bar) where radius(p0) = r0.
        r0: Float
            Reference radius level (in cm) corresponding to p0.

    .. py:method:: setup_star_sed(inputs, wn)
    .. code-block:: pycon

        Read stellar spectrum model: starspec, kurucz, or blackbody

        Returns
        Input stellar flux spectrum (erg s-1 cm-2 cm)

    .. py:method:: validate_species(var, molec=None, elements=[])
    .. code-block:: pycon

        Validate composition variables.

.. py:function:: run(cfile, run_step=None, with_log=True)
.. code-block:: pycon

    Pyrat Bay initialization driver.

    Parameters
    ----------
    cfile: String
        A Pyrat Bay configuration file.
    run_step: String
        DEPRECATED
    with_log: Bool
        Flag to save screen outputs to file (True) or not (False)
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

  5.6703744191844314e-05

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

.. py:data:: parsec
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

  1.66053906892e-24

.. py:data:: me
.. code-block:: pycon

  9.1093837139e-28

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

  2.6867801117984436e+19

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

  1129583354672.9111

.. py:data:: C2
.. code-block:: pycon

  1.4387768775039338

.. py:data:: C3
.. code-block:: pycon

  8.852821669717024e-13

.. py:data:: ROOT
.. code-block:: pycon

  'os.path.realpath(os.path.dirname(__file__) + '/../..') + '/'

.. py:data:: FILTERS
.. code-block:: pycon

  'os.path.realpath(os.path.dirname(__file__) + '/../..') + '/pyratbay/data/filters/'

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

  ['tli', 'atmosphere', 'opacity', 'spectrum', 'radeq', 'retrieval']

.. py:data:: samplers
.. code-block:: pycon

  ['snooker', 'multinest']

.. py:data:: statistics
.. code-block:: pycon

  ['med_central', 'max_central', 'max_like', 'global_max_like']

.. py:data:: transmission_rt
.. code-block:: pycon

  ['transit']

.. py:data:: emission_rt
.. code-block:: pycon

  ['emission', 'emission_two_stream', 'f_lambda']

.. py:data:: eclipse_rt
.. code-block:: pycon

  ['eclipse', 'eclipse_two_stream']

.. py:data:: rt_paths
.. code-block:: pycon

  ['transit', 'eclipse', 'eclipse_two_stream', 'emission', 'emission_two_stream', 'f_lambda']

.. py:data:: retflags
.. code-block:: pycon

  ['temp', 'rad', 'press', 'mol', 'ray', 'cloud', 'patchy', 'mass', 'tstar']

.. py:data:: tmodels
.. code-block:: pycon

  ['isothermal', 'guillot', 'madhu']

.. py:data:: chemmodels
.. code-block:: pycon

  ['free', 'equilibrium']

.. py:data:: radmodels
.. code-block:: pycon

  ['hydro_m', 'hydro_g']

.. py:data:: amodels
.. code-block:: pycon

  ['sodium_vdw', 'potassium_vdw']

.. py:data:: rmodels
.. code-block:: pycon

  ['rayleigh_H', 'rayleigh_H2', 'rayleigh_He', 'rayleigh_e-']

.. py:data:: cmodels
.. code-block:: pycon

  ['deck', 'ccsgray', 'lecavelier']

.. py:data:: h_ion_models
.. code-block:: pycon

  ['h_ion_john1988']


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
        Monotonously decreasing pressure profile (in bar).
    temperature: 1D float ndarray
        Temperature profile for pressure layers (in Kelvin).
    species: 1D string ndarray
        List of atmospheric species.
    abundances: 2D float ndarray
        The species mole mixing ratio (of shape [nlayers,nspecies]).
    radius: 1D float ndarray
        Monotonously increasing radius profile (in cm).
    punits: String
        Pressure units of output.
    runits: String
        Radius units of output.
    header: String
        Header message (comment) to include at the top of the file.

    Examples
    --------
    >>> import numpy as np
    >>> import pyratbay.io as io
    >>> import pyratbay.atmosphere as pa

    >>> atmfile = 'WASP-00b.atm'
    >>> nlayers = 5
    >>> pressure = pa.pressure('1e-8 bar', '1e2 bar', nlayers)
    >>> temperature = pa.tmodels.Isothermal(pressure)(1500.0)
    >>> species = "H2 He H2O".split()
    >>> abundances = [0.8499, 0.15, 1e-4]
    >>> vmr = pa.uniform(abundances, nlayers)
    >>> io.write_atm(atmfile, pressure, temperature, species, vmr,
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
    vmr: 2D float ndarray
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
    ('bar', 'kelvin', 'volume', None)
    ['H2' 'He' 'H2O']
    [1.0000e-08 3.1623e-06 1.0000e-03 3.1623e-01 1.0000e+02]
    [1500. 1500. 1500. 1500. 1500.]
    [[8.499e-01 1.500e-01 1.000e-04]
     [8.499e-01 1.500e-01 1.000e-04]
     [8.499e-01 1.500e-01 1.000e-04]
     [8.499e-01 1.500e-01 1.000e-04]
     [8.499e-01 1.500e-01 1.000e-04]]
    None

.. py:function:: write_spectrum(wl, spectrum, filename, type)
.. code-block:: pycon

    Write a spectrum to file.

    Parameters
    ----------
    wl: 1D float iterable
        Wavelength array in micron units.
    spectrum: 1D float iterable
        Spectrum array. (rp/rs)**2 for transmission (unitless),
        planetary flux for emission (erg s-1 cm-2 cm units).
    filename: String
        Output file name.
    type: String
        Data type (only used for header comments):
        - 'transit' for transit spectra
        - 'eclipse' for secondary eclipse spectra
        - 'emission' for emission flux
        - 'filter' for a instrumental filter transmission

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
    >>> wl = np.linspace(1.1, 1.7, nwave)
    >>> spectrum = np.ones(nwave)
    >>> io.write_spectrum(wl, spectrum, 'sample_spectrum.dat', type='transit')
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

.. py:function:: write_spectra(spectra, wl, temperatures, filename)
.. code-block:: pycon

    Write flux spectra as function of wavelength and temperature to file.

    Parameters
    ----------
    spectra: 1D float iterable
        Flux spectra arrays (erg s-1 cm-2 cm units).
    wl: 1D float iterable
        Wavelength array in microns.
    temperatures: 1D float iterable
        Effective temperature array in Kelvin degrees.
    filename: String
        Output file name.

.. py:function:: read_spectra(filename)
.. code-block:: pycon

    Write flux spectra as function of wavelength and temperature to file.

    Parameters
    ----------
    filename: String
        Input file name.

    Returns
    -------
    spectra: 1D float iterable
        Flux spectra arrays (erg s-1 cm-2 cm units).
    wn: 1D float iterable
        Wavenumber array in cm-1.
    temperatures: 1D float iterable
        Effective temperature array in Kelvin degrees.

.. py:function:: write_opacity(ofile, species, temp, press, wn, opacity)
.. code-block:: pycon

    Write an opacity table as a binary npz file.

    Parameters
    ----------
    ofile: String
        Output filename where to save the opacity data.
        File extension must be .npz
    species: String
        The species name.
    temp: 1D float ndarray
        Temperature array (Kelvin degree).
    press: 1D float ndarray
        Pressure array (bar).
    wn: 1D float ndarray
        Wavenumber array (cm-1).
    opacity: 3D float ndarray
        Tabulated opacities (cm2 molecule-1) of shape [ntemp, nlayers, nwave].

.. py:function:: read_opacity(ofile, extract='all')
.. code-block:: pycon

    Read an opacity table from file.
    Compatible with petitRADTRANS3 opacity files as well.

    Parameters
    ----------
    ofile: String
        Path to a Pyrat Bay opacity file.
    extract: String
        Information to extract, select between:
        - 'arrays' for the species, temperature, pressure, and wavenumber
        - 'opacity' for the opacity grid
        - 'all' to get both and the array dimensions.

    Returns
    -------
    units: dict
        The physical units for the different quantities
    species: String
        The species name.
    temp: 1D float array
        The temperature array (K)
    press: 1D float array
        The pressure array (bar)
    wn: 1D float array
        The wavenumber array (cm-1)
    opacity: 3D float ndarray tuple
        The tabulated opacities (cm2 molecule-1), of shape
        [ntemp, nlayers, nwave].

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

.. py:function:: write_observations(obs_file, inst_names, wl, wl_half_width, depth=None, depth_err=None, depth_units='none')
.. code-block:: pycon

    Write an observation file for use in pyrat bay.
    These can be a combination of tophat pass bands (non-zero wl) or
    paths to files with tabulated passbands (zero wl value).

    Parameters
    ----------
    obs_file: String
        Name of observation file.
    inst_names: 1D string iterable
        Name for the bandpasses.  If this is of type str, use
        the same name for all bands.
    wl: 1D string ndarray
        Central wavelength of the bands (in microns).
        If wl is zero for a data point, assume that the inst_name
        is a file path to a tabulated pass band.
    wl_half_width: 2D float ndarray
        Bandpass half-width (in microns).
    data: 1D float ndarray
        Transit of eclipse depth corresponding to  each band.
        If not None, include this data into the observation file.
    depth_err: 1D float ndarray
        Depth uncertainties.
    depth_units:  String
        Units of input depth data.

    Examples
    --------
    >>> import pyratbay.io as io

    >>> # Observation file with only the pass bands:
    >>> wl = [2.144, 2.333, 2.523]
    >>> half_widths = [0.095, 0.095, 0.095]
    >>> io.write_observations('obs_file.txt', 'HST', wl, half_widths)

    >>> # Observation file with pass bands and data:
    >>> data = np.array([329.6, 344.5, 301.4])
    >>> uncert = np.array([20.4, 21.9, 23.5])
    >>> io.write_observations(
    >>>     'obs_file.txt', 'HST', wl, half_widths, data, uncert, 'ppm',
    >>> )

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
    depth: 1D string list
        The transit or eclipse depths for each filter.
    depth_err: 1D float ndarray
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
    >>> bands, depth, depth_err = io.read_observations(obs_file)

    >>> # File including only the passband info:
    >>> obs_file = 'filters.dat'
    >>> bands = io.read_observations(obs_file)

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
    >>> mol, hit_iso, exo_iso, ratio, mass = \
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
        The cross-section source: 'exomol' or 'taurex' (see note below).
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
        Pressure sample of the opacity file (in bars)
    temperature: 1D float ndarray
        Temperature sample of the opacity file (in Kelvin degrees).
    wavenumber: 1D float ndarray
        Wavenumber sample of the opacity file (in cm-1).
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
    >>> xs, press, temp, wn, species = io.import_xs(filename, 'exomol')

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


pyratbay.tools
______________


.. py:module:: pyratbay.tools

.. py:function:: parse_error_param(var)
.. code-block:: pycon

    Parse error-scaling parameters. There are two options:
    - err_scale_name: Scale uncertainties as a multiplicative factor.
    - err_quad_name: Scale uncertainties by adding in quadrature.

    Parameters
    ----------
    var: String
        Parameter name. Must begin with either "err_scale_" or "err_quad_".

    Returns
    -------
    inst: String
        Instrument name, parsed from the end of var.
    texname: String
        Instrument name as a tex string, for plotting purposes.
    scaling: String
        What type of uncertainty scaling: "scale" or "quadrature".

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> # Valid options:
    >>> print(pt.parse_error_param('err_scale_WFC3'))
    ('WFC3', '$S_\\sigma^{\\rm WFC3}$', 'scale')
    >>> print(pt.parse_error_param('err_quad_NIRSpec_PRISM'))
    ('NIRSpec PRISM', '$\\sigma_{\\rm NIRSpec PRISM}$', 'quadrature')

    >>> # Not valid scaling throws error:
    >>> print(pt.parse_error_param('err_fudging_IRAC1'))
    ValueError: Invalid error scaling parameter 'fudging_IRAC1'. Valid options begin with: ['scale_', 'quad_']

.. py:class:: Data(data, uncert, band_names, offset_models=None, err_models=None)

    .. code-block:: pycon

    .. code-block:: pycon

        Parameters
        ----------
        data: 1D float iterable
            The data values.
        uncert: 1D float iterable
            The uncertainty values.
        band_names: 1D string iterable
            Names for the uncertainties.
        offset_models: String or 1D iterable of strings
            List of data offset models, the strings must match
            a substring of at least one of the band_names, these specific
            data points will be affected by the respective offset model.
        err_models: String or 1D iterable of strings
            List of error inflation model names, must begin with either
            - "err_scale_" to scale uncertainties as a multiplicative factor
            - "err_quad_" to scale uncertainties by adding in quadrature
            followed by a string matching a substring of at least one of
            the band_names, these specific data points will be affected by
            the error scaling model.

        Examples
        --------
        >>> import pyratbay.tools as pt
        >>> import pyratbay.constants as pc
        >>> import pyratbay.io as io
        >>> import matplotlib.pyplot as plt

        >>> # Load a set of uncertainties obtained with JWST instruments:
        >>> obs_file = '/Users/pato/Documents/compendia/ERS/WASP39b/data/synthesis_v02/wasp39b_g395h_lrs.dat'
        >>> bands, depths, uncert = io.read_observations(obs_file)
        >>> band_names = [band.name for band in bands]
        >>> wl = [band.wl0 for band in bands]
        >>> print(set(band_names))
        {'nirspec_g395h_nrs1', 'nirspec_g395h_nrs2', 'miri_lrs'}

        >>> # Offsets for MIRI data only:
        >>> offsets = 'offset_miri'
        >>> obs = pt.Data(depths, uncert, band_names, offset_models=offsets)
        >>> offset_depths = obs.offset_data([400.0], 'ppm')

        >>> fig = plt.figure(0, (8,4))
        >>> plt.clf()
        >>> plt.semilogx(wl, obs.data/pc.ppm, 'o', color='darkorange')
        >>> plt.plot(wl, offset_depths/pc.ppm, '^' , color='blue', mfc='none')
        >>> plt.xlabel('Wavelength (um)')
        >>> plt.ylabel('Depths (ppm)')

        >>> # Multiplicative error scaling (increase by 1.5x):
        >>> # Target NRS1 uncertainties specifically of the NIRSpec instrument
        >>> obs = pt.Data(depths, uncert, band_names, err_models='err_scale_nrs1')
        >>> inflated_err = obs.scale_errors([0.3])

        >>> fig = plt.figure(1, (8,4))
        >>> plt.clf()
        >>> plt.semilogx(wl, obs.uncert/pc.ppm, 'o', color='xkcd:green')
        >>> plt.plot(wl, inflated_err/pc.ppm, '^' , color='blue', mfc='none')
        >>> plt.xlabel('Wavelength (um)')
        >>> plt.ylabel('Uncertainties (ppm)')

        >>> # Add-in-quadrature error scaling:
        >>> # Target all NIRSpec uncertainties (200ppm), for MIRI add 300ppm
        >>> err_models = ['err_quad_nirspec', 'err_quad_miri']
        >>> obs2 = pt.Data(depths, uncert, band_names, err_models=err_models)
        >>> inflated_err2 = obs2.scale_errors([2.3, 2.5], 'ppm')
        >>> plt.plot(wl, inflated_err2/pc.ppm, 's', color='orange', mfc='none')
        >>> plt.tight_layout()

    .. py:method:: offset_data(vals=None, units='none')
    .. code-block:: pycon

        Offset data values

        Parameters
        ----------
        vals: 1D float iterable
            Data offset values for each model.
        units: String
            The units of the scaling value. Options are:
            'none', 'percent', 'ppt', or 'ppm'.

        Returns
        -------
        data: 1D float array
            The offset data array.

    .. py:method:: scale_errors(vals=None, units='none')
    .. code-block:: pycon

        Scale uncertainties according to input scaling values.

        Parameters
        ----------
        vals: 1D float iterable
            Uncertainty scaling value (in log10).
        units: String
            The units of the scaling value for quadrature errors.
            Options are: 'none', 'percent', 'ppt', or 'ppm'.

        Returns
        -------
        uncert: 1D float array
            The scaled uncertainty array.

.. py:function:: check_mpi4py()
.. code-block:: pycon

    Detect when the code was called with MPI and mpi4py module is missing

    Only raise an error when needed (more than one processor required),
    otherwise you might be running multiple runs in parallel but not
    talking to each other.

.. py:function:: check_mpi_is_needed(inputs)
.. code-block:: pycon

    Prevent using parallel processes through MPI when not needed
    (only required for MultiNest runs).

.. py:function:: get_mpi_rank()
.. code-block:: pycon

    Get the MPI rank of the current process (intended for MPI runs).
    If mpi4py is not installed, return zero.

    Returns
    -------
    rank: Interger
        The MPI process rank.

.. py:function:: get_mpi_size()
.. code-block:: pycon

    Get the size of the current group of process (intended for MPI runs).
    If mpi4py is not installed, return one.

    Returns
    -------
    size: Interger
        The size of the MPI group of processes.

.. py:function:: mpi_barrier()
.. code-block:: pycon

    Make an MPI barrier() call. Ignore it if mpi4py is not installed.

.. py:class:: Namespace(args=None, log=None)

    .. code-block:: pycon

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

    .. py:method:: get_path(pname, desc='', exists=False, make_dir=False)
    .. code-block:: pycon

        Extract pname file path (or list of paths) from Namespace,
        return the canonical path.

        Parameters
        ----------
        pname: String
            The parameter name to extract.
        desc: String
            A description to display in case of raising an error
        exists: Bool
            If True, raise an error if the file path does not exists.
        make_dir: Bool
            If True, make directories for the file path if needed.

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

.. py:function:: parse(cfile, with_log=True, mute=False)
.. code-block:: pycon

    Read the command line arguments.

    Parameters
    ----------
    cfile: String
        A Pyrat Bay configuration file.
    with_log: Bool
        Flag to save screen outputs to file (True) or not (False)
        (e.g., to prevent overwritting log of a previous run).
    mute: Bool
        If True, enforce verb to take a value of -1.

    Returns
    -------
    args: Dictionary
        A dictionary containing the input configuration variables.
    log: mc3.Log
        A logging log object.

.. py:function:: parse_bool(args, param, default=False)
.. code-block:: pycon

    Parse a string parameter in args into a Bool.

.. py:function:: parse_str(args, param)
.. code-block:: pycon

    Parse a string parameter in args into a string.

.. py:function:: parse_int(args, param)
.. code-block:: pycon

    Convert a dictionary's parameter from string to integer.
    Set parameter to None if it was not in the dictionary.

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

.. py:function:: parse_var_vals(var_input)
.. code-block:: pycon

    Parse keys that contain variable and values

.. py:class:: Loglike(pyrat)

    .. code-block:: pycon

        Initialize self.  See help(type(self)) for accurate signature.

.. py:function:: weighted_to_equal(posterior_file, get_weighted=False, min_size=15000)
.. code-block:: pycon

    Compute an equally-weighted sample from a weighted-probability sample
    read from a Multinest output.

    Parameters
    ----------
    posterior_file: String
        A MultiNest probability-weighted sample output.
    get_weighted: Bool
        If True, also return the weighted sample.
    min_size: Integer
        Set the minimum sample size for the equally weighted posterior.

    Returns
    -------
    equal_posterior: 2D float array
        An equally-weighted posterior sample with dimensions (nsamples, npars).
    weighted_posterior: 2D float array
        The Multinest probabilty-weighted sample with dimensions
        (nsamples, npars).  This is only returned if get_weighted is True.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> posterior = pt.weighted_to_equal('multinest_output.txt')

    >>> # Bet both equal and weighted samples:
    >>> posterior, weighted = pt.weighted_to_equal(
    >>>     'multinest_output.txt',
    >>>     get_weighted=True,
    >>> )

.. py:function:: posterior_snapshot(retrieval_file, pnames)
.. code-block:: pycon

    Take a snapshot of a retrieval run, plot the histogram and traces
    of the parameters.

.. py:function:: get_multinest_map(stats_file)
.. code-block:: pycon

    Get maximum-a-posteriori (MAP) parameters from a MultiNest output file.

    Parameters
    ----------
    stats_file: String
        Path to a Multinest *stats.dat output file.

    Returns
    -------
    params: 1D float array
        The MAP parameter values.

.. py:function:: multinest_run(pyrat, basename)
.. code-block:: pycon

    A Wrapper of a MultiNest posterior sampling.

    Parameters
    ----------
    pyrat: Pyrat() object
    basename: String
        Basename for output files. May contain path.
        Should not contain a file extension.

    Note
    ----
    For OS X users, it is recommended to set the TMPDIR environment
    variable to "/tmp", e.g., from the command line:
        export TMPDIR=/tmp
    to avoid an MPI error when terminating the execution
    (the call will run to completion in any case)
    https://github.com/open-mpi/ompi/issues/7393#issuecomment-882018321

.. py:function:: posterior_post_processing(cfg_file=None, pyrat=None, suffix='')
.. code-block:: pycon

    Compute quantities of interest from a retrieval posterior distribution.
    The produced data is stored into a pickle file with root name based
    on the logfile.

    Parameters
    ----------
    cfg_file: String
        A pyratbay config file of a retrieval run (already executed,
        so the parameter posterior files must already exist).
    pyrat: a Pyrat instance
        A pyrat object of an already executed retrieval.
        Used if cfg_file is None.

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

.. py:function:: eta(time_seconds, n_completed, n_total, fmt='.2f')
.. code-block:: pycon

    Find most appropriate units to report the remaining time
    (seconds, minutes, hours, days)

    Parameters
    ----------
    time_seconds: Float
        An amount of time in seconds.
    n_completed: Integer
        Number of completed steps.
    n_total: Integer
        Total number of steps to complete.

    Returns
    -------
    delta_time: Float
        The time_seconds in the recalculated units.

.. py:function:: resolve_theme(theme)
.. code-block:: pycon

    Resolve input Theme or color into a mc3.plots.Theme instance.
    Makes sure that input is either None, a mc3.plots.Theme, or
    a value that can be interpreted as a matplotlib color.

    Parameters
    ----------
    theme: Any
        A matplotlib color or a mc3.plots.Theme instance

    Returns
    -------
    theme: mc3.plots.Theme instance
        A Theme computed using the input color.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> import mc3

    >>> # A Theme instance is returned unmodified
    >>> theme = pt.resolve_theme(mc3.plots.THEMES['indigo'])
    >>> # Anything that can be interpreted as matplolib color:
    >>> theme1 = pt.resolve_theme('red')
    >>> theme2 = pt.resolve_theme('xkcd:green')
    >>> theme3 = pt.resolve_theme((0,0,1))

    >>> # If input is None, return None
    >>> theme = pt.resolve_theme(None)

    >>> # Anything else will throw an error:
    >>> theme = pt.resolve_theme('not_a_plt_color')
    ValueError: Invalid color theme: 'not_a_plt_color'

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

.. py:function:: is_number(value)
.. code-block:: pycon

    Check whether a string value can be parsed as a number.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> # These return True
    >>> pt.is_number('1')
    >>> pt.is_number('1.0')
    >>> pt.is_number('-3.14')
    >>> pt.is_number('+3.14')
    >>> pt.is_number('1.0e+02')
    >>> pt.is_number('inf')
    >>> pt.is_number('nan')

    >>> # These return False
    >>> pt.is_number('1.0-3.14')
    >>> pt.is_number('10abcde')
    >>> pt.is_number('1.0e')
    >>> pt.is_number('true')
    >>> pt.is_number('none')

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

.. py:function:: mkdir(file_path)
.. code-block:: pycon

    Create a directory for given file_path if it doesn't exists.
    Creating nested folders is not allowed.

    Parameters
    ----------
    file_path: String
        Path to a file.

    Examples
    --------
    >>> import pyratbay.tools as pt
    >>> log_file = 'NS1/ns_emission_tutorial.log'
    >>> pt.mkdir(log_file)

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

        Initialize self.  See help(type(self)) for accurate signature.

.. py:function:: get_exomol_mol(file)
.. code-block:: pycon

    Parse an exomol file to extract the molecule and isotope name.

    Parameters
    ----------
    file: String
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

.. py:function:: interpolate_opacity(cs_file, temperature=None, pressure=None, wn_mask=None, wl_thinning=1)
.. code-block:: pycon

    Interpolate the cross-section data from an opacity file over a
    desired temperature and pressure array.

    Parameters
    ----------
    cs_file: String
        Path to a cross-section file.
    temperature: 1D float array
        The desired temperature array in K.
        If this is the same as the tabulated temperatures, do not interpolate.
    pressure: 1D float array
        The desired pressure profile in bars.
        If this is the same as the tabulated pressure, do not interpolate.
    wn_mask: 1D bool array
        A mask of wavelength points to take.
    wl_thinning: Integer
        Thinning factor to take every n-th value of the wavenumber array

    Returns
    -------
    interp_cs: 4D float array
        The interpolated cross-section array.

.. py:function:: none_div(a, b)
.. code-block:: pycon

    Non-breaking division when values are None.

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


pyratbay.opacity
________________


.. py:module:: pyratbay.opacity

.. py:class:: Collision_Induced(cia_file, *, wn=None, wl=None, log=None)

    .. code-block:: pycon

        Parameters
        ----------
        cia_file: String
            A CIA cross section file.
        wn: 1D float array
            Wavenumber array (cm-1 units) where to sample the CIA
            (only one of wl or wn should be provided).
        wl: 1D float array
            Wavelength array (micron units) where to sample the CIA
            (only one of wl or wn should be provided).

        Examples
        --------
        >>> import pyratbay.spectrum as ps
        >>> import pyratbay.constants as pc
        >>> import pyratbay.opacity as op

        >>> wl = ps.constant_resolution_spectrum(0.61, 10.01, 15000.0)
        >>> cs_file = f'{pc.ROOT}/pyratbay/data/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat'
        >>> cia = op.Collision_Induced(cs_file, wl=wl)

    .. py:method:: calc_cross_section(temperature)
    .. code-block:: pycon

        Calculate cross-section spectra (cm-1) over temperature
        profiles by interpolating from tabulated values.

        Parameters
        ----------
        temperature: float or 1D float array
            Temperature array (Kelvin) at which to interpolate the
            cross section.

        Returns
        -------
        cross_section: 1D or 2D float array
            Cross section spectra (cm-1 /(molec cm-3)**N)
            where N is the number of species.
            If temperature is scalar, cross_section is 1D of length self.nwave.
            Otherwise, cross_section is a 2D array of shape [nlayers,nwave]

    .. py:method:: calc_extinction_coefficient(temperature, density)
    .. code-block:: pycon

        Calculate extinction-coefficient spectra (cm-1) over temperature
        and density profiles by interpolating in temperature from
        tabulated values.

        Parameters
        ----------
        temperature: 1D float array
            Temperature array (Kelvin) at which to interpolate the
            cross section.
        density: 2D float array
            Number density array (molecules cm-3) of self.species
            over the temperature array.
            If temperature is 1D array, must be of shape [self.nspec, ntemp]
            If temperature is scalar, must be of length self.nspec.

        Returns
        -------
        extinction_coefficient: 1D or 2D float array
            Extinction coefficient spectra (cm-1)
            If temperature is scalar, output is a 1D array of length nwave.
            Otherwise, outout is a 2D array of shape [nlayers,nwave]

.. py:class:: Hydrogen_Ion(wn)

    .. code-block:: pycon

        Parameters
        ----------
        wn: 1D float array
            Wavenumber array where to sample the opacities (cm-1).

    .. py:method:: calc_extinction_coefficient(temperature, density, layer=None)
    .. code-block:: pycon

        Calculate extinction coefficient spectra (cm-1) for given
        temperature and number-density profiles.

        Parameters
        ----------
        temperature: Float or 1D float array
            Temperature profile in Kelvin (of size ntemp)
        density: 1D/2D float array
            Number density profiles (molecules per cm3) of H and electron
            over the given temperature array.
            If temperature is 1D array, must be of shape [ntemp,2]
            If temperature is scalar, density must be of size 2.
        layer: Integer
            If not None, compute extinction coefficient only at selected layer.
            In this case, the output array dimension is reduced by 1.

        Returns
        -------
        extinction_coefficient: 1D/2D float array
            H- extinction coefficient spectrum (cm-1)
            If temperature is scalar, output is a 1D array of length nwave.
            Otherwise, outout is a 2D array of shape [ntemp,nwave]

    .. py:method:: cross_section_bound_free(temperature)
    .. code-block:: pycon

        Compute the bound-free cross section for H- in cm5/H_mol/electron.

        Equation (3) of John 1988, AA, 193, 189.

    .. py:method:: cross_section_free_free(temperature)
    .. code-block:: pycon

        Compute the free-free cross section for H- in cm5/H_mol/electron.

        Equation (6) of John 1988, AA, 193, 189.

.. py:function:: make_tli(dblist, pflist, dbtype, tlifile, wl_low, wl_high, wl_units, log=None)
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
    wl_low: Float
        Lower wavelength boundary to consider, units given by wl_units.
    wl_high: Float
        High wavelength boundary to consider, units given by wl_units.
    wl_units: String
        Wavelength units (when not specified in wl_low nor wl_high).
    log: Log object
        An mc3.utils.Log instance to log screen outputs to file.

.. py:class:: Line_Sample(cs_files, *, pressure=None, temperature=None, min_wl=None, max_wl=None, min_wn=None, max_wn=None, isotope_ratios=None, wl_thinning=1, log=None)

    .. code-block:: pycon

        Read line-sampled cross-section table(s), with units of cm2 molec-1.

        Parameters
        ----------
        cs_files: String or iterable of strings
            Line-sampled cross section file(s) to read.
        pressure: 1D floar array
            Pressure profile where to resample the opacities (bar).
            If None, use the tabulated pressure array from the opacities.
            If not None, it is allowed to extrapolate to lower pressures
            but not to higher pressures.
        min_wl: 1D float ndarray
            Minimum wavelength value to extract from line-sample files (um)
            (only one of min_wl or max_wn should be provided).
        max_wl: 1D float ndarray
            Maximum wavelength value to extract from line-sample files (um)
            (only one of min_wn or max_wl should be provided).
        min_wn: 1D float ndarray
            Minimum wavenumber value to extract from line-sample files (cm-1)
            (only one of min_wn or max_wl should be provided).
        max_wn: 1D float ndarray
            Maximum wavenumber value to extract from line-sample files (cm-1)
            (only one of min_wl or max_wn should be provided).
        wl_thinning: Integer
            Thinning factor to take every n-th value of the wavenumber array
        isotope_ratios: String

        Examples
        --------
        >>> import pyratbay.opacity as op
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt

        >>> # Generate these files from
        >>> # pyratbay.rtfd.io/en/latest/cookbooks/line_list_hitran.html
        >>> cs_files = [
        >>>     'cross_section_R020K_0150-3000K_0.5-5.0um_hitran_H2O.npz',
        >>>     'cross_section_R020K_0150-3000K_0.5-5.0um_hitemp_CO.npz',
        >>> ]
        >>> # Line-sample opacities at wavenumbers > 1.0 cm-1:
        >>> ls = op.Line_Sample(cs_files, min_wl=1.0)

        >>> # Total or per-mol cross sections (cm2 molec-1):
        >>> temp = np.tile(1400.0, ls.nlayers)
        >>> cs = ls.calc_cross_section(temp)
        >>> cs_per_mol = ls.calc_cross_section(temp, per_mol=True)

        >>> # Cross section per species (cm2 molec-1):
        >>> ec = ls.calc_cross_section(temp, per_mol=True)

        >>> # Take a look
        >>> wl = 1e4/ls.wn
        >>> plt.figure(1)
        >>> plt.clf()
        >>> plt.plot(wl, ec[0,35], color='blue', lw=1.0)
        >>> plt.plot(wl, ec[1,35], alpha=0.6, color='orange', lw=1.0)
        >>> plt.yscale('log')

    .. py:method:: calc_cross_section(temperature, layer=None, per_mol=False, pars=None)
    .. code-block:: pycon

        Calculate cross-section spectra (cm2 molec-1) over temperature
        profiles by interpolating from tabulated values.

        Parameters
        ----------
        temperature: 1D float array
            Temperature array (Kelvin) at which to interpolate the
            cross section (must match nlayers size).
        layer: Integer
            If not None, compute cross-sections only at selected layer.
            In this case, the output array dimension is reduced by 1.
        per_mol: bool
            If True, compute cross sections individually per species.
            If False, co-add cross section contributions from all species.
        pars: 1D iterable
            If not None, update the iso_ratio parameters with given values.

        Returns
        -------
        cross_section: 1D, 2D, or 3D float array
            Cross section spectra (cm2 molec-1)
            Output array has shape [nspec, nlayers, nwave].
            If per_mol is False, nspec dimension is removed.
            If layer is not None, nlayers dimension is removed.

    .. py:method:: calc_extinction_coefficient(temperature, density, layer=None, per_mol=False, pars=None)
    .. code-block:: pycon

        Calculate extinction-coefficient spectra (cm-1) for temperature
        and number-density profiles.

        Parameters
        ----------
        temperature: 1D float array
            Temperature array (Kelvin) at which to interpolate the
            cross section (must match nlayers size).
        density: 2D float array
            Number-density profiles (molec cm-3) at each layer.
            Array has shape [nspec,nlayers].
        layer: Integer
            If not None, compute extinction coefficient only at selected layer.
            In this case, the output array dimension is reduced by 1.
        per_mol: bool
            If True, compute extinction coefficients individually per species.
            If False, co-add extinction contributions from all species.
        pars: 1D iterable
            If not None, update the iso_ratio parameters with given values.

        Returns
        -------
        extinction: 1D, 2D, or 3D float array
            Extinction-coefficient spectra (cm-1)
            Output array has shape [nspec, nlayers, nwave].
            If per_mol is False, nspec dimension is removed.
            If layer is not None, nlayers dimension is removed.

    .. py:method:: get_wl(units='um')
    .. code-block:: pycon

        Get the wavelength array in the desired units.

        Parameters
        ----------
        units: String
            Select one from: ['A', 'nm', 'um', 'mm', 'cm', 'm', 'km']

        Returns
        -------
        wl: 1D float ndarray
            The wavelength array from the line-sample grid.

.. py:function:: optical_depth(rt_path, extinction, radius=None, itop=0, ibottom=None, maxdepth=inf, extinction_cloudy=None, raypath=None)
.. code-block:: pycon

    Calculate the optical depth.

    Parameters
    ----------
    rt_path: String
        Radiative-transfer path.
    extinction: 2D float array
        Extinction coefficient (cm-1) at each layer and wavelength channel.
    radius: 1D float array
        Radius (altitude) array of the atmospheric layers (cm), from top
        to bottom.  Used to compute the raypath.
    itop: Integer
        Index at top of the atmosphere, opacity at layer above this
        will be ignored.
    ibottom: Integer
        Index at the bottom of the atmosphere. The optical depth will
        be calculated down to the layer pointed by this index.
    maxdepth: Float
        Maximum optical depth to compute.  The optical depth calculation
        will stop when maxdepth is reached (checked independently in each
        wavelength channel).
    extinction_couldy: 2d float array
        If provided, signals that the optical depth should compute
        separately the depth of a cloudy and a clear atmosphere.
        This array contains the extinction coefficient of the 'cloud'
        absorbers.
    raypath:
        The ray paths, which can be directly provided instead of radius.

    Returns
    -------
    raypath:
        Path followed by the rays.
    depth: 2D float array
        The optical depth at each layer and wavelength channel.
    ideep: 1D integer array
        Indices of each wavelenght channel of the layer where the optical
        depth surpassed maxdepth.
    depth_clear: 2D float array
        If extinction_couldy is provided, the optical depth at each
        layer and wavelength channel for a cloud-less atmosphere.
    ideep_clear: 1D integer array
        Same as ideep but for depth_clear.


pyratbay.opacity.alkali
_______________________


.. py:module:: pyratbay.opacity.alkali

.. py:class:: SodiumVdW(pressure, *, wn=None, wl=None, cutoff=4500.0)

    .. code-block:: pycon

        Parameters
        ----------
        pressure: 1D float array
            Pressure profile (bars) over which the opacities will be
            evalulated.
        wn: 1D float array
            Wavenumber array (cm-1 units) over which the opacities
            will be sampled (only one of wl or wn should be provided).
        wl: 1D float array
            Wavelength array (micron units) over which the opacities
            will be sampled (only one of wl or wn should be provided).
        cutoff: Float
            Maximum wavenumber extent (cm-1) of the line-profiles
            from the center of each line.

    .. py:method:: calc_cross_section(temperature, layer=None)
    .. code-block:: pycon

        Calculate cross section (cm2 molecule-1) at given temperatures
        Saves value into self.cross_section.

        Parameters
        ----------
        temperature: 1D float array
            Temperature profile in Kelvin (must match self.pressure array)
        layer: Interger
            If not None, index at which to calculate the cross section.

        Returns
        -------
        cross_section: 2D/1D float array
            If layer is None, cross-section spectra (cm2 molec-1) at each
            pressure-temperature point (also sets self.cross_section).
            If layer is not None, cross-section spectrum at a single layer.

    .. py:method:: calc_extinction_coefficient(temperature, density, layer=None)
    .. code-block:: pycon

        Calculate extinction coefficient (cm-1) for given temperature
        and number-density profiles.

        Parameters
        ----------
        temperature: 1D float array
            Temperature profile in Kelvin (must match self.pressure array)
        density: 1D float array
            Number density profile (molecules per cm3) of self.species.
            (must match self.pressure array)
        layer: Interger
            If not None, calculate the extinction coefficient at a
            single layer pointed by this index. In this case the output
            is a 1D array.

        Returns
        -------
        extinction_coefficient: 2D float array
            Alkali extinction coefficient spectrum (units of cm-1)
            of shape [nlayers, nwave].

    .. py:method:: voigt_det(temperature)
    .. code-block:: pycon

        Calculate Voigt profile value at the detuning wavelength
        from the center of the lines (at each layer).

        Parameters
        ----------
        temperature: 1D float array
            A temperature profile (kelvin).

        Returns
        -------
        voigt_det: 2D float array
            Voigt-profile values at the detuning wavelength
            of shape [nlayers,nlines].

.. py:class:: PotassiumVdW(pressure, *, wn=None, wl=None, cutoff=4500.0)

    .. code-block:: pycon

        Parameters
        ----------
        pressure: 1D float array
            Pressure profile (bar) over which the opacities will be
            evalulated.
        wn: 1D float array
            Wavenumber array (cm-1 units) over which the opacities
            will be sampled (only one of wl or wn should be provided).
        wl: 1D float array
            Wavelength array (micron units) over which the opacities
            will be sampled (only one of wl or wn should be provided).
        cutoff: Float
            Maximum wavenumber extent (cm-1) of the line-profiles
            from the center of each line.

    .. py:method:: calc_cross_section(temperature, layer=None)
    .. code-block:: pycon

        Calculate cross section (cm2 molecule-1) at given temperatures
        Saves value into self.cross_section.

        Parameters
        ----------
        temperature: 1D float array
            Temperature profile in Kelvin (must match self.pressure array)
        layer: Interger
            If not None, index at which to calculate the cross section.

        Returns
        -------
        cross_section: 2D/1D float array
            If layer is None, cross-section spectra (cm2 molec-1) at each
            pressure-temperature point (also sets self.cross_section).
            If layer is not None, cross-section spectrum at a single layer.

    .. py:method:: calc_extinction_coefficient(temperature, density, layer=None)
    .. code-block:: pycon

        Calculate extinction coefficient (cm-1) for given temperature
        and number-density profiles.

        Parameters
        ----------
        temperature: 1D float array
            Temperature profile in Kelvin (must match self.pressure array)
        density: 1D float array
            Number density profile (molecules per cm3) of self.species.
            (must match self.pressure array)
        layer: Interger
            If not None, calculate the extinction coefficient at a
            single layer pointed by this index. In this case the output
            is a 1D array.

        Returns
        -------
        extinction_coefficient: 2D float array
            Alkali extinction coefficient spectrum (units of cm-1)
            of shape [nlayers, nwave].

    .. py:method:: voigt_det(temperature)
    .. code-block:: pycon

        Calculate Voigt profile value at the detuning wavelength
        from the center of the lines (at each layer).

        Parameters
        ----------
        temperature: 1D float array
            A temperature profile (kelvin).

        Returns
        -------
        voigt_det: 2D float array
            Voigt-profile values at the detuning wavelength
            of shape [nlayers,nlines].

.. py:function:: get_model(name, *args, **kwargs)
.. code-block:: pycon

    Get an alkali model by its name.


pyratbay.opacity.broadening
___________________________


.. py:module:: pyratbay.opacity.broadening

.. py:class:: Lorentz(x0=0.0, hwhm=1.0, scale=1.0)

    .. code-block:: pycon

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
        Atmospheric pressure (bar).
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
    >>> pressure = 1.0  # bar
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
            Minimum atmospheric pressure (bar).

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
        >>> min_press = 1e-5 # bar
        >>> dmin, lmin = b.min_widths(min_temp, max_temp, min_wn, max_mass,
        >>>     min_rad, min_press)
        >>> print('Minimum Doppler half width: {:.2e} cm-1
    '
        >>>       'Minimum Lorentz half width: {:.2e} cm-1'.format(dmin,lmin))
        Minimum Doppler half width: 8.44e-04 cm-1
        Minimum Lorentz half width: 2.21e-07 cm-1


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
            Maximum atmospheric pressure (bar).

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
        >>> max_press = 100.0 # bar
        >>> dmax, lmax = b.max_widths(min_temp, max_temp, max_wn, min_mass,
        >>>     max_rad, max_press)
        >>> print('Maximum Doppler half width: {:.2e} cm-1
    '
        >>>       'Maximum Lorentz half width: {:.2e} cm-1'.format(dmax,lmax))
        Maximum Doppler half width: 4.62e-02 cm-1
        Maximum Lorentz half width: 1.21e+01 cm-1



pyratbay.opacity.clouds
_______________________


.. py:module:: pyratbay.opacity.clouds

.. py:class:: CCSgray(pressure, wn)

    .. code-block:: pycon

        Initialize self.  See help(type(self)) for accurate signature.

    .. py:method:: calc_cross_section()
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

.. py:class:: Deck(pressure, wn)

    .. code-block:: pycon

        Initialize self.  See help(type(self)) for accurate signature.

    .. py:method:: calc_extinction_coefficient(radius, temperature, pars=None, layer=None)
    .. code-block:: pycon

        Calculate gray-cloud deck that's transparent above ptop,
        and becomes instantly opaque at ptop, with
        ptop (bar) = 10**pars[0].

        Parameters
        ----------
        radius: 1D float ndarray
            Atmospheric radius profile (in cm).
        temperature: 1D float ndarray
            Atmospheric temperature profile (in Kelvin degree).
        pars: 1D float iterable
            If not None, update the model parameters with input values
        layer: integer
            If not None, check whether the cloud top is above or below
            the given atmospheric layer at this index.

.. py:class:: Lecavelier(pressure, wl=None, wn=None)

    .. code-block:: pycon

        Initialize self.  See help(type(self)) for accurate signature.

    .. py:method:: calc_cross_section(pars=None)
    .. code-block:: pycon

        Calculate the Rayleigh cross section in cm2 molec-1:
            cross section = k_ray * s0 * (lambda/l0)**alpha_ray,
        parameterized as params = [log10(k_ray), alpha_ray].

        Parameters
        ----------
        pars: 1D iterable
            If not None, update the model parameters with given values.

        Returns
        -------
        cross_section: 1D float array
            Cross section (cm2 molecule-1) as function of wavelength

    .. py:method:: calc_extinction_coefficient(temperature, pars=None, layer=None)
    .. code-block:: pycon

        Calculate extinction-coefficient (cm-1 units) over wavelength
        and pressure arrays.
        The nominal density profile of the absorber is assumed as
        density = pressure / (k*temperature)

        temperature: 1D float array
            Temperature profile (K)
        pars: 1D iterable
            If not None, update the model parameters with given values
            and re-calculate the cross sections.
        layer: Integer
            If not None, compute the extinction coefficient only
            at the given index in density array.

        Returns
        -------
        extinction_coefficient: 2D float array
            The Rayleigh extinction coefficient (cm-1 units).


pyratbay.opacity.linelist
_________________________


.. py:module:: pyratbay.opacity.linelist

.. py:class:: Hitran(dbfile, pffile, log)

    .. code-block:: pycon

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

    .. py:method:: get_iso(molname)
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

    .. py:method:: get_iso(molname)
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

    .. py:method:: get_iso(molname)
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

    .. py:method:: get_iso(molname)
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

    .. py:method:: get_iso(molname)
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

    .. py:method:: get_iso(molname)
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

.. py:class:: Vald(dbfile, ion, pffile, log)

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

    .. py:method:: dbread(wn_init, wn_end, verb)
    .. code-block:: pycon

        Read a VALD database.

        Parameters
        ----------
        wn_init: Scalar
            Initial wavenumber limit (in cm-1).
        wn_end: Scalar
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
        iso_id: 2D integer ndarray
          Isotope index.

    .. py:method:: get_iso(molname)
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

.. py:function:: get_tips_molecules()
.. code-block:: pycon

    Get a list of all TIPS molecules.

    Note tha this list does not contain 'O', hence it's not the same
    list as the HITRAN mol ID list: https://hitran.org/docs/molec-meta/

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

.. py:function:: check_exomol_files(files)
.. code-block:: pycon

    Check that all input exomol files are of the same type.
    Check that all refer to a same molecule.
    Collect molecule and isotopes names.

    Parameters
    ----------
    files: List of strings
        A list of Exomol files.

    Returns
    -------
    file_type: String
        Whether all input files are .pf files (return 'pf'),
        all input files are .states or .states.bz2 files (return 'states'),
        or else return ''.
    molecule: String
        Molecule's name.
    isotopes: List of strings
        List of isotope names.

.. py:function:: tips(molecule, isotopes=None, outfile=None, db_type='as_exomol')
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
    for molecule H2O, with isotopes ['116', '118', '117', '126', '128', '127', '226', '228', '227'],
    and temperature range 1--6000 K.

.. py:function:: exomol_pf(files, outfile=None)
.. code-block:: pycon

    Extract ExoMol partition-function values from input files.
    If requested, write the partition-function into a file for use
    with Pyrat Bay.

    Parameters
    ----------
    files: String or List of strings
        Input Exomol ilenames.  Files must either all correspond to .pf
        files or all correspond to .states files.
        For multiple isotopes, all files must correspond to a same molecule.
    outfile: String
         If not None, save output to file.  If outfile == 'default',
         save output to file named as PF_exomol_molecule.dat

    Returns
    -------
    pf: 2D float ndarray
        TIPS partition function for input molecule.
    isotopes: 1D string list
        List of isotopes.
    temps: 1D float ndarray
        Partition-function temperature samples (K).

    Examples
    --------
    >>> import pyratbay.opacity.partitions as pf
    >>>
    >>> # Extract data from Exomol .pf files
    >>> # wget https://www.exomol.com/db/HCN/1H-12C-14N/Harris/1H-12C-14N__Harris.pf
    >>> # wget https://www.exomol.com/db/HCN/1H-13C-14N/Larner/1H-13C-14N__Larner.pf
    >>> files = ['1H-12C-14N__Harris.pf', '1H-13C-14N__Larner.pf']
    >>> pf_data, isotopes, temps = pf.exomol_pf(files)

.. py:function:: exomol_states(files, tmin, tmax, tstep, outfile=None)
.. code-block:: pycon

    Extract ExoMol partition-function values from input files.
    If requested, write the partition-function into a file for use
    with Pyrat Bay.

    Parameters
    ----------
    files: String or List of strings
        Input Exomol ilenames.  Files must either all correspond to .pf
        files or all correspond to .states files.
        For multiple isotopes, all files must correspond to a same molecule.
    tmin: Float
        Mimimum temperature to sample the partitions.
        Required to sample from .state files only.
    tmax: Float
        Maximum temperature to sample the partitions.
        Required to sample from .state files only.
    tstep: Float
        Temperature step at which to sample the temperature array
        Required to sample from .state files only.
    outfile: String
         If not None, save output to file.  If outfile == 'default',
         save output to file named as PF_exomol_molecule.dat

    Returns
    -------
    pf: 2D float ndarray
        TIPS partition function for input molecule.
    isotopes: 1D string list
        List of isotopes.
    temps: 1D float ndarray
        Partition-function temperature samples (K).

    Examples
    --------
    >>> import pyratbay.opacity.partitions as pf
    >>>
    >>> # Extract data from Exomol .states files
    >>> # wget https://www.exomol.com/db/HCN/1H-12C-14N/Harris/1H-12C-14N__Harris.states.bz2
    >>> # wget https://www.exomol.com/db/HCN/1H-13C-14N/Larner/1H-13C-14N__Larner.states.bz2
    >>> files = [
    >>>     '1H-12C-14N__Harris.states.bz2',
    >>>     '1H-13C-14N__Larner.states.bz2',
    >>> ]
    >>> pf, isotopes, temps = pf.exomol_states(
    >>>     files, tmin=5.0, tmax=5000.0, tstep=5.0,
    >>> )

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


pyratbay.opacity.rayleigh
_________________________


.. py:module:: pyratbay.opacity.rayleigh

.. py:class:: Kurucz(wn, species)

    .. code-block:: pycon

        Parameters
        ----------
        wn: 1D float ndarray
           Wavenumber in cm-1.
        species: String
           The species, which can be H, He, H2, or e-.

    .. py:method:: calc_extinction_coefficient(density, layer=None)
    .. code-block:: pycon

        Calculate extinction-coefficient (cm-1 units) over wavelength
        and layers grid.

        Parameters
        ----------
        density: 1D float array
            Number density of this species over an atmospheric profile
            (molecules cm-3 units)
        layer: Integer
            If not None, compute the extinction coefficient only
            at the given index in density array.

        Returns
        -------
        extinction_coefficient: 2D float array
            The Rayleigh extinction for this species for the given
            atmosphere (cm-1 units).

    .. py:method:: set_wn(wn)
    .. code-block:: pycon

        When wn is updated the cross-sections must be re-calculated.


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

.. py:function:: spectrum(spectrum, wavelength, rt_path, data=None, uncert=None, bands_wl0=None, bands_flux=None, bands_response=None, bands_wl=None, label='model', bounds=None, logxticks=None, resolution=150.0, yran=None, filename=None, fignum=501, axis=None, marker='o', ms=5.0, lw=1.25, fs=14, data_front=True, units=None, dpi=300, theme=None, data_color='black')
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
        Observing geometry: transit, eclipse, or emission.
    data: 1D float ndarray
        Observing data points at each bands_wl0.
    uncert: 1D float ndarray
        Uncertainties of the data points.
    bands_wl0: 1D float ndarray
        The mean wavelength for each band/data point.
    bands_flux: 1D float ndarray
        Band-integrated model spectrum at each bandwl.
    bands_response: Iterable of 1D float ndarrays
        Transmission response curve for each band.
    bands_wl: Iterable of 1D float ndarrays.
        The wavelength arrasy for each bands_response curve.
    label: String
        Label for spectrum curve.
    bounds: Tuple
        Tuple with -2, -1, +1, and, +2 sigma boundaries of spectrum.
        If not None, plot shaded area between +/-1sigma and +/-2sigma
        boundaries.
    logxticks: 1D float ndarray
        If not None, switch the X-axis scale from linear to log, and set
        the X-axis ticks at the locations given by logxticks.
    resolution: Float
        Binning resolution to display the spectra.
    yran: 1D float ndarray
        Figure's Y-axis boundaries.
    filename: String
        If not None, save figure to filename.
    fignum: Integer
        Figure number.
    axis: AxesSubplot instance
        The matplotlib Axes of the figure.
    ms: Float
        Marker sizes.
    lw: Float
        Line widths.
    fs: Float
        Font sizes.
    data_front: Bool
        display the data in front of models
    units: String
        Flux units. Select from: 'percent', 'ppt', 'ppm', 'none'.
    dpi: Integer
        The resolution in dots per inch for saved files.
    theme: string or mc3.plots.Theme object
        A color theme for the models.

    Returns
    -------
    ax: AxesSubplot instance
        The matplotlib Axes of the figure.

.. py:function:: contribution(contrib_func, wl, rt_path, pressure, filename=None, filters=None, fignum=-21, dpi=300)
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
        Layer's pressure array (bars).
    filename: String
        Filename of the output figure.
    filters: 1D string ndarray
        Name of the filter bands (optional).
    fignum: Integer
        Figure number.
    dpi: Integer
        The resolution in dots per inch for saved files.

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

.. py:function:: temperature(pressure, profiles=None, labels=None, colors=None, bounds=None, ax=None, filename=None, theme='blue', alpha=[0.75, 0.5], fs=13, lw=2.0, fignum=504, dpi=300)
.. code-block:: pycon

    Plot temperature profiles.

    Parameters
    ----------
    pressure: 1D float ndarray
        The atmospheric pressure profile in bars.
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
    ax: AxesSubplot instance
        If not None, plot into the given axis.
    filename: String
        If not None, save plot to given file name.
    theme: string or mc3.plots.Theme object
        A color theme for the profiles and credible regions.
    alpha: 2-element float iterable
        Alpha transparency for bounds regions.
    fs: Float
        Labels font sizes.
    lw: Float
        Lines width.
    fignum: Integer
        Figure's number (ignored if axis is not None).
    dpi: Integer
        The resolution in dots per inch for saved files.

    Returns
    -------
    ax: AxesSubplot instance
        The matplotlib Axes of the figure.

.. py:function:: abundance(vol_mix_ratios, pressure, species, highlight=None, xlim=None, colors=None, dashes=None, filename=None, lw=2.0, fignum=505, fs=13, legend_fs=None, ax=None, dpi=300)
.. code-block:: pycon

    Plot atmospheric volume-mixing-ratio abundances.

    Parameters
    ----------
    vol_mix_ratios: 2D float ndarray
        Atmospheric volume mixing ratios to plot [nlayers,nspecies].
    pressure: 1D float ndarray
        Atmospheric pressure [nlayers], in bars.
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
    dpi: Integer
        The resolution in dots per inch for saved files.

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
    >>> vmr = pa.chemistry('equilibrium', pressure, temperature, species).vmr
    >>> ax = pp.abundance(
    >>>     vmr, pressure, species, colors='default',
    >>>     highlight='H2O CH4 CO CO2 NH3 HCN H2 H He'.split())

.. py:data:: default_colors
.. code-block:: pycon

  {'H2O': 'blue', 'CO': 'xkcd:green', 'CO2': 'red', 'CH4': 'gold', 'H2': 'indigo', 'SO2': 'deepskyblue', 'HCN': '0.6', 'NH3': 'darkorange', 'N2': 'darkkhaki', 'H': 'magenta', 'TiO': 'black', 'VO': 'peru', 'Na': 'darkviolet', 'K': 'olive', 'C2H2': 'green', 'C2H4': 'pink', 'He': 'dodgerblue', 'H2S': 'cornflowerblue'}

.. py:function:: posteriors(post_file, theme='blue', data_color='black', plot_species=None, vmr_lims=None, logxticks=None, dpi=300)
.. code-block:: pycon

    Plot contribution functions, temperature profiles, VMRs, and spectra
    derived from a retrieval posterior file.

    Parameters
    ----------
    post_file: String
        A posterior pickle file produced by pt.posterior_post_processing()
        containing post-processed medians and quantiles for atmospheric
        values of interest.
    theme: string or mc3.plots.Theme object
        A color theme for the models.
    data_color: a valid matplotlib color
        The color for the spectrum data points.
    plot_species: 1D string iterable
        List of species to plot in VMR figure.
        If None, default to the species in post_data['active_species'],
        which includes the species that actively contribute to the opacity.
    vmr_limits: 2-element float iterable
        Plotting boundaries for the volume mixing ratio.
    logxticks: 1D float ndarray
        If not None, switch the X-axis scale from linear to log, and set
        the X-axis ticks at the locations given by logxticks.
    dpi: Integer
        The resolution in dots per inch for saved files.

    Examples
    --------
    >>> import pyratbay.plots as pp
    >>> post_file = 'ns_emission_tutorial_posteriors_info.pickle'
    >>> theme = 'red'
    >>> pp.posteriors(post_file, theme='red')

    >>> vmr_lims = 1e-5, 1.0
    >>> pp.posteriors(post_file, theme='red', vmr_lims=vmr_lims)

    >>> plot_species = 'H2O CO H2 He H CH4 CO2 C N O'.split()
    >>> pp.posteriors(post_file, theme='red', plot_species=plot_species)


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
    >>> s = np.trapezoid(flux, wn) * (pc.rsun/pc.au)**2
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
    >>> s = np.trapezoid(flux, wn) * (pc.rsun/pc.au)**2
    >>> print("Solar constant [T={:.0f} K, logg={:.1f}]:  S = {:.1f} W m-2".
    >>>       format(ktemp, klogg, s * 1e-3))
    Solar constant [T=5750 K, logg=4.5]:  S = 1340.0 W m-2
    >>> # Pretty close to the solar constant: ~1361 W m-2

    >>> # Read the whole set of models in file:
    >>> # (in this case, ktemp and klogg are 1D arrays)
    >>> fluxes, wn, ktemp, klogg, continua = ps.read_kurucz(kfile)

.. py:class:: PassBand(filter_file, wl=None, wn=None, counting_type='photon')

    .. code-block:: pycon

        Parameters
        ----------
        filter_file: String
            Path to filter file containing wavelength (um) and passband
            response function in two columns.
        wl: 1D float array
            Wavelength array at which evaluate the passband response's
            in micron units.
            (only one of wl or wn should be provided on call)
        wn: 1D float array
            Wavenumber (cm-1) at which the filter is intended to be evalulated.
        counting_type: String
            Detector counting type (i.e., the response function units),
            choose between 'photon' (default) or 'energy'.

        Examples
        --------
        >>> import pyratbay.spectrum as ps
        >>> import pyratbay.constants as pc
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>>
        >>> # Test on a blackbody spectrum
        >>> wl = np.linspace(1.0, 10.0, 10000)
        >>> hj_spectrum = ps.bbflux(1e4/wl, 2100.0)
        >>> # Create a Spitzer/IRAC2 band and integrate spectrum over it
        >>> band = ps.PassBand(f'{pc.FILTERS}spitzer_irac2.dat', wl=wl)
        >>> band_flux = band(hj_spectrum)
        >>>
        >>> plt.figure(0)
        >>> plt.clf()
        >>> plt.plot(wl, hj_spectrum, c='black')
        >>> plt.plot(band.wl0, band_flux, 'o', c='royalblue')
        >>> plt.plot(band.wl, band.response*5e7, c='0.7')
        >>> plt.ylim(bottom=0.0)
        >>>
        >>> # Now, we can re-evaluate for a bunch of spectra
        >>> temps = np.arange(900, 2500.0, 150.0)
        >>> spectra = [ps.bbflux(1e4/wl, temp) for temp in temps]
        >>> fluxes = [band(spectrum) for spectrum in spectra]
        >>>
        >>> plt.figure(1)
        >>> plt.clf()
        >>> for i,temp in enumerate(temps):
        >>>     color = plt.cm.viridis(i/11)
        >>>     plt.plot(wl, spectra[i], c=color)
        >>>     plt.plot(band.wl0, fluxes[i], 'o', c=color)
        >>> plt.plot(band.wl, band.response*2e7, c='0.7', zorder=-1)
        >>> plt.xscale('log')
        >>> plt.ylim(bottom=0.0)
        >>> plt.xlim(1.0, 10.0)

    .. py:method:: integrate(spectrum, wl=None, wn=None)
    .. code-block:: pycon

        Integrate a spectral function over the passband.

        Parameters
        ----------
        spectrum: 1D float array
            Spectral function to be band-integrated. The spectral sampling
            of this must be set at initialization or with the set_spectrum()
            method.  Otherwise, it can be set at run time with the wn or wl
            arguments.
        wn: 1D float array
            (optional) wavenumber array (cm-1) over which spectrum is sampled.
        wl: 1D float array
            (optional) wavelength array (um) over which spectrum is sampled.

        Returns
        -------
        bandflux: Float
            Band-integrated value of the spectral  function.

    .. py:method:: save_filter(save_file)
    .. code-block:: pycon

        Write filter response function data to file, into two columns
        the wavelength (um) and the response function.

        Parameters
        ----------
        save_file: String
            File where to save the filter data.

    .. py:method:: set_sampling(wl=None, wn=None)
    .. code-block:: pycon

        Interpolate filter response function at specified spectral array.
        The response funciton is normalized such that the integral over
        wavenumber equals one.

        Parameters
        ----------
        wl: 1D float array
            Wavelength array at which evaluate the passband response's
            in micron units (only one of wl or wn should be provided on call)
        wn: 1D float array
            Wavenumber array at which evaluate the passband response's
            in cm-1 units (only one of wl or wn should be provided on call)

        Defines
        -------
        self.response  Normalized interpolated response function
        self.idx       IndicesWavenumber indices
        self.wn        Passband's wavenumber array
        self.wl        Passband's wavelength array

        Returns
        -------
        out_wave: 1D float array
            Same as self.wl or self.wn depending on the input argument.
        out_response: 1D float array
            Same as self.response

        Examples
        --------
        >>> # See examples in help(ps.PassBand.__init__)

.. py:class:: Tophat(wl0, half_width, name='tophat', wl=None, wn=None, counting_type='photon', ignore_gaps=False)

    .. code-block:: pycon

        Parameters
        ----------
        wl0: Float
            The passband's central wavelength (um units).
        half_width: Float
            The passband's half-width (um units).
        name: Str
            A user-defined name for the filter when calling str(self),
            e.g., to identify the instrument provenance of this filter.
        wl: 1D float array
            Wavelength array at which evaluate the passband response's
            in micron units (only one of wl or wn should be provided on call)
        wn: 1D float array
            Wavenumber (cm-1) at which the filter is intended to be evalulated.
        counting_type: String
            Detector counting type (i.e., the response function units),
            choose between 'photon' (default) or 'energy'.
        ignore_gaps: Bool
            If True and there are no points inside the band,
            set the idx, wn, wl, and response variables to None.
            A ps.bin_spectrum() call on such band will return np.nan.
            Otherwise the code will throw a ValueError.

        Examples
        --------
        >>> import pyratbay.spectrum as ps
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>>
        >>> # Test on a blackbody spectrum
        >>> wl = np.linspace(1.0, 10.0, 10000)
        >>> hj_spectrum = ps.bbflux(1e4/wl, 2100.0)
        >>> # Create a top-hat passband and integrate spectrum
        >>> band = ps.Tophat(wl0=4.5, half_width=0.1, wl=wl)
        >>> band_flux = band(hj_spectrum)
        >>>
        >>> plt.figure(0)
        >>> plt.clf()
        >>> plt.plot(wl, hj_spectrum, c='black')
        >>> plt.plot(band.wl0, band_flux, 'o', c='royalblue')
        >>> plt.plot(band.wl, band.response*5e6, c='0.7')
        >>> plt.ylim(bottom=0.0)
        >>>
        >>> # Now, we can re-evaluate for a bunch of spectra
        >>> temps = np.arange(900, 2500.0, 150.0)
        >>> spectra = [ps.bbflux(1e4/wl, temp) for temp in temps]
        >>> fluxes = [band(spectrum) for spectrum in spectra]
        >>>
        >>> plt.figure(1)
        >>> plt.clf()
        >>> for i,temp in enumerate(temps):
        >>>     color = plt.cm.viridis(i/11)
        >>>     plt.plot(wl, spectra[i], c=color)
        >>>     plt.plot(band.wl0, fluxes[i], 'o', c=color)
        >>> plt.plot(band.wl, band.response*2e6, c='0.7', zorder=-1)
        >>> plt.xscale('log')
        >>> plt.ylim(bottom=0.0)
        >>> plt.xlim(1.0, 10.0)

    .. py:method:: integrate(spectrum, wl=None, wn=None)
    .. code-block:: pycon

        Integrate a spectral function over the passband.

        Parameters
        ----------
        spectrum: 1D float array
            Spectral function to be band-integrated. The spectral sampling
            of this must be set at initialization or with the set_spectrum()
            method.  Otherwise, it can be set at run time with the wn or wl
            arguments.
        wn: 1D float array
            (optional) wavenumber array (cm-1) over which spectrum is sampled.
        wl: 1D float array
            (optional) wavelength array (um) over which spectrum is sampled.

        Returns
        -------
        bandflux: Float
            Band-integrated value of the spectral  function.

    .. py:method:: save_filter(save_file)
    .. code-block:: pycon

        Write filter response function data to file, into two columns
        the wavelength (um) and the response function.

        Parameters
        ----------
        save_file: String
            File where to save the filter data.

    .. py:method:: set_sampling(wl=None, wn=None)
    .. code-block:: pycon

        Interpolate filter response function at specified spectral array.
        The response function is normalized such that the integral over
        wavenumber equals one.

        Parameters
        ----------
        wl: 1D float array
            Wavelength array at which evaluate the passband response's
            in um units (only one of wl or wn should be provided on call)
        wn: 1D float array
            Wavenumber array at which evaluate the passband response's
            in cm-1 units (only one of wl or wn should be provided on call)

        Defines
        -------
        self.idx  Wavenumber indices
        self.wn  Passband's wavenumber array
        self.wl  Passband's wavelength array
        self.response  Normalized interpolated response function

        Returns
        -------
        out_wave: 1D float array
            Same as self.wl or self.wn depending on the input argument.
        out_response: 1D float array
            Same as self.response

        Examples
        --------
        >>> # See examples in help(ps.Tophat.__init__)

.. py:function:: constant_resolution_spectrum(wave_min, wave_max, resolution)
.. code-block:: pycon

    Compute a constant resolving-power sampling array.

    Parameters
    ----------
    wave_min: Float
        Lower spectral boundary.  This could be either a wavelength
        or a wavenumber. This is agnositc of units.
    wave_max: Float
        Upper spectral boundary.  This could be either a wavelength
        or a wavenumber. This is agnositc of units.
    resolution: Float
        The sampling resolving power: R = wave / delta_wave.

    Returns
    -------
    wave: 1D float array
        A spectrum array with the given resolving power.

    Examples
    --------
    >>> import numpy as np
    >>> import pyratbay.spectrum as ps

    >>> # A low-resolution wavelength sampling:
    >>> wl_min = 0.5
    >>> wl_max = 4.0
    >>> resolution = 5.5
    >>> wl = ps.constant_resolution_spectrum(wl_min, wl_max, resolution)
    >>> print(wl)
    [0.5        0.6        0.72       0.864      1.0368     1.24416
     1.492992   1.7915904  2.14990848 2.57989018 3.09586821 3.71504185]
    >>> # The actual resolution matches the input:
    >>> wl_mean = 0.5*(wl[1:]+wl[:-1])
    >>> print(wl_mean/np.ediff1d(wl))
    [5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5]

.. py:function:: bin_spectrum(bin_wl, wl, spectrum, half_widths=None, gaps=None)
.. code-block:: pycon

    Bin down a spectrum.

    Parameters
    ----------
    bin_wl: 1D float array
        Central wavelength (um) of the desired binned down spectra.
    wl: 1D float array
        Wavelength samples of the original spectrum.
    spectrum: 1D float array
        Spectral values to be binned down.
    half_widths: 1D float array
        The bin half widths (um).
        If None, assume that the bin edges are at the mid-points
        of the bin_wl array.
    gaps: String
        If None (default) and there are bins that do not cover any value,
        (e.g., when the resolution of wl is similar to bin_wl's), raise error.
        If gaps=='ignore', patch those bins with np.nan values.
        If gaps=='interpolate', patch those bins linearly interpolating.
        Use with care.

    Returns
    -------
    bin_spectrum: 1D float array
        The binned spectrum.

    Notes
    -----
    Probably bad things will happen if bin_wl has a similar
    or coarser resolution than wl.

    Examples
    --------
    >>> import pyratbay.spectrum as ps
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt

    >>> # Make a noisy high-resolution signal
    >>> wl = ps.constant_resolution_spectrum(1.0, 3.0, resolution=5000)
    >>> spectrum = np.sin(3.14*wl) + np.random.normal(1.5, 0.1, len(wl))
    >>> # Bin it down:
    >>> bin_wl = ps.constant_resolution_spectrum(1.0, 3.0, resolution=125)
    >>> bin_spectrum = ps.bin_spectrum(bin_wl, wl, spectrum)

    >>> # Compare original and binned signals
    >>> plt.figure(0)
    >>> plt.clf()
    >>> plt.plot(wl, spectrum, '.', ms=2, color='gray')
    >>> plt.plot(bin_wl, bin_spectrum, color='red')

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
    >>> wn1, irac1 = io.read_spectrum(pc.FILTERS+'spitzer_irac1.dat')
    >>> wn2, irac2 = io.read_spectrum(pc.FILTERS+'spitzer_irac2.dat')
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

.. py:function:: wn_mask(wn, wn_min, wn_max, tol=1e-08)
.. code-block:: pycon

    Get mask of wn values withing given ranges with a extra tolerance
    to account for floating-point (in)precision.

    Parameters
    ----------
    wn: 1D float array
        Wavenumber array to mask.
    wn_min: float
        Minumum wavenumber in mask.
    wn_max: float
        Maximum wavenumber in mask.
    tol: float
        Tolerance factor at mask edges, calculated as delta_wn*tol,
        where delta_wn is the sampling stepsize at the edges.

    Returns
    -------
    wn_mask: 1D bool array
        Mask of wavenumber values within ranges.

.. py:function:: inst_convolution(wl, spectrum, resolution, sampling_res=None)
.. code-block:: pycon

    Convolve a spectrum according to an instrumental resolving power

    Parameters
    ----------
    wl: 1D float array
        Spectral array (can be either wavelength or wavenumber).
    spectrum: 1D float array
        Full-resolution spectrum to be convolved.
    resolution: float
        Instrumental resolving power R = lambda/delta_lambda
        Where delta_lambda is the FHWM of the gaussian to be applied.
    sampling_red: float
        Sampling resolution of the input wl spectrum.

    Examples
    --------
    >>> import pyratbay.spectrum as ps
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt

    >>> wl_min = 1.499
    >>> wl_max = 1.501
    >>> samp_resolution = 100_000
    >>> wl = ps.constant_resolution_spectrum(wl_min, wl_max, samp_resolution)
    >>> nwave = len(wl)

    >>> # A delta at wl0 ~ 1.5
    >>> spectrum = np.zeros(nwave)
    >>> spectrum[np.where(wl>1.5)[0][0]] = 1.0
    >>> wl0 = wl[np.where(wl>1.5)[0][0]]

    >>> resolution = 5_000
    >>> conv = ps.inst_convolution(wl, spectrum, resolution)

    >>> # Plot convolved line and expected FWHM
    >>> half_max = np.amax(conv)/2
    >>> hwhm = 0.5 * wl0 / resolution
    >>> plt.figure(0)
    >>> plt.clf()
    >>> plt.plot(wl, conv, color='salmon', lw=2)
    >>> plt.plot([wl0-hwhm, wl0+hwhm], [half_max,half_max], color='xkcd:blue', lw=2)

.. py:function:: rv_shift(vel_km, wn=None, wl=None)
.. code-block:: pycon

    Apply a radial velocity Doppler shift to a 1D wavelength array.

    Parameters
    ----------
    vel_km: Float
        Radial velocity in km/s.
    wn: 1D float array
        Wavenumber array.
    wl: 1D float array
        Wavelength array.

    Returns
    -------
    wave: 1D float array
        Doppler-shifted wavenumber of wavelebngth array

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

.. py:function:: band_cf(cf, bands_response, wn, bands_idx)
.. code-block:: pycon

    Compute band-averaged contribution functions or transmittances.

    Parameters
    ----------
    cf: 2D float ndarray
        The contribution function or transmittance of
        shape [nlayers, nwave].
    bands_response: List of 1D ndarrays
        List of band transmission response curves.
    wn: 1D float ndarray
        The wavenumber sampling (in cm-1).
    bands_idx: List of 1D ndarrays
        List of wavenumber-indices in wn sampled by bands_response.

    Returns
    -------
    bands_cf: 2D float ndarray
        The band-integrated contribution functions of
        shape [nlayers, nbands].

.. py:function:: convective_flux(pressure, temperature, cp, gravity, mu, rho, alpha=1.5, beta=0.5)
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

.. py:function:: transmission(depth, radius, rstar, ideep=None, atm_itop=0, deck_rsurf=None, deck_itop=None)
.. code-block:: pycon

    Compute a transmission spectrum for transit geometry

    Parameters
    ----------
    depth: 2D float array
        Optical depth at each layer and wavelength channel [nlayers,nwave].
        Atmospheric layers are sorted from top to bottom (i.e.,
        depth[0] is the top-most layer)
    radius: 1D float array
        Radius profile of atmospheric layers (cm)
    rstar: Float
        Stellar radius (cm).
    ideep: 1D integer array
        Index of the 'bottom' of the atmosphere at each wavelength,
        from which to start the integration (e.g., at optically thick regime)
    atm_itop: Integer
        Index of the top of the atmosphere (to integrate only up to the layer)
    cloud_rsurf: Float
        If not None, the radius (cm) of an opaque cloud deck, which
        becomes the bottom integration boundary.
    cloud_itop: Integer
        If not None, index of the atmosphere layer right below the
        opaque cloud deck.

    Returns
    -------
    spectrum: 1D float array
        The transmission spectrum.

.. py:function:: plane_parallel_rt(depth, blackbody, wn, quadrature_mu, quadrature_weights=None, ideep=None, atm_itop=0, cloud_tsurf=None, cloud_itop=None)
.. code-block:: pycon

    Compute the intensity (and optionally flux) spectra under
    plane-parallel geometry for a range of slant angles

    Parameters
    ----------
    depth: 2D float array
        Optical depth at each layer and wavelength channel [nlayers,nwave].
        Atmospheric layers are sorted from top to bottom (i.e.,
        depth[0] is the top-most layer)
    blackbody: 2D float array
        Plank fuction evaluated at each layer and wavelength channel
        (same shape as depth).
    wn: 1D float array
        Wavenumber array (cm-1).
    quadrature_mu: 1D float array
        Cosine of the slant-path angles repect to the normal.
    quadrature_weights: 1D float array
        Weights for each quadrature_mu.  If given, compute and return the
        Gaussian-quadrature integral of the intensity over mu (i.e., flux)
    ideep: 1D integer array
        Index of the 'bottom' of the atmosphere at each wavelength,
        from which to start the integration.
    atm_itop: Integer
        Index of the top of the atmosphere (to integrate only up to the layer)
    cloud_tsurf: Float
        If not None, the temperature at cloud_itop, an opaque cloud
        deck that becomes the bottom integration boundary.
    cloud_itop: Integer
        If not None, index of the atmosphere layer right below an
        opaque cloud deck.

    Returns
    -------
    intensity: 2D float array
        The intensity spectra at each angle mu (erg s-1 cm-2 cm sr-1).
    flux: 1D float array [optional]
        The flux spectrum of the integrated intensities (erg s-1 cm-2 cm).
        (typically, a day-side integrated via the Gaussian-quadrature)

.. py:function:: radiative_equilibrium(pressure, radeq_temps, nsamples, chem_model, two_stream_rt, wavenumber, spec, atm, convection=False, tmin=0.0, tmax=6000.0)
.. code-block:: pycon

    Compute radiative-thermochemical equilibrium atmosphere.
    Currently there is no convergence criteria implemented,
    some 100--300 iterations are typically sufficient to converge
    to a stable temperature-profile solution.

    Parameters
    ----------
    pressure: 1D float array
        Pressure profile in barye units.
    nsamples: Integer
        Number of radiative-equilibrium iterations to run.
    convection: Bool
        If True, include convective-flux transport in the radiative
        equilibrium calculation.

    Returns
    -------
    radeq_temps: 2D float array
        Temperature profiles at each iteration.


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
       The pressure profile (in bars).

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> import pyratbay.constants as pc

    >>> nlayers = 9
    >>> # These are all equivalent:
    >>> p1 = pa.pressure(ptop=1e-6, pbottom=1e2, nlayers=nlayers)
    >>> p2 = pa.pressure(1e-6, 1e2, nlayers, units='bar')
    >>> p3 = pa.pressure('1e-6 bar', '1e2 bar', nlayers)
    >>> p4 = pa.pressure(1e-6*pc.bar, 1e2*pc.bar, nlayers, units='barye')
    >>> print(p1)
    [1.e-06 1.e-05 1.e-04 1.e-03 1.e-02 1.e-01 1.e+00 1.e+01 1.e+02]

.. py:function:: temperature(tmodel, pressure=None, nlayers=None, log=None, params=None)
.. code-block:: pycon

    Temperature profile wrapper.

    Parameters
    ----------
    tmodel: String
        Name of the temperature model.
    pressure: 1D float ndarray
        Atmospheric pressure profile in bars.
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
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, units="bar")
    >>> # Isothermal profile:
    >>> temp_iso = pa.temperature("isothermal", pressure, params=1500.0)
    >>> print(temp_iso)
    [1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500. 1500.]

    >>> # Guillot (2010) temperature profile:
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
        Name of chemistry model, select from: 'free' or 'equilibrium'
    pressure: 1D float ndarray
        Atmospheric pressure profile (bars).
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
    chem_network: Callable
        The atmospheric chemistry-network model.
    species: 1D string array
        Output species, 'tea' species might differ from input species
        if there is no thermochemical data for a given one (the species
        is then removed)
    vmr: 2D float array
        Output volume mixing ratios of shape [nlayers, nspecies]

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
    >>> chem_model = 'equilibrium'
    >>> network, out_species, vmr_tea = pa.chemistry(chem_model, pressure, temperature, species)

    >>> abundances = np.array([
    >>>     5e-4, 3e-5, 2e-4, 1e-8,  1e-6, 1e-14, 1e-13, 5e-10,
    >>>     1e-4, 0.85, 5e-3, 0.14,  3e-23, 1e-23])
    >>> chem_model = 'free'
    >>> network, out_species, vmr_uni = pa.chemistry(
    >>>     chem_model, pressure, temperature, species, q_uniform=abundances)

    >>> # Plot the results:
    >>> ax1 = pp.abundance(
    >>>     vmr_tea, pressure, species,
    >>>     colors='default', xlim=[1e-30, 3.0])

    >>> ax2 = pp.abundance(
    >>>     vmr_uni, pressure, species,
    >>>     colors='default', xlim=[1e-30, 3.0], fignum=506)

.. py:function:: uniform(abundances, nlayers)
.. code-block:: pycon

    Generate an atmospheric file with uniform abundances.
    Save it into atmfile.

    Parameters
    ----------
    abundances: 1D float ndarray
        The species mole mixing ratio.
    nlayers: Integer
       Number of pressure layers.

    Returns
    -------
    vmr: 2D Float ndarray
        Abundance profiles of shape [nlayers,nspecies]

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> nlayers = 11
    >>> species = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    >>> abundances = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    >>> vmr = pa.uniform(abundances, nlayers)
    >>> print(vmr)
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

.. py:function:: hydro_g(pressure, temperature, mu, g, p0=None, r0=None)
.. code-block:: pycon

    Calculate radii using the hydrostatic-equilibrium equation considering
    a constant gravity.

    Parameters
    ----------
    pressure: 1D float ndarray
        Atmospheric pressure for each layer (in bar).
    temperature: 1D float ndarray
        Atmospheric temperature for each layer (in K).
    mu: 1D float ndarray
        Mean molecular mass for each layer (in g mol-1).
    g: Float
        Atmospheric gravity (in cm s-2).
    p0: Float
        Reference pressure level (in bar) where radius(p0) = r0.
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
    >>> temperature = pa.tmodels.Isothermal(pressure)(1500.0)
    >>> mu = np.tile(2.3, nlayers)
    >>> g = pc.G * pc.mjup / pc.rjup**2
    >>> r0 = 1.0 * pc.rjup
    >>> p0 = 1.0  # bar
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
        Atmospheric pressure for each layer (in bar).
    temperature: 1D float ndarray
        Atmospheric temperature for each layer (in K).
    mu: 1D float ndarray
        Mean molecular mass for each layer (in g mol-1).
    mass: Float
        Object's mass (in g).
    p0: Float
        Reference pressure level (in bar) where radius(p0) = r0.
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
    >>> temperature = pa.tmodels.Isothermal(pressure)(1500.0)
    >>> mu = np.tile(2.3, nlayers)
    >>> mplanet = 1.0 * pc.mjup
    >>> r0 = 1.0 * pc.rjup
    >>> p0 = 1.0  # bar
    >>> # Radius profile in Jupiter radii:
    >>> radius = pa.hydro_m(pressure, temperature, mu, mplanet, p0, r0)/pc.rjup
    >>> print(radius)
    [1.05973436 1.05188019 1.04414158 1.036516   1.029001   1.02159419
     1.01429324 1.00709591 1.         0.99300339 0.986104  ]

.. py:function:: hill_radius(smaxis, mplanet, mstar)
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

    Use the ideal gas law to calculate number density in molecules cm-3.

    Parameters
    ----------
    abundances: 2D float array
        Species volume mixing fraction.
        Can have a 2D shape of [nlayers,nmol] or a 1D shape of [nlayers]
    pressure: 1D array
        Atmospheric pressure profile (in bars).
    temperature: 1D array
        Atmospheric temperature profile (in kelvin).

    Returns
    -------
    density: 2D float array
        Atmospheric density in molecules cm-3.

    Examples
    --------
    >>> import pyratbay.atmosphere as pa
    >>> nlayers = 11
    >>> pressure = pa.pressure(1e-8, 1e2, nlayers, units='bar')
    >>> temperature = np.tile(1500.0, nlayers)
    >>> species = ["H2", "He", "H2O", "CO", "CO2", "CH4"]
    >>> abundances = [0.8496, 0.15, 1e-4, 1e-4, 1e-8, 1e-4]
    >>> vmr = pa.uniform(abundances, nlayers)
    >>> dens = pa.ideal_gas_density(vmr, pressure, temperature)
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

.. py:function:: temperature_posterior(posterior, tmodel)
.. code-block:: pycon

    Compute the median and inter-quantiles regions (68% and 95%)
    of a temperature profile posterior.

    Parameters
    ----------
    posterior: 2D float ndarray
        A posterior distribution for the parameters of tmodel
        of shape [nsamples, npars].
    tmodel: Callable
        Temperature-profile model.

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
    >>> import pyratbay.plots as pp
    >>> import numpy as np

    >>> # Non-inverted temperature profile:
    >>> pressure = pa.pressure('1e-6 bar', '1e2 bar', nlayers=100)
    >>> tmodel = pa.tmodels.Guillot(pressure)

    >>> # Simulate posterior where log_kappa' and log_gamma1 vary:
    >>> nsamples = 5000
    >>> posterior = np.tile([-4.0, -1.0, 0.0, 0.0, 1000.0, 0.0], (nsamples,1))
    >>> posterior[:,0] = np.random.normal(-4.0, 0.5, nsamples)
    >>> posterior[:,1] = np.random.normal(-1.0, 0.5, nsamples)

    >>> # Compute posterior and take a look at it:
    >>> tpost = pa.temperature_posterior(posterior, tmodel)
    >>> ax = pp.temperature(pressure, profiles=tpost[0], bounds=tpost[1:])

.. py:function:: qcapcheck(vmr, qcap, ibulk)
.. code-block:: pycon

    Check if the cummulative abundance of traces exceeds qcap.

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

.. py:function:: vmr_scale(vmr, species, vmr_models, vmr_pars, bulk, qsat=None, iscale=None, ibulk=None, bratio=None, invsrat=None)
.. code-block:: pycon

    Scale specified species abundances and balance bulk abundances to
    conserve sum(vmr)=1 in each layer.

    Parameters
    ----------
    vmr: 2D float ndarray
        Volume mixing ratio of atmospheric species [nlayers, nspecies].
    species: 1D string ndarray
        Names of the species in the atmosphere.
    vmr_models: iterable of pyratbay.atmosphere.vmr_models instances
        List of VMR models.  It can also be an individual model.
    vmr_pars: 1D float ndarray
        List of parameters for each model in vmr_models.
    bulk: 1D string ndarray
        Names of the bulk (dominant) species.
    qsat: Float
        Maximum allowed combined abundance for trace species.
    iscale: 1D integer ndarray
        Indices of mol_model species in vmr.
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

    >>> nlayers = 51
    >>> pressure = pa.pressure('1e-7 bar', '100 bar', nlayers)
    >>> vmr = np.tile([0.85, 0.15, 1e-4, 1e-4, 1e-4, 1e-4], (nlayers,1))
    >>> species = ['H2', 'He' ,'H2O', 'CH4', 'CO', 'CO2']
    >>> vmr_models = [
    >>>     pa.vmr_models.IsoVMR('H2O', pressure),
    >>>     pa.vmr_models.IsoVMR('CO', pressure),
    >>> ]
    >>> # VMR with updated H2O and CO abundances:
    >>> vmr_pars = [-3.5, -3.3]
    >>> bulk = ['H2', 'He']
    >>> scaled_vmr = pa.vmr_scale(vmr, species, vmr_models, vmr_pars, bulk)

    >>> # Show abundances at a layer:
    >>> for i,mol in enumerate(species):
    >>>     print(f'VMR_{mol:4s} = {scaled_vmr[0,i]:.5f}')
    VMR_H2   = 0.84914
    VMR_He   = 0.14985
    VMR_H2O  = 0.00032
    VMR_CH4  = 0.00010
    VMR_CO   = 0.00050
    VMR_CO2  = 0.00010


pyratbay.atmosphere.tmodels
___________________________


.. py:module:: pyratbay.atmosphere.tmodels

.. py:class:: Isothermal(pressure)

    .. code-block:: pycon

        Parameters
        ----------
        pressure: 1D float iterable
            Pressure array (bar) where to evaluate the temperature profile.

.. py:class:: Guillot(pressure, gravity=None)

    .. code-block:: pycon

        Parameters
        ----------
        pressure: 1D float ndarray
            Atmospheric pressure profile (bar).
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
        be solved without knowing the temperature a priori, thus making
        this a circular problem (shrug emoji).

.. py:class:: Madhu(pressure)

    .. code-block:: pycon

        Parameters
        ----------
        pressure: 1D float ndarray
            Pressure array in bar.

.. py:function:: get_model(name, *args, **kwargs)
.. code-block:: pycon

    Get a temperature-profile model by its name.


pyratbay.atmosphere.vmr_models
______________________________


.. py:module:: pyratbay.atmosphere.vmr_models

.. py:function:: hybrid_vmr(model, val, net)
.. code-block:: pycon

    Set a free VMR for a molecule on top of a thermochemical-equilibrium VMR.
    Care not to exceed the number of available elements.

.. py:class:: MetalEquil()

    .. code-block:: pycon

        Initialize self.  See help(type(self)) for accurate signature.

.. py:class:: ScaleEquil(name)

    .. code-block:: pycon

        Parameters
        ----------
        name: String
            The element's name, the format must be: [X/H]
            where 'X' is the element name (e.g.: [C/H], [Na/H], ...).

.. py:class:: RatioEquil(name)

    .. code-block:: pycon

        Parameters
        ----------
        name: String
            The elements name ratio, the format must be: X/Y
            where 'X' and 'Y' are the element names (e.g.: C/O, N/O, ...).

.. py:class:: IsoVMR(species, pressure)

    .. code-block:: pycon

        Parameters
        ----------
        species: String
            The atmospheric species name for this VMR profile model.
        pressure: 1D float iterable
            Pressure array (bar) where to evaluate the temperature profile.

.. py:class:: ScaleVMR(species, pressure, vmr0)

    .. code-block:: pycon

        Parameters
        ----------
        species: String
            The atmospheric species name for this VMR profile model.
        pressure: 1D float iterable
            Pressure array (bar) where to evaluate the temperature profile.
        vmr0: 1D float array

.. py:class:: SlantVMR(species, pressure)

    .. code-block:: pycon

        Parameters
        ----------
        species: String
            The atmospheric species name for this VMR profile model.
        pressure: 1D float iterable
            Pressure array (bar) where to evaluate the temperature profile.


pyratbay.pyrat.spectrum
_______________________


.. py:module:: pyratbay.pyrat.spectrum

.. py:class:: Spectrum(inputs, log)

    .. code-block:: pycon

    .. code-block:: pycon

        Make the wavenumber sample from user inputs.

.. py:function:: spectrum(pyrat)
.. code-block:: pycon

    Spectrum calculation driver.

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

