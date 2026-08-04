.. include:: ../../_substitutions.rst

Eclipse retrieval: WASP-18b SOSS
================================

This tutorial shows how perform an atmospheric retrieval of the
secondary-eclipse spectra of WASP-18b, constrained by the JWST/SOSS
observations.  We will replicate the analysis presented in
[Deline2025]_, that is, a retrieval of:

- the combined JWST, Spitzer, TESS, and CHEOPS observations
- 1D atmosphere assuming thermochemical-equilibrium VMRs
- a Madhu temperature profile, *T(p)*
- a Rayleigh opacity model
- a gray, patchy cloud-deck model


We can break the analysis into the following steps:

1. :ref:`wasp18b_obs`

   - :ref:`wasp18b_obs_file`
   - :ref:`wasp18b_sed`
   - :ref:`wasp18b_cross_sec`
   - :ref:`wasp18b_config`
2. :ref:`wasp18b_run`
3. :ref:`wasp18b_stats`

----------------------------------------------------------------------

.. _wasp18b_obs:

Setup
-----

For the setup we will need three ingredients:

#. A **configuration file** to define the system parameters, atmospheric model, posterior sampling, etc.

#. An **observation file** defining the data points: depths, uncertainties, and bin wavelengths

#. **Cross-section files** for the atmospheric species


Lets start with the required input files, and then go over the
configuration file.


.. _wasp18b_obs_file:

Observation file
~~~~~~~~~~~~~~~~

``Pyrat Bay`` observation files tell the code the data that is being
fit.  These are a plain text files containing the transit depth,
uncertainty, and the band (which can be either a tophat or a broadband
passband).  There is always one data point per row.

Here we will use WASP-18b eclipse observations by CHEOPS, TESS, JWST,
and Spitzer from [Coulombe2023]_ and [Deline2025]_.

.. tab-set::

  .. tab-item:: Download observation file
     :selected:

     Here's the WASP-18b observation files ready to use. Click the link
     to see/download the entire file.  Some notes about observations
     files:

     - Lines starting with ``#`` are comments.

     - An optional line with ``@DEPTH_UNITS`` (1) indicates that the
       file contains depths and uncertainties, and (2) sets the units
       for the depths in the following line.

     - A mandatory line with ``@DATA`` indicates where the data starts

     - JWST data is modeled as tophat bands defined by the central
       wavelength, the bin half-width, and (optionally) a name for the
       instrument

     - CHEOPS, TESS, and Spitzer data come from broad-band photometry.
       Their response functions are specify as paths to plain files
       that tabulate the wavelength and response.  CHEOPS, TESS,
       Kepler, and Spitzer passbands are provided in ``Pyrat Bay``,
       and can be accessed via the ``{FILTERS}`` flag. for other bands
       type the full path to the files.

     .. literalinclude:: ../../_static/data/obs_wasp18b_eclipse_all.dat
        :caption: Extract from file: `obs_wasp18b_eclipse_all.dat <../../_static/data/obs_wasp18b_eclipse_all.dat>`__
        :language: ini
        :lines: 1-26


     .. raw:: html

        <details>
        <summary>Click here to show/hide: <a href="../../_static/data/obs_wasp18b_eclipse_jwst.dat">obs_wasp18b_eclipse_jwst.dat</a></summary>

     .. literalinclude:: ../../_static/data/obs_wasp18b_eclipse_jwst.dat
         :caption: File:  `obs_wasp18b_eclipse_all.dat <../../_static/data/obs_wasp18b_eclipse_jwst.dat>`__
         :language: ini

     .. raw:: html

        </details>


  .. tab-item:: Compute observation file

     We will constrain this retrieval to the JWST, Spitzer, CHEOPS,
     and TESS eclipse observations. So we need to collect that
     data. For the JWST spectroscopic observations we will use the
     NAMELESS spectral reduction (available on Zenodo
     https://zenodo.org/records/7907569), which we will model as a
     series of top-hat narrow passbands.

     The CHEOPS, TESS, and Spitzer observations [Deline2025]_ consist
     of broad photometric passbands. For these we will use
     passband filter files  (that are included in ``Pyrat Bay``).

     Simulations with ``Pyrat Bay`` load this information from an
     observation file input. The script below creates observation
     files for (1) the JWST observations, and (2) all photometric and
     spectroscopic observations combined.

     .. code-block:: python

         # Save JWST data
         import numpy as np
         import pyratbay.io as io

         jwst_data = np.loadtxt('NAMELESS_W18b_spectrum.txt', unpack=True)
         jwst_wl, jwst_depths, jwst_uncerts, jwst_half_widths = jwst_data
         njwst = len(jwst_wl)

         # Save JWST data:
         obs_file = 'obs_wasp18b_eclipse_jwst.dat'
         jwst_inst_names = ['NIRISS' for _ in jwst_wl]
         io.write_observations(
             obs_file,
             jwst_inst_names,
             jwst_wl, jwst_half_widths,
             jwst_depths, jwst_uncerts, depth_units='ppm',
         )


         # Save CHEOPS + TESS + JWST + Spitzer data:
         deline2025_depths = [211.9, 340.4, 3098, 3925, 4080, 4350]
         deline2025_uncerts = [ 7.6,   7.0,  112,   25,  230,  205]
         n_photo = len(deline2025_depths)

         depths = np.concatenate([deline2025_depths, jwst_depths])
         depth_uncerts = np.concatenate([deline2025_uncerts, jwst_uncerts])
         # zero wavelength => inst_name is path to passband file (photometry)
         wl = np.concatenate([np.zeros(n_photo), jwst_wl])
         half_widths = np.concatenate([np.zeros(n_photo), jwst_half_widths])
         photo_names = [
             '{FILTERS}cheops.dat',
             '{FILTERS}tess.dat',
             '{FILTERS}spitzer_irac1.dat',
             '{FILTERS}spitzer_irac2.dat',
             '{FILTERS}spitzer_irac3.dat',
             '{FILTERS}spitzer_irac4.dat',
         ]
         inst_names = photo_names + jwst_inst_names

         obs_file = 'obs_wasp18b_eclipse_all.dat'
         io.write_observations(
             obs_file,
             inst_names,
             wl, half_widths, depths, depth_uncerts,
             depth_units='ppm',
         )


.. _wasp18b_cross_sec:

Cross sections
~~~~~~~~~~~~~~

Following the analysis of [Deline2025]_, we will include
line-sampled cross sections for these molecules: |H2O|, CO, |CO2|,
|CH4|, TiO, VO, HCN, |NH3|, and |C2H2|.  Here we will work with
the latests opacity sources for these species from ExoMol and HITEMP.

The current recommendation for sampled cross sections for JWST
retrievals is to adopt a resolution :math:`R>20.000`.  So, here we
will use a cross section grid at :math:`R=25.000`, sampling from
:math:`0.35-10.5` μm in wavelength (to cover the spectral range of the
data), :math:`500-4000` K in temperature, and from
:math:`100-1.0^{-9}` bar in pressure.

Now, be aware that cross sections always have assumptions baked into
them.  Below you can find ready-to-use cross sections and their
assumptions.

.. Alternatively, if you need to adjust the ranges or sampling resolution for a specific project, compute your own cross sections starting from the line-list files (where you can customize at will).

.. tab-set::

  .. tab-item:: Download cross sections
     :selected:

     The Zenodo repository `doi.org/10.5281/zenodo.16965391
     <https://zenodo.org/records/16965391>`__ contains the
     cross-section files that we will use for this JWST atmospheric
     retrieval. See the list below for direct links to the files for
     each molecule.

     These cross sections have been computed assuming an
     |H2|/He-dominated atmosphere, and terrestrial isotopic
     ratios. The lines have Voigt profiles with a wing cut-off at 300
     HWHM and at 25 |kayser|.  The grids sampling are:

     - Wavelength: :math:`0.15-33` μm, at a constant resolution of :math:`R=25.000`
     - Temperature: :math:`200-5000` K, with :math:`\Delta T = 150` K
     - Pressure: :math:`1.0^{-9}-1.0^{3}` bar, equally sampled in log(`p`) with 4 samples per dex.

     .. list-table:: Tabulated cross section files
       :header-rows: 1

       * - Species (source)
         - References
       * - `H2O <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_H2O_exomol_pokazatel.npz>`__ (exomol, pokazatel)
         -  [Polyansky2018]_
       * - `CO <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_CO_hitemp_2019.npz>`__ (HITEMP, li)
         - [Li2015]_
       * - `CO2 <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_CO2_ames_ai3000k.npz>`__ (ames, ai3000k)
         - [Huang2023]_
       * - `CH4 <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_CH4_exomol_mm.npz>`__ (exomol, mm)
         - [Yurchenko2024a]_

       * - `TiO <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_TiO_exomol_toto.npz>`__ (exomol, toto)
         - [McKemmish2019]_
       * - `VO <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_VO_exomol_hyvo.npz>`__ (exomol, hyvo)
         - [Bowesman2024]_

       * - `HCN <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_HCN_exomol_harris_larner.npz>`__ (exomol, harris larner)
         - [Harris2008]_ [Barber2014]_
       * - `NH3 <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_NH3_exomol_coyute.npz>`__ (exomol, coyute)
         - [Coles2019]_ [Yurchenko2024b]_
       * - `C2H2 <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_C2H2_exomol_acety.npz>`__ (exomol, acety)
         - [Chubb2020]_

     .. Note:: If you want to see the source script or need to
         customize the cross sections (e.g., broader temperature
         ranges, finer resolution, different line profiles), follow
         the steps in the `'Compute cross sections'` tab.

  .. tab-item:: Compute cross sections

     TBD

.. _wasp18b_sed:


Stellar SED spectrum
~~~~~~~~~~~~~~~~~~~~

The stellar SED is required to calculate the planet-to-star flux
ratio. Here we will use a stellar SED for WASP-18 from the PHOENIX
models [Husser2013]_, which we will get using the ``Gen TSO``
package. If you haven’t already, install this package with this
shell command:

.. code-block:: shell

    pip install gen_tso

Now, we can create the stellar SED spectrum with this Python script:

.. code-block:: python

    import gen_tso.pandeia_io as pandeia
    import pyratbay.constants as pc
    import pyratbay.spectrum as ps
    import matplotlib.pyplot as plt


    # Use the Gen TSO package to get a PHOENIX SED model for WASP-18 (teff=6430.0, logg=4.31)
    # Closest SED to WASP-18 is an F5V model (Teff=6500K, logg=4.0)
    scene = pandeia.make_scene(
        sed_type='phoenix',
        sed_model='f5v',
    )
    sed_wl, flux = pandeia.extract_sed(scene, wl_range=(0.35,12.0))
    # Convert flux from mJy to erg s-1 cm-2 cm-1
    sed_flux = flux * pc.c / 1e26

    # Lower the resolution to something closer to NIRISS
    bin_wl = ps.constant_resolution_spectrum(0.35, 12.0, resolution=1500.0)
    bin_sed_flux = ps.bin_spectrum(bin_wl, sed_wl, sed_flux, gaps='interpolate')

    # Save to file
    starspec_file = 'phoenix_F5V_6500K_WASP18.dat'
    io.write_spectrum(
        bin_wl,
        bin_sed_flux,
        starspec_file,
        type='emission',
    )

    # Take a look
    plt.figure(0, (7.5,4.5))
    plt.clf()
    plt.subplots_adjust(0.07, 0.11, 0.98, 0.95)
    ax = plt.subplot(111)
    ax.plot(bin_wl, bin_sed_flux, color='royalblue', alpha=0.85)
    ax.set_title('WASP-18 SED Spectrum')
    ax.set_xscale('log')
    ax.set_xlim(0.35, 10.5)
    ax.set_xlabel(r'Wavelength ($\mathrm{\mu}$m)', fontsize=12)
    ax.set_ylabel(r'$F_{\rm p}$ (erg s$^{-1}$ cm$^{-2}$ cm)', fontsize=12)
    ax.tick_params(direction='in', which='both', labelsize=11)

.. plt.savefig('../../figures/phoenix_sed_wasp18.png')

.. image:: ./../../figures/phoenix_sed_wasp18.png
   :width: 70%
   :align: center



----------------------------------------------------------------------

.. _wasp18b_config:

Configuration file
~~~~~~~~~~~~~~~~~~

Lastly, the configuration file will put together the inputs, define
the atmospheric model, and configure the retrieval options.  Here
below is the file we will use for the JWST observation of WASP-39b.

.. wasp18b_retrieval_eclipse_jwst.cfg

.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/wasp18b_retrieval_eclipse_jwst.cfg">wasp18b_retrieval_eclipse_jwst.cfg</a></summary>

.. literalinclude:: ../../_static/data/wasp18b_retrieval_eclipse_jwst.cfg
    :caption: File: wasp18b_retrieval_eclipse_jwst.cfg
    :language: ini

.. raw:: html

   </details>

Lets break this down:


.. tab-set::

  .. tab-item:: General
     :selected:

     .. literalinclude:: ../../_static/data/wasp18b_retrieval_eclipse_jwst.cfg
        :language: ini
        :lines: 3-10

     This first section defines what we want to run. ``runmode``
     indicates that we want a retrieval.  ``logfile`` sets the path to
     the output files.  Note that ``logfile`` can contain a folder,
     which will be created if needed.  Finally, ``verb`` sets the
     screen-output verbosity.


  .. tab-item:: Target

     .. literalinclude:: ../../_static/data/wasp18b_retrieval_eclipse_jwst.cfg
        :language: ini
        :lines: 12-21

     Here we define the observing path of the observation (in this
     case we have a secondary eclipse), the path to the observation file
     discussed above (and the desired output units for plots)

     ``wl_low`` and ``wl_high`` set the and the spectral range to model.
     Note that the wavelenght sampling is partly set by the
     line-sampled opacity files (resolution and maximum wavelength
     coverage).  One can trim the wavelength ranges (as shown here) to
     extract only the region covered by the observations.  One can
     also lower the resolution via a ``wl_thinning = n`` parameter,
     which will take every n-th sample of the opacity files (with
     ``n`` an integer).


     .. literalinclude:: ../../_static/data/wasp18b_retrieval_eclipse_jwst.cfg
        :language: ini
        :lines: 24-35

     And this section defines the system parameters. For an eclipse
     run, the relevant properties will be the stellar radius and SED,
     as well as the planetary mass and radius.

     Note that this ``rplanet`` value is the reference altitute
     situated at the ``ref_pressure`` pressure (this is the constrain
     to compute the layer's :math:`r(p)` profile under hydrostatic
     equilibrium).  Also note that ``ref_pressure`` does not need to be
     at one of the sampled layers (it can be anywhere in between the
     atmosphere pressure range).

     ``starspec`` defines the stellar SED: :math:`F_{\rm
     star}(\lambda)`.  This is a plain-text file with two columns: the
     wavelength (μm) and the surface flux (erg s-1 cm-2 cm).  The
     eclipse depth is ultimately computed as:

     .. math::

         d(\lambda) = F_{\rm star}(\lambda)/F_{\rm planet}(\lambda) *
         (R_{\rm star}/R_{\rm planet})^2,

     with the ``rstar`` and ``rplanet`` parameters defining the stellar
     and planet radii.


  .. tab-item:: Atmosphere

     .. literalinclude:: ../../_static/data/wasp18b_retrieval_eclipse_jwst.cfg
        :language: ini
        :lines: 38-56

     These parameters define the atmospheric-profile models.  The pressure
     parameters are clear, the only constraint is that the bottom pressure
     must be covered by the opacity files.  That is, it's only possible to
     extrapolate to lower pressures (because then the opacities are in the
     Doppler broadening regime, i.e., not dependent on pressure).

     For the temperature profile we will use the [Madhusudhan2009]_
     model.  The parameters will be set below when discussing the
     retrieval parameters.

     For the composition we will adopt VMR profiles in thermochemical
     equilibrium.  We must then define the species to include in the
     atmosphere.  It's very *important* to include not only the
     species that are expected to show up in the spectrum, but also
     the species that are expected to chemically interact with our
     species of interest.  For example, if we expect to model H-
     opacity, we must include electron donors like Na+, K+, Fe+, etc.,
     since the electron density is fundamental to estimate the
     contribution from H-.

     Finally we set the radius-profile model, this is a
     hydrostatic-equilibrium model assuming a variable gravity depending on
     the mass of the planet :math:`g(r) = GM/r^2`.


  .. tab-item:: Absorbers

     .. literalinclude:: ../../_static/data/wasp18b_retrieval_eclipse_jwst.cfg
        :language: ini
        :lines: 59-88

     Now we define the atmospheric absorbers. Make sure that all
     absorber species are included in the atmospheric composition.
     Note that these files impose constraints on the domain that can
     be explored.  The ``sampled_cross_sec`` files determine the maximum
     resolution, spectral range, temperature range, and pressure
     range.  The ``continuum_cross_sec`` files define temperature range
     constraints, but their wavelength ranges can be exceeded.

     In addition to the line-sampled opacities, we add Na and K
     opacity models from [Burrows2000]_, CIA, Rayleigh opacities for
     |H2| and He, and the H- continuum opacity.


  .. tab-item:: Parameters

     .. literalinclude:: ../../_static/data/wasp18b_retrieval_eclipse_jwst.cfg
        :language: ini
        :lines: 90-106

     Here we define the free parameters to modify the elemental
     abundances.  In this case we fit for specific elemental
     metallicities for C an O.  Then we define a catch-all parameter
     for all other metals: [M/H].  These metallicity factors are in
     log10 scale, relative to solar. Thus, values of ``[X/H] = 0`` or
     ``[X/H] = 1`` correspond to 1x and 10x solar respectively.

     .. note:: It is also possible to directly fit ratios between
               species.  This can be set by a parameter named ``X/Y``
               with ``X`` and ``Y`` the elements of interest.  For
               example, a common parameterization is to have a pair of
               parameters for the metallicity and the carbon to oxygen
               ratio: ``[M/H]`` a ``C/O``.

     And then we define the retrievals parameters, their initial values,
     boundaries, and priors.  Since here we will sample the posterior using
     pymultinest [Feroz2009]_ [Buchner2014]_, the most important values are
     the lower and upper boundaries (the initial value is irrelevant for
     the retrieval). The ``step`` value determine which parameters are left
     free to fit (``step>0``) and which are kept fixed at their initial
     value (``step=0``, thus making it trivial to try runs with different
     configurations).

     If desired, one can also set **Gaussian priors** by specifiying the prior
     value and uncertainty after the parameter's ``step``.

     Thus, in summary, this retrieval will fit for the:

     - temperature profile: ``log_p1`` to ``T0`` parameters.
     - composition: ``[M/H]``, ``[C/H]``, and ``[O/H]`` metallicity
       scale factors.


  .. tab-item:: Sampler

     .. literalinclude:: ../../_static/data/wasp18b_retrieval_eclipse_jwst.cfg
        :language: ini
        :lines: 108-121

     Finally, we configure the posterior sampler. In this case we use
     pymultinest [Feroz2009]_ [Buchner2014]_, with 1500 live points.
     ``resume=True`` allows you to pick up a previous run and continue
     from there.

     ``tlow`` and ``thigh`` allow the code to set additional
     temperature-range constraints (beyond those set by the
     temperature-model parameters).

     ``theme`` and ``data_color`` allow you to customize the color of
     the models and data points, respectively, in the output plots.
     Any valid `matplotlib color
     <https://matplotlib.org/stable/users/explain/colors/colors.html#colors-def>`_
     is a valid color.

     The ``wl_ticks`` parameter has two effects: if set, it indicates
     the code to plot wavelengths axes in log scale with the given
     ticks (otherwise defaults to a linear scale).

     The ``post_processing = True`` parameter indicates to compute median
     +/-1sigma, and +/-2sigma statistics out of the posterior distribution.
     Note that this is a post-process step done *after* the posterior
     sampling is finished.  These statistics are computed for the spectra,
     the temperature profiles, contribution functions, and VMRs (along with
     plots of them).  All these data will be neatly packed into a picke
     file.


.. _wasp18b_run:

Retrieval run
-------------

.. Note:: Before running these step, see :ref:`Multinest retrievals
          <ret_multi>` section about MPI / Multinest installation.

To launch the retrieval run, we use the following command from the
prompt.  Since we are using multinest, we will make use of its MPI
parallel-computing capability (thus, the prefix ``mpirun -n 64``):

.. code-block:: shell

    # Launch the retrieval with 64 parallel CPUs
    mpirun -n 64 pbay -c wasp18b_retrieval_eclipse_jwst.cfg

You can adjust the number of CPUs according to your machine/cluster
limitations.  ``Pyrat Bay`` internally uses shared memory to optimize
the memory demand.

That's it. Now we wait until the run is over. This should take from
one to a few days depending on your machine.


Retrieval outputs
~~~~~~~~~~~~~~~~~

TBD

.. _wasp18b_stats:

Detection statistics
--------------------

TBD

