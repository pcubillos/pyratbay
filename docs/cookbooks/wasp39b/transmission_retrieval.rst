.. include:: ../../_substitutions.rst

.. _wasp39b:

Transmission retrieval: WASP-39b JWST
=====================================


This tutorial shows how perform an atmospheric retrieval of the
transmission spectra of WASP-39b, constrained by the JWST
observations.  We will replicate the analysis presented in Welbanks et
al., that is, a retrieval of:

- the combined NIRISS, NIRCam, NIRSpec, and MIRI observations, with
  fixed offsets
- 1D atmosphere with free VMRs (constant-with-altitude)
- a Madhu temperature profile, *T(p)*
- a Rayleigh opacity model
- a gray, patchy cloud-deck model


We can break the analysis into the following steps:

- :ref:`wasp39b_obs`
- :ref:`wasp39b_cross_sec`
- :ref:`wasp39b_config`
- :ref:`wasp39b_run`
- :ref:`wasp39b_stats`

----------------------------------------------------------------------

.. _wasp39b_obs:

Setup
-----

For the setup we will need three ingredients:

#. A **configuration file** to define the system parameters, atmospheric model, posterior sampling, etc.

#. An **observation file** defining the data points: depths, uncertainties, and bin wavelengths

#. **Cross-section files** for the atmospheric species


Lets start with the required input files, and then go over the
configuration file.


.. _wasp39b_obs_file:

Observation file
~~~~~~~~~~~~~~~~

``Pyrat Bay`` observation files tell the code what data is being fit.
These are a plain text files containing the transit depth,
uncertainty, central wavelength (um), bin half-width (um), and label.
There is always one data point per row.  Each data point is being
assumed as a top-hat passband within the respective wavelength range.

Below you can find an extract of the WASP-39b JWST transit data.
Click the link to see/download the entire file.  Note in the header
that comments are allowed, and there are two special flags that let
users define the depth units and where the data starts.


.. literalinclude:: ../../_static/data/obs_wasp39b_transit_jwst.dat
   :caption: File: `obs_wasp39b_transit_jwst.dat <../../_static/data/obs_wasp39b_transit_jwst.dat>`__
   :language: ini
   :lines: 1-20


.. _wasp39b_cross_sec:

Cross sections
~~~~~~~~~~~~~~

Following the analysis of Welbanks et al., we will include
line-sampled cross sections for these molecules: |H2O|, CO, |CO2|,
|CH4|, |SO2|, |H2S|, HCN, |NH3|, and |C2H2|.  Here we will work with
the latests opacity sources for these species from ExoMol and HITEMP.

The current recommendations for sampled cross sections for JWST
retrievals are a resolution :math:`R>20.000`.  So, here we will use a
cross section grid at :math:`R=25.000`, sampling from :math:`0.5-12`
μm in wavelength (to cover the spectral range of the data), from
:math:`500-2500` K in temperature, and from :math:`100-1.0^{-9}` bar
in pressure.

Now, beware that cross section files have many assumptions baked into
them.  In addition, one might need to adjust the ranges or sampling
resolution of the grid for specific project.  Thus, below there are
two options, (a) download and use already made the cross-section
files (b) compute your own cross sections starting from the line-list
files (where you can customize at will).

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
       * - `SO2 <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_SO2_exomol_exoames.npz>`__ (exomol, exoames)
         - [Underwood2016]_
       * - `H2S <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_H2S_exomol_ayt2.npz>`__ (exomol, ayt2)
         - [Azzam2016]_ [Chubb2018]_
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


.. _wasp39b_config:

Configuration file
~~~~~~~~~~~~~~~~~~

Lastly, the configuration file will put together the inputs, define
the atmospheric model, and configure the retrieval options.  Here
below is the file we will use for the JWST observation of WASP-39b.

.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/wasp39b_retrieval_transit_jwst.cfg">wasp39b_retrieval_transit_jwst.cfg</a></summary>

.. literalinclude:: ../../_static/data/wasp39b_retrieval_transit_jwst.cfg
    :caption: File: wasp39b_retrieval_transit_jwst.cfg
    :language: ini

.. raw:: html

   </details>

Lets break this down:


.. tab-set::

  .. tab-item:: General
     :selected:

     .. literalinclude:: ../../_static/data/wasp39b_retrieval_transit_jwst.cfg
        :language: ini
        :lines: 3-10

     This first section defines what we want to run. ``runmode``
     indicates that we want a retrieval.  ``logfile`` sets the path to
     the output files. Note that ``logfile`` can contain a folder,
     which will be created if needed.  Finally, ``verb`` sets the
     screen-output verbosity.


  .. tab-item:: Target

     .. literalinclude:: ../../_static/data/wasp39b_retrieval_transit_jwst.cfg
        :language: ini
        :lines: 12-21

     Here we define the observing path of the observation (in this
     case we have a transit), the location of the observation file
     discussed above (and the desired output units)

     ``wl_low`` and ``wl_high`` set the and the spectral range to model.
     Note that the wavelenght sampling is partly set by the opacity
     files (resolution and maximum wavelength coverage).  One can trim
     the wavelength ranges (as shown here) to extract only the region
     covered by the observations.  One can also lower the resolution
     via a ``wl_thinning = n`` parameter, which will take every n-th
     sample of the opacity files (with ``n`` an integer).


     .. literalinclude:: ../../_static/data/wasp39b_retrieval_transit_jwst.cfg
        :language: ini
        :lines: 24-33

     And this next section defines the system parameters. For a
     transit run, the relevant properties will be the stellar radius
     and the planetary mass, radius, and reference pressure.

     Note that this ``rplanet`` value is the reference altitute
     situated at the ``refpressure`` pressure (this is the constrain
     to compute the layer's :math:`r(p)` profile under hydrostatic
     equilibrium).  Also note that ``refpressure`` does not need to be
     at one of the sampled layers (it can be anywhere in between the
     atmosphere pressure range).


  .. tab-item:: Atmosphere

     .. literalinclude:: ../../_static/data/wasp39b_retrieval_transit_jwst.cfg
        :language: ini
        :lines: 36-55

     These parameters define the atmospheric-profile models.  The pressure
     parameters are clear, the only constraint is that the bottom pressure
     must be covered by the opacity files.  That is, it's only possible to
     extrapolate to lower pressures (because then the opacities are in the
     Doppler broadening regime).

     For the temperature profile we will use the [Madhusudhan2009]_
     model. The parameters will be set below when discussing the retrieval
     parameters.

     For the composition we will adopt constant-with-altitude VMR profiles,
     so we define the species to include and their initial VMR values.

     Finally we set the radius-profile model, this is a
     hydrostatic-equilibrium model assuming a variable gravity depending on
     the mass of the planet :math:`g(r) = GM/r^2`.

  .. tab-item:: Absorbers

     .. literalinclude:: ../../_static/data/wasp39b_retrieval_transit_jwst.cfg
        :language: ini
        :lines: 58-89


     Now we define the atmospheric absorbers.  Make sure that all absorber
     species are defined in the atmospheric composition.  Note that some of
     these set constraints to the domain to be expored.  The
     ``sampled_cross_sec`` files set the maximum resolution, spectral
     range, temperature range, and pressure.  The ``continuum_cross_sec``
     files set temperature range constraints, but one can span beyond their
     wavelength ranges.

     On top of the line sampled opacities, we add Na and K alkali models from
     [Burrows2000]_, CIA, and |H2| and He Rayleigh opacities.  For clouds
     we add a cloud-deck model and a parametric [Lecavelier2008]_ Rayleigh
     model.

  .. tab-item:: Parameters

     .. literalinclude:: ../../_static/data/wasp39b_retrieval_transit_jwst.cfg
        :language: ini
        :lines: 91-132

     Here we define which VMR abundances are considered free
     parameters and which are *'bulk'* (i.e., filler) abundances.  All
     other species' VMRs will be kept fixed at their initial values as
     defined above.

     And then we define the retrievals parameters, their initial values,
     boundaries, and priors.  Since here we will sample the posterior using
     pymultinest [Feroz2009]_ [Buchner2014]_, the most important values are
     the lower and upper boundaries (the initial value is irrelevant for
     the retrieval). The ``step`` value determine which parameters are left
     free to fit (``step>0``) and which are kept fixed at their initial
     value (``step=0``, thus making it trivial to try runs with different
     configurations).

     If desired, one can also set **Gaussian priors** by specifiying the prior
     value and uncertainty after the parameter's ``step``, as done here for
     the planet mass.

     Thus, in summary, this retrieval will fit:

     - temperature profile: ``log_p1`` to ``T0`` parameters.
     - mass and reference pressure: ``M_planet`` (units as in ``mplanet``)
       and ``log_p_ref`` in :math:`\log10(p/{\rm bar})`
     - aundances: ``log_H2O`` to ``log_K`` values set the
       :math:`\log10({\rm VMR})` for the given species.
     - Rayleigh: ``log_k_ray`` and ``alpha_ray`` set the strength and slope
       of the [Lecavelier2008]_ model.
     - cloud deck: ``log_p_cl`` sets the cloud-top pressure in
       :math:`\log10(p/{\rm bar})` units.
     - cloud patchy factor: ``f_patchy`` a value between 0 and 1, which
       affects both cloud deck and [Lecavelier2008]_ models.


  .. tab-item:: Sampler

     .. literalinclude:: ../../_static/data/wasp39b_retrieval_transit_jwst.cfg
        :language: ini
        :lines: 134-147

     Finally, we configure the posterior sampler. In this case we use
     pymultinest [Feroz2009]_ [Buchner2014]_, with 1000 live points.
     ``resume=True`` allows you to pick up a previous run and continue from
     there.

     ``tlow`` and ``thigh`` allow the code to set additional temperature
     range constraints (beyond those set by the temperature-model
     parameters).

     ``theme`` and ``data_color`` allow you to customize the color of the
     models and data points in the output plots.  Any valid `matplotlib
     color
     <https://matplotlib.org/stable/users/explain/colors/colors.html#colors-def>`_
     is a valid color.

     The ``wl_ticks`` parameter has two effects: if set, it makes the code
     to plot wavelengths axes in log scale with the given ticks (otherwise
     defaults to a linear scale).

     The ``post_processing = True`` parameter indicates to compute median
     +/-1sigma, and +/-2sigma statistics out of the posterior distribution.
     Note that this is a post-process step done *after* the posterior
     sampling is finished.  These statistics are computed for the spectra,
     the temperature profiles, contribution functions, and VMRs (along with
     plots of them).  All these data will be neatly packed into a picke
     file.


.. _wasp39b_run:

Retrieval run
-------------

To launch the retrieval run, we use the following command from the
prompt.  Since we are using multinest, we will make use of its MPI
parallel-computing capability (thus, the prefix ``mpirun -n 64``):

.. code-block:: shell

    # Launch the retrieval with 64 parallel CPUs
    mpirun -n 64 pbay -c wasp39b_retrieval_transit_jwst.cfg

You can adjust the number of CPUs according to your machine/cluster
limitations.  ``Pyrat Bay`` internally uses shared memory to optimize
the memory demand.

That's it. Now we wait until the run is over. This should take from
one to a few days depending on your machine.

Retrieval outputs
~~~~~~~~~~~~~~~~~

TBD

.. _wasp39b_stats:

Detection statistics
--------------------

OK, we have now a posterior distribution for species on WASP-39b based
on the JWST observations, for some there are well constrained VMRs,
for others there are upper limits.  We want now to assess the
siginificance of each detection. For this we will run a series of
leave-one-out retrievals with the same configuration as before, but
removing opacity from one species at a time.

Here's an extract of what changed in the cofiguration for the run without |H2O|:

.. code-block:: ini
    :emphasize-lines: 5,19

    # Pyrat Bay run mode [tli pt atmosphere spectrum radeq opacity mcmc]
    runmode = retrieval

    # Output file names
    logfile = ret_wasp39b_no_H2O/WASP39b_jwst_no_H2O.log

    ...

    # Param name    value  lo_bound  hi_bound  step   prior  prior_sigma
    retrieval_params =
        log_p1        -4.0     -9.0       2.0   0.3
        log_p2        -7.2     -9.0       2.0   0.3
        log_p3        -1.0     -2.0       2.0   0.3
        a1            1.50     0.02       2.0   0.02
        a2            0.35     0.02       2.0   0.02
        T0           850.0    800.0    1300.0   30.0
        M_planet     0.266      0.1       0.43  0.05  0.266  0.033
        log_p_ref    -1.0      -9.0       2.0   0.3
        log_H2O      -10.0    -12.0      -0.3   0.0
        log_CO2      -3.00    -12.0      -0.3   0.3

    ...

As you see, we only change two lines:

- edit ``logfile`` to set the proper name for the outputs
- fix the |H2O| free parameter (``step=0``) and set the initial value to a
  negligible VMR (``-10.0``)

And then, there are two more optimizations for a better efficiency:

- to remove the |H2O| file from ``sampled_cross_sec`` to consume less
  resources
- set ``post_processing = False``, this will prevent the code from
  running the post processing after the retrieval.  We will do that in
  a separate call in a background process.

Here are sample config files for leave-one-out runs for |H2O|, |CO2|, and |SO2|:

- `wasp39b_retrieval_transit_jwst_no_H2O.cfg <../../_static/data/wasp39b_retrieval_transit_jwst_no_H2O.cfg>`__
- `wasp39b_retrieval_transit_jwst_no_CO2.cfg <../../_static/data/wasp39b_retrieval_transit_jwst_no_CO2.cfg>`__
- `wasp39b_retrieval_transit_jwst_no_SO2.cfg <../../_static/data/wasp39b_retrieval_transit_jwst_no_SO2.cfg>`__


The recommendation is, since we want to run a series of retrievals, we
write a script like the one below to concatenate one run after the
other.

Note that we added the ``pbay --post ...`` calls after each
posterior sampling, with an ampersand at the end to rnu it in the
background. This is useful since the post-processing uses only a
single CPU and might take a few hours to complete.  Putting it in
background allow us to launch each retrival right after the previous
one.

.. code-block:: shell
   :caption: File: wasp39b_loo_retrievals.sh

    # Launch the retrievals, and then the post-processing in the background
    mpirun -n 64 pbay -c wasp39b_retrieval_transit_jwst_no_H2O.cfg
    pbay --post wasp39b_retrieval_transit_jwst_no_H2O.cfg &

    mpirun -n 64 pbay -c wasp39b_retrieval_transit_jwst_no_CO2.cfg
    pbay --post wasp39b_retrieval_transit_jwst_no_CO2.cfg &

    mpirun -n 64 pbay -c wasp39b_retrieval_transit_jwst_no_SO2.cfg
    pbay --post wasp39b_retrieval_transit_jwst_no_SO2.cfg &


Then to start the retrievals, run this command from the prompt:

.. code-block:: shell

    # Launch leave-one-out retrievals
    sh wasp39b_loo_retrievals.sh
