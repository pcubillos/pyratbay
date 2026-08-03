.. include:: ../../_substitutions.rst

.. _yses1b:

Isotopic-ratio retrieval: YSES-1b
=================================

This tutorial shows how perform an atmospheric retrieval of direct
imaging observations of YSES-1b with SINFONI.  Here we will learn
about atmospheric analyses including:

- ground-based high-resolution data
- direct imaging flux measurements
- isotopic-ratio modeling

We can break the analysis into the following steps:

- :ref:`yses1b_data`
    - :ref:`yses1b_obs`
    - :ref:`yses1b_cross_sec`

- :ref:`yses1b_retrievals`
    - :ref:`yses1b_config`
    - :ref:`yses1b_run`
    - :ref:`yses1b_post`

----------------------------------------------------------------------

.. _yses1b_data:

File inputs
-----------

For the setup we will need three ingredients:

#. A **configuration file** to define the system parameters, atmospheric model, posterior sampling, etc.

#. An **observation file** defining the data points: depths, uncertainties, and bin wavelengths

#. **Cross-section files** for the atmospheric species


Lets start with the required input files, and then go over the
configuration file.


.. _yses1b_obs:

Observation file
~~~~~~~~~~~~~~~~

``Pyrat Bay`` observation files tell the code what data is being fit.
These are a plain text files containing the spectrum, its
uncertainties, central wavelength (μm), bin half-width (μm), and (an
optional) label.  Here, the data point are assumed to be at the
instrumental pixel-level sampling (zero bin half-width).

Below you can find an extract of the YSES-1b emission data from SINFONI.
Click the link to see/download the entire file.


.. literalinclude:: ../../_static/data/obs_YSES1b_sinfoni.dat
   :caption: File: `obs_YSES1b_sinfoni.dat <../../_static/data/obs_YSES1b_sinfoni.dat>`__
   :language: ini
   :lines: 1-14


----------------------------------------------------------------------

.. _yses1b_cross_sec:

Cross section files
~~~~~~~~~~~~~~~~~~~

For this example we will focus on the molecular absorption from |H2O|
and CO, which dominate the emission spectrum in the K band probed by
the SINFONI observation.  We thus must generate line-sampled cross
sections for them at the required resolution to probe individual
isotopic lines. We can split this into two steps:

- Download and format line-by-line data
- Compute cross sections

**Line-by-line data**

For the line-by-line data we will use the HITEMP database for CO and
the Exomol database for |H2O|, pre-processed with ``repack``
[Cubillos2017b]_ to work with only the dominant line transitions.
The table below shows then the data to use:

.. list-table::
  :header-rows: 1

  * - Species
    - Source
    - References
  * - `H2O <https://zenodo.org/records/14266247/files/H2O_exomol_pokazatel_0.24-500.0um_100-3500K_threshold_0.01_lbl.dat>`__
    - Exomol / pokazatel
    - [Polyansky2018]_
  * - `CO <https://hitran.org/files/HITEMP/bzip2format/05_HITEMP2019.par.bz2>`__
    - HITEMP
    - [Li2015]_


Here are the configuration files and script to format the LBL data for
use in ``Pyrat Bay``.  Noteworthy of these files is that the data will
use the TIPS partition function info from [Gamache2021]_ and that they
define the maximum wavelength range to cover.  The script below shows
how to download the LBL data and format it.

.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/lbl_format_H2O_repack_pokazatel.cfg">lbl_format_H2O_repack_pokazatel.cfg</a></summary>

.. literalinclude:: ../../_static/data/lbl_format_H2O_repack_pokazatel.cfg
    :caption: File: lbl_format_H2O_repack_pokazatel.cfg
    :language: ini

.. raw:: html

   </details>

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/lbl_format_CO_hitemp_li.cfg">lbl_format_CO_hitemp_li.cfg</a></summary>

.. literalinclude:: ../../_static/data/lbl_format_CO_hitemp_li.cfg
    :caption: File: lbl_format_CO_hitemp_li.cfg
    :language: ini

.. raw:: html

   </details>


.. code-block:: shell

    # Fetch line-by-line data
    wget https://zenodo.org/records/14266247/files/H2O_exomol_pokazatel_0.24-500.0um_100-3500K_threshold_0.01_lbl.dat
    wget https://hitran.org/files/HITEMP/bzip2format/05_HITEMP2019.par.bz2
    bzip2 -d 05_HITEMP2019.par.bz2

    # Parse line-by-line data into pyratbay format
    pbay -c lbl_format_CO_hitemp_li.cfg
    pbay -c lbl_format_H2O_repack_pokazatel.cfg


**Compute cross-section files**

``Pyrat Bay`` provides the code to sample line-by-line data from the
latests opacity sources into *custom* cross section files at desired
resolution, ranges, and isotopes.  For the SINFONI K-band data, we
will adopt a wavelength sampling resolution of R=200,000 over
:math:`1.9-2.4` μm.  For temperature and pressure we sample between
:math:`150-3000` K and :math:`100-1.0^{-9}` bar, respectively.
These cross sections have been computed assuming an |H2|/He-dominated
atmosphere, and terrestrial isotopic ratios.  The lines have Voigt
profiles with a wing cut-off at 500 HWHM.

Here are the configuration files to compute the cross section tables:

.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/opacity_exomol_H2O_Kband_hires.cfg">opacity_exomol_H2O_Kband_hires.cfg</a></summary>

.. literalinclude:: ../../_static/data/opacity_exomol_H2O_Kband_hires.cfg
    :caption: File: opacity_exomol_H2O_Kband_hires.cfg
    :language: ini

.. raw:: html

   </details>

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/opacity_hitemp_CO_Kband_hires_iso26.cfg">opacity_hitemp_CO_Kband_hires_iso26.cfg</a></summary>

.. literalinclude:: ../../_static/data/opacity_hitemp_CO_Kband_hires_iso26.cfg
    :caption: File: opacity_hitemp_CO_Kband_hires_iso26.cfg
    :language: ini

.. raw:: html

   </details>

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/opacity_hitemp_CO_Kband_hires_iso36.cfg">opacity_hitemp_CO_Kband_hires_iso36.cfg</a></summary>

.. literalinclude:: ../../_static/data/opacity_hitemp_CO_Kband_hires_iso36.cfg
    :caption: File: opacity_hitemp_CO_Kband_hires_iso36.cfg
    :language: ini

.. raw:: html

   </details>

**Per-isotope cross sections**

Note that for |H2O| the cross sections are evaluated assuming the
terrestrial isotopic ratios. For CO instead, we compute *one
cross-section table per isotope* we are interested in the two most
abundant ones ¹²CO and ¹³CO.  This is achived via the
``single_isotope`` key, in which users must specify the isotopolog to
extract, for example:


.. literalinclude:: ../../_static/data/opacity_hitemp_CO_Kband_hires_iso26.cfg
    :caption: Extract from `opacity_hitemp_CO_Kband_hires_iso26.cfg <../../_static/data/opacity_hitemp_CO_Kband_hires_iso26.cfg>`__
    :language: ini
    :lines: 35-37

Note that when the ``single_isotope`` key is set, internally the cross
section with assume an isotopic abundance ratio of 1.0.  Later during
spectral sampling / atmospheric retrieval users must set the isotopic
ratios for the respective cross section(s).

Use the followin script to compute the cross sections:

.. code-block:: shell

    # Compute cross section table for H2O
    pbay -c opacity_exomol_H2O_Kband_hires.cfg

    # Compute cross section table for CO per isotope
    pbay -c opacity_hitemp_CO_Kband_hires_iso26.cfg
    pbay -c opacity_hitemp_CO_Kband_hires_iso36.cfg




----------------------------------------------------------------------


.. _yses1b_retrievals:

Retrieval analysis
------------------

.. _yses1b_config:

Configuration file
~~~~~~~~~~~~~~~~~~


The configuration file will put together the inputs, define the
atmospheric model, and configure the retrieval options.  Here is the
configuration file for the YSES-1b retrieval analysis:

.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/yses1b_retrieval_emission_sinfoni.cfg">yses1b_retrieval_emission_sinfoni.cfg</a></summary>

.. literalinclude:: ../../_static/data/yses1b_retrieval_emission_sinfoni.cfg
    :caption: File: yses1b_retrieval_emission_sinfoni.cfg
    :language: ini

.. raw:: html

   </details>


.. _yses1b_isotopic:

Isotopic ratio modeling
.......................

Here are the relevant keys to model isotopic ratios:

.. literalinclude:: ../../_static/data/yses1b_retrieval_emission_sinfoni.cfg
   :language: ini
   :lines: 54-59, 79-86, 98

- ``sampled_cross_sec`` defines the input line-sampled cross section
  data for the isotopes of interest (in this case, ¹²CO and ¹³CO)
- ``isotope_ratios`` determines which isotopic ratios are being
  modeled.  Each row configures an isotope via these 3 fields:

  - ``cso_iso`` identifies the files in ``sampled_cross_sec`` to
    assign to these isotopes. E.g., here the files containing
    *'iso26'* and *'iso36'* in their names.
  - ``label`` sets a label for the isotope (see below for
    ``iso_ratio`` and ``retrieval_params``)
  - ``iso_ratio`` sets the value for the isotope. There are two
    options: one is to directly set the value in log scale (e.g., here
    the ratio for for ¹³CO is set to :math:`\log_{10}(13{\rm CO}) =
    -1.96`).  Alternatively, an isotope ratio can be set as a *filler*
    value via the ``fill_iso`` format.  E.g., here the ¹²CO ratio is
    set to ¹²CO = 1.0 - ¹³CO. This is typically reserved for the
    dominant/most abundant isotope.

    .. note:: Note that a filler isotope ratio can be set in
              combination with multiple other isotopes.  For example
              if ¹²CO, ¹³CO, and C¹⁸O were defined. ¹²CO could be set
              as a filler of the other two (¹²CO = 1.0 - ¹³CO - C¹⁸O)
              by setting its value to ``fill_36_28``

- ``retrieval_params`` enables to use isotopic ratios as free
  parameters. The parameter name syntax is ``iso_label``, where
  ``label`` is the label defined in ``isotope_ratios``.


.. _yses1b_imaging:

Direct imaging
..............

Here are the relevant keys for direct imaging data:

.. literalinclude:: ../../_static/data/yses1b_retrieval_emission_sinfoni.cfg
   :language: ini
   :lines: 12-14, 24-25

- ``rt_path`` indicates that the input data is measured as incoming
  flux in W m⁻² μm⁻¹, as received on Earth
- ``distance`` sets the distance between Earth and the target, required
  to calculate the received flux

.. _yses1b_hires:

High-resolution data
....................

Here are the relevant keys for high-resolution data:

.. literalinclude:: ../../_static/data/yses1b_retrieval_emission_sinfoni.cfg
   :language: ini
   :lines:  15-18, 84-87

- ``obsfile_hires`` (as opposed to ``obsfile``) indicates that the
  input data is high-resolution data evaluated at/near the pixel
  sampling.  In this case the code will compute the model spectrum at
  the resolution defined by the sampled cross sections, then convolve
  to the ``inst_resolution`` resolving power, and then evaluate at the
  wavelength of the data
- ``inst_resolution`` sets the instrumental resolving power
- In ``retrieval_params``, the ``rv_shift`` parameter enables the code
  to apply a Doppler shift to the input data (in km s⁻¹ units)

Overview
........

.. tab-set::

  .. tab-item:: General
     :selected:

     .. literalinclude:: ../../_static/data/yses1b_retrieval_emission_sinfoni.cfg
        :language: ini
        :lines: 3-10

     This first section defines what we want to run. ``runmode``
     indicates that we want a retrieval.  ``logfile`` sets the path to
     the output files. Note that ``logfile`` can contain a folder,
     which will be created if needed.  Finally, ``verb`` sets the
     screen-output verbosity.


  .. tab-item:: Target

     .. literalinclude:: ../../_static/data/yses1b_retrieval_emission_sinfoni.cfg
        :language: ini
        :lines: 12-22

     Here we define the observing path of the observation (in this
     case a direct flux measurement) and the location of the
     observation file.

     ``wl_low`` and ``wl_high`` set the and the spectral range to model.
     Note that the wavelenght sampling is partly set by the opacity
     files (resolution and maximum wavelength coverage).  One can trim
     the wavelength ranges (as shown here) to extract only the region
     covered by the observations.  One can also lower the resolution
     via a ``wl_thinning = n`` parameter, which will take every n-th
     sample of the opacity files (with ``n`` an integer).


     .. literalinclude:: ../../_static/data/yses1b_retrieval_emission_sinfoni.cfg
        :language: ini
        :lines: 24-29

     This section defines the system parameters. For a direct flux
     measurement the relevant properties will be the distance to the
     target, planetary mass and radius, and reference pressure.

     Note that this ``rplanet`` value is the reference altitute
     situated at the ``ref_pressure`` pressure (this is the constrain
     to compute the layer's :math:`r(p)` profile under hydrostatic
     equilibrium).  Also note that ``ref_pressure`` does not need to be
     at one of the sampled layers (it can be anywhere in between the
     atmosphere pressure range).


  .. tab-item:: Atmosphere

     .. literalinclude:: ../../_static/data/yses1b_retrieval_emission_sinfoni.cfg
        :language: ini
        :lines: 31-52

     These parameters define the atmospheric-profile models.
     For the temperature profile here we use the [Madhusudhan2009]_
     model. The parameters will be set below when discussing the retrieval
     parameters.

     For the composition we will adopt thermochemical-equilibrium
     abundances, parameterized by the metallicity and C/O ratio.  Note
     that many more species than those with cross sections are
     defined. This is necessary to perform correct
     thermochemical-equilibrium calculations.

     Finally we set the radius-profile model, this is a
     hydrostatic-equilibrium model assuming a variable gravity depending on
     the mass of the planet :math:`g(r) = GM/r^2`.


  .. tab-item:: Absorbers

     .. literalinclude:: ../../_static/data/yses1b_retrieval_emission_sinfoni.cfg
        :language: ini
        :lines: 54-72

     Now we define the atmospheric absorbers.  Make sure that all
     absorber species are defined in the atmospheric composition.
     Note that some of these set constraints to the domain to be
     expored.  The ``sampled_cross_sec`` files set the maximum
     resolution, spectral range, temperature range, and pressure.  The
     ``continuum_cross_sec`` files set temperature range constraints,
     but one can span beyond their wavelength ranges.  On top of the
     line sampled opacities, we include |H2| and He Rayleigh opacities
     and a cloud-deck model.

  .. tab-item:: Parameters

     .. literalinclude:: ../../_static/data/yses1b_retrieval_emission_sinfoni.cfg
        :language: ini
        :lines: 74-100

     Here we first define how to parameterize the atmospheric
     composition and the isotopic ratios.

     Then we define the retrievals parameters, their initial values,
     boundaries, and priors.  Since here we will sample the posterior using
     pymultinest [Feroz2009]_ [Buchner2014]_, the most important values are
     the lower and upper boundaries (the initial value is irrelevant for
     the retrieval). The ``step`` value determine which parameters are left
     free to fit (``step>0``) and which are kept fixed at their initial
     value (``step=0``, thus making it trivial to try runs with different
     configurations).

     If desired, one can also set **Gaussian priors** by specifiying the prior
     value and uncertainty after the parameter's ``step``, as done here for
     the planet mass (prior from Bohn+2020).


  .. tab-item:: Sampler

     .. literalinclude:: ../../_static/data/yses1b_retrieval_emission_sinfoni.cfg
        :language: ini
        :lines: 102-114

     Finally, we configure the posterior sampler. In this case we use
     pymultinest [Feroz2009]_ [Buchner2014]_, with 1000 live points.
     ``resume=True`` allows you to pick up a previous run and continue from
     there.

     ``tlow`` and ``thigh`` allow the code to set additional temperature
     range constraints.

     ``theme`` and ``data_color`` allow you to customize the color of the
     models and data points in the output plots.  Any valid `matplotlib
     color
     <https://matplotlib.org/stable/users/explain/colors/colors.html#colors-def>`_
     is a valid color.

.. The ``post_processing = True`` parameter indicates to compute
     median +/-1sigma, and +/-2sigma statistics out of the posterior
     distribution.  Note that this is a post-process step done *after*
     the posterior sampling is finished.  These statistics are
     computed for the spectra, the temperature profiles, contribution
     functions, and VMRs (along with plots of them).  All these data
     will be neatly packed into a picke file.


.. _yses1b_run:

Retrieval run
~~~~~~~~~~~~~

To launch the retrieval run, we use the following command from the
prompt.  Since we are using multinest, we will make use of its MPI
parallel-computing capability (thus, the prefix ``mpirun -n 64``):

.. code-block:: shell

    # Launch the retrieval with 64 parallel CPUs
    mpirun -n 64 pbay -c yses1b_retrieval_emission_sinfoni.cfg

You can adjust the number of CPUs according to your machine/cluster
limitations.  ``Pyrat Bay`` internally uses shared memory to optimize
the memory demand.

That's it. Now we wait until the run is over. This should take from
one to a few days depending on your machine.

.. _yses1b_post:

Retrieval outputs
~~~~~~~~~~~~~~~~~

TBD
