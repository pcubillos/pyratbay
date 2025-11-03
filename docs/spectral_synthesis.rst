.. include:: _substitutions.rst

.. _spectral_synthesis:

Spectral Synthesis
==================

This tutorial shows how compute transmission, emission, or eclipse
spectra with ``Pyrat Bay``.

- :ref:`spec_config`
- :ref:`spec_system`
- :ref:`spec_star`
- :ref:`spec_wavelength`
- :ref:`spec_observations`
- :ref:`spec_atmosphere`
- :ref:`spec_opacities`
- :ref:`spec_parameters`
- :ref:`spec_demo`


.. _spec_config:

Configuration File
------------------

To compute spectra, use a configuration file with the ``runmode`` key
set to ``spectrum``.  Spectrum runs also require a ``logfile`` key,
which sets the name of the output log and spectrum files.


**Observing Geometry**

The ``rt_path`` key sets the radiative-transfer scheme and observing
geometry to use.  These are the options:


.. list-table::
   :header-rows: 1
   :widths: 5, 10, 20, 30

   * - Observing geometry
     - ``rt_path``
     - Output spectrum
     - Comments

   * - Transmission
     - ``transit``
     - (|Rp|/|Rs|)\ :sup:`2`
     - Transmission spectrum

   * - Eclipse
     - ``eclipse``
     - |Fp|/|Fs|
     - Occultation spectrum

   * - Eclipse
     - ``eclipse_two_stream``
     - |Fp|/|Fs|
     - Occultation spectrum

   * - Emission
     - ``emission``
     - |Fp| (erg s\ :sup:`-1` cm\ :sup:`-2` cm)
     - Flux at the planet's surface

   * - Emission
     - ``emission_two_stream``
     - |Fp| (erg s\ :sup:`-1` cm\ :sup:`-2` cm)
     - Flux at the planet's surface

   * - Emission
     - ``f_lambda``
     - |Fp| (W m\ :sup:`-2` μm\ :sup:`-1`)
     - Flux measured at Earth


Here is a sample configuration file to compute a
transmission spectrum:


.. literalinclude:: ../examples/tutorial/spectrum_transmission.cfg
   :language: ini
   :caption: File `spectrum_transmission.cfg <../../_static/data/tutorial/spectrum_transmission.cfg>`__


The output spectrum can be set with the ``specfile`` key, otherwise it
is taken from the ``logfile`` name (replacing the `.log` file
extension with `.dat`)


.. _spec_system:

System parameters
-----------------

The system parameters have multiple uses.

Hill radius
~~~~~~~~~~~

The ``mstar``, ``mplanet``, and ``smaxis`` keys set the stellar mass,
planetary mass, and orbital semi-major axis.  If these keys are set in
the configuration file, the code will compute the planetary Hill
radius (:math:`R_{\rm H} = a \sqrt[3]{M_{\rm p}/3M_{\rm s}}`).  In
such case, ``Pyrat Bay`` will neglect atmospheric layers at altitudes
larger than :math:`R_{\rm H}`, since they should not be
gravitationally bound to the planet.



Radius ratio
~~~~~~~~~~~~

To compute eclipse depths from the emission spectra (|Fp|), the user
needs to set the ``rstar`` and ``rplanet`` keys, which define the
stellar and planetary radius.  The eclipse depths can then be computed as:

.. math::
    {\rm Eclipse\ depth} = \frac{F_{\rm p}}{F_{\rm s}}
                  \left(\frac{R_{\rm p}}{R_{\rm s}}\right)^2

The ``tstar`` key sets the stellar effective temperature, which can be
used to define a stellar blackbody spectrum (|Fs|, see :ref:`starspec`).


.. _spec_star:

Stellar Spectrum
----------------

The stellar spectrum is required to compute eclipse depths as the
planet-to-star flux spectrum (for transit calculations, a stellar
spectrum is not required).  ``Pyrat Bay`` provides several options to
set a stellar spectrum.

.. tab-set::

  .. tab-item:: Custom spectrum
     :selected:

     Users can use their own custom stellar spectra via the
     ``starspec`` argument.  This must point to a plain file file
     containing a spectrum in two columns: the first column has the
     wavelength array in microns, the second column has the flux
     spectrum in erg s\ :sup:`-1` cm\ :sup:`-2` cm units.

     .. code-block:: ini

         # Custom stellar spectrum file
         starspec = inputs/WASP18_spectrum.dat

  .. tab-item:: Kurucz model

     Users can use a Kurucz stellar model [Castelli2003]_ via the
     ``kurucz`` argument of the configuration file, pointing to a
     Kurucz model.  These models can be downloaded from `this link
     <http://kurucz.harvard.edu/grids/>`__.  The code selects the
     correct Kurucz model based on the stellar temperature and surface
     gravity values:

     .. code-block:: ini

         # Kurucz stellar spectrum
         tstar = 5700
         log_gstar = 4.5
         kurucz = inputs/fp00k2odfnew.pck

  .. tab-item:: Black body

     By defining the stellar effective temperature ``tstar``, the code
     will adopt a blackbody spectrum for the star (unless the
     ``starspec`` or ``kurucz`` arguments have been set).

     .. code-block:: ini

         # Stellar effective temperature (K)
         tstar = 5700


.. _spec_atmosphere:

Atmosphere Model
----------------

Users can choose to compute the
The ``atmfile`` key sets the input atmospheric model from which to
compute the spectrum.  If the file pointed by ``atmfile`` does not
exist, the codel will attempt to produce it (provided all necessary
input parameters are set in the configuration file).  The atmospheric
model can be produced with ``Pyrat Bay`` or be a custom input from the
user.


The user can re-compute the temperature profile of the atmosphere by
specifying the ``tmodel`` and ``tpars`` keys (see :ref:`temp_profile`).

Radius profile
~~~~~~~~~~~~~~

The ``radmodel`` key sets the model to compute the atmospheric
layers's altitude assuming hydrostatic equilibrium.  This table shows
the currently available models:

=====================  =========================
Models (``radmodel``)  Comments
=====================  =========================
hydro_m                Hydrostatic equilibrium with :math:`g(r)=GM/r^2`
hydro_g                Hydrostatic equilibrium with constant gravity
[undefined]            Take radius profile from input atmospheric file if exists
=====================  =========================

See the :ref:`altitude_profile` section for details.

The ``refpressure``, ``rplanet``, ``mplanet`` and ``gplanet`` keys set
the planetary reference pressure and radius level (:math:`p_0` and
:math:`R_0`), the planetary mass (:math:`M_p`) and planetary surface
gravity (:math:`g`), respectively.

.. Note:: Note that the user can supply its own atmospheric altitude
          profile (through the input atmospheric model), possibly not
          in hydrostatic equilibrium.  In this case, do not set the
          ``radmodel`` key.


.. _spec_wavelength:

Spectrum sampling
-----------------

The ``wllow`` and ``wlhigh`` keys set the wavelength boundaries for
the output spectrum (values must contain units; otherwise, set the
units with the ``wlunits`` key).


The ``wnstep`` sets the sampling rate in |kayser|.  Note that this
will be the output sampling rate.  Internally, ``Pyrat Bay`` must
compute line profiles at a higher resolution to ensure not to
undersample the line profiles.  The ``wnosamp`` key (an integer) sets
the oversampling factor of the high-resolution sampling relative to
``wnstep`` (that is, the high-resolution sampling rate is
``wnstep/wnosamp``).  Typical values for the optical/IR are ``wnstep =
1.0`` and ``wnosamp = 2000``.

.. https://en.wikipedia.org/wiki/Highly_composite_number

Alternatively, the user can request a constant-resolution output by
setting the ``resolution`` key (where the resolution is
:math:`R=\lambda/\Delta\lambda`).


.. _spec_opacities:

Opacities
---------

  - cross sections
  - clouds


Line-list data
~~~~~~~~~~~~~~

Use the ``tlifile`` key to include TLI file(s) containing LBL
opacities (to create a TLI file, see :ref:`tlitutorial`).  The user
can include zero, one, or multiple TLI files if desired.

Note that the ``tlifile`` opacities will be neglected if the
configuration file sets input LBL opacities through the ``extfile``
(see the rules in :ref:`opacity_io`).

.. _cia_opacity:

Collision Induced Absorption
----------------------------

Use the ``csfile`` key to include opacities from cross-section files.
A cross-section file contains opacity (in |kayser| amagat\ :sup:`-2`
units) tabulated as a function of temperature and wavenumber.  Since
this format is tabulated in wavenumber, it is best suited for
opacities that vary smoothly with wavenumber, like collision-induced
absorption (CIA).  However, the code can also process LBL opacities,
as long as the files follow the right format (more on this later).

The following table list the most-commonly used CIA opacity sources:

========== ========== ========== ===================== =====================
Sources    Species    T range    |nu| range (|kayser|) References
========== ========== ========== ===================== =====================
`HITRAN`_  |H2|--|H2|  200--3000 1.0--500.0            [Richard2012]_ [Karman2019]_
`HITRAN`_  |H2|--H    1000--2500 1.0--100.0            [Richard2012]_ [Karman2019]_
`HITRAN`_  |H2|--He    200--9900 0.5--500.0            [Richard2012]_ [Karman2019]_
`Borysow`_ |H2|--|H2|   60--7000 0.6--500.0            [Borysow2001]_ [Borysow2002]_
`Borysow`_ |H2|--He     50--7000 0.5--31.25            [Borysow1988]_ [Borysow1989a]_ [Borysow1989b]_ [Jorgensen2000]_
========== ========== ========== ===================== =====================


.. _Borysow: https://www.astro.ku.dk/~aborysow/programs/index.html
.. _HITRAN: https://hitran.org/cia


For the **HITRAN** CIA database, ``Pyrat Bay`` provides these shell commands to re-format the downloaded CIA files

.. code-block:: shell

    # Download and format HITRAN H2-H2 CIA file for Pyrat Bay:
    $ wget https://hitran.org/data/CIA/H2-H2_2011.cia
    $ pbay -cs hitran H2-H2_2011.cia

    # And for HITRAN H2-He CIA
    $ wget https://hitran.org/data/CIA/H2-H2_2011.cia
    $ pbay -cs hitran H2-He_2011.cia

For the **Borysow** CIA database, the code provides already-formatted
files for |H2|-|H2| in the 60--7000K and 0.6--500 um range \[`here
<https://github.com/pcubillos/pyratbay/blob/master/pyratbay/data/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat>`__\]
(this file pieces together the tabulated |H2|-|H2| files described in
the references above); and for |H2|-He in the 50--7000K and 0.5--31.25
um range \[`here
<https://github.com/pcubillos/pyratbay/blob/master/pyratbay/data/CIA/CIA_Borysow_H2He_0050-7000K_0.5-031um.dat>`__\]
(this file was created using a re-implementation of the code described
in the references above).  The user can access these files via the
``{ROOT}`` shortcut, as in the example below:

.. code-block:: python

    csfile =
        {ROOT}pyratbay/data/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat
        {ROOT}pyratbay/data/CIA/CIA_Borysow_H2He_0050-7000K_0.5-031um.dat


Alkali Opacity
~~~~~~~~~~~~~~

Use the ``alkali`` key to include opacities from alkali species.
Currently, the code provides the [Burrows2000]_ models for the sodium
and potassium resonant lines, based on van der Waals and statistical
theory.  The following table lists the available alkali model names:

====================  ========= =========================
Models (``alkali``)   Species   References
====================  ========= =========================
sodium_vdw            Na        [Burrows2000]_
potassium_vdw         K         [Burrows2000]_
====================  ========= =========================

This implementation adopts the line parameters from the VALD database
[Piskunov1995]_ and collisional-broadening half-width from [Iro2005]_.

.. _rayleigh_opacity:

Rayleigh Opacity
~~~~~~~~~~~~~~~~

The ``rayleigh`` key sets Rayleigh opacity models.  The following
table lists the available Rayleigh model names:

===================== ======= ============================  ===
Models (``rayleigh``) Species Parameter names               References
===================== ======= ============================  ===
lecavelier            ---     ``log_k_ray``, ``alpha_ray``  [Lecavelier2008]_
dalgarno_H            H       ---                           [Dalgarno1962]_
dalgarno_He           He      ---                           [Kurucz1970]_
dalgarno_H2           |H2|    ---                           [Kurucz1970]_
===================== ======= ============================  ===

The Dalgarno Rayleigh models are tailored for H, He, and |H2| species,
and thus are not parametric.  The Lecavelier Rayleigh model is more flexible
and allows the user to modify the absorption strength and wavelength
dependency according to:

.. math::
    k(\lambda) = \kappa_{\rm ray} \kappa_0 \left(\frac{\lambda}{\lambda_0}\right)^{\alpha_{\rm ray}},

where :math:`\lambda_0=0.35` um and :math:`\kappa_0=5.31 \times
10^{-27}` cm\ :sup:`2` molecule\ :sup:`-1` are constants, and
:math:`\log(\kappa_{\rm ray})` and :math:`\alpha_{\rm ray}` are
fitting parameters that can be set through the ``rpars`` key.
Adopting values of :math:`\log(\kappa_{\rm ray})=0.0` and
:math:`\alpha_{\rm ray}=-4` reduces the Rayleigh opacity to that
expected for the |H2| molecule.

.. note:: Be aware that the implementation of the Lecavelier model
          uses the |H2| number-density profile (:math:`n_{\rm H2}`, in
          molecules cm\ :sup:`-3`) to compute the extinction
          coefficient (in |kayser| units) for a given atmospheric
          model: :math:`e(\lambda) = k(\lambda)\ n_{\rm H2}`.  We do
          this, because we are mostly interested in |H2|-dominated
          atmospheres, and most people consider a nearly constant |H2|
          profile.  Obviously, this needs to be fixed at some point in
          the future for a more general use.

.. _cloud_opacity:


Cloud Models
~~~~~~~~~~~~

Use the ``clouds`` key to include aerosol/haze/cloud opacities.
Currently, the code provides simple gray cloud models (listed below),
but soon we will include more complex Mie-scattering clouds for use in
forward- and retrieval modeling.  The following table lists the
currently available cloud model names:


And these are the available haze/cloud models (``clouds`` parameter):

=================== ============================================ =============================
Models (``clouds``) Parameter names                              Description
=================== ============================================ =============================
deck                ``log_p_cl``                                 Opaque gray cloud deck
ccsgray             ``log_k_gray``, ``log_p_top``, ``log_p_bot`` Constant gray cross-section
=================== ============================================ =============================

Use the ``cpars`` key to set the cloud model parameters.  The '*deck*'
model makes the atmosphere instantly opaque at the :math:`\log(p_{\rm cl})` pressure
(in bar units).

The '*ccsgray*' model creates a constant cross-section opacity between
the :math:`p_{\rm t}` and :math:`p_{\rm b}` pressures (in bar units),
and the |f| parameter scaling the opacity as: :math:`k = \kappa_{\rm gray}\ \kappa_0`, with
:math:`\kappa_0=5.31 \times 10^{-27}` cm\ :sup:`2` molecule\
:sup:`-1`.  This model uses the |H2| number-density profile to compute
the extinction coefficient as: :math:`e(\lambda) = k\ n_{\rm H2}` (in
|kayser| units).  Since the |H2| mixing ratio is generally constant,
the '*ccsgray*' opacity will scale linearly with pressure over the
atmosphere.

.. For any of these type of models, the user can include multiple
   models, simply by concatenating multiple models (and parameters)
   one after the other in the config file.

Patchy Cloud/Hazes
~~~~~~~~~~~~~~~~~~

Set the ``fpatchy`` argument to compute transmission spectra from a
linear combination of a clear and cloudy/hazy spectra.  The
cloudy/hazy component will include the opacity defined by the
:ref:`cloud_opacity` and the ``lecavelier`` :ref:`rayleigh_opacity`.
For example, for a 45% cloudy / 55% clear atmosphere, set:

.. code-block:: python

  # Patchy fraction, value between [0--1]:
  fpatchy = 0.45


Flux dilution factor
~~~~~~~~~~~~~~~~~~~~

Set the ``f_dilution`` argument to set an flux dilution factor
[Taylor2020]_, with values between 0--1, which compensates for
emission from an inhomogeneous atmosphere.  The dilution factor
represents the fractional area of the hottest region on the planet
(assuming that the colder regions flux is negligible in comparison).

.. code-block:: python

  # Flux dilution factor, value between [0--1]:
  f_dilution = 0.85



Abundances Scaling
------------------

Use the ``molvars`` and ``molpars`` keys to modify the abundance of
certain atmospheric species.  There are two sets of options depending
on whether the atmosphere is modeled in chemical equilibrium or not.
The following tables describe the available models (when the
corresponding ``molpars`` is :math:`x`).

For runs with free abundances:

============= =========================================== =====
``molvars``   Scaling                                     Description
============= =========================================== =====
``log_mol``   :math:`\log_{10}{\rm VMR} = x`              Set log(VMR) of species 'mol' to given value (constant with altitude)
``scale_mol`` :math:`{\rm VMR}(p) = {\rm VMR}_0(p)\ 10^x` Scale existing VMR of species 'mol' abundance by given value
============= =========================================== =====

In this case the parameters modify directly the VMR of specific
species.  To preserve the sum of the VMR at 1.0 at each layer, the
code will adjust the values of custom '*bulk*' species defined using
the ``bulk`` key.  A good practice is to set here the dominant species
in an atmosphere (e.g., ``bulk = H2 He`` for primary atmospheres).  If there
is more than one '*bulk*' species, the code preserves the relative
VMRs ratios between the bulk species.

For example, the following configuration will set uniform mole mixing
fractions for |H2O| and CO of :math:`10^{-3}` and :math:`10^{-4}`,
respectively; and adjust the abundances of |H2| and He to
preserve a total mixing fraction of 1.0 at each layer:

.. code-block:: python

  molvars = log_H2O log_CO
  molpars = -3.0    -4.0
  bulk = H2 He


For runs in thermochemical equilibrium (``chemistry = tea``):

================ ================================= =====
``molvars``      Scaling                           Description
================ ================================= =====
``metal``        :math:`{\rm [M/H]} = x`           Set metallicity (dex units) of all metal species (everything except H and He)
``[X/H]``        :math:`{\rm [X/H]} = x`           Set metallicity (dex units) of element 'X' relative to solar (overrided metal)
``X/Y``          ...                               Set abundance of element 'X' relative to that of element 'Y' (note not in dex units)
================ ================================= =====


key, whereas :math:`Q_0` is the abundance of the given species taken
from the atmospheric file ``atmfile``.
Note that the user can specify as many scaling parameters as wished,
as long as there are corresponding values for these three keys
(``molvars``, ``molpars``).


For example, the following configuration will compute abundances in
thermochemical equilibrium assuming 30x solar abundances for carbon
and oxygen, and 10x solar for all other metals:


.. code-block:: python

  chemistry = tea
  molvars = metal [C/H] [O/H]
  molpars = 1.0   1.5    1.5


.. _spec_observations:

Observations
------------

Filter Pass-bands
-----------------

Use the ``filters`` key to set the path to instrument filter
pass-bands (see [link to formats?]).  These can be used to compute
band-integrated values for the transmission or eclipse-depth spectra.


Observed Data
-------------

Use the ``data`` and ``uncert`` keys to set values for observed
transit- or eclipse-depth values and their uncertainties,
respectively.  Logically, you want to set data values corresponding to
the ``filters`` pass-bands.

Use the ``dunits`` key to specify the units of the ``data`` and
``uncert`` values (default: ``dunits = none``).  Typical values are:
'*none*', '*percent*', or '*ppm*' (see :ref:`units` section).

.. note:: Note that the ``filters``, ``data``, and ``uncert`` keys are
          not strictly required for a spectrum run, but they will
          allow the code to plot these information if requested (see
          [link to plots]).

.. _spec_parameters:

Fitting parameters
------------------


Number of CPUs
~~~~~~~~~~~~~~

The ``ncpu`` key sets the number of CPUs to use when computing LBL
opacities or when running retrievals (default: ``ncpu = 1``).

Verbosity
~~~~~~~~~

The ``verb`` key sets the screen-output and logfile verbosity level.
Higher ``verb`` values will display increasingly levels of detail
according to the following table:

========  =====================
``verb``  Screen Outputs
========  =====================
<0        Errors
0         Warnings
1         Headlines
2         Details
3         Debug
========  =====================


Plane-parallel Hemispheric Integration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For eclipse geometry, the code computes the emergent intensity under
the plane-parallel approximation, and then it integrates (sums)
intensity spectra at different angles with respect to the normal to
model the emitted flux spectrum.  The ``raygrid`` sets the angles (in
degrees) where to evaluate these intensities (default: ``raygrid = 0
20 40 60 80``).  The user can set custom values for these angles as
long as: (1) the first value is zero (normal to the planet's
'surface'), (2) they lie in the [0,90) range, and (3) they are
increasing order.

Alternatively, the user can set the ``quadrature`` key to perform a
Gaussian-quadrature integration, where the ``quadrature`` value sets
number of Gaussian-quadrature points (in which case, ``raygrid`` will
be ignored).

.. Plots: logxticks
          yran

----------------------------------------------------------------------

.. _spec_demo:

Examples
--------

.. note:: Before running this example, make sure that you have
   generated the TLI file from the :ref:`tli_tutorial_example`,
   generated the atmospheric profiles from the
   :ref:`abundance_tutorial_example`, and download the configuration
   file shown above, e.g., with these shell commands:

   .. code-block:: shell

       tutorial_path=https://raw.githubusercontent.com/pcubillos/pyratbay/master/examples/tutorial
       wget $tutorial_path/spectrum_transmission.cfg


In an interactive run, a spectrum run returns a '*pyrat*' object that
contains all input, intermediate, and output variables used to compute
the spectrum.  The following Python script computes and plots a
transmission spectrum using the configuration file found at the top of
this tutorial:

.. code-block:: python

    import matplotlib.pyplot as plt
    plt.ion()

    import pyratbay as pb
    import pyratbay.constants as pc

    pyrat = pb.run('spectrum_transmission.cfg')

    # Plot the resulting spectrum:
    wl = 1.0 / (pyrat.spec.wn*pc.um)
    depth = pyrat.spec.spectrum / pc.percent
    wl_ticks = [0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 5.0]

    plt.figure(-3, (7,4))
    plt.clf()
    ax = plt.subplot(111)
    plt.semilogx(wl, depth, "-", color='orange', lw=1.0)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticks(wl_ticks)
    plt.xlim(0.3, 5.0)
    plt.ylabel("Transit depth (Rp/Rs)$^2$ (%)")
    plt.xlabel("Wavelength (um)")

    # Or, alternatively:
    ax = pyrat.plot_spectrum()

And the results should look like this:

.. image:: ./figures/pyrat_transmission-spectrum_tutorial.png
    :width: 70%
    :align: center


.. note:: Note that although the user can define most input units,
          nearly all variables are stored in CGS units in the
          '*pyrat*' object.

The '*pyrat*' object is modular, and implements several convenience
methods to plot and display its content, as in the following example:

.. code-block:: python

    # pyrat object's string representation:
    print(pyrat)

    # String representation of the spectral variables:
    print(pyrat.spec)
