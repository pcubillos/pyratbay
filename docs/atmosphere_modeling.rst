.. |H2O| replace:: H\ :sub:`2`\ O
.. |CO2| replace:: CO\ :sub:`2`
.. |CH4| replace:: CH\ :sub:`4`
.. |H2|  replace:: H\ :sub:`2`

.. |kappa|  replace:: :math:`\log_{10}(\kappa')`
.. |gamma1| replace:: :math:`\log_{10}(\gamma_1)`
.. |gamma2| replace:: :math:`\log_{10}(\gamma_2)`
.. |alpha|  replace:: :math:`\alpha`
.. |beta|   replace:: :math:`\beta`
.. |Tirr|   replace:: :math:`T_{\rm irr}`
.. |Tint|   replace:: :math:`T_{\rm int}`

.. |logp1| replace:: :math:`\log_{10}(p_1)`
.. |logp2| replace:: :math:`\log_{10}(p_2)`
.. |logp3| replace:: :math:`\log_{10}(p_3)`
.. |a1|    replace:: :math:`a_1`
.. |a2|    replace:: :math:`a_2`
.. |T0|    replace:: :math:`T_0`


.. _atmospheretutorial:

Atmosphere Modeling
===================

This documentation shows how to model 1D planetary atmospheres.  There
are four properties that can be modeled:

1. :ref:`pressure`
2. `Temperature profile <#temperature-profiles>`_
3. `Abundance profiles (volume mixing ratios) <#abundance-profiles>`_
4. :ref:`radius`

Users can model these profiles 

--------------------------------------------------------

Regardless of which profiles are computed, in an interactive run the
code returns a five-element tuple containing the pressure profile
(bar), the temperature profile (Kelvin), the abundance profiles
(volume mixing fraction), the species names, and the altitude profile
(cm).  The outputs that were not calculated are set to ``None``.
Also, regardless of the input units, the output variables will always
be in CGS units.

In the config file, the user can set the ``atmfile`` argument to
specify an input atmospheric file from where to read pressure,
temperature, volume mixing ratios, and/or altitude profiles.  The
``output_atmfile`` argument instead can be set to specify a file name
where to store the outputs.

.. _pressure:

Pressure profile
----------------

``pyratbay`` offers the ``pa.pressure()`` function to compute pressure
profiles equi-spaced in log-pressure.  Users need to provide the the
pressure at the top of the atmosphere ``ptop``, at the bottom
``pbottom``, the number of layers ``nlayers``, and (optionally) the
units ``units``.  See :ref:`units` for a list of
available pressure units.  


.. tab-set::

  .. tab-item:: default units
     :selected:
  
     .. code-block:: python
    
       import pyratbay.atmosphere as pa
      
       # Generate pressure profile (default units are bars):
       press = pa.pressure(ptop=1e-8, pbottom=1e2, nlayers=61)

       print(press)

     Expected output:

     .. code-block:: none

      [1.00000000e-08 1.46779927e-08 2.15443469e-08 3.16227766e-08
       4.64158883e-08 6.81292069e-08 1.00000000e-07 1.46779927e-07
       2.15443469e-07 3.16227766e-07 4.64158883e-07 6.81292069e-07
       1.00000000e-06 1.46779927e-06 2.15443469e-06 3.16227766e-06
       4.64158883e-06 6.81292069e-06 1.00000000e-05 1.46779927e-05
       2.15443469e-05 3.16227766e-05 4.64158883e-05 6.81292069e-05
       1.00000000e-04 1.46779927e-04 2.15443469e-04 3.16227766e-04
       4.64158883e-04 6.81292069e-04 1.00000000e-03 1.46779927e-03
       2.15443469e-03 3.16227766e-03 4.64158883e-03 6.81292069e-03
       1.00000000e-02 1.46779927e-02 2.15443469e-02 3.16227766e-02
       4.64158883e-02 6.81292069e-02 1.00000000e-01 1.46779927e-01
       2.15443469e-01 3.16227766e-01 4.64158883e-01 6.81292069e-01
       1.00000000e+00 1.46779927e+00 2.15443469e+00 3.16227766e+00
       4.64158883e+00 6.81292069e+00 1.00000000e+01 1.46779927e+01
       2.15443469e+01 3.16227766e+01 4.64158883e+01 6.81292069e+01
       1.00000000e+02]


  .. tab-item:: string inputs

     .. code-block:: python
    
       import pyratbay.atmosphere as pa
      
       # Generate pressure profile (specify units in pressure boundaries):
       press = pa.pressure(ptop='1e-8 bar', pbottom='1e2 bar', nlayers=61)

       print(press)

     Expected output:

     .. code-block:: none

      [1.00000000e-08 1.46779927e-08 2.15443469e-08 3.16227766e-08
       4.64158883e-08 6.81292069e-08 1.00000000e-07 1.46779927e-07
       2.15443469e-07 3.16227766e-07 4.64158883e-07 6.81292069e-07
       1.00000000e-06 1.46779927e-06 2.15443469e-06 3.16227766e-06
       4.64158883e-06 6.81292069e-06 1.00000000e-05 1.46779927e-05
       2.15443469e-05 3.16227766e-05 4.64158883e-05 6.81292069e-05
       1.00000000e-04 1.46779927e-04 2.15443469e-04 3.16227766e-04
       4.64158883e-04 6.81292069e-04 1.00000000e-03 1.46779927e-03
       2.15443469e-03 3.16227766e-03 4.64158883e-03 6.81292069e-03
       1.00000000e-02 1.46779927e-02 2.15443469e-02 3.16227766e-02
       4.64158883e-02 6.81292069e-02 1.00000000e-01 1.46779927e-01
       2.15443469e-01 3.16227766e-01 4.64158883e-01 6.81292069e-01
       1.00000000e+00 1.46779927e+00 2.15443469e+00 3.16227766e+00
       4.64158883e+00 6.81292069e+00 1.00000000e+01 1.46779927e+01
       2.15443469e+01 3.16227766e+01 4.64158883e+01 6.81292069e+01
       1.00000000e+02]


  .. tab-item:: units argument

     .. code-block:: python
          
       import pyratbay.atmosphere as pa
      
       # Generate pressure profile (specify units):
       press = pa.pressure(ptop=1e-8, pbottom=1e2, units='bar', nlayers=61)

       print(press)

     Expected output:

     .. code-block:: none

      [1.00000000e-08 1.46779927e-08 2.15443469e-08 3.16227766e-08
       4.64158883e-08 6.81292069e-08 1.00000000e-07 1.46779927e-07
       2.15443469e-07 3.16227766e-07 4.64158883e-07 6.81292069e-07
       1.00000000e-06 1.46779927e-06 2.15443469e-06 3.16227766e-06
       4.64158883e-06 6.81292069e-06 1.00000000e-05 1.46779927e-05
       2.15443469e-05 3.16227766e-05 4.64158883e-05 6.81292069e-05
       1.00000000e-04 1.46779927e-04 2.15443469e-04 3.16227766e-04
       4.64158883e-04 6.81292069e-04 1.00000000e-03 1.46779927e-03
       2.15443469e-03 3.16227766e-03 4.64158883e-03 6.81292069e-03
       1.00000000e-02 1.46779927e-02 2.15443469e-02 3.16227766e-02
       4.64158883e-02 6.81292069e-02 1.00000000e-01 1.46779927e-01
       2.15443469e-01 3.16227766e-01 4.64158883e-01 6.81292069e-01
       1.00000000e+00 1.46779927e+00 2.15443469e+00 3.16227766e+00
       4.64158883e+00 6.81292069e+00 1.00000000e+01 1.46779927e+01
       2.15443469e+01 3.16227766e+01 4.64158883e+01 6.81292069e+01
       1.00000000e+02]



.. _temp_profile:

Temperature profiles
--------------------

Currently, there are three available temperature models:

..   :widths: 7, 20, 25

.. list-table:: Temperature profile models
   :header-rows: 1

   * - Model name
     - Parameter names
     - References
   * - ``isothermal``
     - ``T_iso``
     - ---
   * - ``guillot``
     - ``log_kappa'``, ``log_gamma1``, ``log_gamma2``, ``alpha``, ``T_irr``, ``T_int``
     - [Line2013]_
   * - ``madhu``
     - ``log_p1``, ``log_p2``, ``log_p3``, ``a1``, ``a2``, ``T0``
     - [Madhusudhan2009]_


Any of these models can be used either as stand-alone functions or via
the ``pb.run()`` function with a configuration file.


Interactive notebooks
^^^^^^^^^^^^^^^^^^^^^

This Notebook explains the model parameters and shows how to use the
temperature models in a Python script:

- `Temperature profiles tutorial <cookbooks/temperature_profiles.ipynb>`__


Configuration files
^^^^^^^^^^^^^^^^^^^

Alternatively, temperature profiles can be generated from
configuration files, which can be run from the command line or from
interactive Python sessions. Here are examples for each of the models:

.. tab-set::

  .. tab-item:: isothermal
     :selected:

     .. raw:: html
     
        <details>
        <summary>Click here to show/hide: temperature_profile_isothermal.cfg</summary>
     
     .. literalinclude:: ./_static/data/temperature_profile_isothermal.cfg
         :caption: File: `temperature_profile_isothermal.cfg <./_static/data/temperature_profile_isothermal.cfg>`__

     .. raw:: html
     
        </details>

     Copy this configuration file to your local folder.  Then users
     can generate temperature profiles, either from an interactive
     python session, as in the following script:

     .. code-block:: python
          
       import pyratbay as pb

       # Generate an atmosphere object with the profiles:
       atm = pb.run("temperature_profile_isothermal.cfg")

       # The atm object contains the temperature profile, among other properties:
       # (e.g., see also atm.press for the pressure array)
       print(atm.temp)

     .. code-block:: none

       [990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990.
        990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990.
        990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990.
        990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990. 990.
        990. 990. 990. 990. 990.]

     Which will create an output .atm file with the
     pressure-temperature profile (same root file name as ``logfile``
     in the configuration file).  Alternatively, users can execute
     this script from the command line:

     .. code-block:: shell

       pbay -c temperature_profile_isothermal.cfg


  .. tab-item:: guillot

     .. raw:: html
     
        <details>
        <summary>Click here to show/hide: temperature_profile_guillot.cfg</summary>
     
     .. literalinclude:: ./_static/data/temperature_profile_guillot.cfg
         :caption: File: `temperature_profile_madhu.cfg <./_static/data/temperature_profile_guillot.cfg>`__
     
     .. raw:: html
     
        </details>

     Copy this configuration file to your local folder.  Then users
     can generate temperature profiles, either from an interactive
     python session, as in the following script:

     .. code-block:: python
          
       import pyratbay as pb

       # Generate an atmosphere object with the profiles:
       atm = pb.run("temperature_profile_guillot.cfg")

       # The atm object contains the temperature profile, among other properties:
       # (e.g., see also atm.press for the pressure array)
       print(atm.temp)

     .. code-block:: none

       [ 893.13809757  893.13809472  893.13809066  893.13808487  893.13807664
         893.13806493  893.13804829  893.13802468  893.13799121  893.13794384
         893.13787688  893.13778235  893.13764913  893.13746171  893.13719851
         893.13682967  893.13631394  893.13559461  893.13459405  893.13320653
         893.13128901  893.1286492   893.12503098  893.12009655  893.11340618
         893.10439669  893.092362    893.07644255  893.05563579  893.02884952
         892.99503549  892.95346546  892.90425101  892.84926945  892.79374973
         892.74890619  892.7361951   892.79400726  892.98786822  893.42537591
         894.27686781  895.80153311  898.37522069  902.5089626   908.83533653
         918.02632402  930.60377954  946.63515057  965.39004504  985.14468916
        1003.38444149 1017.57191204 1026.32163644 1030.23462361 1031.37402834
        1031.6122709  1031.73878338 1031.91095189 1032.16328573 1032.53332577
        1033.07575077]

     Which will create an output .atm file with the
     pressure-temperature profile (same root file name as ``logfile``
     in the configuration file).  Alternatively, users can execute
     this script from the command line:

     .. code-block:: shell

       pbay -c temperature_profile_guillot.cfg


  .. tab-item:: madhu

     .. raw:: html
     
        <details>
        <summary>Click here to show/hide: temperature_profile_madhu.cfg</summary>
     
     .. literalinclude:: ./_static/data/temperature_profile_madhu.cfg
         :caption: File: `temperature_profile_madhu.cfg <./_static/data/temperature_profile_madhu.cfg>`__
     
     .. raw:: html
     
        </details>

     Copy this configuration file to your local folder.  Then users
     can generate temperature profiles, either from an interactive
     python session, as in the following script:

     .. code-block:: python
          
       import pyratbay as pb

       # Generate an atmosphere object with the profiles:
       atm = pb.run("temperature_profile_madhu.cfg")

       # The atm object contains the temperature profile, among other properties:
       # (e.g., see also atm.press for the pressure array)
       print(atm.temp)

     .. code-block:: none

       [ 850.31978367  850.67004908  851.24532311  852.09429012  853.24712497
         854.71853396  856.51414567  858.63564929  861.08344173  863.85759587
         866.95812108  870.38501736  874.13828472  878.21792316  882.62393267
         887.35631326  892.41506491  897.80018765  903.51168146  909.54954634
         915.9137823   922.60438934  929.62136744  936.96471663  944.63443689
         952.63052822  960.95299063  969.60182411  978.57702867  987.8786043
         997.50655101 1007.46086879 1017.74155765 1028.34861758 1039.28204858
        1050.54166067 1062.12637581 1074.0308543  1086.2350518  1098.68260322
        1111.25730855 1123.79165051 1136.13907285 1148.27829152 1160.35458503
        1172.60739618 1185.23773307 1198.28698302 1211.52304248 1224.33040902
        1235.73387588 1244.72660511 1250.79507726 1254.20443381 1255.76698116
        1256.34273652 1256.51152735 1256.5505651  1256.55759081 1256.55849942
        1256.55849942]

     Which will create an output .atm file with the
     pressure-temperature profile (same root file name as ``logfile``
     in the configuration file).  Alternatively, users can execute
     this script from the command line:

     .. code-block:: shell

       pbay -c temperature_profile_madhu.cfg


.. _abundance_profile:

Abundance Profiles
------------------

Currently, there are two models to set the base volume-mixing-ratio
chemistry (``chemistry`` argument): ``uniform`` or ``tea``.  Each one
requires a different set of arguments, which is described in the table
and sections below:

====================== ==================================================== ====
Models (``chemistry``) Required arguments [optional arguments]              References
====================== ==================================================== ====
uniform                ``species``, ``uniform``, [``vmr_vars``]             ---
tea                    ``species``, [``vmr_vars``]                          [Blecic2016]_
====================== ==================================================== ====


Uniform abundances
^^^^^^^^^^^^^^^^^^

To produce a uniform-abundance model, the configuration file must
contain the ``species`` key specifying a list of the name of the
species to include in the atmosphere, and the ``uniform`` key
specifying the mole mixing fraction for each of the species listed in
``species``.

Here is an example of a uniform atmosphere configuration file (`tutorial_atmosphere_uniform.cfg
<https://github.com/pcubillos/pyratbay/blob/master/examples/tutorial_atmosphere_uniform.cfg>`_):

.. literalinclude:: ../examples/tutorial_atmosphere_uniform.cfg


Thermochemical-equilibrium Abundances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Pyrat Bay`` computes thermochemical equilibrium abundances (TEA) via
the chemcat package, by minimizing the Gibbs free energy at each
layer.  To produce a TEA model, the configuration file must set
``chemistry=tea``.  The ``species`` argument sets the species to
include in the atmosphere.

The TEA run assumes a solar elemental composition from [Asplund2021]_
as the base for the thermochemical equilibrium model; however, the
user can customize the elemental abundances using the ``vmr_vars``
argument.  The table below shows the available options.

================= ===
``vmr_vars``      Notes
================= ===
``[M/H]``         Metallicity of all elemental species (dex units, with respect to solar)
``[X/H]``         Metallicity of element ``X`` (dex units, with respect to solar).  Overrides ``[M/H]``
``X/Y``           Abundance of element ``X`` relative to element ``Y``.  Overrides ``[M/H]`` for ``X``, but ``Y`` can be previously modified by ``[M/H]`` or ``[Y/H]``.
================= ===


.. ``log_mol``       log10(VMR) of species ``mol`` (constant with altitude)


Here is an example of a thermochemical-equilibrium atmosphere
configuration file (`tutorial_atmosphere_tea.cfg
<https://github.com/pcubillos/pyratbay/blob/master/examples/tutorial_atmosphere_tea.cfg>`_):

.. literalinclude:: ../examples/tutorial_atmosphere_tea.cfg

----------------------------------------------------------------------

.. _abundance_tutorial_example:

Abundance-profile Examples
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Before running this example, download the configuration
          files shown above, e.g., with these shell commands:

          .. code-block:: shell

              tutorial_path=https://raw.githubusercontent.com/pcubillos/pyratbay/master/examples
              wget $tutorial_path/tutorial_atmosphere_tea.cfg
              wget $tutorial_path/tutorial_atmosphere_uniform.cfg

The following Python script creates and plots the abundance
Aprofiles for the configuration files shown above:

.. code-block:: python

    import matplotlib.pyplot as plt
    plt.ion()

    import pyratbay as pb
    import pyratbay.plots as pp

    # Generate a uniform and a thermochemical-equilibrium atmospheric model:
    atm_tea = pb.run("tutorial_atmosphere_tea.cfg")
    atm_uni = pb.run("tutorial_atmosphere_uniform.cfg")

    # Plot the results:
    plt.figure(12, (6,5))
    plt.clf()
    ax1 = plt.subplot(211)
    ax1 = pp.abundance(
        atm_tea.vmr, atm_tea.press, atm_tea.species,
        colors='default', xlim=[1e-12, 3.0], legend_fs=8, ax=ax1,
    )
    ax2 = plt.subplot(212)
    ax2 = pp.abundance(
        atm_uni.vmr, atm_uni.press, atm_uni.species,
        colors='default', xlim=[1e-12, 3.0], legend_fs=0, ax=ax2,
    )
    plt.tight_layout()

And the results should look like this:

.. image:: ./figures/pyrat_atmosphere_tutorial.png
   :width: 70%
   :align: center

----------------------------------------------------------------------

.. _radius:

Radius Profile
--------------

If the user sets the ``radmodel`` key, the code will to compute the
atmospheric altitude profile (radius profile).  The currently
available models solve for the hydrostatic-equilibrium equation,
combined with the ideal gas law with a pressure-dependent gravity
(``radmodel=hydro_m``, recommended):

.. math::
   \frac{dr}{r^2} = -\frac{k_{\rm B}T}{\mu G M_p} \frac{dp}{p},

or a constant surface gravity (``radmodel=hydro_g``):

.. math::
   dr = -\frac{k_{\rm B}T}{\mu g} \frac{dp}{p},

where :math:`M_{\rm p}` is the mass of the planet, :math:`T(p)` is the
atmospheric temperature, :math:`\mu(p)` is the atmospheric mean
molecular mass, :math:`k_{\rm B}` is the Boltzmann constant, and
:math:`G` is the gravitational constant.  Note that :math:`T(p)` and
:math:`\mu(p)` are computed from the models of the :ref:`temp_profile`
and :ref:`abundance_profile`, respectively.



To obtain the particular solution of these differential equations,
the user needs to supply a pair of radius--pressure reference values
to define the boundary condition :math:`r(p_0) = R_0`.  The
``rplanet`` and ``refpressure`` keys set :math:`R_0` and :math:`p_0`,
respectively.
The ``mplanet`` and ``gplanet`` keys set the planetary mass (:math:`M_p`)
and surface gravity (:math:`g`) corresponding to ``rplanet``, respectively.

.. Note:: Note that the user needs only to define one variable between
    ``mplanet`` and ``gplanet`` since they are related
    through the equation: :math:`g(R_0) = G M_p / R_0^2`.  Whenever
    ``mplanet`` is set to a given value, ``gplanet`` is updated
    consistently, and vice versa.

Note that the selection of the :math:`\{p_0,R_0\}` pair is arbitrary.
A good practice is to choose values close to the transit radius of
the planet.  Although the pressure at the transit radius is a priori
unknown for a give particular case [Griffith2014]_, its value lies
at around 0.1 bar.

Here is an example of a hydrostatic-equilibrium atmosphere
configuration file (`atmosphere_hydro_m.cfg
<https://github.com/pcubillos/pyratbay/blob/master/examples/tutorial/atmosphere_hydro_m.cfg>`_):

.. literalinclude:: ../examples/tutorial/atmosphere_hydro_m.cfg


Altitude-profile Examples
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Before running this example, download the configuration
          files shown above, e.g., with these shell commands:

          .. code-block:: shell

              tutorial_path=https://raw.githubusercontent.com/pcubillos/pyratbay/master/examples/tutorial
              wget $tutorial_path/atmosphere_hydro_m.cfg
              wget $tutorial_path/atmosphere_hydro_g.cfg

The following Python script creates and plots the profiles
for the configuration file shown above:

.. code-block:: python

    import matplotlib.pyplot as plt
    plt.ion()

    import pyratbay as pb
    import pyratbay.constants as pc

    # Kepler-11c mass and radius:
    pressure, temp, q, species, radius = pb.run("atmosphere_hydro_m.cfg")
    pressure, temp, q, species, radius_g = pb.run("atmosphere_hydro_g.cfg")

    # Plot the results:
    plt.figure(12, (6,5))
    plt.clf()
    ax = plt.subplot(111)
    ax.semilogy(radius_g/pc.rearth, pressure/pc.bar, lw=2, c='navy', label='constant g')
    ax.semilogy(radius/pc.rearth, pressure/pc.bar, lw=2, c='orange', label='g = g(p)')
    ax.set_ylim(1e2, 1e-6)
    ax.set_xlabel(r'Radius $(R_{\oplus})$', fontsize=12)
    ax.set_ylabel('Pressure (bar)', fontsize=12)
    ax.tick_params(labelsize=11)
    ax.legend(loc='upper left', fontsize=12)
    plt.tight_layout()

And the results should look like this:

.. image:: ./figures/pyrat_hydrostatic_tutorial.png
   :width: 70%
   :align: center

