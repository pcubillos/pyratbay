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

Atmosphere Tutorial
===================

This run mode generates a 1D atmospheric model (pressure, temperature,
abundances, and altitude profiles), and saves it to a file.  At minimum,
the user needs to provide the arguments required to compute the
pressure and temperature profiles.  Further, ``Pyrat Bay`` will
compute volume-mixing ratio (abundance) profiles only if the user sets the
``chemistry`` provides the respective
arguments.  Likewise, the code will compute the altitude profiles only
if the user provides the respective arguments (which also require that
the abundance profiles are defined).

Regardless of which profiles are computed, in an interactive run, the
code returns a five-element tuple containing the pressure profile
(barye), the temperature profile (Kelvin degree), the abundance
profiles (volume mixing fraction), the species names, and the altitude
profile (cm).  The variables not calculated profiles will be ``None``.
Also, regardless of the input units, the output variables will always
be in CGS units.


Pressure Profile
----------------

The pressure profile model is an equi-spaced array in log-pressure,
determined by the pressure at the top of the atmosphere ``ptop``, at
the bottom of the atmosphere ``pbottom``, and the number of layers
``nlayers``.

The units for the ``ptop`` and ``pbottom`` pressures may be defined
in-place along with the variable value (e.g., ``pbottom = 100 bar``)
or through the ``punits`` key (in-place units take precedence over
``punits``).  See :ref:`units` for a list of available units.



Temperature Profiles
--------------------

Currently, there are three temperature models (``tmodel`` argument) to
compute temperature profiles.  Each one requires a different set of
parameters (``tpars``), which is described in the table and sections below:


=================== ==================================================== ====
Models (``tmodel``) Parameters (``tpars``)                               References
=================== ==================================================== ====
isothermal          :math:`T_0`                                          ---
tcea                |kappa|, |gamma1|, |gamma2|, |alpha|, |Tirr|, |Tint| [Line2013]_
madhu               |logp1|, |logp2|, |logp3|, |a1|, |a2|, |T0|          [Madhusudhan2009]_
=================== ==================================================== ====


Isothemal
^^^^^^^^^

The isothermal model is the simplest model, having just one free
parameter (``tpars``): the temperature (:math:`T_0`) at all layers.

Here is an example of an isothermal atmosphere configuration file:

.. literalinclude:: ../examples/tutorial/pt_isothermal.cfg


Three-channel Eddington Approximation (TCEA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The tcea model has six parameters: |kappa|, |gamma1|, |gamma2|,
|alpha|, |Tirr|, and |Tint| as defined in [Line2013]_, except for
:math:`\kappa'`, which corresponds to :math:`\kappa'
\equiv \kappa/g`.

Here is an example of a tcea atmosphere configuration file:

.. literalinclude:: ../examples/tutorial/pt_tcea.cfg


Madhu profiles
^^^^^^^^^^^^^^

The madhu model has six parameters: |logp1|, |logp2|, |logp3|, |a1|,
|a2|, and |T0|, as described in [Madhusudhan2009]_, where the pressure
values must be given in bars.  A thermally inverted
profile will result when :math:`p_1 < p_2`; a non-inverted profile
will result when :math:`p_2 < p_1`.  The pressure parameters must also
satisfy: :math:`p_1 < p_3`.

.. literalinclude:: ../examples/tutorial/pt_madhu.cfg


-------------------------------------------------------------------

Temperature-profile Examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following Python script creates and plots the pressure-temperature
profiles for the configuration files shown above:


.. code-block:: python

  import matplotlib.pyplot as plt
  plt.ion()

  import pyratbay as pb
  import pyratbay.constants as pc

  # Before running the rest of this script, you'll need to copy the
  # config files to your current folder, you can find them here, e.g.:
  # print(f'{pc.ROOT}examples/tutorial/pt_isothermal.cfg')

  # Generate PT profiles:
  press, t_iso   = pb.run("pt_isothermal.cfg")[0:2]
  press, t_tcea  = pb.run("pt_tcea.cfg")[0:2]
  press, t_madhu = pb.run("pt_madhu.cfg")[0:2]

  # Plot the PT profiles:
  plt.figure(11)
  plt.clf()
  plt.semilogy(t_iso,   press/pc.bar, color='b', lw=2, label='isothermal')
  plt.semilogy(t_tcea,  press/pc.bar, color='r', lw=2, label='tcea')
  plt.semilogy(t_madhu, press/pc.bar, color='g', lw=2, label='madhu')
  plt.ylim(100, 1e-5)
  plt.xlim(800, 1200)
  plt.legend(loc="best", fontsize=12)
  plt.xlabel("Temperature  (K)", fontsize=12)
  plt.ylabel("Pressure  (bar)", fontsize=12)
  plt.tight_layout()


And the results should look like this:

.. image:: ./figures/pyrat_PT_tutorial.png
   :width: 70%
   :align: center


.. note:: If any of the required variables is missing form the
          configuration file, ``Pyrat Bay`` will throw an error
          indicating the missing value, and **stop executing the
          run.** Similarly, ``Pyrat Bay`` will throw a warning for a
          missing variable that was defaulted, and **continue
          executing the run.**


----------------------------------------------------------------------

Abundance Profiles
------------------

Currently, there are two chemistry models (``chemistry`` argument) to
compute volume-mixing-ratio abundances: uniform or thermochemical
equilibrium.  Each one requires a different set of arguments, which is
described in the table and sections below:

====================== ==================================================== ====
Models (``chemistry``) Required arguments [optional arguments]              References
====================== ==================================================== ====
uniform                ``species``, ``uniform``                                     ---
tea                    ``species``, ``elements``, [``xsolar``, ``escale`` ] [Blecic2016]_
====================== ==================================================== ====


Uniform abundances
^^^^^^^^^^^^^^^^^^

To produce a uniform-abundance model, the configuration file must
contain the ``species`` key specifying a list of the name of the
species to include in the atmosphere, and the ``uniform`` key
specifying the mole mixing fraction for each of the species listed in
``species``.  An atmosphere config file must also set the ``atmfile``
key specifying the output atmospheric file name.

Here is an example of a uniform atmosphere configuration file:

.. literalinclude:: ../examples/tutorial/atmosphere_uniform.cfg


Thermochemical-equilibrium Abundances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Pyrat Bay`` computes abundances in thermochemical equilibrium via
the TEA package [Blecic2016]_, by minimizing the Gibbs free energy at
each layer.  To produce a TEA model, the configuration file must
contain the ``species`` key specifying the species to include in the
atmosphere, the ``elements`` key specifying the elemental composition.
An atmosphere config file must also set the ``atmfile`` key specifying
the output atmospheric file name.

The TEA run assumes a solar elemental composition from [Asplund2009]_;
however, the user can scale the abundance of metals by setting the
``xsolar`` key, or can scale individual elemental abundances by setting
the ``escale`` key element and scale factor pairs:

Here is an example of a thermochemical-equilibrium atmosphere
configuration file:

.. literalinclude:: ../examples/tutorial/atmosphere_tea.cfg

----------------------------------------------------------------------

Abundance-profile Examples
^^^^^^^^^^^^^^^^^^^^^^^^^^

The following Python script creates and plots the abundance 
Aprofiles for the configuration files shown above:

.. code-block:: python

    import matplotlib.pyplot as plt
    plt.ion()

    import pyratbay as pb
    import pyratbay.plots as pp

    # Before running the rest of this script, you'll need to copy the
    # config files to your current folder, you can find them here, e.g.:
    # print(f'{pb.constants.ROOT}examples/tutorial/atmosphere_tea.cfg')

    # Generate a uniform and a thermochemical-equilibrium atmospheric model:
    pressure, temp, q_tea, species, radius = pb.run("atmosphere_tea.cfg")
    pressure, temp, q_uniform, species, radius = pb.run("atmosphere_uniform.cfg")

    # Plot the results:
    plt.figure(12, (6,5))
    plt.clf()
    ax1 = pp.abundance(
        q_tea, pressure, species, colors='default', xlim=[1e-11, 10.0],
        legend_fs=0, ax=plt.subplot(211))
    ax2 = pp.abundance(
        q_uniform, pressure, species, colors='default', xlim=[1e-11, 10.0],
        legend_fs=8, ax=plt.subplot(212))
    plt.tight_layout()

And the results should look like this:

.. image:: ./figures/pyrat_atmosphere_tutorial.png
   :width: 70%
   :align: center

