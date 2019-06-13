.. |H2O| replace:: H\ :sub:`2`\ O
.. |CO2| replace:: CO\ :sub:`2`
.. |CH4| replace:: CH\ :sub:`4`
.. |H2|  replace:: H\ :sub:`2`

.. |kappa|  replace:: :math:`\log_{10}(\kappa)`
.. |gamma1| replace:: :math:`\log_{10}(\gamma_1)`
.. |gamma2| replace:: :math:`\log_{10}(\gamma_2)`
.. |alpha|  replace:: :math:`\alpha`
.. |beta|   replace:: :math:`\beta`


.. _pttutorial:

Pressure-Temperature Tutorial
=============================

This run mode creates a 1D profile of pressure-temperature layers.
Currently, there are four available temperature models that can be set
with the ``tmodel`` key. Each one of these require a different set of
parameters (``tpars``).  The models, parameters, and references are
listed in the following table:

=================== ============================================ ==== 
Models (``tmodel``) Parameters (``tpars``)                       References
=================== ============================================ ====
isothermal          :math:`T_0`                                  ---
tcea                |kappa|, |gamma1|, |gamma2|, |alpha|, |beta| [Line2013]_
madhu_inv           :math:`a_1, a_2, p_1, p_2, p_3, T_3`         [Madhusudhan2009]_
madhu_noinv         :math:`a_1, a_2, p_1, p_2, p_3, T_3`         [Madhusudhan2009]_
=================== ============================================ ====


Pressure Profile
----------------

The pressure profile is an equi-spaced array in log-pressure,
determined by the pressure at the top of the atmosphere ``ptop``, at
the bottom of the atmosphere ``pbottom``, and the number of layers
``nlayers``.

The units for the ``ptop`` and ``pbottom`` pressures may
be defined in place (as in the sample above) or may be defined with
the ``punits`` key (in-place units take precedence over ``punits``).
See :ref:`units` for a list of available units.


Temperature Profiles
--------------------

Isothemal
^^^^^^^^^

The isothermal model is the simplest model, having just one free
parameter (``tpars``): the temperature (:math:`T_0`) at all layers.

Here is an example of a PT configuration file:

.. literalinclude:: ../examples/tutorial/tutorial_pt-isothermal.cfg


Three-channel Eddington Approximation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The tcea model has five parameters: |kappa|,
|gamma1|, |gamma2|, |alpha|, and |beta| as defined in
[Line2013]_.  This model also requires the stellar radius
(``rstar``), the orbital semi-major axis (``smaxis``), the planetary
surface gravity (``gplanet``), the stellar effective temperature
(``tstar``), and the planetary internal temperature (``tint``).

.. literalinclude:: ../examples/tutorial/tutorial_pt-tcea.cfg

Note that the units for ``gplanet`` are cm s\ :sup:`-2`, and the units
for temperature keys (like ``tstar`` and ``tint``) are Kelvin.

.. note:: ``Pyrat Bay`` can compute the planetary surface gravity
          (``gplanet``) from the planetary mass (``mplanet``) and
          radius (``rplanet``).

Madhu profiles
^^^^^^^^^^^^^^

The madhu_inv model has six parameters: :math:`a_1, a_2, p_1, p_2, p_3,
T_3`, whereas the madhu_noinv model has five
parameters: :math:`a_1, a_2, p_1, p_3, T_3` as defined in
[Madhusudhan2009]_.  **[To be reviewed]**

-------------------------------------------------------------------

Examples
--------

The PT run mode returns a two-element tuple with the pressure and
temperature arrays (in CGS units).  The following Python script
creates (and plots) an isothermal and a TCEA pressure-temperature
profile, using the parameters shown in the previous sections:


.. code-block:: python

  import matplotlib.pyplot as plt
  plt.ion()

  import pyratbay as pb
  import pyratbay.constants as pc

  # Generate PT profiles:
  press, T_iso  = pb.run("tutorial_pt-isothermal.cfg")
  press, T_tcea = pb.run("tutorial_pt-tcea.cfg")

  # Plot the PT profiles:
  plt.figure(11)
  plt.clf()
  plt.semilogy(T_iso,  press/pc.bar, color='b', lw=2, label='Isothermal')
  plt.semilogy(T_tcea, press/pc.bar, color='r', lw=2, label='TCEA')
  plt.ylim(100, 1e-5)
  plt.xlim(1200, 1800)
  plt.legend(loc="best")
  plt.xlabel("Temperature  (K)")
  plt.ylabel("Pressure  (bar)")

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

