.. |H2O| replace:: H\ :sub:`2`\ O
.. |CO2| replace:: CO\ :sub:`2`
.. |CH4| replace:: CH\ :sub:`4`
.. |H2|  replace:: H\ :sub:`2`

.. _pttutorial:

Pressure-Temperature Tutorial
=============================

This mode creates a 1D set of pressure-temperature layers.  The
pressure array is equi-spaced in log-pressure.  This mode produces a
pdf image of the pressure-temperature profile and it returns the
pressure and temperature arrays.

The temperature model (``tmodel``) can be isothermal,
three-channel Eddington approximation (TCEA) [Line2013]_, or the Madhusudhan parameterized model for thermally inverted (MadhuInv) or non-inverted (MadhuNoInv) atmospheres [Madhusudhan2009]_.  The
number of model parameter (``tparams``) and other system parameters
depend on the temperature model.

Here is an example of a PT configuration file:

.. code-block:: python

  [pyrat]

  # Pyrat Bay run mode, select from: [tli pt atmosphere spectrum opacity mcmc]
  runmode = pt

  # Pressure array:
  punits  = bar    ; Default pressure units
  pbottom = 100.0  ; Bottom-layer pressure  (default units: punits)
  ptop    = 1e-5   ; Top-layer pressure (default units: punits)
  nlayers = 100    ; Number of atmospheric layers

  # Temperature-profile model, select from: [isothermal TCEA MadhuInv MadhuNoInv]
  tmodel  = isothermal
  tparams = 1500.0
  #    log10(kappa) log10(g1) log10(g2) alpha beta
  tparams = -3.0    -0.25     0.0       0.0   1.0

  # System parameters:
  radunits = km
  rstar    = 1.27 rsun  ; Stellar radius (default units: radunits)
  tstar    = 5800.0     ; Stellar effective temperature in K
  smaxis   = 0.045 au   ; Semi-major axis (default units: radunits)
  gplanet  = 800.0      ; Planetary surface gravity in cm s-2
  tint     = 100.0      ; Planetary internal temperature in K

  # Verbosity level [1--5]:
  verb = 4

The isothermal model has one free parameter: the temperature.
The TCEA model has five parameters: :math:`\log_{10}(\kappa),
\log_{10}(\gamma1), \log_{10}(\gamma2), \alpha, \beta` as defined in
[Line2013]_.  The TCEA model also requires the stellar radius
(``rstar``), the orbital semi-major axis (``smaxis``), the planetary
surface gravity (``gplanet``), the stellar effective temperature
(``tstar``), and the planetary internal temperature (``tint``).
The MadhuInv model has six parameters: :math:`a_1, a_2, p_1, p_2, p_3,
T_3` as defined in [Madhusudhan2009]_. The MadhuNoInv model has five
parameters: :math:`a_1, a_2, p_1, p_3, T_3` as defined in
[Madhusudhan2009]_.

To create an isothermal pressure-temperature profile run from the
Python interpreter:

.. code-block:: python

  # Generate an isothermal PT profile (output values in CGS units):
  pressure, T_isothermal = pb.pbay.run("tutorial_pt-isothermal.cfg")
  # Generate a TCEA PT profile:
  pressure, T_tcea = pb.pbay.run("tutorial_pt-tcea.cfg")

Note that the only difference between these configuration files are the
``tmodel`` and ``tparams`` varables.

Plot the profiles:

.. code-block:: python

  # Plot the PT profiles:
  plt.figure(-1)
  plt.clf()
  plt.semilogy(T_isothermal, pressure/pb.constants.bar, color="b",
               lw=2, label='Isothermal')
  plt.semilogy(T_tcea, pressure/pb.constants.bar, color="r",
               lw=2, label='TCEA')
  plt.ylim(100, 1e-5)
  plt.legend(loc="best")
  plt.xlabel("Temperature  (K)")
  plt.ylabel("Pressure  (bar)")
  plt.savefig("pyrat_PT_tutorial.pdf")

.. image:: ./figures/pyrat_PT_tutorial.png
   :width: 70%
   :align: center

.. note:: If any of the required variables is missing form the
          configuration file, ``Pyrat Bay`` will throw an error
          indicating the missing value, and **stop executing the
          run.**

.. note:: Similarly, ``Pyrat Bay`` will throw a warning for a missing
          variable that was defaulted, and **continue executing the run.**

