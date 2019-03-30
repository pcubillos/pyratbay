.. |H2O| replace:: H\ :sub:`2`\ O
.. |CO2| replace:: CO\ :sub:`2`
.. |CH4| replace:: CH\ :sub:`4`
.. |H2|  replace:: H\ :sub:`2`

.. _atmapi:

Atmosphere API
==============

This mode generates a 1D atmospheric model (pressure, temperature,
abundances).  So far, ``Pyrat Bay`` implements uniform- and
thermochemical-equilibrium-abundance profiles (through the ``TEA`` sub
module).  In the interactive run, the code returns the pressure,
temperature, and the 2D array of abundances.

The configuration file for this mode only has a few extra parameters
in addition of the PT mode:

.. code-block:: python

  [pyrat]

  # Pyrat Bay run mode, select from: [tli pt atmosphere spectrum opacity mcmc]
  runmode = atmosphere
  ...
  # Atmospheric model:
  atmfile  = WASP-00b.atm            ; Input/output atmospheric file
  elements = H He C N O Na           ; Input elemental composition
  species  = H2 He Na H2O CH4 CO CO2 ; Output species composition
  xsolar   = 1.0                     ; Solar-metallicity scaling factor
  #uniform  = 0.85 0.149 3e-6 4e-4 1e-4 4e-4 1e-7 ; Uniform abundances

``atmfile`` sets the output atmospheric file. ``species`` determines
the species present in the atmosphere.

To decide between a uniform or a TEA model, include or exclude the
``uniform`` variable, respectively.  The ``uniform`` values set the
abundances of each species in the ``species`` list, respectively.

A TEA run computes the abundances from a given elemental
solar-abundances list (``elements``).  The ``xsolar`` variable allows
the user to scale the elemental metallic abundances (everything but H
and He).

To generate the atmospheric model, run from the Python interpreter:

.. code-block:: python

  # Generate a TEA atmospheric model:
  pressure, temperature, abundances = pb.pbay.run("tutorial_atmosphere-tea.cfg")
  # Generate a uniform-abundance atmospheric model:
  pressure, temperature, abundances = pb.pbay.run("tutorial_atmosphere-uniform.cfg")

The ``atmosphere`` subpackage offers the ``readatm`` function to read an
atmospheric model.

.. code-block:: python

  # Read the atmospheric files:
  spec, press, temp, q_tea     = pb.atmosphere.readatm("WASP-00b.atm")
  spec, press, temp, q_uniform = pb.atmosphere.readatm("WASP-00c.atm")

  # Plot the results:
  plt.figure(-2)
  plt.clf()
  ax = plt.subplot(211)
  for i in np.arange(len(spec)):
    plt.loglog(q_tea[:,i], press, label=spec[i], lw=2)

  plt.ylim(np.amax(press), np.amin(press))
  plt.xlim(1e-10, 1.0)
  plt.legend(loc='best', fontsize=11)
  plt.ylabel("Pressure  (bar)")
  ax = plt.subplot(212)
  for i in np.arange(len(spec)):
    plt.loglog(q_uniform[:,i], press, label=spec[i], lw=2)

  plt.ylim(np.amax(press), np.amin(press))
  plt.xlim(1e-10, 1.0)
  plt.xlabel("Mole mixing fraction")
  plt.ylabel("Pressure  (bar)")
  plt.savefig("pyrat_atmosphere_tutorial.pdf")

.. image:: ./figures/pyrat_atmosphere_tutorial.png
   :width: 70%
   :align: center


