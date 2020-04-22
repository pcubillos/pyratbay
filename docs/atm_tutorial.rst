
.. _atmosphereapi:

Atmosphere Tutorial
===================

This run mode generates a 1D atmospheric model (pressure, temperature,
abundances), and saves it to file.  Currently, ``Pyrat Bay``
implements two abundance models: uniform- and
thermochemical-equilibrium-abundance (TEA) profiles.


Uniform Abundances
------------------

To produce a uniform-abundance model, the configuration file must
contain the ``species`` key specifying a list of the name of the
species to include in the atmosphere, and the ``uniform`` key
specifying the mole mixing fraction for each of the species listed in
``species``.  An atmosphere config file must also set the ``atmfile``
key specifying the output atmospheric file name.

Here is an example of a uniform atmosphere configuration file:

.. literalinclude:: ../examples/tutorial/atmosphere_uniform.cfg


Thermochemical Equilibrium Abundances
-------------------------------------

``Pyrat Bay`` computes abundances in thermochemical equilibrium via
the TEA package ([Blecic2016]_), by minimizing the Gibbs free energy at
each layer.  To produce a TEA model, the configuration file must
contain the ``species`` key specifying the species to include in the
atmosphere, the ``elements`` key specifying the elemental composition.
An atmosphere config file must also set the ``atmfile`` key specifying
the output atmospheric file name.

The TEA run assumes a solar elemental composition from [Asplund2009]_;
however, the user can enhance the metallicity of metals by setting the
``xsolar`` key, or can input custom elemental abundances by setting
the ``solar`` key with the path to a file containing the elemental
compositions (must follow the format of `this file
<https://github.com/pcubillos/pyratbay/blob/master/inputs/AsplundEtal2009.txt>`_).

Here is an example of a thermochemical-equilibrium atmosphere
configuration file:

.. literalinclude:: ../examples/tutorial/atmosphere_tea.cfg

----------------------------------------------------------------------

Examples
--------

In an interactive run, an atmosphere run returns the 1D pressure and
temperature profiles (CGS units), and the 2D abundance profiles (mole
mixing fractions).  The following Python script creates (and plots) a
uniform and a TEA atmospheric model (using the parameters shown in the
previous sections):

.. code-block:: python

    import matplotlib.pyplot as plt
    plt.ion()

    import pyratbay as pb
    import pyratbay.atmosphere as pa
    import pyratbay.io as io

    # Generate a uniform and a thermochemical-equilibrium atmospheric model:
    pressure, temperature, abundances = pb.run("atmosphere_tea.cfg")
    pressure, temperature, abundances = pb.run("atmosphere_uniform.cfg")

    # Read the atmospheric files:
    units, species, press, temp, q_tea, rad = io.read_atm("WASP-00b.atm")
    units, species, press, temp, q_uniform, rad = io.read_atm("WASP-00c.atm")

    # Plot the results:
    plt.figure(12)
    plt.clf()
    ax = plt.subplot(211)
    for q, spec in zip(q_tea.T, species):
        plt.loglog(q, press, label=spec, lw=2)

    plt.ylim(np.amax(press), np.amin(press))
    plt.xlim(1e-9, 1.0)
    plt.ylabel("Pressure (bar)", fontsize=11)
    ax = plt.subplot(212)
    for q, spec in zip(q_uniform.T, species):
        plt.loglog(q, press, label=spec, lw=2)

    plt.ylim(np.amax(press), np.amin(press))
    plt.xlim(1e-9, 1.0)
    plt.xlabel("Volume mixing fraction", fontsize=11)
    plt.ylabel("Pressure (bar)", fontsize=11)
    plt.legend(loc='best', fontsize=10)
    plt.tight_layout()

And the results should look like this:

.. image:: ./figures/pyrat_atmosphere_tutorial.png
   :width: 70%
   :align: center
