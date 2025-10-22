.. include:: _substitutions.rst

.. _getstarted:

Getting Started
===============

``Pyrat Bay`` is a multi-purpose package that enables the modeling of
exoplanet atmospheres and their spectra.  These can involve several
different physical processes.  The following table summarizes the
modeling capabilities enabled by ``Pyrat Bay``:

.. list-table:: 
   :header-rows: 1
   :widths: 10, 35, 15

   * - Calculation
     - Description
     - Output
   * - :doc:`line_sampling`
     - Sample line-transition data (Exomol, HITEMP) into cross-section spectra at a
       fixed grid of pressures, temperatures, and wavenumbers
     - Cross sections tables
   * - :doc:`atmosphere_modeling`
     - Compute 1D temperature, volume-mixing ratios, and radius
       profiles of a planetary atmosphere as function of pressure
     - Atmospheric :math:`T(p)`, :math:`{\rm VMR}(p)`, and :math:`r(p)` profiles
   * - :doc:`spectral_synthesis`
     - Radiative-transfer calculations given an input exoplanet atmosphere
     - Transit-depth, eclipse-depth, and/or emission spectra
   * - :doc:`atmospheric_retrievals`
     - Given an exoplanet parametric model and a spectroscopic observation,
       infer the exoplanet atmospheric properties
     - Posterior distribution of planetary model parameters
   * - :doc:`radiative_equilibrium`
     - Radiative-transfer calculations across an exoplanet atmosphere
     - Equilibrium :math:`T(p)` and :math:`{\rm VMR}(p)`, emission spectra


Any of these steps can be run from the command line or interactively
(though the Python Interpreter or a Jupyter Notebook).  To streamline
execution, ``Pyrat Bay`` provides a single command to execute any of
these runs.  To set up any of these runs, ``Pyrat Bay`` uses
configuration files following the standard `INI format
<https://docs.python.org/3/library/configparser.html#supported-ini-file-structure>`_.

The :ref:`qexample` section below demonstrates a simple
forward-modeling spectrum run.  The following sections provide a 
detailed for each of the running modes.  Finally, most of the low- and
mid-level routines used for these calculations are available
through the ``Pyrat Bay`` sub modules (see :ref:`API`).




---------------------------------------------------------------------

.. _install:

Installation
------------

To install ``Pyrat Bay`` run the following command from the terminal:

.. code-block:: shell

    pip install pyratbay

.. Or if you prefer conda:
   code-block:: shell
   conda install -c conda-forge "pyratbay>=2.0.0b4"


Alternatively (e.g., for developers), clone the repository to your local machine with the following terminal commands:

.. code-block:: shell

    git clone https://github.com/pcubillos/pyratbay
    cd pyratbay
    pip install -e .


``Pyrat Bay`` (version 2.0+) has been extensively tested to work on
Unix/Linux and OS X machines and is available for Python 3.9+.

---------------------------------------------------------------------

.. _qexample:

Quick Example
-------------

The following command-line scripts show how to calculate transmission
and eclipse spectra for an exoplanet atmosphere between 0.4 and 8.0
um.  First, create a directory to place input and output files, e.g.:

.. code-block:: shell

   mkdir run_demo
   cd run_demo

Download the H2O line-list database from the HITRAN server and unzip it:

.. code-block:: shell

   # Download HITRAN H2O line list
   wget https://www.cfa.harvard.edu/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Compressed-files/01_hit12.zip
   unzip 01_hit12.zip


Now download the ``Pyrat Bay`` configuration files here below:

.. raw:: html

   <details>
   <summary>Click here to show/hide: tutorial_tli_hitran_H2O.cfg</summary>

.. literalinclude:: ./_static/data/tutorial_tli_hitran_H2O.cfg
    :caption: File: `tutorial_tli_hitran_H2O.cfg <./_static/data/tutorial_tli_hitran_H2O.cfg>`__

.. raw:: html

   </details>

   <details>
   <summary>Click here to show/hide: tutorial_spectrum_transmission.cfg</summary>

.. literalinclude:: ./_static/data/tutorial_spectrum_transmission.cfg
    :caption: File: `tutorial_spectrum_transmission.cfg <./_static/data/tutorial_spectrum_transmission.cfg>`__

.. raw:: html

   </details>

   <details>
   <summary>Click here to show/hide: tutorial_spectrum_eclipse.cfg</summary>

.. literalinclude:: ./_static/data/tutorial_spectrum_eclipse.cfg
    :caption: File: `tutorial_spectrum_eclipse.cfg <./_static/data/tutorial_spectrum_eclipse.cfg>`__

.. raw:: html

   </details>


Execute these commands from the shell to create a
Transition-Line-Information file for |H2O|, and then to use it to
compute transmission and emission spectra:

.. code-block:: shell

   # Format line-by-line opacity:
   pbay -c tutorial_tli_hitran_H2O.cfg

   # Compute transmission and emission spectra:
   pbay -c tutorial_spectrum_transmission.cfg
   pbay -c tutorial_spectrum_eclipse.cfg


Outputs
^^^^^^^

That's it, now let's see the results.  The screen outputs and any
warnings raised are saved into log files.  The output spectrum is
saved to a separate file, to see it, run this Python script (on
interactive mode, I suggest starting the session with ``ipython
--pylab``):

.. code-block:: python

  import pyratbay as pb
  import pyratbay.constants as pc
  import pyratbay.spectrum as ps
  import pyratbay.io as io
  import matplotlib
  import matplotlib.pyplot as plt
  plt.ion()


  wl, transmission = io.read_spectrum("transmission_spectrum_tutorial.dat", wn=False)
  wl, eclipse = io.read_spectrum("eclipse_spectrum_tutorial.dat", wn=False)

  bin_wl = ps.constant_resolution_spectrum(0.4, 8.0, resolution=200)
  bin_transit = ps.bin_spectrum(bin_wl, wl, transmission)
  bin_eclipse = ps.bin_spectrum(bin_wl, wl, eclipse)

  fig = plt.figure(0)
  plt.clf()
  fig.set_size_inches(7,5)
  plt.subplots_adjust(0.12, 0.1, 0.98, 0.95, hspace=0.15)
  ax = plt.subplot(211)
  plt.plot(wl, transmission/pc.percent, color="royalblue", label="transmission model", lw=1.0)
  plt.plot(bin_wl, bin_transit/pc.percent, "salmon", lw=1.5, label='binned model')
  plt.xscale('log')
  plt.ylabel('Transit depth (%)')
  ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
  ax.set_xticks([0.5, 0.7, 1.0, 2.0, 3.0, 4.0, 6.0])
  ax.tick_params(which='both', direction='in')
  plt.xlim(0.4, 8.0)
  plt.ylim(1.88, 2.18)
  plt.legend(loc="upper left")

  ax = plt.subplot(212)
  plt.plot(wl, eclipse/pc.ppm, "royalblue", label="eclipse model", lw=1.0)
  plt.plot(bin_wl, bin_eclipse/pc.ppm, "salmon", lw=1.5, label='binned model')
  plt.xscale('log')
  plt.xlabel(r"Wavelength  (um)")
  plt.ylabel(r"$F_{\rm p}/F_{\rm s}$ (ppm)")
  ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
  ax.set_xticks([0.5, 0.7, 1.0, 2.0, 3.0, 4.0, 6.0])
  ax.tick_params(which='both', direction='in')
  plt.xlim(0.4, 8.0)
  plt.ylim(0, 3200)
  plt.legend(loc="upper left")
  plt.draw()
  plt.savefig("pyrat_spectrum_demo.png", dpi=300)

The output figure should look like this:

.. image:: ./figures/pyrat_spectrum_demo.png
   :width: 70%
   :align: center

---------------------------------------------------------------------

Command-line runs
-----------------

As shown above, ``Pyrat Bay`` enables a command-line entry point to
execute any of the runs listed above:

.. code-block:: shell

    pbay -c config_file.cfg

The configuration file determines what run mode to execute by setting
the ``runmode`` key.  Each of these modes have different
required/optional keys, which are detailed in further sections.

This same entry point offers a couple of secondary processes (display
version, re-format files). To display these options, run:

.. code-block:: shell

    pbay -h


Interactive runs
----------------

The same process can be executed from the Python Interpreter or in a
Jupyter Notebook:

.. code-block:: python

    import pyratbay as pb
    pyrat = pb.run('tutorial_spectrum_transmission.cfg')
    ax = pyrat.plot_spectrum()

The output vary depending on the selected run mode.  Additional low-
and mid-level routines are also available through this package (see
the :ref:`API`).

------------------------------------------------------------------------

In the following sections you can find a more detailed description and
examples of how to run ``Pyrat Bay`` for each available configuration.
