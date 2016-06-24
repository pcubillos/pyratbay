.. _getstarted:

Getting Started
===============

The Pyrat Bay package provides multiple sub packages to compute
radiative-transfer calculations for exoplanets:

The ``pyrat`` is the main package that provides the radiative-transfer
code that computes an emission or transmission spectrum for a given
atmospheric model.  The ``lineread`` package formats online-available
line-by-line opacity databases, used later by ``pyrat``.  The ``pbay``
package provides the retrieval framework (using a Markov-chain Monte
Carlo algorithm, MCMC) to model and constrain exoplanet atmospheres.

Additional packages provide specific function to read stellar spectra
(``starspec``); generate, read, and write 1D atmospheric models
(``atmosphere``), read and process waveband filters (``wine``),
provide universal and astrophysical constants (``constants``),
plotting (``plots``) and additional tools (``tools``).

Each one of these can be run either from the shell prompt (through the
executable file ``pyratbay/pbay.py``) or in an interactive session
through the Python interpreter.


System Requirements
-------------------

Pyrat-Bay (version 0.0) is known to work (at least) on Unix/Linux
(Ubuntu, Fedora) and OSX (10.9+) machines, with the following
software:

* Python (version 2.7)
* Numpy (version 1.8.2+)
* Scipy (version 0.13.3+)
* Matplotlib (version 1.3.1+)

MC3 may work with previous versions of these software.
However we do not guarantee it nor provide support for that.

.. _install:

Install
-------

To obtain the current stable ``Pyrat Bay`` code, clone the repository
to your local machine with the following terminal commands:

.. code-block:: shell

  topdir=`pwd`
  git clone --recursive https://github.com/pcubillos/pyratbay

Compile
-------

Compile the C programs:

.. code-block:: shell

  cd $topdir/pyratbay
  make

.. To remove the program binaries, execute (from the respective directories):
   code-block:: shell
   make clean


Quick Example: pyrat forward modeling
-------------------------------------

The following script quickly you calculate a water transmission
spectrum between 0.5 and 5.5 um.  These instructions are meant to be
executed from a Shell terminal.  After you installed and compiled the
package, create a working directory to place the files and execute the
programs, e.g.:

.. code-block:: shell

   cd $topdir
   mkdir run01
   cd run01

Download the water line-transition database from the HITRAN server:

.. code-block:: shell

   # Using wget:
   wget --user=HITRAN --password=getdata -N https://www.cfa.harvard.edu/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Compressed-files/01_hit12.zip
   # Or alternatively, curl:
   curl -u HITRAN:getdata https://www.cfa.harvard.edu/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Compressed-files/01_hit12.zip -o 01_hit12.zip

Unzip the file:

.. code-block:: shell

   unzip 01_hit12.zip


Copy and execute the demo configuration file for ``lineread`` to make a
Transition-Line-Information (TLI) file:

.. code-block:: shell

   cp $topdir/pyratbay/examples/lineread_demo/demo_hitran.cfg .
   $topdir/pyratbay/pbay.py -c demo_hitran.cfg

Copy the ``pyrat`` configuration file and run it to compute the
transmission and emission spectra:

.. code-block:: shell

   cp $topdir/pyratbay/examples/pyrat_demo/* .
   $topdir/pyratbay/pyrat.py -c demo_transmission.cfg
   $topdir/pyratbay/pyrat.py -c demo_emission.cfg

Outputs
^^^^^^^

That's it, now let's see the results.  The screen outputs and any
warnings raisedare are saved into log files.  The output spectrum is
saved to a separate file, to see it, run this Python script:

.. code-block:: python

  import matplotlib.pyplot as plt
  import sys
  sys.path.append("../pyratbay/")
  import pyratbay as pb
  wl, transmission = pb.starspec.readpyrat("./transmission_spectrum_demo.dat", wn=False)
  wl, emission     = pb.starspec.readpyrat("./emission_spectrum_demo.dat", wn=False)
  
  plt.figure(0)
  plt.clf()
  ax = plt.subplot(211)
  plt.semilogx(wl, transmission, "b", label="Pyrat transmission model")
  plt.xlabel(r"${\rm Wavelength\ \ (um)}$")
  plt.ylabel(r"${\rm Modulation\ \ (R_p/R_s)^2}$")
  ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
  ax.set_xticks([0.5, 0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0])
  plt.xlim(0.5, 5.5)
  plt.ylim(0.018, 0.0205)
  plt.legend(loc="upper left")

  ax = plt.subplot(212)
  plt.semilogx(wl, emission, "b", label="Pyrat emission model")
  plt.xlabel(r"${\rm Wavelength\ \ (um)}$")
  plt.ylabel(r"${\rm Emission\ \ (erg\ s^{-1} cm^{-2} cm)}$")
  ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
  ax.set_xticks([0.5, 0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0])
  plt.ylim(0, 60000)
  plt.xlim(0.5, 5.5)
  plt.legend(loc="upper left")
  plt.show()
