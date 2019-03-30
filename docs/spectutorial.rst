.. |H2O| replace:: H\ :sub:`2`\ O
.. |CO2| replace:: CO\ :sub:`2`
.. |CH4| replace:: CH\ :sub:`4`
.. |H2|  replace:: H\ :sub:`2`

.. _spectutorial:

Spectrum Tutorial
=================

``Pyrat Bay`` offers a unique driver command that allow multiple
running modes.  A configuration file defines all of the run settings,
inputs, and outputs.

As said, all of these calls can either be run from the shell or from
the Python interpreter.  However, the interpreter has the advantage of
that you can interact with the outputs.


Running Modes
-------------

``Pyrat Bay`` offers a sequence of running modes (``runmode``):

+----------------+------------------------------------------------------------+
|  Run mode      | Description                                                |
+================+============================================================+
| ``tli``        | Generate a transition-line-information file (used for      |
|                | spectral computation)                                      |
+----------------+------------------------------------------------------------+
| ``pt``         | Compute a temperature-pressure profile                     |
+----------------+------------------------------------------------------------+
| ``atmosphere`` | Generate a 1D atmospheric model (pressure, temperature,    |
|                | and abundances)                                            |
+----------------+------------------------------------------------------------+
| ``spectrum``   | Compute forwad-modeling spectra (transmission or emission) |
+----------------+------------------------------------------------------------+
| ``opacity``    | Generate an extinction-coefficient table (to speed up      |
|                | spectra computation)                                       |
+----------------+------------------------------------------------------------+
| ``mcmc``       | Run an atmospheric retrieval                               |
+----------------+------------------------------------------------------------+


Running modes that require a previous step (e.g., a spectrum
computation requires an atmospheric model), can do all required
calculations in a single run, as long as the user provides the
necessary parameters for each step in the configuration file.
Depending on the selected running mode, the returned outputs will
differ.

The following examples show how to run each of these modes from the
Python interpreter.

.. note:: All of the tutorial commands are briefly summarized into
          this file:  `/pyratbay/examples/tutorial/run.py <https://github.com/pcubillos/pyratbay/blob/master/examples/tutorial/run.py>`_.

For a better organization, we recommend to work from a newly created
folder, e.g., from a 'run_tutorial/' folder at the same level than the
repository, i.e.:

.. code-block:: shell

   #   parentdir/
   #   |-- pyratbay/
   #   `-- run_tutorial/

Start by copying the configuration files into the working directory
and start a Python session:

.. code-block:: shell

   cp ../pyratbay/examples/tutorial/tutorial_*.cfg .
   ipython --pylab

Also you will need to download a Kurucz stellar model file.  You can
use the following shell command:

.. code-block:: shell

  curl http://kurucz.harvard.edu/grids/gridp00odfnew/fp00k2odfnew.pck -o fp00k2odfnew.pck

Be sure to include this script each time you open a Python session:

.. code-block:: python

  # This script assumes that you are in a folder at the same level than
  # the repository, i.e.:
  #    parentdir/
  #    |-- pyratbay/
  #    `-- run_tutorial/
  #  Alternatively, set the appropriate path in sys.path.append().

  import os
  import sys
  import matplotlib
  import numpy as np
  import matplotlib.pyplot as plt
  plt.ion()

  # Edit the path to the Pyrat-Bay package if necessary:
  sys.path.append("../pyratbay")
  import pyratbay as pb



spectrum Mode
.............

This mode computes a transmission or emission spectrum.  Since this
mode requires an atmospheric model, the ``atmfile`` variable can work
as either input or output.  If the atmospheric file already exists, it
will take it as input; if the atmospheric file doesn't exist, the code
will generate it (provided the configuration file contains the
required arguments).

Here is an example configuration file for this mode:

.. code-block:: python

  [pyrat]

  # Pyrat Bay run mode, select from: [tli pt atmosphere spectrum opacity mcmc]
  runmode = spectrum

  # Atmospheric model:
  atmfile  = WASP-00b.atm   ; Input/output atmospheric file

  # TLI opacity files:
  linedb  = ./HITRAN_H2O_0.3-5.0um.tli

  # Cross-section opacity files:
  csfile  = ../pyratbay/inputs/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat
            ../pyratbay/inputs/CIA/CIA_Borysow_H2He_1000-7000K_0.5-400um.dat

  # Wavelength sampling options:
  wlunits = um
  wllow   =  0.3 um ; Spectrum lower boundary (default units: wlunits)
  wlhigh  =  5.0 um ; Spectrum higher boundary (default units: wlunits)

  # Wavenumber options:
  wnunits = cm
  wnstep  = 1.0   ; Sampling rate (default units: wnunits)
  wnosamp = 2160  ; Wavenumber over-sampling rate

  # System parameters:
  radunits = km         ; Default distance units
  punits   = bar        ; Default pressure units
  rstar    = 1.27 rsun  ; Stellar radius (default units: radunits)
  rplanet  = 1.0 rjup   ; Planetary radius (default units: radunits)
  mplanet  = 0.31 mjup  ; Planetary mass
  refpressure = 0.1     ; Reference pressure at rplanet (default units: punits)

  # Maximum optical depth to calculate:
  maxdepth = 10.0

  # Observing geometry, select between: [transit eclipse]
  path = transit

  # Rayleigh models, select from: [lecavelier dalgarno_H dalgarno_He dalgarno_H2]
  rayleigh = lecavelier
  # Lecavelier parameters (log10(H2-cross-section), alpha):
  rpars    = 0.0 -4.0

  # Haze models, select from: [deck ccsgray]
  hazes = deck   ; Opaque gray cloud deck model
  hpars = -0.5   ; log10(cloud top pressure[bar])

  # Alkali opacity, select from: [SodiumVdWst PotassiumVdWst]
  alkali = SodiumVdWst

  # Number of CPUs to use for parallel processing:
  nproc = 7

  # Verbosity level [1--5]:
  verb  = 4

  # Output file names:
  logfile    = ./transmisison_tutorial.log
  outspec    = ./transmisison_spectrum_tutorial.dat

.. note:: Pro-tip: By specifying the planetary mass (``mplanet``) and
          radius (``rplanet``), ``Pyrat Bay`` will automatically
          compute ``gplanet`` using Newton's gravitational law.


For a transmission-spectrum configuration (``path=transit``) ``Pyrat
Bay`` computes the modulation spectrum, a unitless quantity
proportional to the squared planet-to-star radius ratio
(:ref:`spectrum`).  For an emission-spectrum configuration
(``path=eclipse``) ``Pyrat Bay`` computes the day-side hemisphere
integrated flux-emission spectrum (evaluated at the surface of the
planet) in erg s\ :sup:`-1` cm\ :sup:`-2` cm (:ref:`spectrum`).

Besides the cross-section and line-by-line opacities, ``Pyrat Bay``
provides the following alkali resonant-line models (``alkali`` parameter):

====================  ========= =========================
Alkali Models         Species   Reference      
====================  ========= =========================
SodiumVdWst           Na        [Burrows2000]_
PotassiumVdWst        K         [Burrows2000]_
====================  ========= =========================

These are the available Rayleigh models (``rayleigh`` parameter):

====================  ======== =======================  ===
Rayleigh Models       Species  Parameters               Reference
====================  ======== =======================  ===
lecavelier            Any      :math:`\log(f), \alpha`  [Lecavelier2008]_
dalgarno_H            H        None                     [DalgarnoWilliams1962]_
dalgarno_He           He       None                     [Kurucz1970]_
dalgarno_H2           |H2|     None                     [Kurucz1970]_
====================  ======== =======================  ===

And these are the available haze/cloud models (``hazes`` parameter):

====================  ================================================= ==
Haze/Cloud Models       Parameters                                      Comments
====================  ================================================= ==
deck                  :math:`\log(p_{\rm top})`                         Opaque gray cloud deck at :math:`p_{\rm top}` pressure. 
ccsgray               :math:`\log(f), \log(p_{\rm t}), \log(p_{\rm b}`) Constant cross-section (:math:`\log(f)`) gray cloud between :math:`p_{\rm t}` and :math:`p_{\rm b}`.
====================  ================================================= ==

The Rayleigh and haze parameters are input throught the ``rpars`` and
``hpars``, respectively.  For any of these type of models, the user
can include multiple models, simply by concatenating multiple models
(and parameters) one after the other in the config file.

.. TBD: Add details about the models' parameters in a different page.

To compute a ``Pyrat`` model spectrum run the following script:

.. code-block:: python

  pyrat = pb.pbay.run("tutorial_spectrum.cfg")

This returns a ``pyrat`` object that contains all the input,
intermediate, and output variables used.  Until I got a decent
documentation working, take a look at `objects.py
<https://github.com/pcubillos/pyratbay/blob/master/pyratbay/pyrat/objects.py>`_
to see the object's structure.

.. note:: Note that although the user can define the input units,
          (nearly) all variables are stored in CGS units in the Pyrat
          object.

To plot the resulting spectrum you can use this script:

.. code-block:: python

  plt.figure(-3)
  plt.clf()
  ax = plt.subplot(111)
  plt.semilogx(1e4/pyrat.spec.wn, pyrat.spec.spectrum, "b-")
  ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
  ax.set_xticks([0.3, 0.4, 0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0])
  plt.xlim(0.3, 5.0)
  plt.ylabel("Modulation  (Rp/Rs)^2")
  plt.xlabel("Wavelength  (um)")

.. image:: ./figures/pyrat_transmission-spectrum_tutorial.png
   :width: 70%
   :align: center

Or alternatively, use this ``plots`` subpackage's routine:

.. code-block:: python

  ax = pb.plots.spectrum(pyrat=pyrat, gaussbin=2)

  ax.set_xscale('log')
  ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
  ax.set_xticks([0.3, 0.4, 0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0])
  plt.show()
  plt.savefig("pyrat_transmission-spectrum_tutorial.pdf")

If you want to compute emission spectra, all you need to do is to
change ``path`` to ``eclipse`` and re run:

.. code-block:: python

  pyrat = pb.pbay.run("tutorial_spectrum.cfg")

  

References
----------

.. [Burrows2000] `Burrows et al. (2000): The Near-Infrared and Optical Spectra of Methane Dwarfs and Brown Dwarfs <http://adsabs.harvard.edu/abs/2000ApJ...531..438B>`_
.. [Cubillos2017] `An Algorithm to Compress Line-transition Data for Radiative-transfer Calculations <http://adsabs.harvard.edu/abs/2017ApJ...850...32C>`_
.. [DalgarnoWilliams1962] `Dalgarno & Williams (1962): Rayleigh Scattering by Molecular Hydrogen <http://adsabs.harvard.edu/abs/1962ApJ...136..690D>`_
.. [Irwin1981] `Irwin (1981): Polynomial partition function approximations of 344 atomic and molecular species <http://adsabs.harvard.edu/abs/1981ApJS...45..621I>`_
.. [Kurucz1970] `Atlas: a Computer Program for Calculating Model Stellar Atmospheres <http://adsabs.harvard.edu/abs/1970SAOSR.309.....K>`_
.. [Laraia2011] `Laraia et al. (2011): Total internal partition sums to support planetary remote sensing <http://adsabs.harvard.edu/abs/2011Icar..215..391L>`_
.. [Lecavelier2008] `Lecavelier des Etangs et al. (2008): Rayleigh Scattering in the Transit Spectrum of HD 189733b <http://adsabs.harvard.edu/abs/2008A%26A...481L..83L>`_
.. [Line2013] `A Systematic Retrieval Analysis of Secondary Eclipse Spectra. I. A Comparison of Atmospheric Retrieval Techniques <http://adsabs.harvard.edu/abs/2013ApJ...775..137L>`_
.. [Madhusudhan2009] `Madhusudhan & Seager (2009): A Temperature and Abundance Retrieval Method for Exoplanet Atmospheres. <http://adsabs.harvard.edu/abs/2009ApJ...707...24M>`_
.. [PS1997] `Partridge & Schwenke (1997): The determination of an accurate isotope dependent potential energy surface for water from extensive ab initio calculations and experimental data <http://adsabs.harvard.edu/abs/1997JChPh.106.4618P>`_
.. [Plez1998] `Plez (1998): A new TiO line list <http://adsabs.harvard.edu/abs/1998A%26A...337..495P>`_
.. [Richard2012] `New section of the HITRAN database: Collision-induced absorption (CIA) <http://adsabs.harvard.edu/abs/2012JQSRT.113.1276R>`_
.. [Rothman2010] `Rothman et al. (2010): HITEMP, the high-temperature molecular spectroscopic database <http://adsabs.harvard.edu/abs/2010JQSRT.111.2139R>`_
.. [Rothman2013] `Rothman et al. (2013): The HITRAN2012 molecular spectroscopic database <http://adsabs.harvard.edu/abs/2013JQSRT.130....4R>`_
.. [Schwenke1998] `Schwenke (19988): Opacity of TiO from a coupled electronic state calculation parametrized by AB initio and experimental data <http://adsabs.harvard.edu/abs/1998FaDi..109..321S>`_
.. [Tennyson2016] `Tennyson et al. (2016): The ExoMol database: Molecular line lists for exoplanet and other hot atmospheres <http://adsabs.harvard.edu/abs/2016JMoSp.327...73T>`_
.. [terBraak2006] `ter Braak (2006): A Markov Chain Monte Carlo version of the genetic algorithm Differential Evolution <http://dx.doi.org/10.1007/s11222-006-8769-1>`_
.. [BraakVrugt2008] `ter Braak & Vrugt (2008): Differential Evolution Markov Chain with snooker updater and fewer chains <http://dx.doi.org/10.1007/s11222-008-9104-9>`_
