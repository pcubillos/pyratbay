.. |H2O| replace:: H\ :sub:`2`\ O
.. |CO2| replace:: CO\ :sub:`2`
.. |CH4| replace:: CH\ :sub:`4`
.. |H2|  replace:: H\ :sub:`2`

.. _api:

API
===

There is :ref:`atmapi`.

``Pyrat Bay`` offers a unique driver command that allow multiple
running modes.  A configuration file defines all of the run settings,
inputs, and outputs.

All of these calls can either be run from the shell or from
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


TLI Mode
........

This mode formats the Line-by-line (LBL) line-transition information
into a TLI file, used by ``Pyrat Bay`` to compute opacities.  The
following table list the available data bases (Note that cross-section
opacity files (CS) are not processed into TLI files):

==================== ============================= ==== ====== =========
Source               Species                       Type Format Reference
==================== ============================= ==== ====== =========
HITRAN               |H2O|, CO, |CO2|, |CH4| (+43) LT   LBL    [Rothman2013]_
HITEMP               |H2O|, CO, |CO2|, NO, OH      LT   LBL    [Rothman2010]_
ExoMol               |H2O|, CO, |CO2|, |CH4| (+)   LT   LBL/CS [Tennyson2016]_
Partridge & Schwenke |H2O|                         LT   LBL    [PS1997]_
Schwenke             TiO                           LT   LBL    [Schwenke1998]_
Plez                 VO                            LT   LBL    [Plez1998]_
Borysow              |H2|-|H2|, |H2|-He            CIA  CS     TBD
HITRAN               |H2|-|H2|, |H2|-He (+12)      CIA  CS     [Richard2012]_
==================== ============================= ==== ====== =========

``Pyrat Bay`` is also compatible with the ``repack`` code to compress
Exomol/HITEMP databases [Cubillos2017]_.

.. ExoMol               |H2O|, CO, |CO2|, |CH4| (+)   LT   LBL    Coming Soon


Here is an example of a TLI configuration file:

.. code-block:: python

   [pyrat]
   # For syntax see:  https://docs.python.org/2/library/configparser.html

   # Pyrat Bay run mode, select from: [tli pt atmosphere spectrum opacity mcmc]
   runmode = tli

   # List of line-transtion databases:
   dblist = ./01_hit12.par
   # Type of line-transition database:
   dbtype  = hit
   # List of partition functions for each database:
   pflist = ctips

   # Initial wavelength (microns):
   iwl =  0.3
   # Final wavelength (microns):
   fwl =  5.0

   # Output TLI file:
   outfile = ./HITRAN_H2O_0.3-5.0um.tli

   # Verbosity level [1--5]:
   verb  = 4

A TLI run requires as input the set of LBL database files
(``dblist``), DB type (``dbtype``), and partition function file
(``pflist``).  Multiple DB files from multiple species can be set in a
same configuration file, as long as one sets the corresponding list of
DB types and partition-function files.  The following table shows the
available DBs and source URLs:

====================  =============================   ====== ===
Database              Species                         dbtype URL
====================  =============================   ====== ===
Partridge & Schwenke  |H2O|                           ps     http://kurucz.harvard.edu/molecules/h2o/h2ofastfix.bin
HITRAN                |H2O|, CO, |CO2|, |CH4| (+43)   hit    http://cfa.harvard.edu/hitran
HITEMP                |H2O|, CO, |CO2|, NO, OH        hit    http://cfa.harvard.edu/hitran
Exomol                |H2O|, CO, |CO2|, |CH4| (+43)   emol   http://www.exomol.com/
Schwenke              TiO                             ts     http://kurucz.harvard.edu/molecules/tio/tioschwenke.bin
Plez                  VO                              vo     http://www.pages-perso-bertrand-plez.univ-montp2.fr
repack                Exomol/HITEMP/schwenke-TiO      repack https://github.com/pcubillos/repack
====================  =============================   ====== ===

.. VALD                  TBD                             vald   TBD

The following table lists the available partition-function files and
source URLs.  See the :ref:`sscripts` section to format the online
partition-function files into the ``Pyrat Bay`` format.

====================  =====================  ===
Database              Temperature range (K)  URL
====================  =====================  ===
Partridge & Schwenke  10-6000                http://kurucz.harvard.edu/molecules/h2o/h2opartfn.dat
HITRAN and HITEMP     70-3000                ctips*
Schwenke TiO          10-6000                http://kurucz.harvard.edu/molecules/tio/tiopart.dat
Plez VO               1000-7000              poly**
====================  =====================  ===

\* For the HITRAN and HITEMP databases, ``Pyrat Bay``
provides a modified version of the Total Internal Partition Sums
(TIPS) code [Laraia2011]_ to calculate the partition functions.

\** The VO database uses a polynomial formula from [Irwin1981]_.

.. note:: Before running the tli tutorial, download the HITRAN |H2O|
          file as in :ref:`qexample`.

To create the TLI file, run from the Python interpreter:

.. code-block:: python

   # Make a TLI file with opacity line-transition info:
   pb.pbay.run("tutorial_tli.cfg")

The output TLI file will include only the lines within the specified
wavelength ranges (``iwl`` and ``fwl``).  The screen output will be
stored to an ASCII log file with the same name as the TLI file.

PT Mode
.......

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


atmosphere Mode
...............

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
.  For an emission-spectrum configuration
(``path=eclipse``) ``Pyrat Bay`` computes the day-side hemisphere
integrated flux-emission spectrum (evaluated at the surface of the
planet) in erg s\ :sup:`-1` cm\ :sup:`-2` cm.

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


opacity Mode
............

If you plan to generate multiple spectra for a same planet with
different atmospheric models, ``Pyrat Bay`` allows you to generate an
opacity table (as function of wavelength, temperature, and pressure)
for each species with LBL opacity data (i.e., 4D in total).

Once this grid is created, the code will interpolate the extinction
coefficient from the grid instead of repeating the line-by-line
calculations, significantly speeding up the spectrum calculation.

To create/use an extinction-coefficient grid, the configuration file
just need the following variables (in addition to a spectrum run):

.. code-block:: python

  [pyrat]

  # Pyrat Bay run mode, select from: [tli pt atmosphere spectrum opacity mcmc]
  runmode = opacity
  ...
  # Opacity file name and temperature range and step
  extfile = ./opacity_100-3000K_1.0-5.0um.dat
  tmin    =  100   ; Minimum temperature for grid
  tmax    = 3000   ; Maximum temperature for grid
  tstep   =  100   ; Temperature step for grid
  nproc   =    7   ; Number of parallel processors

The ``extfile`` variable sets the file name of the input/output
extinction-coefficient file.  The ``tmin``, ``tmax``, and ``tstep``
variables set the temperature sampling rate of the grid.  The
``nproc`` variable (default ``nproc=1``) sets the number of parallel
processors used to compute the extinction-coefficient.

The following table describes what ``Pyrat Bay`` outputs depending on
the ``runmode``, whether ``extfile`` was set in the configuration
file, and whether the extinction file already exists:

+-----------+-----------+-------------+---------------------------------------+
| Run mode  |``extfile``| File exists | Output                                |
+===========+===========+=============+=======================================+
| opacity   | defined   | No          | Generate new grid file                |
+           +-----------+-------------+---------------------------------------+
|           | defined   | Yes         | Overwrite grid file                   |
+           +-----------+-------------+---------------------------------------+
|           | undefined | ---         | Error                                 |
+-----------+-----------+-------------+---------------------------------------+
| spectrum  | defined   | No          | Generate grid and compute spectrum    |
+           +-----------+-------------+---------------------------------------+
|           | defined   | Yes         | Use existing grid to compute spectrum |
+           +-----------+-------------+---------------------------------------+
|           | undefined | ---         | LBL-opacity spectrum calculation      |
+-----------+-----------+-------------+---------------------------------------+
| mcmc      | defined   | ---         | Same as ``runmode=spectrum``          |
+           +-----------+-------------+---------------------------------------+
|           | undefined | ---         | Error                                 |
+-----------+-----------+-------------+---------------------------------------+

As always, to generate an extinction-coefficient grid, run the
following script:

.. code-block:: python

  pyrat = pb.pbay.run("tutorial_opacity.cfg")


mcmc Mode
.........

This mode allows you to fit spectra to observed exoplanet data.
``Pyrat Bay`` incorporates the ``MC3`` package
(`github.com/pcubillos/MCcubed
<https://github.com/pcubillos/MCcubed>`_) to retrieve best-fitting
parameters and credible regions for the atmospheric parameters in a
Bayesian (MCMC) framework.

Here is an extract of an mcmc configuration file, showing the new
required variables:

.. code-block:: python

  [pyrat]

  # Pyrat Bay run mode, select from: [tli pt atmosphere spectrum opacity mcmc]
  runmode = mcmc
  ...
  # Filter bandpasses:
  filter = ../pyratbay/inputs/filters/tutorial/tutorial_band01.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band02.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band03.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band04.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band05.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band06.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band07.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band08.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band09.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band10.dat

  # Eclipse data:
  data =   0.000072  0.000066  0.000078  0.000120  0.000135
           0.000160  0.000196  0.000232  0.000312  0.000344
  uncert = 0.000023  0.000021  0.000020  0.000019  0.000018
           0.000017  0.000016  0.000015  0.000014  0.000014

  # Kurucz stellar model:
  kurucz = ./fp00k2odfnew.pck

  # Retrieval abundances:
  bulk     = H2 He    ; Bulk (dominant) abundance species
  molscale = H2O      ; Variable-abundance species

  # Temperature-profile model, select from: [isothermal TCEA MadhuInv MadhuNoInv]
  tmodel = TCEA

  # Retrieval models, select from: [pt rad mol ray haze]
  retflag = pt mol

  # Fitting parameters:
  #         log(kappa) log(g1) log(g2)  alpha  beta  log(fH2O)
  params   = -0.6      -0.4     0.0     0.0    1.0      0.0
  pmin     = -3.0      -0.9    -1.3     0.0    0.5     -6.0
  pmax     =  1.0       1.0     0.7     1.0    1.1      3.0
  stepsize =  0.01      0.01    0.0     0.0    0.01     0.01

  # MCMC parameters:
  walk     = snooker   ; MCMC algorithm, select from: mrw, demc, snooker
  nsamples = 50000     ; Total number of MCMC samples
  nchains  =     7     ; Number of parallel MCMC chains
  burnin   =  1000     ; Burn-in iterations per chain
  thinning =     1     ; Chains thinning factor

  # MCMC temperature boundaries  (TBD: merge with tmin/tmax)
  tlow  =  100
  thigh = 3000


.. note:: Note that an ``mcmc`` run requires the user to set an
          extinction-coefficient grid (``extfile``) to allow the code
          to finish within a Hubble time.


The observational data is input through the ``filter``, ``data``, and
``uncert`` variables, which correspond to the filter transmission
files, the eclipse or transit values (corresponding to each filter),
and the data uncertainties, respectively.

For eclipse geometry, the user needs to input a stellar flux model.
``Pyrat Bay`` currently incorporates `Kurucz models
<http://kurucz.harvard.edu/grids.html>`_ Through the ``kurucz``
variable (marcs and Phoenix TBI).  The code selects the correct Kurucz
model based on the stellar temperature (``tstar``) and surface gravity
(``gstar``).

Through the ``retflag`` parameter, the user defines which models will
vary in the retrieval.  Currently the available options are ``pt`` for
the temperature model, ``rad`` for the planetary radius at
``refpressure``, ``mol`` for the species abundances, ``ray`` for the
Rayleigh models, and ``haze`` for the haze/cloud models.

The ``params`` variable encapsulates **all** of the MCMC model
parameter into a single array.  The user is responsible for listing
the MCMC parameters in the right order, and the right number of
parameters.  The order is always: ``pt, rad, mol, ray, haze``.
The number of ``pt``, ``ray``, and ``haze`` parameters depends on the
``tmodel``, ``rayleigh``, and ``haze`` models, respectively.
The number of ``rad`` parameters is always one, planetary radius in
**kilometers**.
The number of ``mol`` parameters is the number of ``molscale``
species.

The ``molscale`` variable set the species with variable abundance.  To
do so, the code scales the whole initial species abundance profile
(:math:`q_X^0(p)`) with the abundance free parameter (:math:`f_X`) as:

.. math::   q_X(p) = q_X^0(p) \times 10^{f_X}
   :label: eq:fabundance

To preserve the sum of the mixing ratios at each layer, the code
implements the ``bulk`` variable, which sets the species used to
balance the abundances such that the mixing ratio equals one at each
layer.

The ``pmin`` and ``pmax`` variables set the boundaries for each
parameter.  The ``stepsize`` variable sets the initial random jump of
the parameters.  If ``stepsize=0`` for a given parameter, the
parameter will remain fixed at its initial value.

Finally, ``walk`` defines the MCMC sampling algorithm: Set
``walk=snooker`` (default, recommended), for the DEMC-z algorithm with
snooker propsals [BraakVrugt2008]_; or ``walk=demc`` for the
Differential-Evolution MCMC algorithm [terBraak2006]_.  ``nsamples``
sets the total number of MCMC samples, ``nchains`` sets the number of
parallel MCMC chains, ``burnin`` sets the number of removed iterations
at the beginning of each chain, and ``thinning`` the thinning factor.

Just like before, to run the MCMC modeling, simply execute this command:

.. code-block:: python

  pyrat = pb.pbay.run("tutorial_mcmc.cfg")

  
.. _sscripts:

Scripts
-------

The `scripts
<https://github.com/pcubillos/pyratbay/tree/master/scripts>`_ folder
provide Python executable files (from shell) that reformat
cross-section data from the given online format (Borysow, EXOMOL,
HITRAN) into the ``Pyrat Bay`` format.

Additionally, there are executable files that reformat the
partition-function files from the given online format (Partridge &
Schwenke's |H2O|, Schwenke's TiO, and Barklem's) into the ``Pyrat
Bay`` format.

More explicit details are TBD. For the moment read the file's
docstrings for use.

