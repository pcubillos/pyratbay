.. _line_list_hitran:

HITRAN/HITEMP line lists
========================

This tutorial shows how to fetch HITRAN/HITEMP line lists, and sample
them into cross-section files for use in ``Pyrat Bay``
radiative-transfer calculations.

   .. Note::
    You can also find this tutorial as a `jupyter notebook here
    <https://github.com/pcubillos/pyratbay/blob/master/docs/cookbooks/line_list_hitran.ipynb>`_.


``Pyrat Bay`` has a two-step process to process line lists:

1. **Convert line lists** from their original format (e.g., HITRAN
   ``.par`` files, ExoMol ``.states/.trans`` files) **into
   transition-line information files (TLI files)**. This is simple a
   re-formatting step, the data is still kept as the info per
   line-transition (wavelengths, *gf*, *Elow*, isotope). TLI files can
   readily be used for ``Pyrat Bay`` radiative-transfer calculations,
   but such runs are slow as the code computes the lines shape and
   strength *on the fly* to obtain the cross sections.

2. **Conver TLI files into cross-section tables** (saved as Numpy
   ``.npz`` files). This step evaluates (i.e., *samples*) the
   line-transition information over a grid of [wavelength, temperature,
   pressure], which involves computing the line shape and strength of
   all lines at each given wl, pressure, and temperature value of the
   grid. Cross-section tables are ideal for radiative-transfer
   calculations, since the code simply interpolates from them (and
   therefore, these calculations are fast).

The main issue with cross-section is that they are not too flexible (one
might want to change, e.g., the wavelength resolution or line broadening
parameters, for which the user would need to re-generate cross-sections
from the TLI files). For this reason ``Pyrat Bay`` was designed with
this two-step approach.

Download HITRAN/HITEMP data
---------------------------

You can find HITRAN and HITEMP line lists from their website:

-  https://hitran.org/lbl
-  https://hitran.org/hitemp

For this demo, we will get the HITRAN/H2O and HITEMP/CO line lists. We
can do this with the following prompt commands:

.. code:: shell

   # Download the data
   $ wget https://www.cfa.harvard.edu/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Compressed-files/01_hit12.zip
   $ wget https://hitran.org/hitemp/data/bzip2format/05_HITEMP2019.par.bz2

   # Unzip the data
   $ bzip2 -d 05_HITEMP2019.par.bz2
   $ unzip 01_hit12.zip

Compute TLI files
-----------------

The easiest way to generate TLI files is via configuration files and the
command line. The config file below
(`tli_hitran_H2O_cookbook.cfg <https://github.com/pcubillos/pyratbay/blob/master/docs/cookbooks/tli_hitran_H2O_cookbook.cfg>`__)
converts the HITRAN/H2O line-list (see ``dblist``) into a TLI file (see
``tlifile`` or ``logfile``).

Partition-function information must also be provided (see ``pflist``).
If you set ``pflist=tips``, ``Pyrat Bay`` will use the partition
functions from `Gamache et
al. (2017) <https://ui.adsabs.harvard.edu/abs/2017JQSRT.203...70G>`__.
Alternatively, one can input the path to a partition-function file [TBD:
Explain how].

Lastly, the user can specify the wavelength range of the extracted data
(see ``wllow`` and ``wlhigh``). Normally one set the widest know range
(to avoid needing to re-calculating TLI files if a future calculation
needs it), but for sake of this demo, we will extract just over the
region that we need:

.. code:: ini

   [pyrat]

   # Select Pyrat Bay run mode: [tli atmosphere spectrum opacity retrieval radeq]
   runmode = tli

   # Output log and TLI file (if you ommit `tlifile`, it will be automatically generated from the logfile):
   logfile = HITRAN_H2O_0.5-5.0um.log
   tlifile = HITRAN_H2O_0.5-5.0um.tli


   # List of line-transtion databases:
   dblist = 01_hit12.par

   # Type of line-transition database, select from:
   # [hitran exomol repack]
   dbtype = hitran

   # List of partition functions for each database:
   pflist = tips

   # Initial and final wavelength:
   wllow = 0.5 um
   wlhigh = 5.0 um

   # Verbosity level (<0:errors, 0:warnings, 1:headlines, 2:details, 3:debug):
   verb = 2

To generate the TLI file for HITEMP/CO we use a similar config file
(`tli_hitran_CO_cookbook.cfg <https://github.com/pcubillos/pyratbay/blob/master/docs/cookbooks/tli_hitran_CO_cookbook.cfg>`__):

.. code:: ini

   [pyrat]

   # Select Pyrat Bay run mode: [tli atmosphere spectrum opacity retrieval radeq]
   runmode = tli

   # Output log and TLI file (if you ommit `tlifile`, it will be automatically generated from the logfile):
   logfile = HITRAN_CO_0.5-5.0um.log


   # List of line-transtion databases:
   dblist = 01_hit12.par

   # Type of line-transition database, select from:
   # [hitran exomol repack]
   dbtype = hitran

   # List of partition functions for each database:
   pflist = tips

   # Initial and final wavelength:
   wllow = 0.5 um
   wlhigh = 5.0 um

   # Verbosity level (<0:errors, 0:warnings, 1:headlines, 2:details, 3:debug):
   verb = 2

To generate the tli files, we run these ``Pyrat Bay`` prompt commands:

.. code:: shell

   $ pbay -c tli_hitran_H2O_cookbook.cfg
   $ pbay -c tli_hitran_CO_cookbook.cfg

Compute cross-section tables
----------------------------

As with TLI files, cross-section files can be generated via
configuration files and the command line. The config file below
(`opacity_hitran_H2O_cookbook.cfg <https://github.com/pcubillos/pyratbay/blob/master/docs/cookbooks/opacity_hitran_H2O_cookbook.cfg>`__)
computes a cross-section table (output name ``extfile``).

These parameters define each array of the cross-section table:

-  The ``pbottom``, ``ptop``, and ``nlayers`` parameters define the
   pressure sampling array
-  The ``tmin``, ``tmax``, and ``tstep`` parameters define the
   temperature sampling array
-  The ``wllow``, ``wlhigh``, and ``resolution`` parameters define the
   spectral array at a constant resolution (alternatively, one can
   replace ``resolution`` with ``wnstep`` to sample at constant
   :math:`\Delta`\ wavenumber, units in cm\ :math:`^{-1}`)

For the composition (``species``), make sure to include the molecule for
which we are computing the cross-sections. Also, include the
*background* gas, which is relevant for the pressure broadening (here,
we assume a H2/He-dominated atmosphere). Only the VMR values of the
background gasses are important, trace-gas VMRs are irrelevant (see
``chemistry`` or ``uniform``. ``tmodel`` and ``tpars`` are needed to
define the atmosphere’s temperature profile, but for an opacity run,
these do not impact the calculations.

Lastly, the user can set ``ncpu`` (recommended) to speed up the
calculations using parallel computing.

.. code:: ini

   [pyrat]

   # Select Pyrat Bay run mode: [tli atmosphere spectrum opacity retrieval radeq]
   runmode = opacity

   # Output log and cross-section file:
   # (if you ommit extfile it will be automatically generated from logfile name)
   logfile = cross_section_R020K_0150-3000K_0.5-5.0um_hitran_H2O.log
   extfile = cross_section_R020K_0150-3000K_0.5-5.0um_hitran_H2O.npz

   # Pressure sampling:
   pbottom = 100 bar
   ptop = 1e-8 bar
   nlayers = 51

   # Temperature profile (needed, but not relevant for cross-section generation)
   tmodel = isothermal
   tpars = 1000.0

   # A simplified H2/He-dominated composition
   chemistry = uniform
   species = H2  He  H2O  CO
   uniform = 0.85 0.15 1e-4 1e-4


   # Wavelength sampling
   wllow = 0.5 um
   wlhigh = 5.0 um
   resolution = 20000.0
   # Line-profile wings extent (in HWHM from center):
   vextent = 300.0

   # Input TLI file:
   tlifile = HITRAN_H2O_0.5-5.0um.tli

   # Cross-section temperature sampling:
   tmin =  150
   tmax = 3000
   tstep = 150

   # Number of CPUs for parallel processing:
   ncpu = 16

   # Verbosity level (<0:errors, 0:warnings, 1:headlines, 2:details, 3:debug):
   verb = 2

The configuration file for the CO cross section is similar
(`opacity_hitran_CO_cookbook.cfg <https://github.com/pcubillos/pyratbay/blob/master/docs/cookbooks/opacity_hitran_CO_cookbook.cfg>`__):

.. code:: ini

   [pyrat]

   # Select Pyrat Bay run mode: [tli atmosphere spectrum opacity retrieval radeq]
   runmode = opacity

   # Output log and cross-section file:
   logfile = cross_section_R020K_0150-3000K_0.5-5.0um_hitemp_CO.log

   # Pressure sampling:
   pbottom = 100 bar
   ptop = 1e-8 bar
   nlayers = 51

   # Temperature profile (needed, but not relevant for cross-section generation)
   tmodel = isothermal
   tpars = 1000.0

   # A simplified H2/He-dominated composition
   chemistry = uniform
   species = H2  He  H2O  CO
   uniform = 0.85 0.15 1e-4 1e-4


   # Wavelength sampling
   wllow = 0.5 um
   wlhigh = 5.0 um
   resolution = 20000.0
   # Line-profile wings extent (in HWHM from center):
   vextent = 500.0

   # Input TLI file:
   tlifile = HITRAN_CO_0.5-5.0um.tli

   # Cross-section temperature sampling:
   tmin =  150
   tmax = 3000
   tstep = 150

   # Number of CPUs for parallel processing:
   ncpu = 16

   # Verbosity level (<0:errors, 0:warnings, 1:headlines, 2:details, 3:debug):
   verb = 2

To generate the cross-section files, we run these ``Pyrat Bay`` prompt
commands:

.. code:: shell

   $ pbay -c opacity_hitran_H2O_cookbook.cfg
   $ pbay -c opacity_hitran_CO_cookbook.cfg
