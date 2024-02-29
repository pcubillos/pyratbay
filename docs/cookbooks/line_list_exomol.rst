.. _line_list_exomol:

ExoMol line lists
=================

This tutorial shows how to fetch ExoMol line lists, and sample them into
cross-section files for use in ``Pyrat Bay`` radiative-transfer
calculations.

   .. Note::
    You can also find this tutorial as a `jupyter notebook here
    <https://github.com/pcubillos/pyratbay/blob/master/docs/cookbooks/line_list_exomol.ipynb>`_.


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

Download ExoMol data
--------------------

You can find HITRAN and HITEMP line lists from their website:

-  https://www.exomol.com/data/molecules

There are 3 types of files to fetch, the ``.trans`` files, the
``.states`` files, and the ``.pf`` files. For this demo, we will get
those fot the HCN line lists. We can do this with the following prompt
commands:

.. code:: shell

   # Download the data
   $ wget https://www.exomol.com/db/HCN/1H-12C-14N/Harris/1H-12C-14N__Harris.states.bz2
   $ wget https://www.exomol.com/db/HCN/1H-12C-14N/Harris/1H-12C-14N__Harris.trans.bz2
   $ wget https://www.exomol.com/db/HCN/1H-12C-14N/Harris/1H-12C-14N__Harris.pf

   $ wget https://www.exomol.com/db/HCN/1H-13C-14N/Larner/1H-13C-14N__Larner.states.bz2
   $ wget https://www.exomol.com/db/HCN/1H-13C-14N/Larner/1H-13C-14N__Larner.trans.bz2
   $ wget https://www.exomol.com/db/HCN/1H-13C-14N/Larner/1H-13C-14N__Larner.pf

   # Unzip the data
   $ bzip2 -d *.bz2

Format the partition function
-----------------------------

Before generating the TLI file, we will format the partition-function
files from ExoMol for use in ``Pyrat Bay``. We can do this with the
following prompt command where we first specify the source (``exomol``)
and then list all *‘.pf’* files of interest (one can combine multiple
isotopologues of a species into a single file):

.. code:: shell

   $ pbay -pf exomol 1H-12C-14N__Harris.pf 1H-13C-14N__Larner.pf

This will produce the *PF_exomol_HCN.dat* file, which can be passed as
input for the TLI config file.

Compute TLI files
-----------------

The easiest way to generate TLI files is via configuration files and the
command line. The config file below
(`tli_exomol_HCN_cookbook.cfg <https://github.com/pcubillos/pyratbay/blob/master/docs/cookbooks/tli_exomol_HCN_cookbook.cfg>`__)
converts the ExoMol/HCN line-lists (see ``dblist``) into a TLI file (see
``tlifile`` or ``logfile``).

Partition-function information must also be provided (see ``pflist``).
As in this demo (see above), this is the path to a partition-function
file (either a unique PF file for all ``dblist`` files, or one PF file
for each ``dblist`` file). Alternatively, one can set ``pflist = tips``
to use the partition functions from `Gamache et
al. (2017) <https://ui.adsabs.harvard.edu/abs/2017JQSRT.203...70G>`__.

Lastly, the user can specify the wavelength range of the extracted data
(see ``wllow`` and ``wlhigh``). Normally one want to the widest possible
range (to avoid needing to re-calculating TLI files if a future
calculation needs it), but for sake of this demo, we will extract just
over a narrow region:

.. code:: ini

   [pyrat]

   # Select Pyrat Bay run mode: [tli atmosphere spectrum opacity retrieval radeq]
   runmode = tli

   # Output log and TLI file (if you ommit `tlifile`, it will be automatically generated from the logfile):
   logfile = Exomol_HCN_1.0-3.0um.log

   # List of line-transtion databases (.trans files):
   # (make sure the corresponding .states files are in the same folder)
   dblist =
       1H-12C-14N__Harris.trans  
       1H-13C-14N__Larner.trans

   # Type of line-transition database, select from:
   # [hitran exomol repack pands tioschwenke voplez]
   dbtype = exomol

   # List of partition functions for each database:
   pflist = PF_exomol_HCN.dat

   # Initial and final wavelength:
   wllow = 1.0 um
   wlhigh = 3.0 um

   # Verbosity level (<0:errors, 0:warnings, 1:headlines, 2:details, 3:debug):
   verb = 2

To generate the tli files, we run these ``Pyrat Bay`` prompt commands:

.. code:: shell

   $ pbay -c tli_exomol_HCN_cookbook.cfg

Compute cross-section tables
----------------------------

As with TLI files, cross-section files can be generated via
configuration files and the command line. The config file below
(`opacity_exomol_HCN_cookbook.cfg <https://github.com/pcubillos/pyratbay/blob/master/docs/cookbooks/opacity_exomol_HCN_cookbook.cfg>`__)
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
   logfile = cross_section_R020K_0150-3000K_1.0-3.0um_exomol_HCN_harris-larner.log

   # Pressure sampling:
   pbottom = 100 bar
   ptop = 1e-8 bar
   nlayers = 51

   # Temperature profile (needed, but not relevant for cross-section generation)
   tmodel = isothermal
   tpars = 1000.0

   # A simplified H2/He-dominated composition
   chemistry = uniform
   species = H2  He  HCN
   uniform = 0.85 0.15 1e-4


   # Wavelength sampling
   wllow = 1.0 um
   wlhigh = 3.0 um
   resolution = 20000.0
   # Line-profile wings extent (in HWHM from center):
   vextent = 300.0

   # Input TLI file:
   tlifile = Exomol_HCN_1.0-3.0um.tli

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

   $ pbay -c opacity_exomol_HCN_cookbook.cfg
