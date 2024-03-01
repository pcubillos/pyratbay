.. _line_list_repack:

Repack line lists
=================


Sometimes, computing cross-section spectra from billions of line
transitions becomes unfeasible. For such cases, the ``repack`` tool
[Cubillos2017b]_ helps
to identify and retain the strong line transitions that dominate the
spectrum. ``repack`` effectively the line list down to millions, without
significantly impacting the cross section spectra. This tutorial shows
how to fetch line lists that have been pre-processed with ``repack``,
and sample them into cross-section files for use in ``Pyrat Bay``
radiative-transfer calculations.

   .. Note::
    You can also find this tutorial as a `jupyter notebook here
    <https://github.com/pcubillos/pyratbay/blob/master/docs/cookbooks/line_list_repack.ipynb>`_.


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

Download Repack data
--------------------

You can find ExoMol-repacked line lists from this Zenodo repository
(DOI: 10.5281/zenodo.3768503):

-  https://zenodo.org/records/3768504

The following table shows species already processed with ``repack`` and
are ready for use:

+-----------------+-----------------+-----------------+-------------------+
| Species         | Source          | Isotopologues   | References        |
+=================+=================+=================+===================+
| `H2O <https:    | pokazatel       | `116 <http      | [Polyansky2018]_  |
| //zenodo.org/re |                 | s://www.exomol. |                   |
| cords/3768504/f |                 | com/data/molecu |                   |
| iles/H2O_exomol |                 | les/H2O/1H2-16O |                   |
| _pokazatel_0.24 |                 | /POKAZATEL/>`__ |                   |
| -500.0um_100-35 |                 |                 |                   |
| 00K_threshold_0 |                 |                 |                   |
| .01_lbl.dat>`__ |                 |                 |                   |
+-----------------+-----------------+-----------------+-------------------+
| `NH3 <https:/   | coyute & byte   | `4111 <h        | [Yurchenko2015]_, |
| /zenodo.org/rec |                 | ttps://www.exom | [Coles2019]_      |
| ords/3768504/fi |                 | ol.com/data/mol |                   |
| les/NH3_exomol_ |                 | ecules/NH3/14N- |                   |
| coyute-byte_0.5 |                 | 1H3/CoYuTe/>`__ |                   |
| -500.0um_100-35 |                 | &               |                   |
| 00K_threshold_0 |                 | `5111 <ht       |                   |
| .03_lbl.dat>`__ |                 | tps://www.exomo |                   |
|                 |                 | l.com/data/mole |                   |
|                 |                 | cules/NH3/15N-1 |                   |
|                 |                 | H3/BYTe-15/>`__ |                   |
|                 |                 |                 |                   |
+-----------------+-----------------+-----------------+-------------------+
| `TiO            | toto            | `66 <h          | [McKemmish2019]_  |
| <https://zenodo |                 | ttps://www.exom |                   |
| .org/records/37 |                 | ol.com/data/mol |                   |
| 68504/files/TiO |                 | ecules/TiO/46Ti |                   |
| _exomol_toto_0. |                 | -16O/Toto/>`__, |                   |
| 33-500um_100-35 |                 | `76 <h          |                   |
| 00K_threshold_0 |                 | ttps://www.exom |                   |
| .01_lbl.dat>`__ |                 | ol.com/data/mol |                   |
|                 |                 | ecules/TiO/47Ti |                   |
|                 |                 | -16O/Toto/>`__, |                   |
|                 |                 | `86 <h          |                   |
|                 |                 | ttps://www.exom |                   |
|                 |                 | ol.com/data/mol |                   |
|                 |                 | ecules/TiO/48Ti |                   |
|                 |                 | -16O/Toto/>`__, |                   |
|                 |                 | `96 <h          |                   |
|                 |                 | ttps://www.exom |                   |
|                 |                 | ol.com/data/mol |                   |
|                 |                 | ecules/TiO/49Ti |                   |
|                 |                 | -16O/Toto/>`__, |                   |
|                 |                 | `06 <h          |                   |
|                 |                 | ttps://www.exo  |                   |
|                 |                 | mol.com/data/mo |                   |
|                 |                 | lecules/TiO/50T |                   |
|                 |                 | i-16O/Toto/>`__ |                   |
+-----------------+-----------------+-----------------+-------------------+
| `VO             | vomyt           | `16             | [McKemmish2016]_  |
| <https://zenodo |                 | <https://www.ex |                   |
| .org/records/37 |                 | omol.com/data/m |                   |
| 68504/files/VO_ |                 | olecules/VO/51V |                   |
| exomol_vomyt_0. |                 | -16O/VOMYT/>`__ |                   |
| 29-500um_100-35 |                 |                 |                   |
| 00K_threshold_0 |                 |                 |                   |
| .01_lbl.dat>`__ |                 |                 |                   |
+-----------------+-----------------+-----------------+-------------------+
| `HCN <https://  | harris & larner | `124 <http      | [Harris2006]_,    |
| zenodo.org/reco |                 | s://www.exomol. | [Harris2008]_     |
| rds/3768504/fil |                 | com/data/molecu |                   |
| es/HCN_exomol_h |                 | les/HCN/1H-12C- |                   |
| arris-larner_0. |                 | 14N/Harris/>`__ |                   |
| 56-500um_100-35 |                 | &               |                   |
| 00K_threshold_0 |                 | `134 <http      |                   |
| .01_lbl.dat>`__ |                 | s://www.exomol. |                   |
|                 |                 | com/data/molecu |                   |
|                 |                 | les/HCN/1H-13C- |                   |
|                 |                 | 14N/Larner/>`__ |                   |
|                 |                 |                 |                   |
+-----------------+-----------------+-----------------+-------------------+
| `CO2 <htt       | ucl4000         | `266 <http      | [Yurchenko2020]_  |
| ps://zenodo.org |                 | s://www.exomol. |                   |
| /records/376850 |                 | com/data/molecu |                   |
| 4/files/CO2_exo |                 | les/CO2/12C-16O |                   |
| mol_ucl4000_0.5 |                 | 2/UCL-4000/>`__ |                   |
| -500.0um_100-35 |                 |                 |                   |
| 00K_threshold_0 |                 |                 |                   |
| .01_lbl.dat>`__ |                 |                 |                   |
+-----------------+-----------------+-----------------+-------------------+
| `CH4 <https     | yt10to10        | `2111 <htt      | [Yurchenko2013]_, |
| ://zenodo.org/r |                 | ps://www.exomol | [Yurchenko2014]_  |
| ecords/3768504/ |                 | .com/data/molec |                   |
| files/CH4_exomo |                 | ules/CH4/12C-1H |                   |
| l_yt10to10_0.82 |                 | 4/YT10to10/>`__ |                   |
| -500.0um_100-35 |                 |                 |                   |
| 00K_threshold_0 |                 |                 |                   |
| .03_lbl.dat>`__ |                 |                 |                   |
+-----------------+-----------------+-----------------+-------------------+
| `CH4 <https     | yt34to10        | `2111 <htt      | [Yurchenko2014]_, |
| ://zenodo.org/r |                 | ps://www.exomol | [Yurchenko2017]_  |
| ecords/3768504/ |                 | .com/data/molec |                   |
| files/CH4_exomo |                 | ules/CH4/12C-1H |                   |
| l_yt34to10_0.83 |                 | 4/YT34to10/>`__ |                   |
| -500.0um_100-35 |                 |                 |                   |
| 00K_threshold_0 |                 |                 |                   |
| .03_lbl.dat>`__ |                 |                 |                   |
+-----------------+-----------------+-----------------+-------------------+
| `SO2 <http      | exoames         | `266 <ht        | [Underwood2016]_  |
| s://zenodo.org/ |                 | tps://www.exomo |                   |
| records/3768504 |                 | l.com/data/mole |                   |
| /files/SO2_exom |                 | cules/SO2/32S-1 |                   |
| ol_exoames_1.25 |                 | 6O2/ExoAmes>`__ |                   |
| -100.0um_100-35 |                 |                 |                   |
| 00K_threshold_0 |                 |                 |                   |
| .03_lbl.dat>`__ |                 |                 |                   |
+-----------------+-----------------+-----------------+-------------------+
| `H2S <h         | ayt2            | `112            | [Azzam2016]_,     |
| ttps://zenodo.o |                 | <https://www.ex | [Chubb2018]_      |
| rg/records/3768 |                 | omol.com/data/m |                   |
| 504/files/H2S_e |                 | olecules/H2S/1H |                   |
| xomol_ayt2_0.28 |                 | 2-32S/AYT2/>`__ |                   |
| -500.0um_100-30 |                 |                 |                   |
| 00K_threshold_0 |                 |                 |                   |
| .01_lbl.dat>`__ |                 |                 |                   |
+-----------------+-----------------+-----------------+-------------------+
| `C2H2 <ht       | acety           | `2211 <ht       | [Chubb2020]_      |
| tps://zenodo.or |                 | tps://www.exomo |                   |
| g/records/37685 |                 | l.com/data/mole |                   |
| 04/files/C2H2_e |                 | cules/C2H2/12C2 |                   |
| xomol_acety_1.0 |                 | -1H2/aCeTY/>`__ |                   |
| -500.0um_100-35 |                 |                 |                   |
| 00K_threshold_0 |                 |                 |                   |
| .03_lbl.dat>`__ |                 |                 |                   |
+-----------------+-----------------+-----------------+-------------------+
| `C2H4 <ht       | mayty           | `221111 <ht     | [Mant2018]_       |
| tps://zenodo.or |                 | tps://www.exomo |                   |
| g/records/376   |                 | l.com/data/mole |                   |
| 8504/files/>`__ |                 | cules/C2H4/12C2 |                   |
|                 |                 | -1H4/MaYTY/>`__ |                   |
+-----------------+-----------------+-----------------+-------------------+

For this demo, we will work with the VO repack line lists. We fetch the
data can do this with the following prompt commands:

.. code:: shell

   # Download the data
   $ wget https://zenodo.org/records/3768504/files/VO_exomol_vomyt_0.29-500um_100-3500K_threshold_0.01_lbl.dat
   $ wget https://www.exomol.com/db/VO/51V-16O/VOMYT/51V-16O__VOMYT.pf

Format the partition function
-----------------------------

Before generating the TLI file, we will format the partition-function
files from ExoMol for use in ``Pyrat Bay``. We can do this with the
following prompt command where we first specify the source (``exomol``)
and then list all *‘.pf’* files of interest (one can combine multiple
isotopologues of a species into a single file):

.. code:: shell

   $ pbay -pf exomol 51V-16O__VOMYT.pf

This will produce the *PF_exomol_VO.dat* file, which can be passed as
input for the TLI config file.

Compute TLI files
-----------------

The easiest way to generate TLI files is via configuration files and the
command line. The config file below
(`tli_repack_exomol_VO_cookbook.cfg <https://github.com/pcubillos/pyratbay/blob/master/docs/cookbooks/tli_repack_exomol_VO_cookbook.cfg>`__)
converts the repack ExoMol/VO line-lists (see ``dblist``) into a TLI
file (see ``tlifile`` or ``logfile``).

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
   logfile = Exomol_repack_VO_0.5-3.0um.log

   # List of line-transtion databases:
   dblist = VO_exomol_vomyt_0.29-500um_100-3500K_threshold_0.01_lbl.dat

   # Type of line-transition database, select from:
   # [hitran exomol repack]
   dbtype = repack

   # List of partition functions for each database:
   pflist = PF_exomol_VO.dat

   # Initial and final wavelength:
   wllow = 0.5 um
   wlhigh = 3.0 um

   # Verbosity level (<0:errors, 0:warnings, 1:headlines, 2:details, 3:debug):
   verb = 2

To generate the tli files, we run these ``Pyrat Bay`` prompt commands:

.. code:: shell

   $ pbay -c tli_repack_exomol_VO_cookbook.cfg

Compute cross-section tables
----------------------------

As with TLI files, cross-section files can be generated via
configuration files and the command line. The config file below
(`opacity_exomol_VO_cookbook.cfg <https://github.com/pcubillos/pyratbay/blob/master/docs/cookbooks/opacity_exomol_VO_cookbook.cfg>`__)
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
   logfile = cross_section_R020K_0150-3000K_0.3-3.0um_exomol_HCN_vomyt.log

   # Pressure sampling:
   pbottom = 100 bar
   ptop = 1e-8 bar
   nlayers = 51

   # Temperature profile (needed, but not relevant for cross-section generation)
   tmodel = isothermal
   tpars = 1000.0

   # A simplified H2/He-dominated composition
   chemistry = uniform
   species = H2  He  VO
   uniform = 0.85 0.15 1e-4


   # Wavelength sampling
   wllow = 0.3 um
   wlhigh = 3.0 um
   resolution = 20000.0
   # Line-profile wings extent (in HWHM from center):
   vextent = 300.0

   # Input TLI file:
   tlifile = Exomol_repack_VO_0.3-3.0um.tli

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

   $ pbay -c opacity_exomol_VO_cookbook.cfg
