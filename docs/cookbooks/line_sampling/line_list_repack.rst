Repack line sampling
====================

Sometimes, computing cross-section spectra from billions of line
transitions becomes unfeasible. For such cases, the ``repack`` tool [Cubillos2017b]_ helps
to identify and retain the strong line transitions that dominate the
spectrum. ``repack`` effectively the line list down to millions, without
significantly impacting the cross section spectra. This tutorial shows
how to fetch line lists that have been pre-processed with ``repack``,
and sample them into cross-section files for use in ``Pyrat Bay``
radiative-transfer calculations.

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

-------------------------------------------------

Download Repack data
--------------------

You can find ExoMol-repacked line lists from this Zenodo repository
(DOI: 10.5281/zenodo.3768503):

-  https://doi.org/10.5281/zenodo.3768503


The following table shows species already processed with ``repack`` and
are ready for use:

.. csv-table:: Available linelists
   :header: "Species/link", "Database name", "Isotopologues", "References"
   :widths: 8, 10, 13, 20

   `H2O <https://zenodo.org/records/14046762/files/H2O_exomol_pokazatel_0.24-500.0um_100-3500K_threshold_0.01_lbl.dat>`__, pokazatel, 116, [Polyansky2018]_

   `NH3 <https://zenodo.org/records/14046762/files/NH3_exomol_coyute-byte_0.5-500.0um_100-3500K_threshold_0.03_lbl.dat>`__, coyute & byte, "4111, 5111", [Yurchenko2015]_ [Coles2019]_

   `TiO <https://zenodo.org/records/14046762/files/TiO_exomol_toto_0.33-500um_100-3500K_threshold_0.01_lbl.dat>`__, toto, "66, 76, 86, 96, 06", [McKemmish2019]_

   `VO <https://zenodo.org/records/14046762/files/VO_exomol_vomyt_0.29-500um_100-3500K_threshold_0.01_lbl.dat>`__, vomyt, "16", [McKemmish2016]_

   `HCN <https://zenodo.org/records/14046762/files/HCN_exomol_harris-larner_0.56-500um_100-3500K_threshold_0.01_lbl.dat>`__, harris & larner, "124, 134", [Harris2006]_ [Harris2008]_

   `CO2 <https://zenodo.org/records/14046762/files/CO2_exomol_ucl4000_0.5-500.0um_100-3500K_threshold_0.01_lbl.dat>`__, ucl4000, 266, [Yurchenko2020]_

   CO2, ai3000k, "266, 366, 628, 627", [Huang2023]_

   `CH4 <https://zenodo.org/records/14046762/files/CH4_exomol_yt10to10_0.82-500.0um_100-3500K_threshold_0.03_lbl.dat>`__, yt10to10, 2111, [Yurchenko2013]_ [Yurchenko2014]_

   `CH4 <https://zenodo.org/records/14046762/files/CH4_exomol_yt34to10_0.83-500.0um_100-3500K_threshold_0.03_lbl.dat>`__, yt34to10, 2111, [Yurchenko2014]_ [Yurchenko2017]_

   `SO2 <https://zenodo.org/records/14046762/files/SO2_exomol_exoames_1.25-100.0um_100-3500K_threshold_0.03_lbl.dat>`__, exoames, 266, [Underwood2016]_

   H2S, ayt2, 112, [Azzam2016]_ [Chubb2018]_

   `C2H2 <https://zenodo.org/records/14046762/files/C2H2_exomol_acety_1.0-500.0um_100-3500K_threshold_0.03_lbl.dat>`__, acety, 2211, [Chubb2020]_

   `C2H4 <https://zenodo.org/records/14046762/files/C2H4_exomol_mayty_1.4-500um_100-3500K_threshold_0.03_lbl.dat>`__, mayty, 221111, [Mant2018]_



For this demo, we will work with the VO repack line lists. We fetch the
data can do this with the following prompt commands:

.. code:: shell

   # Download the data
   wget https://zenodo.org/records/14046762/files/VO_exomol_vomyt_0.29-500um_100-3500K_threshold_0.01_lbl.dat
   wget https://www.exomol.com/db/VO/51V-16O/VOMYT/51V-16O__VOMYT.pf

-------------------------------------------------

Format partition functions
--------------------------

Before generating the TLI file, we will format the partition-function
files from ExoMol for use in ``Pyrat Bay``. We can do this with the
following prompt command where we first specify the source (``exomol``)
and then list all *‘.pf’* files of interest (one can combine multiple
isotopologues of a species into a single file):

.. code:: shell

   pbay -pf exomol 51V-16O__VOMYT.pf

This will produce the *PF_exomol_VO.dat* file, which can be passed as
input for the TLI config file.

-------------------------------------------------

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


.. literalinclude:: ../../_static/data/tli_repack_exomol_VO_cookbook.cfg
    :caption: File: `tli_repack_exomol_VO_cookbook.cfg <../../_static/data/tli_repack_exomol_VO_cookbook.cfg>`_
    :language: ini


To generate the tli files, we run these ``Pyrat Bay`` prompt commands:

.. code:: shell

   pbay -c tli_repack_exomol_VO_cookbook.cfg

-------------------------------------------------

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


.. literalinclude:: ../../_static/data/opacity_exomol_VO_cookbook.cfg
    :caption: File: `opacity_exomol_VO_cookbook.cfg <../../_static/data/opacity_exomol_VO_cookbook.cfg>`_
    :language: ini


To generate the cross-section files, run this ``Pyrat Bay`` prompt
command:

.. code:: shell

   pbay -c opacity_exomol_VO_cookbook.cfg

-------------------------------------------------

Here's a Python script to take a look at the output cross section:

.. code:: python

   import pyratbay.io as io
   import matplotlib
   import matplotlib.pyplot as plt


   sizes, units, arrays, cross_section = io.read_opacity('cross_section_R020K_0150-3000K_0.3-3.0um_exomol_VO_vomyt.npz')
   molecs, temps, pressure, wn = arrays

   p = 35
   wl = 1e4/wn
   colors = 'royalblue', 'salmon'

   fig = plt.figure(0)
   plt.clf()
   fig.set_size_inches(7, 3)
   plt.subplots_adjust(0.1, 0.145, 0.98, 0.9)
   ax = plt.subplot(111)
   for i,t in enumerate([1,12]):
       label = f'T = {temps[t]:.0f} K'
       plt.plot(
           wl, cross_section[0,t,p], lw=1.0,
           color=colors[i], alpha=0.9, label=label,
       )
   plt.xscale('log')
   plt.yscale('log')
   ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
   ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
   ax.set_xticks([0.3, 0.5, 1.0, 2.0, 3.0])
   plt.xlim(0.3, 3.0)
   plt.ylim(1e-26, 1e-14)
   plt.title('Exomol VO (vomyt)')
   plt.xlabel('Wavelength (um)')
   plt.ylabel(r'Cross section (cm$^{2}$ / molecule)')
   plt.legend(loc='upper right')
   ax.tick_params(which='both', direction='in')


.. image:: ../../figures/VO_cross_section.png
   :width: 90%
   :align: center
