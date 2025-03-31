.. |H2O| replace:: H\ :sub:`2`\ O
.. |CO2| replace:: CO\ :sub:`2`
.. |CH4| replace:: CH\ :sub:`4`
.. |H2|  replace:: H\ :sub:`2`
.. |N2O| replace:: N\ :sub:`2`\ O
.. |NO2| replace:: NO\ :sub:`2`


.. _line_sampling:

Line Sampling
=============


`Pyrat Bay` enable users to generate their own line-sampled
cross-section files out of line-by-line data.  In this process, users
can make all the customization that they deem necessary (e.g.,
line-wing cutoffs, sampling rates, wavelength ranges, temperature
ranges, pressure ranges).

There are two steps needed to compute cross sections:

1. Convert line lists from their **original format** (e.g., HITRAN or
   ExoMol) into transition-line information files (**TLI**).

2. Conver **TLI** files into **cross-section tables** (saved as Numpy
   ``.npz`` files).



Available Databases
-------------------

``Pyrat Bay`` can process line-lists from the two main sources of
line-transition data of interest for exoplanet atmospheres: HITRAN [Gordon2022]_ and
Exomol [Tennyson2016]_.  Additionally, ``Pyrat Bay`` is also compatible with
``repack`` [Cubillos2017b]_, a package extracting only the strong line
transitions from large HITEMP, ExoMol, or AMES databases (reducing
from billions to millions of lines).

The tables below show the main species of interest to model exoplanet
atmospheres from each database (non-exhaustive, there are move species
available in each database).

.. list-table:: Available linelists from HITRAN
   :header-rows: 1
   :widths: 7, 20, 25

   * - Source
     - Species (label)
     - References
   * - HITEMP
     - `H2O <https://hitran.org/hitemp/>`__ (2010)
     - [Rothman2010]_
   * - 
     - `CO2 <https://hitran.org/hitemp/>`__ (2024)
     - [Hargreaves2025]_
   * - 
     - `CO <https://hitran.org/hitemp/>`__ (2019)
     - [Li2015]_
   * - 
     - `CH4 <https://hitran.org/hitemp/>`__ (2020)
     - [Hargreaves2020]_
   * - 
     - `N2O <https://hitran.org/hitemp/>`__ (2019)
     - [Hargreaves2019]_
   * - 
     - `NO <https://hitran.org/hitemp/>`__ (2019)
     - [Hargreaves2019]_
   * - HITRAN
     - `H2O, NH3, and many others <https://hitran.org/lbl/>`__
     - [Gordon2022]_



.. list-table:: Available linelists from ExoMol
   :header-rows: 1
   :widths: 7, 20, 25

   * - Source
     - Species (label)
     - References
   * - ExoMol
     - `H2O <https://www.exomol.com/data/molecules/H2O/>`__ (pokazatel)
     - [Polyansky2018]_
   * - 
     - `CO2 <https://www.exomol.com/data/molecules/CO2/12C-16O2/UCL-4000/>`__ (ucl4000)
     - [Yurchenko2020]_
   * - 
     - `CH4 <https://www.exomol.com/data/molecules/CH4/12C-1H4/MM/>`__ (mm)
     - [Yurchenko2024a]_
   * - 
     - `NH3 <https://www.exomol.com/data/molecules/NH3/>`__ (coyute)
     - [Coles2019]_ [Yurchenko2024b]_
   * - 
     - `TiO <https://www.exomol.com/data/molecules/TiO/>`__ (toto)
     - [McKemmish2019]_
   * - 
     - `VO <https://www.exomol.com/data/molecules/VO/51V-16O/HyVO/>`__ (hyvo)
     - [Bowesman2024]_
   * - 
     - `HCN <https://www.exomol.com/data/molecules/HCN/>`__ (harris & larner)
     - [Harris2008]_ [Barber2014]_
   * - 
     - `SO2 <https://www.exomol.com/data/molecules/SO2/32S-16O2/ExoAmes/>`__ (exoames)
     - [Underwood2016]_
   * - 
     - `H2S <https://www.exomol.com/data/molecules/H2S/1H2-32S/AYT2/>`__ (ayt2)
     - [Azzam2016]_ [Chubb2018]_
   * - 
     - `C2H2 <https://www.exomol.com/data/molecules/C2H2/12C2-1H2/aCeTY/>`__ (acety)
     - [Chubb2020]_



.. list-table:: Available linelists from repack
   :header-rows: 1
   :widths: 7, 20, 25

   * - Source
     - Species (label)
     - References
   * - repack
     - `H2O <https://zenodo.org/api/records/14266247/draft/files/H2O_exomol_pokazatel_0.24-500.0um_100-3500K_threshold_0.01_lbl.dat>`__ (exomol, pokazatel)
     -  [Cubillos2017b]_ [Polyansky2018]_
   * - 
     - `CO2 <https://zenodo.org/api/records/14266247/draft/files/CO2_exomol_ucl4000_0.5-500.0um_100-3500K_threshold_0.01_lbl.dat>`__ (exomol, ucl4000)
     - [Cubillos2017b]_ [Yurchenko2020]_
   * - 
     - `CO2 <https://zenodo.org/api/records/14266247/draft/files/CO2_ames_ai3000k_0.5-50.0um_100-3500K_threshold_0.01_lbl.dat>`__ (ames, ai3000k) 
     - [Cubillos2017b]_ 
   * - 
     - `CH4 <https://zenodo.org/api/records/14266247/draft/files/CH4_exomol_mm_0.83-50.0um_100-3000K_threshold_0.03_lbl.dat>`__ (exomol, mm)
     - [Cubillos2017b]_ [Yurchenko2024a]_
   * - 
     - `NH3 <https://zenodo.org/api/records/14266247/draft/files/NH3_exomol_coyute_0.5-500.0um_100-3000K_threshold_0.01_lbl.dat>`__ (exomol, coyute)
     - [Cubillos2017b]_ [Coles2019]_ [Yurchenko2024b]_
   * - 
     - `TiO <https://zenodo.org/api/records/14266247/draft/files/TiO_exomol_toto_0.33-500um_100-3500K_threshold_0.01_lbl.dat>`__ (exomol, toto)
     - [Cubillos2017b]_ [McKemmish2019]_
   * - 
     - `VO <https://zenodo.org/api/records/14266247/draft/files/VO_exomol_hyvo_0.22-50um_100-3500K_threshold_0.01_lbl.dat>`__ (exomol, hyvo)
     - [Cubillos2017b]_ [Bowesman2024]_
   * - 
     - `HCN <https://zenodo.org/api/records/14266247/draft/files/HCN_exomol_harris-larner_0.56-500um_100-3500K_threshold_0.01_lbl.dat>`__ (exomol, harris larner)
     - [Cubillos2017b]_ [Harris2008]_ [Barber2014]_
   * - 
     - `SO2 <https://zenodo.org/api/records/14266247/draft/files/SO2_exomol_exoames_1.25-100.0um_100-3500K_threshold_0.03_lbl.dat>`__ (exomol, exoames)
     - [Cubillos2017b]_ [Underwood2016]_
   * - 
     - `H2S <https://zenodo.org/api/records/14266247/draft/files/H2S_exomol_ayt2_0.28-500.0um_100-3500K_threshold_0.01_lbl.dat>`__ (exomol, ayt2) 
     - [Cubillos2017b]_ [Azzam2016]_ [Chubb2018]_
   * - 
     - `C2H2 <https://zenodo.org/api/records/14266247/draft/files/C2H2_exomol_acety_1.0-500.0um_100-3500K_threshold_0.03_lbl.dat>`__ (exomol, acety)
     - [Cubillos2017b]_ [Chubb2020]_
   * - 
     - `C2H4 <https://zenodo.org/api/records/14266247/draft/files/C2H4_exomol_mayty_1.4-500um_100-3500K_threshold_0.03_lbl.dat>`__ (exomol, mayty)
     - [Cubillos2017b]_ [Mant2018]_



.. _sample_tli_cfg:

Line-sampling Tutorials
-----------------------

The following links show step-by-step tutorials to compute
cross-section tables (i.e., line sampling) from each database:

- :doc:`cookbooks/line_sampling/line_list_hitran`
- :doc:`cookbooks/line_sampling/line_list_exomol`
- :doc:`cookbooks/line_sampling/line_list_repack`


Computing cross sections often require to sample millions-to-billions
of line transitions. Thus, it can both become a computer-intensive
calculation and require/produce large files. ``Pyrat Bay`` was
designed with this two-step approach to optimize these calculations.

Below there are a few useful notes to keep in mind when computing
cross sections.


Partition functions
~~~~~~~~~~~~~~~~~~~

While line-list databases do not change often with time, thus, TLI
files containing the line transitions can typically be computed only
once, and being used for many projects.  For this reason, it is
recommended to compute TLI files with the widest possible available
wavelength range.

TLI calculations also requires the input of partition functions, which
depend on the temperature.  Ideally, you also want to compute TLI
files with the widest temperature range.  The tutorials above show how
to obtain partition functions from each database.  In addition, this
tutorial below shows more in general how to handle partition functions
(e.g., how to compute PFs at high temperatures):

- :doc:`cookbooks/partition_functions`


Cross sections
~~~~~~~~~~~~~~

In contrast, cross sections might need to be more frequently
re-computed for specific project to adjust, for example, to the
sampling resolution or grid boundaries (in `T`, `p`, or
:math:`\lambda`), or need to separate between different isotopes of a
same species.  Use the tutorials above as a template, and modifies
them to adjust to the scientific requirements/machine capabilities
that each project requires.


If in a hurry and want to immediately start computing spectra, ``Pyrat
Bay`` is also compatible with the `petitRADTRANS
<https://petitradtrans.readthedocs.io/en/latest/content/available_opacities.html#high-resolution-opacities-lbl-lambda-delta-lambda-10-6>`_
cross section files [Molliere2019]_.  These files can be directly used
as input in spectrum or retrieval calculations (and can be used in
combination with the ``wn_thinning`` argument of the configuration
files to reduce the sampling resolution).  More documentation on this
is coming `soon`.
