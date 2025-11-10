.. include:: _substitutions.rst


.. _cross_sections:

Cross Sections
==============

This page shows the full list of cross sections (also named
opacities) that can be included in ``Pyrat Bay`` runs.

- :ref:`cs_sampled`
- :ref:`cs_cia`
- :ref:`cs_alkali`
- :ref:`cs_rayleigh`
- :ref:`cs_h_ion`
- :ref:`cs_clouds`

----------------------------------------------------------------------


.. _cs_sampled:

Sampled cross sections
----------------------

Sampled cross sections are 3D arrays that provide the cross section
for a given species (in cm\ :sup:`2` molecule\ :sup:`-1` units),
sampled over a grid of temperature (K), pressure (bar), and wavenumber
(|kayser|) arrays.  Sampled cross sections in ``Pyrat Bay`` are Numpy
``.npz`` files, see :py:func:`write_opacity()
<pyratbay.io.write_opacity>` in the API for more details.


The Zenodo repositories `doi.org/10.5281/zenodo.16965390
<https://zenodo.org/records/16965390>`__ and
`doi.org/10.5281/zenodo.17060936
<https://zenodo.org/records/17060936>`__ provide pre-computed
cross-section files, intended for a broad range of applications.
See the list below for direct links to the files
for each molecule.

These cross sections have been computed assuming an
|H2|/He-dominated atmosphere, and terrestrial isotopic
ratios. The lines have Voigt profiles with a wing cut-off at 300
HWHM and at 25 |kayser|.  The grids sampling are:

- Wavelength: :math:`0.15-33` μm, at a constant resolution of :math:`R=25.000`
- Temperature: :math:`200-5000` K, with :math:`\Delta T = 150` K
- Pressure: :math:`1.0^{-9}-1.0^{3}` bar, equally sampled in log(`p`) with 4 samples per dex.

.. list-table:: Cross section files
  :header-rows: 1

  * - Species
    - Source
    - References

  * - `AlF <https://zenodo.org/records/17060937/files/cross_section_0.15-33.0um_0200-5000K_R025K_AlF_exomol_mollist.npz>`__
    - Exomol, mollist
    - [Bernath2020]_

  * - `C2H2 <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_C2H2_exomol_acety.npz>`__
    - Exomol, acety
    - [Chubb2020]_

  * - `CH4 <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_CH4_exomol_mm.npz>`__
    - Exomol, mm
    - [Yurchenko2024a]_

  * - `CaH <https://zenodo.org/records/17060937/files/cross_section_0.15-33.0um_0200-5000K_R025K_CaH_exomol_xab.npz>`__
    - Exomol, xab
    - [Owens2022]_

  * - `CO <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_CO_hitemp_2019.npz>`__
    - HITEMP, li
    - [Li2015]_

  * - `CO2 <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_CO2_ames_ai3000k.npz>`__
    - ames, ai3000k
    - [Huang2023]_

  * - `CS <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_CS_exomol_jnk.npz>`__
    - Exomol, jnk
    - [Paulose2015]_

  * - `CS2 <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_CS2_hitran_2020.npz>`__
    - HITRAN, 2020
    - [Gordon2022]_

  * - `FeH <https://zenodo.org/records/17060937/files/cross_section_0.15-33.0um_0200-5000K_R025K_FeH_exomol_mollist.npz>`__
    - Exomol, mollist
    - [Dulic2003]_ [Bernath2020]_

  * - `HCl <https://zenodo.org/records/17060937/files/cross_section_0.15-33.0um_0200-5000K_R025K_HCl_exomol_hitran.npz>`__
    - Exomol, hitran
    - [Dulic2003]_ [Bernath2020]_

  * - `H2O <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_H2O_exomol_pokazatel.npz>`__
    - Exomol, pokazatel
    -  [Polyansky2018]_

  * - `H2S <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_H2S_exomol_ayt2.npz>`__
    - Exomol, ayt2
    - [Azzam2016]_ [Chubb2018]_

  * - `HCN <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_HCN_exomol_harris_larner.npz>`__
    - Exomol, harris larner
    - [Harris2008]_ [Barber2014]_

  * - `HF <https://zenodo.org/records/17060937/files/cross_section_0.15-33.0um_0200-5000K_R025K_HF_exomol_coxon_hajig.npz>`__
    - Exomol, coxon-hajig
    - [Li2013]_ [Coxon2015]_ [Somogyi2021]_

  * - `KCl <https://zenodo.org/records/17060937/files/cross_section_0.15-33.0um_0200-5000K_R025K_KCl_exomol_barton.npz>`__
    - Exomol, barton
    - [Barton2014]_

  * - `KOH <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_KOH_exomol_oyt4.npz>`__
    - Exomol, oyt4
    - [Owens2021]_

  * - `NaCl <https://zenodo.org/records/17060937/files/cross_section_0.15-33.0um_0200-5000K_R025K_NaCl_exomol_barton.npz>`__
    - Exomol, barton
    - [Barton2014]_

  * - `NaH <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_NaH_exomol_rivlin.npz>`__
    - Exomol, rivlin
    - [Rivlin2015]_

  * - `NH3 <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_NH3_exomol_coyute.npz>`__
    - Exomol, coyute
    - [Coles2019]_ [Yurchenko2024b]_

  * - `OCS <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_OCS_exomol_oyt8.npz>`__
    - Exomol, oyt8
    - [Owens2024]_

  * - `OH <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_OH_hitemp_2022.npz>`__
    - HITEMP, 2022
    - [Gordon2022]_

  * - `PH <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_PH_exomol_laty.npz>`__
    - Exomol, laty
    - [Langleben2019]_

  * - `PH3 <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_PH3_exomol_salty.npz>`__
    - Exomol, salty
    - [Sousa-Silva2014]_

  * - `PN <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_PN_exomol_pain.npz>`__
    - Exomol, pain
    - [Semenov2025]_

  * - `PO <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_PO_exomol_pops.npz>`__
    - Exomol, pops
    - [Prajapat2017]_

  * - `SH <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_SH_exomol_gyt.npz>`__
    - Exomol, gyt
    - [Gorman2019]_

  * - `SiH <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_SiH_exomol_sightly.npz>`__
    - Exomol, sightly
    - [Yurchenko2018]_

  * - `SiH4 <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_SiH4_exomol_oy2t.npz>`__
    - Exomol, oy2t
    - [Owens2017]_

  * - `SiO <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_SiO_exomol_siouvenir.npz>`__
    - Exomol, siouvenir
    - [Yurchenko2022]_

  * - `SiS <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_SiS_exomol_oy2t.npz>`__
    - Exomol, ucty
    - [Upadhyay2018]_

  * - `SO <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_SO_exomol_solis.npz>`__
    - Exomol, solis
    - [Brady2023]_

  * - `SO2 <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_SO2_exomol_exoames.npz>`__
    - Exomol, exoames
    - [Underwood2016]_

  * - `TiO <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_TiO_exomol_toto.npz>`__
    - Exomol, toto
    - [McKemmish2019]_

  * - `VO <https://zenodo.org/records/16965391/files/cross_section_0.15-33.0um_0200-5000K_R025K_VO_exomol_hyvo.npz>`__
    - Exomol, hyvo
    - [Bowesman2024]_


Other sources
~~~~~~~~~~~~~

- If this sampling (e.g., different ranges or higher resolution) is
  not suitable for a project or if cross sections for other species
  are needed, users can compute their own cross sections from line
  lists from Exomol or HITRAN/HITEMP. More details in this section:
  :doc:`line sampling <line_sampling>`.

- Alternatively, ``Pyrat Bay`` is also compatible with `petitRADTRANS <https://petitradtrans.readthedocs.io/en/latest/content/available_opacities.html#high-resolution-opacities-lbl-lambda-delta-lambda-10-6>`_ cross-section files [Molliere2019]_.

Usage
~~~~~

In a configuration file, use the ``sampled_cross_sec`` key to list the
cross sections to be included in a run:

.. code-block:: ini

  # Line-sampled cross sections
  sampled_cross_sec =
      inputs/cross_section_0.15-33.0um_0200-5000K_R025K_H2O_exomol_pokazatel.npz
      inputs/cross_section_0.15-33.0um_0200-5000K_R025K_CO_hitemp_2019.npz
      inputs/cross_section_0.15-33.0um_0200-5000K_R025K_CO2_ames_ai3000k.npz
      inputs/cross_section_0.15-33.0um_0200-5000K_R025K_CH4_exomol_mm.npz
      inputs/cross_section_0.15-33.0um_0200-5000K_R025K_SO2_exomol_exoames.npz
      inputs/cross_section_0.15-33.0um_0200-5000K_R025K_H2S_exomol_ayt2.npz


Whenever sampled cross sections are used in ``Pyrat Bay``, output
spectra will be computed at the given wavelength sampling of the cross
sections.  This samping can be trimmed down or down sampled (if
desired) with the following keys:

.. code-block:: ini

    # Wavelength sampling (keep every second point between boundaries)
    wl_low = 1.0 um
    wl_high = 12.0 um
    wl_thinning = 2


Lastly, cross sections can also be used in stand-alone scripts. See
the following notebook for details:

- `Sampled cross sections notebook <cookbooks/opacity_line_sample.ipynb>`__

----------------------------------------------------------------------

.. _cs_cia:

Continuum cross sections
------------------------



This input is intended for continuum cross sections that vary slowly
with wavelength and temperature (and no dependency on pressure).  In
practice, continuum cross sections are used for **collision induced
absorptions (CIA)**.

``Pyrat Bay`` provides CIA continuum opacities for the |H2|--|H2| and
|H2|--He pairs, the two most relevant CIA sources for primary
atmospheres:

.. _Borysow: https://www.astro.ku.dk/~aborysow/programs/index.html


.. list-table::
   :header-rows: 1
   :widths: 10 13 15 30

   * - CIA
     - T range (K)
     - |lambda| range (μm)
     - References

   * - `H₂--H₂ <https://github.com/pcubillos/pyratbay/blob/master/pyratbay/data/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat>`__
     - 60 -- 7000
     - 0.6 -- 500.0
     - [Borysow2001]_ [Borysow2002]_

   * - `H₂--He <https://github.com/pcubillos/pyratbay/blob/master/pyratbay/data/CIA/CIA_Borysow_H2He_0050-7000K_0.5-031um.dat>`__
     - 50 -- 7000
     - 0.5 -- 31.0
     - [Borysow1988]_ [Borysow1989a]_ [Borysow1989b]_ [Jorgensen2000]_

.. raw:: html

    <h3>Usage</h3>

In a configuration file, use the ``continuum_cross_sec`` key to list
the CIA files to be included:

.. code-block:: ini

    # CIA cross sections
    continuum_cross_sec =
        {ROOT}pyratbay/data/CIA/CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat
        {ROOT}pyratbay/data/CIA/CIA_Borysow_H2He_0050-7000K_0.5-031um.dat

Note the ``{ROOT}`` in the path indicates the code to find the CIA file
provided by ``Pyrat Bay`` (so, no need to ever edit this path).


For other CIA pairs, e.g., from `HITRAN <https://hitran.org/cia>`__
[Karman2019]_, ``Pyrat Bay`` provides these shell commands to
re-format the downloaded CIA files.  Here's one example (to run from
the prompt):

.. code-block:: shell

    # Download and format HITRAN H2-H2 CIA file for Pyrat Bay:
    wget https://hitran.org/data/CIA/main/H2-H2_2011.cia
    pbay -cs hitran H2-H2_2011.cia


Lastly, cross sections can also be created and used in stand-alone
scripts.  See the following notebook for details: `CIA cross sections notebook <cookbooks/opacity_cia.ipynb>`__

----------------------------------------------------------------------

.. _cs_alkali:

Alkali
------

For the sodium and potasium alkali resonant doublets, ``Pyrat Bay``
provides the line-profile models from [Burrows2000]_.

=================  ========= =========================
Models             Species   References
=================  ========= =========================
``sodium_vdw``     Na        [Burrows2000]_
``potassium_vdw``  K         [Burrows2000]_
=================  ========= =========================

These profiles are based on van der Waals and statistical theory,
adopting the line parameters from VALD [Piskunov1995]_ and
collisional-broadening half-width from [Iro2005]_.

.. raw:: html

    <h3>Usage</h3>

In a configuration file, use the ``alkali`` key to list the cross
sections to be included in a run.  the optional ``alkali_cutoff`` key
sets a cutoff from the line centers (in |kayser| units) at which to
stop computing the profiles:

.. code-block:: ini

  # Alkali cross sections [sodium_vdw potassium_vdw]
  alkali =
      sodium_vdw
      potassium_vdw

  # Alkali profile cutoff (defaulted to 4500 cm-1)
  alkali_cutoff = 4500.0


Lastly, alkali cross sections can also be used in stand-alone
scripts. See the following notebook for details: `Alkali cross
sections notebook <cookbooks/opacity_alkali.ipynb>`__.

----------------------------------------------------------------------


.. _cs_rayleigh:

Rayleigh
--------

The ``rayleigh`` key sets Rayleigh opacity models.  The following
table lists the available Rayleigh model names:

===============  =======  ===
Models           Species  References
===============  =======  ===
``rayleigh_H``   H        [Kurucz1970]_ [Dalgarno1962]_
``rayleigh_He``  He       [Kurucz1970]_ [Dalgarno1962]_
``rayleigh_H2``  |H2|     [Kurucz1970]_ [Dalgarno1962]_
===============  =======  ===

.. ``rayleigh_e-``  |e-|     ---                           [Kurucz1970]_

These models are tailored for H, He, and |H2| species,
and thus req

.. raw:: html

    <h3>Usage</h3>


- In a configuration file, use the ``rayleigh`` key to list the
  cross sections to be included in a run:

.. code-block:: ini

  # Rayleigh cross sections
  rayleigh =
      rayleigh_H2
      rayleigh_He
      rayleigh_H


- Lastly, Rayleigh cross sections can also be used in stand-alone
  scripts. See the following notebook for details: `Rayleigh cross
  sections notebook <cookbooks/opacity_rayleigh.ipynb>`__

----------------------------------------------------------------------


.. _cs_h_ion:

|H-| opacity
------------

|H-| absorption becomes significant at the high temperatures expected
for ultra Hot Jupiters, where molecular hydrogen dissociates to give
way to atomic and ionic hydrogen as the most abundant species.
``Pyrat Bay`` implements the |H-| cross-section model from
[John1988]_, which accounts for bound-free photo-ionization and
free-free scattering.

.. raw:: html

    <h3>Usage</h3>

- In a configuration file, use the ``h_ion`` key to include the |H-|
  cross section:

.. code-block:: python

    # H- bound-free and free-free opacity
    h_ion = h_ion_john1988

- Lastly, |H-| cross sections can also be used in stand-alone
  scripts. See the following notebook for details: `H bound-free and
  free-free cross sections notebook <cookbooks/opacity_h_ion.ipynb>`__


----------------------------------------------------------------------

.. _cs_clouds:

Clouds
------

The ``clouds`` key sets cloud opacity models.  The table below lists
the available models:

==============  ============================================ ===
Model           Parameter names                              Comments
==============  ============================================ ===
``lecavelier``  ``log_k_ray``, ``alpha_ray``                 [Lecavelier2008]_
``deck``        ``log_p_cl``                                 Opaque gray cloud deck
``ccsgray``     ``log_k_gray``, ``log_p_top``, ``log_p_bot`` Constant gray cross-section
==============  ============================================ ===

.. tab-set::

  .. tab-item:: lecavelier
     :selected:

     The ``lecavelier`` model implements a parametric power-law cross
     section (i.e., non-gray), allowing users to simulate
     Rayleigh-like absorption:

     .. math::
         k(\lambda) = \kappa_{\rm ray}\ \kappa_0 \left(\frac{\lambda}{\lambda_0}\right)^{\alpha_{\rm ray}}.

     Two model parameters modify the strength (``log_k_ray`` =
     :math:`\log_{10}(\kappa_{\rm ray})`) and slope (``alpha_ray`` =
     :math:`\alpha_{\rm ray}`) of the cross section.  Given the
     constants :math:`\lambda_0=0.35` um and :math:`\kappa_0=5.31
     \times 10^{-27}` cm\ :sup:`2` molecule\ :sup:`-1`, evaluating at
     ``log_k_ray = 0.0`` and ``alpha_ray = -4.0`` results in a cross
     section similar to |H2| Rayleigh for a primary atmosphere.


  .. tab-item:: deck

      The ``deck`` model imposes an opaque gray cloud deck at a
      pressure defined by the model parameter ``log_p_cl`` =
      :math:`\log(p_{\rm cl}/{\rm bar})`.

      .. note:: A technical note. To avoid computing spectra that
          depend on the pressure sampling, this cloud model does not
          simply apply a large opacity at the atmospheric layer
          closest to ``log_p_cl``.  Rather, the code interpolates the
          atmospheric profile, to define a '*surface*' located exactly
          at ``log_p_cl``.

  .. tab-item:: ccsgray

      The ``ccsgray`` model creates a constant cross-section opacity
      (:math:`k = \kappa_{\rm gray}\ \kappa_0`) between two pressures,
      with :math:`\kappa_0=5.31 \times 10^{-27}` cm\ :sup:`2`
      molecule\ :sup:`-1` a constant.

      Three model parametes define the cross section, ``log_k_gray`` =
      :math:`\log10 (\kappa_{\rm gray})`, and the log-pressure ranges:
      ``log_p_top`` and ``log_p_bot`` (in bar units).


.. raw:: html

    <h3>Usage</h3>

In a configuration file, use the ``clouds`` key to include
cloud-opacity models. Each row includes the name of the model,
followed by their parameters.


.. code-block:: ini

  # Cloud cross sections and parameters [lecavelier deck ccsgray]
  clouds =
      lecavelier  1.0 -4.0
      deck       -2.0


- Some cloud cross sections can also be used in stand-alone scripts.  See the
following notebook for details: `Lecavelier cross sections notebook
<cookbooks/opacity_lecavelier.ipynb>`__.


Patchy cloud fraction
~~~~~~~~~~~~~~~~~~~~~

The optional ``fpatchy`` key will produce spectra from a linear
combination of a clear and cloudy spectrum (i.e., spectra including
and excluding the cloud opacities listed above).  For example, for a
40% cloudy / 60% clear atmosphere, set:

.. code-block:: python

  # Patchy fraction, value between [0--1]:
  fpatchy = 0.4

