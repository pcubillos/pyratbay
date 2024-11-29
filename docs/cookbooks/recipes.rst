.. _cookbook:

Cookbooks
=========

This section contains relatively short code snippets (mostly
interactive Python interpreter scripts) to accomplish specific tasks
using ``Pyrat Bay``.

-------------------------------------------------------------------

Recipes list:

- **Atmospheric models**

  - `Temperature profiles tutorial <temperature_profiles.ipynb>`_
  - `VMR free profiles tutorial <vmr_free_profiles.ipynb>`_
  - :ref:`vmr_equilibrium_profiles` (TBD)
  - :ref:`radius_profiles` (TBD)

- **Line-by-line opacity sampling**

  - :doc:`line_sampling/line_list_hitran`
  - `ExoMol line-list tutorial <line_list_exomol.ipynb>`_
  - `repack line-list tutorial <line_list_repack.ipynb>`_

- **Opacity models**

  - `Line-sample opacity tutorial <opacity_line_sample.ipynb>`_
  - `Alkali opacity tutorial <opacity_alkali.ipynb>`_
  - `Collision-induced absorption tutorial <opacity_cia.ipynb>`_
  - `Rayleigh opacity tutorial <opacity_rayleigh.ipynb>`_
  - `H- bound-free/free-free opacity <opacity_h_ion.ipynb>`_

- **Miscelaneous**

  - `Passbands <passbands.ipynb>`_ shows how to create and use instrumental response passbands
  - :ref:`partition_functions` (TBD)
  - :ref:`radiative_equilibrium` (TBD)
  - :ref:`transmission_simulation` (TBD)
  - :ref:`emission_simulation` (TBD)

- **End-to-end analyses**

  - :ref:`transmission_retrieval` (TBD)
  - `Emission retrieval <wasp18b/notebook_emission_retrieval.ipynb>`_
  - :doc:`cross_sections_uhj/cross_sections_uhj`
  - :ref:`radiative_equilibrium_simulation` (TBD)
  - :ref:`JWST_proposal_simulation` (TBD)

- :ref:`compendia` lists compendia of peer-reviewed articles with scripts that reproduced the published material

.. toctree::
   :caption: These are the available recipes
   :maxdepth: 2
   :hidden:

   temperature_profiles
   vmr_free_profiles
   line_sampling/line_list_hitran
   line_list_exomol
   line_list_repack
   opacity_line_sample
   opacity_alkali
   opacity_cia
   opacity_rayleigh
   opacity_h_ion
   passbands
   wasp18b/notebook_emission_retrieval
   cross_sections_uhj/cross_sections_uhj
   compendia

