.. _cookbook:

Cookbooks
=========

This section contains relatively short code snippets (mostly
interactive Python interpreter scripts) to accomplish specific tasks
using ``Pyrat Bay``.

-------------------------------------------------------------------

Recipes list:

- **Atmospheric models**

  - `Temperature profiles tutorial <temperature_profiles.ipynb>`__
  - `VMR free profiles tutorial <vmr_free_profiles.ipynb>`__
  - :ref:`vmr_equilibrium_profiles` (TBD)
  - :ref:`radius_profiles` (TBD)

- **Line-by-line opacity sampling**

  - :doc:`line_sampling/line_list_hitran`
  - :doc:`line_sampling/line_list_exomol`
  - :doc:`line_sampling/line_list_repack`

- **Opacity models**

  - `Line-sample opacity tutorial <opacity_line_sample.ipynb>`__
  - `Alkali opacity tutorial <opacity_alkali.ipynb>`__
  - `Collision-induced absorption tutorial <opacity_cia.ipynb>`__
  - `Rayleigh opacity tutorial <opacity_rayleigh.ipynb>`__
  - `H- bound-free/free-free opacity <opacity_h_ion.ipynb>`__

- **Miscelaneous data manipulation**

  - `Instrumental Passbands <passbands.ipynb>`__
  - :doc:`partition_functions`
  - :ref:`radiative_equilibrium` (TBD)
  - :ref:`transmission_simulation` (TBD)
  - :ref:`emission_simulation` (TBD)

- **End-to-end analyses**

  - :ref:`transmission_retrieval` (TBD)
  - `Emission retrieval <wasp18b/notebook_emission_retrieval.ipynb>`__
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
   line_sampling/line_list_exomol
   line_sampling/line_list_repack
   opacity_line_sample
   opacity_alkali
   opacity_cia
   opacity_rayleigh
   opacity_h_ion

   passbands
   partition_functions
   wasp18b/notebook_emission_retrieval
   cross_sections_uhj/cross_sections_uhj
   compendia

