.. _cookbook:

Cookbooks
=========

This section contains relatively short code snippets (mostly
interactive Python interpreter scripts) to accomplish specific tasks
using ``Pyrat Bay``.

-------------------------------------------------------------------

Recipes list:

- **Atmospheric models**

  - `Temperature profile <./temperature_profiles.ipynb>`_
  - `VMR free profiles <vmr_free_profiles.ipynb>`_
  - :ref:`vmr_equilibrium_profiles` (TBD)
  - :ref:`radius_profiles` (TBD)

- **Line-by-line opacity sampling**

  - :ref:`line_list_hitran`
  - :ref:`line_list_exomol`
  - :ref:`line_list_repack`

- **Opacity models**

  - :ref:`opacity_line_sample`
  - :ref:`opacity_alkali`
  - :ref:`opacity_cia`
  - :ref:`opacity_rayleigh`
  - :ref:`opacity_h_ion` (TBD)

- **Miscelaneous**

  - :ref:`passbands` shows how to create and use instrumental response passbands
  - :ref:`partition_functions` (TBD)
  - :ref:`radiative_equilibrium` (TBD)
  - :ref:`transmission_simulation` (TBD)
  - :ref:`emission_simulation` (TBD)

- **End-to-end analyses**

  - :ref:`transmission_retrieval` (TBD)
  - :ref:`emission_retrieval` (TBD)
  - :ref:`radiative_equilibrium_simulation` (TBD)
  - :ref:`JWST_proposal_simulation` (TBD)

- :ref:`compendia` lists compendia of peer-reviewed articles with scripts that reproduced the published material

.. toctree::
   :caption: These are the available recipes
   :maxdepth: 2
   :hidden:

   temperature_profiles
   vmr_free_profiles
   line_list_hitran
   line_list_exomol
   line_list_repack
   opacity_line_sample
   passbands
   opacity_alkali
   opacity_cia
   opacity_rayleigh
   compendia

