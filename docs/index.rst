.. Pyrat-Bay documentation master file, created by
   sphinx-quickstart on Fri Jan  8 16:23:24 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pyrat Bay:
==========

**Python Radiative Transfer in a Bayesian framework**
-----------------------------------------------------

|Build Status|
|docs|
|PyPI|
|conda|
|License|


.. raw:: html

    <embed>
    <span class="__dimensions_badge_embed__"
        data-doi="10.1093/mnras/stab1405"
        data-style="small_circle"
        data-legend="always">
    </span>
    <script async src="https://badge.dimensions.ai/badge.js" charset="utf-8">
    </script>
    </embed>


-------------------------------------------------------------------


:Author:       Patricio Cubillos and contributors (see :ref:`team`)
:Contact:      `patricio.cubillos[at]oeaw.ac.at`_
:Organizations: `Space Research Institute (IWF)`_
:Web Site:     https://github.com/pcubillos/pyratbay
:Date:         |today|

Features
========

``Pyrat Bay`` is an efficient, user-friendly Python tool to compute
radiative-transfer spectra, and fit exoplanet atmospheric properties.
This package offers:

- Transmission or emission spectra of exoplanet transit or eclipses,
  respectively.
- Forward-model or retrieval calculations.

The radiative-transfer include opacity sources from:

- Line-by-line molecular absorption
- Collision-induced absorption
- Rayleigh scattering absorption
- Na and K alkali resonant lines
- Gray and Mie (soon) aerosol opacities

Bayesian (MCMC) posterior sampling of atmospheric parameters:

- Molecular abundances
- Temperature profile
- Pressure-radius
- Rayleigh and cloud properties

.. _team:

Contributors
============

- `Patricio Cubillos`_ (IWF, Austria) `patricio.cubillos[at]oeaw.ac.at`_
- Jasmina Blecic (NYU, Abu Dhabi)


Documentation
=============

.. toctree::
   :maxdepth: 2
   :includehidden:

   getstarted
   tli_tutorial
   atm_tutorial
   spec_tutorial
   opac_tutorial
   retrieval_tutorial
   cookbooks/recipes
   api
   units
   references
   contributing
   license


Be Kind
=======

If you found ``Pyrat Bay`` useful for your research, please cite this article:
  `Cubillos & Blecic (2021): The Pyrat Bay Framework for Exoplanet Atmospheric Modeling: A Population Study of Hubble/WFC3 Transmission Spectra <https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.2675C>`_, *MNRAS, 505, 2675.*

Please prefer to channel your feedback or inquiries through the Github issue tracker: `<https://github.com/pcubillos/pyratbay>`_, or alternatively through this email: `patricio.cubillos[at]oeaw.ac.at`_.

``Pyrat Bay`` is open-source software under the GNU GPL v2 license (see
:ref:`license`) and is compatible with Python>=3.6.


.. _Patricio Cubillos: https://github.com/pcubillos/
.. _patricio.cubillos[at]oeaw.ac.at: patricio.cubillos@oeaw.ac.at
.. _Space Research Institute (IWF): http://iwf.oeaw.ac.at/


.. |Build Status| image:: https://github.com/pcubillos/pyratbay/actions/workflows/python-package.yml/badge.svg?branch=master
    :target: https://github.com/pcubillos/pyratbay/actions/workflows/python-package.yml?query=branch%3Amaster

.. |docs| image:: https://readthedocs.org/projects/pyratbay/badge/?version=latest
    :target: https://pyratbay.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |PyPI| image:: https://img.shields.io/pypi/v/pyratbay.svg
    :target:      https://pypi.org/project/pyratbay/
    :alt: Latest Version

.. |conda| image:: https://img.shields.io/conda/vn/conda-forge/pyratbay.svg
    :target: https://anaconda.org/conda-forge/pyratbay

.. |License| image:: https://img.shields.io/github/license/pcubillos/pyratbay.svg?color=blue
    :target: https://pyratbay.readthedocs.io/en/latest/license.html

