.. Pyrat-Bay documentation master file, created by
   sphinx-quickstart on Fri Jan  8 16:23:24 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pyrat Bay:
==========

**Python Radiative Transfer in a Bayesian framework**
-----------------------------------------------------

|Build Status|
|PyPI|
|docs|
|License|

.. TBD: This is MC3's badge! Replace with Pyrat Bay's when published.
   raw:: html

        <embed>
        <span class="__dimensions_badge_embed__"
            data-doi="10.3847/1538-3881/153/1/3"
            data-style="small_circle"
            data-legend="always">
        </span>
        <script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
        </embed>

.. There will probably wont be a |DOI| (as in Zenodo's DOI)


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

- `Patricio Cubillos`_ (IWF) `patricio.cubillos[at]oeaw.ac.at`_
- Jasmina Blecic (NYU Abu Dhabi)


Documentation
=============

.. toctree::
   :maxdepth: 3

   getstarted
   tli_tutorial
   pt_tutorial
   atm_tutorial
   spec_tutorial
   opac_tutorial
   mcmc_tutorial
   api
   units
   references
   contributing
   license

.. add: references section

Be Kind
=======

Please reference this paper if you found ``Pyrat Bay`` useful for your research:
  `Cubillos & Blecic (2019): An Homogeneous Retrieval Analysis of Exoplanet Transmission Spectra I: The Pyrat-Bay Framework <TBD>`_, in prep.

We welcome your feedback, but do not necessarily guarantee support.
Please send feedback or inquiries to:

  Patricio Cubillos (`patricio.cubillos[at]oeaw.ac.at`_)

Pyrat-Bay is open-source software under the GNU GPL v2 license (see
:ref:`license`), and is compatible with Python 3.6+.


Documentation for Previous Releases
===================================

- `Pyrat Bay version 0.0.50 <http://geco.oeaw.ac.at/patricio/PB_v0.0.50.pdf>`_ (and earlier versions).


.. _Patricio Cubillos: https://github.com/pcubillos/
.. _patricio.cubillos[at]oeaw.ac.at: patricio.cubillos@oeaw.ac.at
.. _Space Research Institute (IWF): http://iwf.oeaw.ac.at/

.. |Build Status| image:: https://travis-ci.com/pcubillos/pyratbay.svg?branch=master
   :target: https://travis-ci.com/pcubillos/pyratbay

.. |docs| image:: https://readthedocs.org/projects/pyratbay/badge/?version=latest
    :target: https://pyratbay.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |PyPI| image:: https://img.shields.io/pypi/v/pyratbay.svg
    :target:      https://pypi.org/project/pyratbay/
    :alt: Latest Version

.. |License| image:: https://img.shields.io/github/license/pcubillos/pyratbay.svg?color=blue
    :target: https://pcubillos.github.io/pyratbay/license.html

.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.0000000.svg
    :target: https://doi.org/10.5281/zenodo.0000000
