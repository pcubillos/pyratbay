.. Pyrat-Bay documentation master file, created by
   sphinx-quickstart on Fri Jan  8 16:23:24 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pyrat Bay: Python Radiative Transfer in a Bayesian framework
============================================================

:Author:       Patricio Cubillos and contributors (see :ref:`team`)
:Contact:      `patricio.cubillos[at]oeaw.ac.at`_
:Organizations: `Space Research Institute (IWF)`_
:Web Site:     https://github.com/pcubillos/pyratbay
:Date:         |today|

Features
========

``Pyrat Bay`` is an efficient, user-friendly tool to compute and fit
radiative-transfer spectra.  This package offers:

- Forward-model radiative-transfer calculation of:

  - Transmission spectra (transit observations)
  - Emission spectra (eclipse observations)
  - Line-by-line calculation
  - Opacity sources:

    - Line-by-line molecular absorption
    - Collision-induced absorption
    - Rayleigh scattering absorption
    - Na and K alkali resonant lines
    - Gray and Mie (soon) scattering opacity

- Bayesian (MCMC) posterior sampling of atmospheric parameters:

  - Molecular abundances
  - Temperature profile
  - Pressure-radius
  - Rayleigh and cloud top levels

.. note:: ``Pyrat Bay`` is temporarily proprietary software.  If you
          want to take a look before we release it, send me an email
          at `patricio.cubillos[at]oeaw.ac.at`_.


.. _team:

Contributors
============

- `Patricio Cubillos`_ (IWF) `patricio.cubillos[at]oeaw.ac.at`_
- Jasmina Blecic (NYU Abu Dhabi)
- Joe Harrington (UCF)

Documentation
=============

.. toctree::
   :maxdepth: 3

   getstarted
   tlitutorial
   pttutorial
   atmtutorial
   spectutorial
   opactutorial
   mcmctutorial
   api
   units
   contributing
   license

.. add: references section

Be Kind
=======

Please reference this paper if you found ``Pyrat Bay`` useful for your research:
  `Cubillos et al. (2019): Yet Another Open-source Radiative-Transifer Code for Exoplanet Modeling`_, in preparation.

We welcome your feedback, but do not necessarily guarantee support.
Please send feedback or inquiries to:

  Patricio Cubillos (`patricio.cubillos[at]oeaw.ac.at`_)

Pyrat-Bay is (temporarily) proprietary software (see :ref:`license`). 


Documentation for Previous Releases
===================================

- `Pyrat Bay version 0.0.50 <http://geco.oeaw.ac.at/patricio/PB_v0.0.50.pdf>`_ (and earlier versions).


.. _Patricio Cubillos: https://github.com/pcubillos/
.. _patricio.cubillos[at]oeaw.ac.at: patricio.cubillos@oeaw.ac.at
.. _Space Research Institute (IWF): http://iwf.oeaw.ac.at/
.. _Cubillos et al. (2019)\: Yet Another Open-source Radiative-Transifer Code for Exoplanet Modeling: https://github.com/pcubillos/pyratbay/lalala
