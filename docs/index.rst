.. Pyrat-Bay documentation master file, created by
   sphinx-quickstart on Fri Jan  8 16:23:24 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pyrat-Bay: Python Radiative Transfer in a Bayesian framework
============================================================

:Author:       Patricio Cubillos and collaborators (see :ref:`team`)
:Contact:      `patricio.cubillos[at]oeaw.ac.at`_
:Organizations: `Space Research Institute (IWF)`_
:Web Site:     https://github.com/pcubillos/Pyrat-Bay
:Date:         |today|

Features
--------

Pyrat-Bay is a efficient, user-friendly tool to compute and fit radiative-transfer spectra.  This package offers:

- Forward-model radiative-transfer calculation of:

  - Transmission spectra (transit observations)
  - Emission spectra (eclipse observations)
  - Line-by-line calculation
  - Opacity sources:

    - Line-by-line molecular absorption
    - Collision-induced absorption
    - Rayleigh scattering absorption

- Bayesian (Markov-chain Monte Carlo) posterior sampling of atmospheric parameters:

  - Molecular abundances
  - Temperature profile

.. _team:

Team Members
------------

- `Patricio Cubillos`_ (UCF, IWF) `patricio.cubillos[at]oeaw.ac.at`_

License
-------

Pyrat-Bay is open-source open-development software under the TBD :ref:`license`.

Be Kind
-------

Please reference this paper if you found this module useful for your research:
  `Cubillos et al. 2016: Yet Another Open-source Radiative-Transifer Code for Exoplanet Modeling`_, in preparation.

We welcome your feedback, but do not necessarily guarantee support (I will try though).
Please send feedback or inquiries to:

  Patricio Cubillos (`patricio.cubillos[at]oeaw.ac.at`_)


Contents
========

.. toctree::
   :maxdepth: 3

   getstarted
   tutorial
   lineread
   pyrat
   pyrat-bay
   license

.. _Patricio Cubillos: https://github.com/pcubillos/
.. _patricio.cubillos[at]oeaw.ac.at: patricio.cubillos@oeaw.ac.at
.. _Space Research Institute (IWF): http://iwf.oeaw.ac.at/
.. _Cubillos et al. 2016\: Yet Another Open-source Radiative-Transifer Code for Exoplanet Modeling: https://github.com/pcubillos/MCcubed/lalala
