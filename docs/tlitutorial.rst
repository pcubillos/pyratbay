.. |H2O| replace:: H\ :sub:`2`\ O
.. |CO2| replace:: CO\ :sub:`2`
.. |CH4| replace:: CH\ :sub:`4`
.. |H2|  replace:: H\ :sub:`2`

.. _tlitutorial:

TLI Tutorial
============

This mode formats the Line-by-line (LBL) line-transition information
into a TLI file, used by ``Pyrat Bay`` to compute opacities.  The
following table list the available data bases (Note that cross-section
opacity files (CS) are not processed into TLI files):

==================== ============================= ==== ====== =========
Source               Species                       Type Format Reference
==================== ============================= ==== ====== =========
HITRAN               |H2O|, CO, |CO2|, |CH4| (+43) LT   LBL    [Rothman2013]_
HITEMP               |H2O|, CO, |CO2|, NO, OH      LT   LBL    [Rothman2010]_
ExoMol               |H2O|, CO, |CO2|, |CH4| (+)   LT   LBL/CS [Tennyson2016]_
Partridge & Schwenke |H2O|                         LT   LBL    [PS1997]_
Schwenke             TiO                           LT   LBL    [Schwenke1998]_
Plez                 VO                            LT   LBL    [Plez1998]_
Borysow              |H2|-|H2|, |H2|-He            CIA  CS     TBD
HITRAN               |H2|-|H2|, |H2|-He (+12)      CIA  CS     [Richard2012]_
==================== ============================= ==== ====== =========

``Pyrat Bay`` is also compatible with the ``repack`` code to compress
Exomol/HITEMP databases [Cubillos2017]_.


Here is an example of a TLI configuration file:

.. code-block:: python

   [pyrat]
   # For syntax see:  https://docs.python.org/2/library/configparser.html

   # Pyrat Bay run mode, select from: [tli pt atmosphere spectrum opacity mcmc]
   runmode = tli

   # List of line-transtion databases:
   dblist = ./01_hit12.par
   # Type of line-transition database:
   dbtype  = hit
   # List of partition functions for each database:
   pflist = ctips

   # Initial wavelength (microns):
   iwl =  0.3
   # Final wavelength (microns):
   fwl =  5.0

   # Output TLI file:
   outfile = ./HITRAN_H2O_0.3-5.0um.tli

   # Verbosity level [1--5]:
   verb  = 4

A TLI run requires as input the set of LBL database files
(``dblist``), DB type (``dbtype``), and partition function file
(``pflist``).  Multiple DB files from multiple species can be set in a
same configuration file, as long as one sets the corresponding list of
DB types and partition-function files.  The following table shows the
available DBs and source URLs:

====================  =============================   ====== ===
Database              Species                         dbtype URL
====================  =============================   ====== ===
Partridge & Schwenke  |H2O|                           ps     http://kurucz.harvard.edu/molecules/h2o/h2ofastfix.bin
HITRAN                |H2O|, CO, |CO2|, |CH4| (+43)   hit    http://cfa.harvard.edu/hitran
HITEMP                |H2O|, CO, |CO2|, NO, OH        hit    http://cfa.harvard.edu/hitran
Exomol                |H2O|, CO, |CO2|, |CH4| (+43)   emol   http://www.exomol.com/
Schwenke              TiO                             ts     http://kurucz.harvard.edu/molecules/tio/tioschwenke.bin
Plez                  VO                              vo     http://www.pages-perso-bertrand-plez.univ-montp2.fr
repack                Exomol/HITEMP/schwenke-TiO      repack https://github.com/pcubillos/repack
====================  =============================   ====== ===

.. VALD                  TBD                             vald   TBD

The following table lists the available partition-function files and
source URLs.  See the :ref:`sscripts` section to format the online
partition-function files into the ``Pyrat Bay`` format.

====================  =====================  ===
Database              Temperature range (K)  URL
====================  =====================  ===
Partridge & Schwenke  10-6000                http://kurucz.harvard.edu/molecules/h2o/h2opartfn.dat
HITRAN and HITEMP     70-3000                ctips*
Schwenke TiO          10-6000                http://kurucz.harvard.edu/molecules/tio/tiopart.dat
Plez VO               1000-7000              poly**
====================  =====================  ===

\* For the HITRAN and HITEMP databases, ``Pyrat Bay``
provides a modified version of the Total Internal Partition Sums
(TIPS) code [Laraia2011]_ to calculate the partition functions.

\** The VO database uses a polynomial formula from [Irwin1981]_.

.. note:: Before running the tli tutorial, download the HITRAN |H2O|
          file as in :ref:`qexample`.

To create the TLI file, run from the Python interpreter:

.. code-block:: python

   # Make a TLI file with opacity line-transition info:
   pb.pbay.run("tutorial_tli.cfg")

The output TLI file will include only the lines within the specified
wavelength ranges (``iwl`` and ``fwl``).  The screen output will be
stored to an ASCII log file with the same name as the TLI file.

