.. |H2O| replace:: H\ :sub:`2`\ O
.. |CO2| replace:: CO\ :sub:`2`
.. |CH4| replace:: CH\ :sub:`4`
.. |H2|  replace:: H\ :sub:`2`
.. |N2O| replace:: N\ :sub:`2`\ O
.. |NO2| replace:: NO\ :sub:`2`


.. _tlitutorial:

TLI Tutorial
============


This mode formats the line-by-line (LBL) opacity data into a
transition-line information (TLI) file, the format used by ``Pyrat
Bay`` to compute opacities.

Available Databases
-------------------

The following table lists the available
LBL opacity databases that are compatible with ``Pyrat Bay``:

.. list-table::
   :header-rows: 1

   * - Source
     - Species
     - References
   * - `HITRAN`_/`HITEMP`_
     - |H2O|, CO, |CO2|, |CH4| (+others)
     - [Rothman2010]_ [Rothman2013]_ [Li2015]_ [Gordon2017]_ [Hargreaves2019]_
   * - `ExoMol`_
     - |H2O|, CO, |CO2|, |CH4|, TiO (+others)
     - [Tennyson2016]_
   * - `Partridge & Schwenke`_
     - |H2O|
     - [PS1997]_
   * - `Schwenke`_
     - TiO
     - [Schwenke1998]_
   * - Plez
     - VO
     - [Plez1998]_


``Pyrat Bay`` is also compatible with the ``repack`` code to compress
ExoMol or HITEMP databases (see more details in [Cubillos2017b]_).

.. note:: Note that besides these LBL opacities, there are also
          :ref:`collision-induced absorption (CIA)<cia_opacity>`
          opacities that vary smoothly with wavelength.  These
          opacities are processed as tabulated cross-section (CS)
          files as a function of temperature and wavelength.


.. _sample_tli_cfg:

Sample Configuration File
-------------------------

Here is an example of a TLI configuration file (`tli_hitran_H2O.cfg
<https://github.com/pcubillos/pyratbay/blob/master/examples/tutorial/tli_hitran_H2O.cfg>`_):

.. literalinclude:: ../examples/tutorial/tli_hitran_H2O.cfg

Databases
---------

The ``dblist`` key specifies the opacity database to read (e.g., the
'`.par`' HITRAN files, or the '`.trans`' ExoMol files).  A TLI
configuration file can contain as many ``dblist`` values as desired;
this is particularly useful for molecules whose LBL data is split into
multiple files (e.g., ExoMol or HITEMP data).

The ``dbtype`` key indicates of what type are the databases listed in
``dblist``.  There must be either one ``dbtype`` value for each ``dblist``
value, or a single ``dbtype`` value that applies to all ``dblist`` values.
The following table lists the ``dbtype`` values for each
database type:

============================ =========== 
Database                     dbtype      
============================ =========== 
HITRAN/HITEMP                hitran      
ExoMol                       exomol      
repack                       repack      
Partridge & Schwenke (|H2O|) pands       
Schwenke (TiO)               tioschwenke 
Plez (VO)                    voplez      
============================ =========== 

..
  note:: It is possible as well to combine multiple species in a
         single TLI run. [but I don't recommend it]


Partition Functions
-------------------

The partition function is a temperature-dependent property of each
isotope that is required to compute line intensities.
There must be either one ``pflist`` value for each ``dblist``
value, or a single ``pflist`` value that applies to all
``dblist`` values.

For the hitran databases, ``Pyrat Bay`` provides the Total Internal
Partition Sums (TIPS) data (see [Laraia2011]_ and [Gamache2017]_).  In
this case the user can set the ``pflist`` value to '`tips`'.

The exomol, pands, and tioschwenke databases provide the partition
function as ascii tables.  In this case, the user must download the
partition-function files, and format them into the ``Pyrat Bay``
format.  For example (from the command line):

.. code-block:: shell

    # Download NH3 ExoMol partition-function data to current dictory:
    wget http://www.exomol.com/db/NH3/14N-1H3/BYTe/14N-1H3__BYTe.pf
    wget http://www.exomol.com/db/NH3/15N-1H3/BYTe-15/15N-1H3__BYTe-15.pf

    # Re-format to pyratbay format:
    pbay -pf exomol 14N-1H3__BYTe.pf 15N-1H3__BYTe-15.pf


The user can also generate TIPS partition-function files from the
command line:

.. code-block:: shell

    # pbay -pf tips 
    pbay -pf tips H2O

    # Use the 'as_exomol' flag to generate a TIPS partition-function file
    # to be used for exomol database (note the output isotope naming differs):
    pbay -pf tips H2O as_exomol


Wavelength Boundaries
---------------------

The user can specify any desired wavelength boundaries for a TLI run,
by setting the ``wllow`` and ``wlhigh`` keys.  The values for these
keys must specify the units, along with the numeric value (see
:ref:`units` for a list of available units).


.. _tli_tutorial_example:

TLI Example
-----------

The ``tlifile`` key (a required key) must specify the name of the
output TLI file.  The output TLI file will include only the lines
within the specified wavelength ranges (``wllow`` and ``wlhigh``).

When running, ``Pyrat Bay`` will also create a log file containing the
screen output.  The log file name can be set explicitly through the
``logfile`` key; otherwise, it will be created from the ``tlifile``
name, changing the extension to '`.log`'.


.. note:: Before running the tli tutorial, download the HITRAN |H2O|
          file as in the :ref:`qexample`, and download the
          configuration file shown above, for example, with these shell
          commands:

          .. code-block:: shell

              tutorial_path=https://raw.githubusercontent.com/pcubillos/pyratbay/master/examples/tutorial
              wget $tutorial_path/tli_hitran_H2O.cfg

To create the TLI file with this :ref:`sample_tli_cfg`, run from the
command line:

.. code-block:: shell

   # Make a TLI file with opacity line-transition info:
   pbay -c tli_hitran_H2O.cfg


.. _HITRAN: https://hitran.org/lbl/
.. _HITEMP: https://hitran.org/hitemp/
.. _ExoMol: http://www.exomol.com/data/molecules/
.. _Partridge & Schwenke: http://kurucz.harvard.edu/molecules/h2o/
.. _Schwenke: http://kurucz.harvard.edu/molecules/tio/
