HITRAN line sampling
====================

This tutorial shows how to fetch HITRAN/HITEMP line lists, and sample
them into cross-section files for use in ``Pyrat Bay``
radiative-transfer calculations.

``Pyrat Bay`` has a two-step process to process line lists:

1. **Convert line lists** from their original format (HITRAN's
   ``.par`` files) **into
   transition-line information files (TLI files)**. This is simple a
   re-formatting step, the data is still kept as the info per
   line-transition (wavelengths, *gf*, *Elow*, isotope). TLI files can
   readily be used for ``Pyrat Bay`` radiative-transfer calculations,
   but such runs are slow as the code computes the lines shape and
   strength *on the fly* to obtain the cross sections.

2. **Conver TLI files into cross-section tables** (saved as Numpy
   ``.npz`` files). This step evaluates (i.e., *samples*) the
   line-transition information over a grid of [wavelength, temperature,
   pressure], which involves computing the line shape and strength of
   all lines at each given wl, pressure, and temperature value of the
   grid. Cross-section tables are ideal for radiative-transfer
   calculations, since the code simply interpolates from them (and
   therefore, these calculations are fast).

The main issue with cross-section is that they are not too flexible (one
might want to change, e.g., the wavelength resolution or line broadening
parameters, for which the user would need to re-generate cross-sections
from the TLI files). For this reason ``Pyrat Bay`` was designed with
this two-step approach.

Download data
-------------

You can find HITRAN and HITEMP line lists from their website:

-  https://hitran.org/lbl
-  https://hitran.org/hitemp

For this demo, we will get the HITRAN/H2O and HITEMP/CO line lists.  We
can do this with the following prompt commands:

.. code:: shell

   # Download the data
   wget https://www.cfa.harvard.edu/HITRAN/HITRAN2012/HITRAN2012/By-Molecule/Compressed-files/01_hit12.zip
   wget https://hitran.org/hitemp/data/bzip2format/05_HITEMP2019.par.bz2

   # Unzip the data
   bzip2 -d 05_HITEMP2019.par.bz2
   unzip 01_hit12.zip


Compute TLI files
-----------------

The easiest way to generate TLI files is via configuration files and
the command line. The config file below converts the HITRAN/H2O
line-list.

.. literalinclude:: ../../_static/data/tli_hitran_H2O_cookbook.cfg
    :caption: File: `tli_hitran_H2O_cookbook.cfg <../../_static/data/tli_hitran_H2O_cookbook.cfg>`_
    :language: ini


For the HITEMP/CO line list we use a similar config file:

.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/tli_hitran_CO_cookbook.cfg">tli_hitran_CO_cookbook.cfg</a></summary>

.. literalinclude:: ../../_static/data/tli_hitran_CO_cookbook.cfg
    :caption: File: tli_hitran_CO_cookbook.cfg
    :language: ini

.. raw:: html

   </details>


The ``tlifile`` and ``logfile`` parameters set the name of the output
files.  The ``dblist`` parameter sets the name(s) of the input HITRAN
file(s), along with the ``dbtype`` parameter which specifies the
format of the input data.

In addition to the line lists, the code requuires the partition-function data for the isotopes.  Set ``pflist=tips``  to use the HITRAN partition
function data from [Gamache2017]_ and [Gamache2021]_ (which is readily available in ``Pyrat Bay``).

.. Note:: 

    Note that the partition function is a temperature dependent value,
    and thus the temperature range of the input partition function
    sets the minimum and maximum temperature values at which the cross
    section can be evaluated.

.. Alternatively, one can input the path to a partition-function file [TBD: Explain how].


Lastly, the ``wllow`` and ``wlhigh`` parameters set the wavelength
range of the extracted data.  Normally
one sets the widest know range (to avoid needing to re-calculating TLI
files for future applications), but for sake of this demo, we
will extract just over the region that we need.


To generate the TLI files, we run these ``Pyrat Bay`` prompt commands:

.. code:: shell

   pbay -c tli_hitran_H2O_cookbook.cfg
   pbay -c tli_hitran_CO_cookbook.cfg


----------------------------


Compute cross-section tables
----------------------------

As with TLI files, cross-section files can be generated via
configuration files and the command line.  The config file below
computes a cross-section table (output name ``extfile``):

.. literalinclude:: ../../_static/data/opacity_hitran_H2O_cookbook.cfg
    :caption: File: `opacity_hitran_H2O_cookbook.cfg <../../_static/data/opacity_hitran_H2O_cookbook.cfg>`_
    :language: ini


The configuration file for the CO cross section is similar:

.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/opacity_hitran_CO_cookbook.cfg">opacity_hitran_CO_cookbook.cfg</a></summary>

.. literalinclude:: ../../_static/data/opacity_hitran_CO_cookbook.cfg
    :caption: File: opacity_hitran_CO_cookbook.cfg
    :language: ini

.. raw:: html

   </details>


These parameters define each array of the cross-section table:

-  The ``pbottom``, ``ptop``, and ``nlayers`` parameters define the
   pressure sampling array
-  The ``tmin``, ``tmax``, and ``tstep`` parameters define the
   temperature sampling array
-  The ``wllow``, ``wlhigh``, and ``resolution`` parameters define the
   spectral array at a constant resolution (alternatively, one can
   replace ``resolution`` with ``wnstep`` to sample at constant
   :math:`\Delta`\ wavenumber, units in cm\ :math:`^{-1}`)

For the composition (``species``), make sure to include the molecule for
which we are computing the cross-sections. Also, include the
*background* gas, which is relevant for the pressure broadening (here,
we assume a H2/He-dominated atmosphere). Only the VMR values of the
background gasses are important, trace-gas VMRs are irrelevant (see
``chemistry`` or ``uniform``. ``tmodel`` and ``tpars`` are needed to
define the atmosphere’s temperature profile, but for an opacity run,
these do not impact the calculations.

Lastly, the user can set ``ncpu`` (recommended) to speed up the
calculations using parallel computing.


To generate the cross-section files, run these ``Pyrat Bay`` prompt commands:

.. code:: shell

   pbay -c opacity_hitran_H2O_cookbook.cfg
   pbay -c opacity_hitran_CO_cookbook.cfg
