.. _cross_sections_uhj:


Cross Sections for an Ultra Hot Jupiter
=======================================

This script shows how to cmopute line-sampled molecular cross
sections, to be used for atmospheric modeling and retrievals.  There
are four main steps to compute cross-section files:

1. `Download line-list data <#line-lists>`_
2. `Compute partition functions <#partition-functions>`_
3. `Compute TLI files <#tli-files>`_
4. `Sample cross sections <#cross-sections>`_

Note that the first three steps are typically executed only once,
allowing you to reuse the output files across projects.  If you have
already the tli files that you need, go directly to `Step 4
<#cross-sections>`__.

`Step 4 <#cross-sections>`__ may need to be executed on a per-project
basis, depending on your specific requirements (e.g., different
spectral, temperature, or pressure ranges; or varying resolutions).

--------------------------------------------------------

1. Line Lists
-------------

First we need to download the line lists, typically sourced from the
ExoMol or HITRAN/HITEMP databases.  You will likely only need to
complete this step once unless a new or updated line list becomes
available. It may be better to store this data in a general directory
on your machine. 

For this project, we will focus on the molecular absorbers relevant
for an ultra-hot Jupiter (WASP-18b). The table below lists the line
lists that we will use.


.. list-table:: Molecular line lists
  :header-rows: 1

  * - Species (source)
    - References

  * - `H2O <https://www.exomol.com/data/molecules/H2O/1H2-16O/POKAZATEL>`__ (exomol, pokazatel)
    -  [Polyansky2018]_

  * - `CO <https://hitran.org/hitemp>`__ (HITEMP, li)
    - [Li2015]_

  * - `CO2 <https://data.nas.nasa.gov/ai3000k>`__ (ames, ai3000k) 
    - [Huang2023]_

  * - `CH4 <https://www.exomol.com/data/molecules/CH4/12C-1H4/MM>`__ (exomol, mm)
    - [Yurchenko2024a]_

  * - `TiO <https://www.exomol.com/data/molecules/TiO>`__ (exomol, toto)
    - [McKemmish2019]_

  * - `VO <https://www.exomol.com/data/molecules/VO/51V-16O/HyVO>`__ (exomol, hyvo) 
    - [Bowesman2024]_

  * - `HCN <https://www.exomol.com/data/molecules/HCN/>`__ (exomol, harris larner)
    - [Harris2008]_ [Barber2014]_

  * - `NH3 <https://www.exomol.com/data/molecules/NH3>`__ (exomol, coyute)
    - [Coles2019]_ [Yurchenko2024b]_

  * - `C2H2 <https://www.exomol.com/data/molecules/C2H2/12C2-1H2/aCeTY>`__ (exomol, acety)
    - [Chubb2020]_


These datasets contain around 200 billion transitions, so it would be
impractical to work with the full set of line lists.  Instead, we work
with the datasets processed with ``repack`` [Cubillos2017b]_, which
filtered the dominant transitions at each wavelength, and thus
reduced the number of lines to model.  The file below contains the
actual files to download.


.. literalinclude:: ../../_static/data/wasp18b_line_lists_data.txt
    :caption: File: `wasp18b_line_lists_data.txt <../../_static/data/wasp18b_line_lists_data.txt>`__  
    :language: none


Here's a script to download these files from the command line using
``wget``.  First, make sure to copy the text file above to your current folder:

.. code-block:: shell

    # Download line lists (note there are several GBs of data):
    wget -i wasp18b_line_lists_data.txt

    # Unzip the HITEMP data:
    bzip2 -d 05_HITEMP2019.par.bz2

--------------------------------------------------------


2. Partition Functions
----------------------

In addition to the line-list data, we need the partition functions for
each molecule and isotope.  For this we will use again the data
provided by ExoMol and HITRAN TIPS ([Gamache2017]_ [Gamache2021]_).

What's *really important* is to be aware of the temperature range of
the partition functions.  For an ultra hot Jupiter like WASP-18b, we
want cross sections as hot as ~5000K.  So, the strategy will be to use
the tabulated partition functions files (.pf files) when they cover
the temperatures we need. Otherwise, we will compute the partition
from the .states files.  For CO2 I provide the partitions computed
from the Ames states files.

This file below contains links to the input partition-function data
from Exomol, ames, and HITRAN:


.. literalinclude:: ../../_static/data/wasp18b_partition_function_data.txt
    :caption: File: `wasp18b_partition_function_data.txt <../../_static/data/wasp18b_partition_function_data.txt>`__  
    :language: none


``Pyrat Bay`` provides the commands to parse and compute the
partitions into the required format.  Here's a script to download
these files above and process the partition functions (First, make
sure to copy the text file above to your current folder):

.. code-block:: shell

    # Download partition function files
    wget -i wasp18b_partition_function_data.txt

    # Compute HITRAN partitions from TIPS data:
    pbay -pf tips CO

    # Convert Exomol partitions from .pf files:
    pbay -pf exomol 1H2-16O__POKAZATEL.pf
    pbay -pf exomol 12C-1H4__MM.pf
    pbay -pf exomol 12C2-1H2__aCeTY.pf
    pbay -pf exomol 51V-16O__HyVO.pf
    pbay -pf exomol 46Ti-16O__Toto.pf 47Ti-16O__Toto.pf 48Ti-16O__Toto.pf 49Ti-16O__Toto.pf 50Ti-16O__Toto.pf

    # Compute Exomol partitions from .states files:
    pbay -pf states 5.0 5000.0 5.0  1H-12C-14N__Harris.states.bz2 1H-13C-14N__Larner.states.bz2
    pbay -pf states 5.0 5000.0 5.0  14N-1H3__CoYuTe.states.bz2 15N-1H3__CoYuTe-15.states.bz2


.. note:: More info about calculating partition functions can be found
    in the :doc:`../partition_functions` tutorial.


--------------------------------------------------------

.. _tli:

3. TLI Files
------------

Now we have all the needed inputs.  Lets return to our root directory
(the one containing the ``inputs/`` folder).  The next step is to
format the line-list and partition-function data into the format for
use in ``Pyrat Bay``, these are called transmission line information
(TLI) files.

Here below is the H2O/Exomol configuration files that run this step,
for example:

.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/tli_exomol_H2O_pokazatel.cfg">tli_exomol_H2O_pokazatel.cfg</a></summary>

.. literalinclude:: ../../_static/data/tli_exomol_H2O_pokazatel.cfg
    :caption: File: tli_exomol_H2O_pokazatel.cfg

.. raw:: html

   </details>

There are two main things to configure during this step:

1. Set the wavelength range to consider.  Best practice is to include
   the full wavelength range that is available.  That way there is no
   need to recalculate TLI files for future projects at other
   wavelengths.  In `Step 4 <#cross-sections>`__ you will have the
   option to fine tune the wavelength range for the cross-sections.

2. The partition-function input is the one file determining what is the
   available temperature range.

Here are all the TLI configuration files:

.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/uhj_configs_tli.txt">uhj_configs_tli.txt</a></summary>

.. literalinclude:: ../../_static/data/uhj_configs_tli.txt
    :caption: File: uhj_configs_tli.txt
    :language: none

.. raw:: html

   </details>


Copy this file to your current folder to download the TLI
configuration files with:

.. code-block:: shell

    wget -i uhj_configs_tli.txt


Now you can compute the TLI files using this ``Pyrat Bay`` shell command:

.. code-block:: shell

    pbay -c tli_exomol_C2H2_acety.cfg
    pbay -c tli_exomol_H2O_pokazatel.cfg
    pbay -c tli_exomol_HCN_harris-larner.cfg
    pbay -c tli_exomol_NH3_coyute-byte.cfg
    pbay -c tli_exomol_TiO_toto.cfg
    pbay -c tli_exomol_VO_vomyt.cfg
    pbay -c tli_hitemp_CH4_2020.cfg
    pbay -c tli_hitemp_CO_li2019.cfg
    pbay -c tli_hitemp_CO2_2010.cfg


This may take a while since you are processing several millions of
line transitions, but once you have generated these TLI files, you
wont likely need to run this step again.

-----------------------------------------------------------

.. _cross-section:

4. Cross Sections
-----------------

The final step is to sample the line lists into tabulated data. In
this step, the code computes the strength and Voigt profile of each
transition of a given molecule, and co-adds them into a cross-section
grid (cm\ :math:`^2` molecule\ :math:`^{-1}`), which is sampled over
pressure, temperature, and wavelength.

.. Note:: Customization of the pressure, temperature, and wavelength
    sampling will vary depending on the application.  For example,
    radiative-equilibrium applications typically require a broad
    wavelength coverage (~0.3-30 µm).  In constrast, atmospheric
    retrievals may focus only on the spectral range covered by the
    observations, thus allowing to have higher spectral resolutions
    than you could with radiative-equilibrium run.  It's a
    trade-off between science requirements and computational
    constraints.

Here we will focus on a emission atmospheric retrieval for an
ultra-hot Jupiter constrained by optical and infrared observations.
Here are the configuration files to sample the cross sections:

.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg">opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg</a></summary>

.. literalinclude:: ../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
     :caption: File: opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO.cfg">opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO.cfg</a></summary>

.. literalinclude:: ../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO.cfg
     :caption: File: opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO2.cfg">opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO2.cfg</a></summary>

.. literalinclude:: ../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO2.cfg
     :caption: File: opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO2.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CH4.cfg">opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CH4.cfg</a></summary>

.. literalinclude:: ../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CH4.cfg
     :caption: File: opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CH4.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_HCN.cfg">opacity_0250-4000K_0.35-12.0um_R025K_exomol_HCN.cfg</a></summary>

.. literalinclude:: ../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_HCN.cfg
     :caption: File: opacity_0250-4000K_0.35-12.0um_R025K_exomol_HCN.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_NH3.cfg">opacity_0250-4000K_0.35-12.0um_R025K_exomol_NH3.cfg</a></summary>

.. literalinclude:: ../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_NH3.cfg
     :caption: File: opacity_0250-4000K_0.35-12.0um_R025K_exomol_NH3.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_TiO.cfg">opacity_0250-4000K_0.35-12.0um_R025K_exomol_TiO.cfg</a></summary>

.. literalinclude:: ../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_TiO.cfg
     :caption: File: opacity_0250-4000K_0.35-12.0um_R025K_exomol_TiO.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_VO.cfg">opacity_0250-4000K_0.35-12.0um_R025K_exomol_VO.cfg</a></summary>

.. literalinclude:: ../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_VO.cfg
     :caption: File: opacity_0250-4000K_0.35-12.0um_R025K_exomol_VO.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_C2H2.cfg">opacity_0250-4000K_0.35-12.0um_R025K_exomol_C2H2.cfg</a></summary>

.. literalinclude:: ../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_C2H2.cfg
     :caption: File: opacity_0250-4000K_0.35-12.0um_R025K_exomol_C2H2.cfg

.. raw:: html

   </details>


|

Lets use the H2O cross-section configuration file to walk through the
relevant parameters:

.. literalinclude:: ../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
    :caption: Extract from: opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
    :language: ini
    :lines: 3-12

This is the boilerplate indicating what to run (``runmode``), the
output file names (the output cross section file will have the same
name as ``logfile`` but as a .npz file), and ``ncpu`` sets how many
parallel CPUs you want to use (use as many as you can without crashing
your machine).


.. literalinclude:: ../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
    :caption: Extract from: opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
    :language: ini
    :lines: 14-26

While the configuration file needs to define an atmosphere, the
relevant parameters here are the ones for the pressure profile.  These
will set the layers at which we will sample the cross sections.  The
``tmodel`` and ``tpars`` are just a filler here (the temperature grid
will be defined later).  Similarly, for the composition (``species``)
we only need to take care that the molecule being sample is in the
atmospheric composition.

Note that the number of pressure layers of the cross section table
does not need to be exactly that used later in a radiative-transfer
calculation.  Here you can set a relatively coarser grid if needed
(when you run retrievals, ``Pyrat Bay`` can evaluate over a finer
pressure grid if requested).


.. literalinclude:: ../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
    :caption: Extract from: opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
    :language: ini
    :lines: 31-35

This section defines the wavelength sampling. ``wl_low`` and
``wl_high`` set the ranges (we want to cover the TESS and JWST
observing ranges), whereas ``resolution`` sets the resolving power of
the spectra (we want a resolution >= 25K to avoid having sampling
biases).  Lastly, the ``vextent`` parameter sets the extent of the
Voigt profile when sampling each line transition (this is the distance
in cm\ :sup:`-1` from the line center; for this we want at least something
> ~300--500 cm\ :sup:`-1`).


.. literalinclude:: ../../_static/data/opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
    :caption: Extract from: opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
    :language: ini
    :lines: 37-40

Then we set the temperature grid. This is a linear grid from ``tmin``
to ``tmax`` with a step size of ``tstep``.  Note that you cannot
sample beyond the temperature ranges given in the partition functions
of the inputs.  That would require extrapolation, which is not too
scientific.


|

Now you can compute the TLI files using this ``Pyrat Bay`` shell command:

.. code-block:: shell

    pbay -c opacity_0250-4000K_0.35-12.0um_R025K_exomol_C2H2.cfg
    pbay -c opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
    pbay -c opacity_0250-4000K_0.35-12.0um_R025K_exomol_HCN.cfg
    pbay -c opacity_0250-4000K_0.35-12.0um_R025K_exomol_NH3.cfg
    pbay -c opacity_0250-4000K_0.35-12.0um_R025K_exomol_TiO.cfg
    pbay -c opacity_0250-4000K_0.35-12.0um_R025K_exomol_VO.cfg
    pbay -c opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CH4.cfg
    pbay -c opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO2.cfg
    pbay -c opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO.cfg


This may take a while since you are processing several millions of
line transitions, but once you have generated these TLI files, you
wont likely need to run this step again.

-----

Once you have the cross-section files needed for yout project, go back to the :doc:`WASP-18b <../wasp18b/notebook_emission_retrieval>` retrieval notebook.
