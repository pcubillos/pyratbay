.. _cross_sections_uhj2:


Cross Sections for an Ultra Hot Jupiter
=======================================

This documentation demonstrates how to download and process molecular
opacity line lists into tabulated cross sections. There are four main
steps to compute cross-section files:

1. `Fetch Line Lists <#line-lists>`_
2. `Fetch Partition Functions <#partition-functions>`_
3. `Compute TLI Files <#tli-files>`_
4. `Sample Cross Sections <#cross-sections>`_

Note that the first three steps are typically executed only once,
allowing you to reuse the output files across projects. However, `Step
4 <#cross-sections>`_ may need to be repeated on a per-project
basis, depending on your specific requirements (e.g., different
spectral, temperature, or pressure ranges; or varying resolutions).

Now go back to the :doc:`WASP-18b <../wasp18b/notebook_emission_retrieval>` retrieval notebook.

--------------------------------------------------------

1. Line Lists
-------------

In this section, we'll download molecular line lists, typically
sourced from the ExoMol or HITRAN/HITEMP databases. You will likely
only need to complete this step once unless a new or updated line list
becomes available. It may be better to store this data in a general
directory on your machine. To create a folder for storing line lists,
run:

.. code-block:: shell

   mkdir inputs
   cd inputs

For this project, we will focus on the molecular absorbers relevant
for an ultra-hot Jupiter (WASP-18b). The table below lists the
molecular line-lists to download and their sources.


+----------+--------+--------------------------+
| Molecule | Source | Line List / Reference    |
+==========+========+==========================+
| CH4      | HITEMP | Hargreaves et al. (2020) |
+----------+--------+--------------------------+
| CO       | HITEMP | Li et al. (2019)         |
+----------+--------+--------------------------+
| CO2      | HITEMP | Rothman et al. (2010)    |
+----------+--------+--------------------------+
| H2O      | ExoMol | pokazatel                |
+----------+--------+--------------------------+
| HCN      | ExoMol | larner/harris            |
+----------+--------+--------------------------+
| NH3      | ExoMol | coyute/byte              |
+----------+--------+--------------------------+
| TiO      | ExoMol | toto                     |
+----------+--------+--------------------------+
| VO       | ExoMol | vomyt                    |
+----------+--------+--------------------------+
| C2H2     | ExoMol | acety                    |
+----------+--------+--------------------------+

The file below contains links to download all the required data.


.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="uhj_line_lists_data.txt">uhj_line_lists_data.txt</a></summary>

.. literalinclude:: uhj_line_lists_data.txt
    :caption: File: uhj_line_lists_data.txt
    :language: none

.. raw:: html

   </details>

Note that for ExoMol data, we will fetch line lists processed with
``repack`` (`Cubillos 2017, ApJ 850
<https://ui.adsabs.harvard.edu/abs/2017ApJ...850...32C>`_). This
package identifies the strong lines dominating the spectrum from the
weak ones, which get discarded, speeding up the sampling process by
reducing the line lists from billions of transitions to only a few
hundred million.

On Linux/OSX, you can copy this file and then download the line-list
data using the ``wget`` shell command (note these are several GB of
data):

.. code-block:: shell

    wget -i uhj_line_lists_data.txt


Now, unpack the HITEMP data with these shell commands:

.. code-block:: shell

    bzip2 -d 05_HITEMP2019.par.bz2
    bzip2 -d 06_HITEMP2020.par.bz2
    unzip '*.zip'
    rm -f *.zip

--------------------------------------------------------


2. Partition Functions
----------------------

In addition to the line-list data, to compute cross sections you will
need the partition functions for each molecules.  This file below
contains the links to the partition functions to extract from the
ExoMol database (the rest we will source from HITRAN).


.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="partition_function_data.txt">partition_function_data.txt</a></summary>

.. literalinclude:: partition_function_data.txt
    :caption: File: partition_function_data.txt
    :language: none

.. raw:: html

   </details>


Copy this file to your ``inputs/`` folder and then download the
partition-function files with this shell command:

.. code-block:: shell

    wget -i partition_function_data.txt


Now we need to format the ExoMol partition function files into the
right format for ``Pyrat Bay``.  For that run this shell commands:

.. code-block:: shell

    pbay -pf exomol 1H2-16O__POKAZATEL.pf
    pbay -pf exomol 1H-12C-14N__Harris.pf 1H-13C-14N__Larner.pf
    pbay -pf exomol 12C2-1H2__aCeTY.pf
    pbay -pf exomol 46Ti-16O__Toto.pf 47Ti-16O__Toto.pf 48Ti-16O__Toto.pf 49Ti-16O__Toto.pf 50Ti-16O__Toto.pf
    pbay -pf exomol 51V-16O__VOMYT.pf


For the other molecules, we will use the HITRAN partition-function
data (Gamache et al. `2017
<https://ui.adsabs.harvard.edu/abs/2017JQSRT.203...70G>`_, `2021
<https://ui.adsabs.harvard.edu/abs/2021JQSRT.27107713G>`_), which are
readily availabel in ``Pyrat Bay`` (no need to download files).  To
generate the partition function run the following shell commands:

.. code-block:: shell

    pbay -pf tips CO
    pbay -pf tips CO2
    pbay -pf tips CH4
    pbay -pf tips NH3 as_exomol

Note that for NH3 we are using the HITRAN partition functions for the
ExoMol line list (because this partition function samples up to 6000K,
which we need for atmospheres of ultra-hot Jupiters).  Thus the
``as_exomol`` argument makes the ouput file to label the isotope names
as in the ExoMol format.

--------------------------------------------------------

.. _tli:

3. TLI Files
------------

Now we have all the needed inputs.  Lets return to our root directory
(the one containing the ``inputs/`` folder).

The next step is to format the line-list and partition-function input
data into the format for use in ``Pyrat Bay``, these are called
transmission line information (TLI) files.

Here below is the H2O/Exomol configuration files that run this step,
for example:

.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="tli_exomol_H2O_pokazatel.cfg">tli_exomol_H2O_pokazatel.cfg</a></summary>

.. literalinclude:: tli_exomol_H2O_pokazatel.cfg
    :caption: File: tli_exomol_H2O_pokazatel.cfg

.. raw:: html

   </details>

A couple of things to note:

- The configuration file indicates the wavelength range to
  consider. Best practice is to include the full wavelength range
  available from the line list.  That way you can create a single TLI
  file that you can use for all of your future projects.  In `Step 4
  <#cross-sections>`_ you will have the option to fine tune the
  wavelength range for specific projects.

- The partition-function input is the one file determining what is the
  available temperature range.

Here are all the TLI configuration files:

.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="config_files_tli.txt">config_files_tli.txt</a></summary>

.. literalinclude:: config_files_tli.txt
    :caption: File: config_files_tli.txt
    :language: none

.. raw:: html

   </details>


Copy this file to your current folder to download the TLI
configuration files with:

.. code-block:: shell

    wget -i config_files_tli.txt


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
this step, each line transition is processed to compute its Voigt
profile, which is then sampled over a specified wavelength range, and
coadded with all other lines for the molecule. We do this across a
regular grid of temperatures and pressures, enabling later use in
radiative-transfer calculations.

.. Note::
    Depending on the application for the cross-section data, you
    may need to set specific parameters. For example,
    radiative-equilibrium applications typically require broad
    wavelength coverage (~0.3–30 µm) to capture the spectral regions
    where most of the stellar and planetary flux is concentrated.  In
    constrast, atmospheric retrievals may focus only on the spectral
    range covered by the observations, thus allowing to have higher
    spectral resolutions than you could with radiative-equilibrium
    run.  It's all a trade-off between science requirements and
    computational constraints.

Here we will focus on a emission atmospheric retrieval for an
ultra-hot Jupiter.

.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg">opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg</a></summary>

.. literalinclude:: opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
:caption: File: opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO.cfg">opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO.cfg</a></summary>

.. literalinclude:: opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO.cfg
:caption: File: opacity_0250-4000K_0.35-12.0um_R025K_hitemp_CO.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="opacity_0250-4000K_0.35-12.0um_R020K_hitemp_CO2.cfg">opacity_0250-4000K_0.35-12.0um_R020K_hitemp_CO2.cfg</a></summary>

.. literalinclude:: opacity_0250-4000K_0.35-12.0um_R020K_hitemp_CO2.cfg
:caption: File: opacity_0250-4000K_0.35-12.0um_R020K_hitemp_CO2.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="opacity_0250-4000K_0.35-12.0um_R020K_hitemp_CH4.cfg">opacity_0250-4000K_0.35-12.0um_R020K_hitemp_CH4.cfg</a></summary>

.. literalinclude:: opacity_0250-4000K_0.35-12.0um_R020K_hitemp_CH4.cfg
:caption: File: opacity_0250-4000K_0.35-12.0um_R020K_hitemp_CH4.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="opacity_0250-4000K_0.35-12.0um_R020K_exomol_HCN.cfg">opacity_0250-4000K_0.35-12.0um_R020K_exomol_HCN.cfg</a></summary>

.. literalinclude:: opacity_0250-4000K_0.35-12.0um_R020K_exomol_HCN.cfg
:caption: File: opacity_0250-4000K_0.35-12.0um_R020K_exomol_HCN.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="opacity_0250-4000K_0.35-12.0um_R020K_exomol_NH3.cfg">opacity_0250-4000K_0.35-12.0um_R020K_exomol_NH3.cfg</a></summary>

.. literalinclude:: opacity_0250-4000K_0.35-12.0um_R020K_exomol_NH3.cfg
:caption: File: opacity_0250-4000K_0.35-12.0um_R020K_exomol_NH3.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="opacity_0250-4000K_0.35-12.0um_R020K_exomol_TiO.cfg">opacity_0250-4000K_0.35-12.0um_R020K_exomol_TiO.cfg</a></summary>

.. literalinclude:: opacity_0250-4000K_0.35-12.0um_R020K_exomol_TiO.cfg
:caption: File: opacity_0250-4000K_0.35-12.0um_R020K_exomol_TiO.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="opacity_0250-4000K_0.35-12.0um_R020K_exomol_VO.cfg">opacity_0250-4000K_0.35-12.0um_R020K_exomol_VO.cfg</a></summary>

.. literalinclude:: opacity_0250-4000K_0.35-12.0um_R020K_exomol_VO.cfg
:caption: File: opacity_0250-4000K_0.35-12.0um_R020K_exomol_VO.cfg

.. raw:: html

   </details>
.. raw:: html

   <details>
   <summary>Click here to show/hide: <a href="opacity_0250-4000K_0.35-12.0um_R020K_exomol_C2H2.cfg">opacity_0250-4000K_0.35-12.0um_R020K_exomol_C2H2.cfg</a></summary>

.. literalinclude:: opacity_0250-4000K_0.35-12.0um_R020K_exomol_C2H2.cfg
:caption: File: opacity_0250-4000K_0.35-12.0um_R020K_exomol_C2H2.cfg

.. raw:: html

   </details>

|

Lets use the H2O cross-section configuration file to walk through the
relevant parameters:

.. literalinclude:: opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
    :caption: Extract from: opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
    :language: ini
    :lines: 3-12

This is the boilerplate indicating what to run (``runmode``), the
output file names (the output cross section file will have the same
name as ``logfile`` but as a .npz file), and ``ncpu`` sets how many
parallel CPUs you want to use (use as many as you can without crashing
your machine).


.. literalinclude:: opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
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


.. literalinclude:: opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
    :caption: Extract from: opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
    :language: ini
    :lines: 31-35

This section defines the wavelength sampling. ``wllow`` and
``wlhigh`` set the ranges (we want to cover the TESS and JWST
observing ranges), whereas ``resolution`` sets the resolving power of
the spectra (we want a resolution >= 25K to avoid having sampling
biases).  Lastly, the ``vextent`` parameter sets the extent of the
Voigt profile when sampling each line transition (this is the distance
in cm\ :sup:`-1` from the line center; for this we want at least something
> ~300--500 cm\ :sup:`-1`).


.. literalinclude:: opacity_0250-4000K_0.35-12.0um_R025K_exomol_H2O.cfg
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
