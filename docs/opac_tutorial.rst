.. |H2O| replace:: H\ :sub:`2`\ O
.. |CO2| replace:: CO\ :sub:`2`
.. |CH4| replace:: CH\ :sub:`4`
.. |H2|  replace:: H\ :sub:`2`

.. _opactutorial:

Opacity-grid Tutorial
=====================

To speed up the spectrum calculations, the ``Pyrat Bay`` can compute
cross-section tables that sample the LBL opacities over a pressure,
temperature, and wavenumber grid.  This is particularly useful when
planning to model several spectra, and is required for MCMC runs.

Sample Configuration File
-------------------------

Here is an example of an opacity-table configuration file (`opacity.cfg
<https://github.com/pcubillos/pyratbay/blob/master/examples/tutorial/opacity.cfg>`_):

.. literalinclude:: ../examples/tutorial/opacity.cfg


The ``exttable`` key sets the name of the opacity table (which is stored
as a Numpy npz file, thus, the file must have a .npz extension).  The
``tmin``, ``tmax``, and ``tstep`` keys set the boundaries and sampling
rate of the temperature grid.  The pressure grid will be taken from
the atmospheric model, and the wavenumber grid will be taken from the
wavelength/wavenumber keys in the configuration file.

The ``tlifile`` key determines the LBL opacities to include in the
opacity table.  Note that one can include opacities from multiple
species into a single opacity table.

Once this opacity table is created, runs that generate spectra will
interpolate from this table to compute extinction coefficients, which
is much faster than rather than the LBL calculations.

.. _opacity_io:

Opacity Table as Input/output
-----------------------------

Depending on the ``runmode``, whether ``extfile`` is set in the
configuration file, and whether the extinction file already exists,
``Pyrat Bay`` will operate differently.
The following table describes what the code will output:

+-----------+-----------+-------------+---------------------------------------+
|``runmode``|``extfile``| File exists | Output                                |
+===========+===========+=============+=======================================+
| opacity   | defined   | No          | Generate new opacity table            |
+           +-----------+-------------+---------------------------------------+
|           | defined   | Yes         | Overwrite existing opacity table      |
+           +-----------+-------------+---------------------------------------+
|           | undefined | ---         | Error                                 |
+-----------+-----------+-------------+---------------------------------------+
| spectrum  | defined   | Yes         | Use existing table to compute spectrum|
+           +-----------+-------------+---------------------------------------+
|           | defined   | No          | Generate grid and compute spectrum    |
+           +-----------+-------------+---------------------------------------+
|           | undefined | ---         | LBL-opacity spectrum calculation      |
+-----------+-----------+-------------+---------------------------------------+
| mcmc      | defined   | ---         | Same as ``runmode = spectrum``        |
+           +-----------+-------------+---------------------------------------+
|           | undefined | ---         | Error                                 |
+-----------+-----------+-------------+---------------------------------------+

.. note:: When the opacity table is an **output**, the configuration
   file must specify a single file for ``extfile`` (which may contain
   opacity for multiple species though).  However, when the opacity
   table is an **input**, the ``extfile`` may contain multiple opacity
   tables (as long as all of them use the same pressure, temperature,
   and wavenumber gridding).


Examples
--------


.. note:: Before running this example, make sure that you have
   generated the TLI file from the :ref:`tli_tutorial_example`,
   generated the atmospheric profiles from the
   :ref:`abundance_tutorial_example`, and download the configuration
   file shown above, e.g., with these shell commands:

   .. code-block:: shell

      tutorial_path=https://raw.githubusercontent.com/pcubillos/pyratbay/master/examples/tutorial
      wget $tutorial_path/opacity.cfg


As in a spectrum run, an opacity run returns a '*pyrat*' object in an
interactive run.  The following Python script computes an opacity file
using the configuration file found at the top of this tutorial:

.. code-block:: python

    import pyratbay as pb

    pyrat = pb.run('opacity.cfg')
