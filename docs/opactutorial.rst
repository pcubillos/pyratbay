.. |H2O| replace:: H\ :sub:`2`\ O
.. |CO2| replace:: CO\ :sub:`2`
.. |CH4| replace:: CH\ :sub:`4`
.. |H2|  replace:: H\ :sub:`2`

.. _opactutorial:

Opacity-grid Tutorial
=====================


opacity Mode
............

If you plan to generate multiple spectra for a same planet with
different atmospheric models, ``Pyrat Bay`` allows you to generate an
opacity table (as function of wavelength, temperature, and pressure)
for each species with LBL opacity data (i.e., 4D in total).

Once this grid is created, the code will interpolate the extinction
coefficient from the grid instead of repeating the line-by-line
calculations, significantly speeding up the spectrum calculation.

To create/use an extinction-coefficient grid, the configuration file
just need the following variables (in addition to a spectrum run):

.. code-block:: python

  [pyrat]

  # Pyrat Bay run mode, select from: [tli pt atmosphere spectrum opacity mcmc]
  runmode = opacity
  ...
  # Opacity file name and temperature range and step
  extfile = ./opacity_100-3000K_1.0-5.0um.dat
  tmin    =  100   ; Minimum temperature for grid
  tmax    = 3000   ; Maximum temperature for grid
  tstep   =  100   ; Temperature step for grid
  nproc   =    7   ; Number of parallel processors

The ``extfile`` variable sets the file name of the input/output
extinction-coefficient file.  The ``tmin``, ``tmax``, and ``tstep``
variables set the temperature sampling rate of the grid.  The
``nproc`` variable (default ``nproc=1``) sets the number of parallel
processors used to compute the extinction-coefficient.

The following table describes what ``Pyrat Bay`` outputs depending on
the ``runmode``, whether ``extfile`` was set in the configuration
file, and whether the extinction file already exists:

+-----------+-----------+-------------+---------------------------------------+
| Run mode  |``extfile``| File exists | Output                                |
+===========+===========+=============+=======================================+
| opacity   | defined   | No          | Generate new grid file                |
+           +-----------+-------------+---------------------------------------+
|           | defined   | Yes         | Overwrite grid file                   |
+           +-----------+-------------+---------------------------------------+
|           | undefined | ---         | Error                                 |
+-----------+-----------+-------------+---------------------------------------+
| spectrum  | defined   | No          | Generate grid and compute spectrum    |
+           +-----------+-------------+---------------------------------------+
|           | defined   | Yes         | Use existing grid to compute spectrum |
+           +-----------+-------------+---------------------------------------+
|           | undefined | ---         | LBL-opacity spectrum calculation      |
+-----------+-----------+-------------+---------------------------------------+
| mcmc      | defined   | ---         | Same as ``runmode=spectrum``          |
+           +-----------+-------------+---------------------------------------+
|           | undefined | ---         | Error                                 |
+-----------+-----------+-------------+---------------------------------------+

As always, to generate an extinction-coefficient grid, run the
following script:

.. code-block:: python

  pyrat = pb.pbay.run("tutorial_opacity.cfg")


