.. include:: _substitutions.rst

.. _retrieval_tutorial:

Retrievals
==========

This section shows how to setup atmospheric retrievals with ``Pyrat Bay``.

- :ref:`ret_config`
- :ref:`ret_free_params`
- :ref:`ret_examples`

----------------------------------------------------------------------

.. _ret_config:

Setting up
----------

Multinest
~~~~~~~~~

Since version 2.0, ``Pyrat Bay`` enables atmospheric retrievals using
the nested sampling algorithm [Skilling2004]_ [Skilling2006]_, via the
MultiNest implementation [Feroz2009]_ [Buchner2014]_.  This is the
recommended retrieval algorithm.


make sure to install multinest and MPI on your machine. This can be
quite specific for each machine, so I cannot help much there. Here and
here are some installation guides that may help.

Then install their Python wrappers, e.g., with these commands:

.. code-block:: shell

    pip install pymultinest
    pip install mpi4py



Configuration file
~~~~~~~~~~~~~~~~~~

These are the requirements for the configuration file of a retrieval run:


.. tab-set::

  .. tab-item:: General
     :selected:

     .. code-block:: ini

        # Pyrat Bay run mode: [tli atmosphere spectrum opacity retrieval radeq]
        runmode = retrieval

        # Output log and spectrum file names:
        logfile = ret_wasp18b_all/WASP18b_eclipse_all.log

        # Observing geometry [transit eclipse f_lambda]
        rt_path = transit

     First, we need to set ``runmode = retrieval`` in the configuration
     file. Second, set a folder (and optionally, a path) to store
     the output files.  Lastly, specify the observing geometry.

     .. note:: ``logfile`` can point to a folder that does not exist
               yet, which will be created when launching the run.


  .. tab-item:: Observations

     .. code-block:: ini

        # The observations
        dunits = ppm
        obsfile = inputs/obs_wasp18b_eclipse_all.dat

     A retrieval run must include the ``obsfile`` argument to define
     the observed dataset (a transit, eclipse, or direct-imaging
     spectrum) which will constrain the model parameters. TBD: link to
     obsfile explanation.

     The ``dunits`` sets the units of the data *outputs* (i.e.,
     figures and screen outputs).


  .. tab-item:: Parameters

     .. code-block:: ini

        # Name      value   pmin    pmax    step  prior  prior_sigma
        retrieval_params =
            T_iso       2850.0  500.0  3000.0    1.0
            M_planet    0.27      0.1     0.4    1.0  0.266  0.033
            log_p_ref   -1.0     -9.0     2.0    1.0
            [M/H]        1.5     -2.0     2.5    1.0
            C/O          0.5     -2.0     2.5    1.0

     Define ``retrieval_params`` to set the free parameters to
     retrieve.  The full list of possible parameters is shown in the
     :ref:`ret_free_params` section.


  .. tab-item:: Sampler

     .. code-block:: ini

         # Retrieval setup:
         sampler = multinest
         nlive = 1000
         resume = True
         post_processing = True

         # Retrieval temperature boundaries:
         tlow  =  800
         thigh = 2500

     Finally, set the retrieval algorithm and other specific
     configurations.  For multinest, we need to define the number of
     live points ``nlive``.

     The ``resume`` option indicates whether to resume sampling from a
     previous run, or start from scratch.


.. _ret_launch:

Launch / multi-processing
~~~~~~~~~~~~~~~~~~~~~~~~~

Retrievals runs in single-CPU mode can be started from the command-line as:

.. code-block:: shell

    pbay -c wasp39b_retrieval_transit_jwst.cfg


Normally, we will want to use parallel computing to speed up our
runs.  Multinest runs can be parallelized via MPI.  This is the command
to do so, e.g., with 128 parallel CPUs:

.. code-block:: shell

    mpirun -n 128 pbay -c wasp39b_retrieval_transit_jwst.cfg


----------------------------------------------------------------------


.. _ret_free_params:

Retrieval parameters
--------------------

This is the list of all available free parameters
(``retrieval_params`` argument) to be fit in a retrieval. Note that
there are requirements to enable some of them.


.. raw:: html

   <style>
     table.retrieval-params {
       border-collapse: collapse;
       width: 100%;
       font-size: 0.95em;
     }
     table.retrieval-params th,
     table.retrieval-params td {
       border: 1px solid #bbb;
       padding: 3px 5px;
       vertical-align: top;
     }
     table.retrieval-params th {
       background-color: #f2f2f2;
       text-align: left;
     }
     table.retrieval-params td code {
       background-color: #f9f9f9;
       padding: 1px 4px;
       border: 1px solid #bbb;
       border-radius: 3px;
       font-family: monospace;
     }
   </style>

   <table class="retrieval-params">
     <tr>
       <th>Type</th>
       <th>Parameter</th>
       <th>Model requirement</th>
       <th>Comments</th>
     </tr>

     <tr><th rowspan="15">Temperature</th><td><code>T_iso</code></td><td><code>tmodel = isothermal</code></td><td>see <a href="atmosphere_modeling.html#temperature-profile">Isothermal TP</a> section</td></tr>
     <tr><td colspan="3" style="border-bottom: 2px solid #999; padding: 0px 0px"></td></tr>
     <tr><td><code>log_p1</code></td><td rowspan="6"><code>tmodel = madhu</code></td><td rowspan="6">see <a href="atmosphere_modeling.html#temperature-profile">Madhu TP</a> section</td></tr>
     <tr><td><code>log_p2</code></td></tr>
     <tr><td><code>log_p3</code></td></tr>
     <tr><td><code>a1</code></td>    </tr>
     <tr><td><code>a2</code></td>    </tr>
     <tr><td><code>T0</code></td>    </tr>

     <tr><td colspan="3" style="border-bottom: 2px solid #999; padding: 0px 0px"></td></tr>
     <tr><td><code>log_kappa'</code></td><td rowspan="6"><code>tmodel = guillot</code></td><td rowspan="6">see <a href="atmosphere_modeling.html#temperature-profile">Guillot TP</a> section</td></tr>
     <tr><td><code>log_gamma1</code></td></tr>
     <tr><td><code>log_gamma2</code></td></tr>
     <tr><td><code>alpha</code></td>     </tr>
     <tr><td><code>T_irr</code></td>     </tr>
     <tr><td><code>T_int</code></td>     </tr>

     <tr><td colspan="4" style="border-bottom: 3px solid #999;"></td></tr>
     <tr><th rowspan="11">Composition</th><td><code>log_X</code></td><td><code>log_X</code> in <code>vmr_vars</code></td><td>Constant VMR profile.<br><code>X</code> is a species name</td></tr>
     <tr><td><code>slope_X</code></td><td rowspan="5"><code>slant_X</code> in <code>vmr_vars</code></td><td rowspan="5"><code>X</code> is a species name</td></tr>
     <tr><td><code>log_VMR0_X</code></td></tr>
     <tr><td><code>log_p0_X</code></td>  </tr>
     <tr><td><code>min_log_X</code></td> </tr>
     <tr><td><code>max_log_X</code></td> </tr>
     <tr><td colspan="3" style="border-bottom: 2px solid #999; padding: 0px 0px"></td></tr>
     <tr><td><code>[M/H]</code></td><td rowspan="3"><code>chemistry = tea</code></td><td>metallicity scale factor</td></tr>
     <tr><td><code>[X/H]</code></td><td><code>X</code> is an element name</td></tr>
     <tr><td><code>X/Y</code></td>  <td><code>X</code> and <code>Y</code> are element names</td></tr>
     <tr><td><code>log_X</code></td><td><code>chemistry = tea</code> and<br><code>log_X</code> in <code>vmr_vars</code></td><td>Constant VMR embedded in equilibrium atmosphere. <code>X</code> is a species name</td></tr>

     <tr><td colspan="4" style="border-bottom: 3px solid #999;"></td></tr>
     <tr>
         <th rowspan="1">Isotopic ratios</th>
         <td><code>iso_X</code></td>
         <td><code>X</code> in <code>isotope_ratios</code></td>
         <td><span>\( \log_{10}({f}) \)</span>, where <span>\( f \)</span> is the isotopic fraction for isotope <code>X</code></td></tr>

     <tr><td colspan="4" style="border-bottom: 3px solid #999;"></td></tr>
     <tr><th rowspan="4">Clouds</th><td><code>log_p_cl</code></td><td><code>deck</code> in <code>clouds</code></td><td>pressure at top of opaque cloud deck.<br><span>\( \log_{10}(p/{\rm bar}) \)</span> units</td></tr>
     <tr><td><code>log_k_ray</code></td><td><code>lecavelier</code> in <code>rayleigh</code></td><td></td></tr>
     <tr><td><code>alpha_ray</code></td><td><code>lecavelier</code> in <code>rayleigh</code></td><td></td></tr>
     <tr><td><code>f_patchy</code></td><td></td><td></td></tr>

     <tr><td colspan="4" style="border-bottom: 3px solid #999;"></td></tr>
     <tr><th rowspan="6">Solo</th><td><code>R_planet</code></td><td></td><td>radius at <code>log_p_ref</code></td></tr>
     <tr><td><code>log_p_ref</code></td><td></td><td>pressure at <code>R_planet</code>.<br><span>\( \log_{10}(p/{\rm bar}) \)</span> units</td></tr>
     <tr><td><code>M_planet</code></td><td></td><td></td></tr>
     <tr><td><code>f_dilution</code></td><td></td><td>Emission dilution factor as in <a href="references.html#Taylor2020">TBD</a> </td></tr>
     <tr><td><code>T_eff</code></td><td></td><td>stellar effective temperature</td></tr>
     <tr><td><code>rv_shift</code></td><td></td><td>radial-velocity offset in km s<sup>-1</sup> </td></tr>

     <tr><td colspan="4" style="border-bottom: 3px solid #999;"></td></tr>
     <tr><th rowspan="4">Data</th><td><code>offset_X</code></td><td><code>offset_X</code> in <code>offset_inst</code></td><td>shift depth to data points containing <code>X</code> in name. Units as given in <code>dunits</code> argument</td></tr>
     <tr><td colspan="3" style="border-bottom: 2px solid #999; padding: 0px 0px;"></td></tr>
     <tr><td><code>scale_X</code></td><td><code>scale_X</code> in <code>uncert_scaling</code></td><td>scale uncertainty of data points containing <code>X</code> in name</td></tr>
     <tr><td><code>quadrature_X</code></td><td><code>quadrature_X</code> in <code>uncert_scaling</code></td><td>add noise in quadrature to data points containing <code>X</code> in name. Units as given in <code>dunits</code> argument</td></tr>
   </table>


The ``retrieval_params`` argument defines the free parameters that can
be retrieved.  See for example the extract below:


.. code-block:: ini

   #   Name       value      min     max   step  prior  prior_sigma
   retrieval_params =
       T_iso       2850    500.0  3000.0    1.0
       M_planet    0.27      0.1     0.4    1.0  0.266  0.033
       log_p_ref   -1.0     -9.0     2.0    1.0
       [M/H]        1.5     -2.0     3.0    1.0
       C/O          0.5      0.0     2.0    1.0

In the ``retrieval_params`` variable, users configure one free
parameter per row, where each input field defines (from left to
right):

- the free-parameter ``name`` (see Table above)
- an initial ``value`` (useful for forward-model calculations or when keeping a fixed value)
- the ``minimum`` and ``maximum`` values to be explored for the parameter
- a ``step`` which defines whether the parameter is free or remains fixed
- (optionally) a Gaussian ``prior`` estimation and uncertainty


----------------------------------------------------------------------

Priors
~~~~~~

By default, parameters priors are assumed to be uniform between the
minimum and maximum boundaries defined in ``retrieval_params``.  To
impose Gaussian priors to a free parameter, users can define the prior
value and uncertainties in the fields after ``step``.  In the example
above, we impossed a Gaussian prior on the planetary mass of :math:`M
= 0.266 \pm 0.033` |Mj|, while keeping unfiform priors for all other
parameters.

Similarly, users can impose uneven priors by defining both prior
uncertainties.  For example, to impose a prior mass estimation of the
form :math:`M = 0.266^{+0.3}_{-0.4}` |Mj|:

.. code-block:: ini
   :emphasize-lines: 1,4

   #   Name       value      min     max   step  prior  prior_lo  prior_hi
   retrieval_params =
       T_iso       2850    500.0  3000.0    1.0
       M_planet    0.27      0.1     0.4    1.0  0.266     0.040     0.030
       log_p_ref   -1.0     -9.0     2.0    1.0
       [M/H]        1.5     -2.0     3.0    1.0
       C/O          0.5      0.0     2.0    1.0



.. _ret_shared_params:

Fixed parameters
~~~~~~~~~~~~~~~~

Users can keep a parameter fixed by setting the ``step`` field to
zero.  In such case the parameter will adopt the provided value during
a run.  For example, to explore metallicities while keeping a solar
C/O ratio, we would set:

.. code-block:: ini
   :emphasize-lines: 7

   #   Name       value      min     max   step  prior  prior_sigma
   retrieval_params =
       T_iso       2850    500.0  3000.0    1.0
       M_planet    0.27      0.1     0.4    1.0  0.266  0.033
       log_p_ref   -1.0     -9.0     2.0    1.0
       [M/H]        0.0     -2.0     3.0    1.0
       C/O          0.59     0.0     2.0    0.0

Shared parameters
~~~~~~~~~~~~~~~~~

Users can have multiple parameters sharing the value of a single free
parameter by setting ``step`` to a negative value, which is the index
of the other parameter to share.

For example, if we want to explore abundances for an atmosphere with
two parameters, one for the alkali metals (Na and K) and one for
everything else, we can define a sodium free parameter [Na/H] (the
*fifth* row in ``retrieval_params``, below) and a potassium free
parameter [K/H] that shares the metallicity of sodium (with index
negative 5):

.. code-block:: ini
   :emphasize-lines: 8

   #   Name       value      min     max   step  prior  prior_sigma
   retrieval_params =
       T_iso       2850    500.0  3000.0    1.0
       M_planet    0.27      0.1     0.4    1.0  0.266  0.033
       log_p_ref   -1.0     -9.0     2.0    1.0
       [M/H]        0.0     -2.0     3.0    1.0
       [Na/H]       0.0     -2.0     3.0    1.0
       [K/H]        0.0     -2.0     3.0   -5.0

----------------------------------------------------------------------

.. _ret_examples:

Full examples
-------------

Here are a couple of examples to reproduce retrieval analyses from
peer-reviewed articles:

- :doc:`cookbooks/wasp39b/transmission_retrieval`
- :doc:`cookbooks/wasp18b/eclipse_retrieval`

.. TBD: High-resolution direct-imaging


