.. |H2O| replace:: H\ :sub:`2`\ O
.. |CO2| replace:: CO\ :sub:`2`
.. |CH4| replace:: CH\ :sub:`4`
.. |H2|  replace:: H\ :sub:`2`

.. _mcmctutorial:

MCMC Tutorial
=============


This mode allows you to fit spectra to observed exoplanet data.
``Pyrat Bay`` incorporates the ``MC3`` package
(`github.com/pcubillos/MCcubed
<https://github.com/pcubillos/MCcubed>`_) to retrieve best-fitting
parameters and credible regions for the atmospheric parameters in a
Bayesian (MCMC) framework.

Here is an extract of an mcmc configuration file, showing the new
required variables:

.. code-block:: python

  [pyrat]

  # Pyrat Bay run mode, select from: [tli pt atmosphere spectrum opacity mcmc]
  runmode = mcmc
  ...
  # Filter bandpasses:
  filter = ../pyratbay/inputs/filters/tutorial/tutorial_band01.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band02.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band03.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band04.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band05.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band06.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band07.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band08.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band09.dat
           ../pyratbay/inputs/filters/tutorial/tutorial_band10.dat

  # Eclipse data:
  data =   0.000072  0.000066  0.000078  0.000120  0.000135
           0.000160  0.000196  0.000232  0.000312  0.000344
  uncert = 0.000023  0.000021  0.000020  0.000019  0.000018
           0.000017  0.000016  0.000015  0.000014  0.000014

  # Kurucz stellar model:
  kurucz = ./fp00k2odfnew.pck

  # Retrieval abundances:
  bulk     = H2 He    ; Bulk (dominant) abundance species
  molscale = H2O      ; Variable-abundance species

  # Temperature-profile model, select from [isothermal tcea madhu_inv madhu_noinv]
  tmodel = tcea

  # Retrieval models, select from: [pt rad mol ray haze]
  retflag = pt mol

  # Fitting parameters:
  #         log(kappa) log(g1) log(g2)  alpha  beta  log(fH2O)
  params   = -0.6      -0.4     0.0     0.0    1.0      0.0
  pmin     = -3.0      -0.9    -1.3     0.0    0.5     -6.0
  pmax     =  1.0       1.0     0.7     1.0    1.1      3.0
  stepsize =  0.01      0.01    0.0     0.0    0.01     0.01

  # MCMC parameters:
  walk     = snooker   ; MCMC algorithm, select from: mrw, demc, snooker
  nsamples = 50000     ; Total number of MCMC samples
  nchains  =     7     ; Number of parallel MCMC chains
  burnin   =  1000     ; Burn-in iterations per chain
  thinning =     1     ; Chains thinning factor

  # MCMC temperature boundaries  (TBD: merge with tmin/tmax)
  tlow  =  100
  thigh = 3000


.. note:: Note that an ``mcmc`` run requires the user to set an
          extinction-coefficient grid (``extfile``) to allow the code
          to finish within a Hubble time.


The observational data is input through the ``filter``, ``data``, and
``uncert`` variables, which correspond to the filter transmission
files, the eclipse or transit values (corresponding to each filter),
and the data uncertainties, respectively.

For eclipse geometry, the user needs to input a stellar flux model.
``Pyrat Bay`` currently incorporates `Kurucz models
<http://kurucz.harvard.edu/grids.html>`_ Through the ``kurucz``
variable (marcs and Phoenix TBI).  The code selects the correct Kurucz
model based on the stellar temperature (``tstar``) and surface gravity
(``gstar``).

Through the ``retflag`` parameter, the user defines which models will
vary in the retrieval.  Currently the available options are ``pt`` for
the temperature model, ``rad`` for the planetary radius at
``refpressure``, ``mol`` for the species abundances, ``ray`` for the
Rayleigh models, and ``haze`` for the haze/cloud models.

The ``params`` variable encapsulates **all** of the MCMC model
parameter into a single array.  The user is responsible for listing
the MCMC parameters in the right order, and the right number of
parameters.  The order is always: ``pt, rad, mol, ray, haze``.
The number of ``pt``, ``ray``, and ``haze`` parameters depends on the
``tmodel``, ``rayleigh``, and ``haze`` models, respectively.
The number of ``rad`` parameters is always one, planetary radius in
**kilometers**.
The number of ``mol`` parameters is the number of ``molscale``
species.

The ``molscale`` variable set the species with variable abundance.  To
do so, the code scales the whole initial species abundance profile
(:math:`q_X^0(p)`) with the abundance free parameter (:math:`f_X`) as:

.. math::   q_X(p) = q_X^0(p) \times 10^{f_X}
   :label: eq:fabundance

To preserve the sum of the mixing ratios at each layer, the code
implements the ``bulk`` variable, which sets the species used to
balance the abundances such that the mixing ratio equals one at each
layer.

The ``pmin`` and ``pmax`` variables set the boundaries for each
parameter.  The ``stepsize`` variable sets the initial random jump of
the parameters.  If ``stepsize=0`` for a given parameter, the
parameter will remain fixed at its initial value.

Finally, ``walk`` defines the MCMC sampling algorithm: Set
``walk=snooker`` (default, recommended), for the DEMC-z algorithm with
snooker propsals [BraakVrugt2008]_; or ``walk=demc`` for the
Differential-Evolution MCMC algorithm [terBraak2006]_.  ``nsamples``
sets the total number of MCMC samples, ``nchains`` sets the number of
parallel MCMC chains, ``burnin`` sets the number of removed iterations
at the beginning of each chain, and ``thinning`` the thinning factor.

Just like before, to run the MCMC modeling, simply execute this command:

.. code-block:: python

  pyrat = pb.pbay.run("tutorial_mcmc.cfg")


