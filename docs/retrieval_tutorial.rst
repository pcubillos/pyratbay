.. |H2O| replace:: H\ :sub:`2`\ O
.. |CO2| replace:: CO\ :sub:`2`
.. |CH4| replace:: CH\ :sub:`4`
.. |H2|  replace:: H\ :sub:`2`

.. _retrievaltutorial:

Retrieval tutorial
==================


This mode allows you to fit spectra to observed exoplanet data.
``Pyrat Bay`` uses the ``MC3`` package (`https://github.com/pcubillos/mc3
<https://github.com/pcubillos/mc3>`_) to retrieve best-fitting
parameters and credible regions for the atmospheric parameters in a
Bayesian (MCMC) framework.


Sample Configuration File
-------------------------

Here is an example of an opacity-table configuration file (`retrieval_eclipse.cfg
<https://github.com/pcubillos/pyratbay/blob/master/examples/tutorial/retrieval_eclipse.cfg>`_):

.. literalinclude:: ../examples/tutorial/retrieval_eclipse.cfg

.. note:: Note that a '*retrieval*' run requires the user to define
          opacity table(s) (``extfile``) to allow the code to finish
          within a Hubble time.


Observed Data
-------------

A retrieval run must include the ``data``, ``uncert``, and ``filters``
keys to define the observed data that will guide the retrieval
sampler.

The ``data`` key sets the observed transit depths for transmission
observations, or eclipse depths for emission observations.  The
``uncert`` key sets the data uncertainties, and the ``filters`` key
set the path to the filter pass-band files corresponding to each data
value.  Use the ``dunits`` key to specify the units of the ``data`` and
``uncert`` values (default: ``dunits = none``).


For eclipse geometry, the user needs to input a stellar flux model to
compute eclipse depths (see :ref:`starspec`).


Retrieval Models and Parameters
-------------------------------

Use the ``retflag`` key to set which model parameters will vary in the
retrieval.  The following table list the available options:

=========== ===========================
``retflag``  Description
=========== ===========================
temp         Include the ``tpars`` temperature parameters in the retrieval
rad          Include the ``rplanet`` value as a retrieval parameter
press        Include log ``refpressure`` value as a retrieval parameter
mol          Include the ``molpars`` abundance parameters in the retrieval
ray          Include the ``rpars`` Rayleigh parameters in the retrieval
cloud        Include the ``cpars`` clouds parameters in the retrieval
patchy       Include the ``fpatchy`` value as a retrieval parameter
mass         Include the ``mplanet`` value as a retrieval parameter
=========== ===========================

Use the ``params`` key to set initial values for the retrieval
parameters.  Note that the user need to set these parameters in the
same order as listed above, i.e.: first the temperature parameters,
then radius, then the mass, the abundance parameters, and so on.

Details on the available models and their parameters are described in
the :ref:`atmospheretutorial` and :ref:`spectutorial`.  The number of
``rad``, ``press``, and ``mass`` parameters is always one each.
The units for log ``refpressure`` are always bars.
The units for ``rplanet`` and ``mplanet`` are set by the ``runits``
and ``mpunits`` keys, or else, will adopted from the units
specified in the ``rplanet`` and ``mplanet`` inputs.

.. note:: Pro-tip: You can use the ``params`` key to run spectrum
          forward models as well.

Use the ``pmin`` and ``pmax`` keys to set minimum and maximum
retrieval boundaries for each parameter listed in ``params``,
respectively.

Use the ``pstep`` key to handle the behavior of the retrieval
parameters.  A ``pstep`` value of zero keeps a parameter fixed during
the retrieval.  A negative-integer ``pstep`` value indicates that such
parameter is shared with the indexed value (i.e., a parameter with a
``pstep`` value of -2 will take the value from the second parameter).


Parameter Priors
----------------

Currently, ``Pyrat Bay`` enables uniform (default) or Gaussian priors
for the retrieval parameters.  Use the ``prior``, ``priorlow``, and
``priorup`` keys to set the mean, upper- and lower-standard deviation
value of the Gaussian priors, respectively.

If both the ``priorlow`` and ``priorup`` values are zero for a given
parameter, the prior will remain uniform between the ``pmin`` and
``pmax`` values.

For example, for the configuration file above, the following config
sets uniform priors for the temperature parameters, and a Gaussian
prior for the |H2O| abundance parameter of :math:`\log({\rm H2O}) = -3.5 \pm 0.3`:

.. code-block:: python

    #          log(kappa) log(g1) log(g2) alpha beta log(H2O)
    prior    = 0.0        0.0     0.0     0.0   0.0  -3.5
    priorlow = 0.0        0.0     0.0     0.0   0.0   0.3
    priorup  = 0.0        0.0     0.0     0.0   0.0   0.3



Retrieval Setup
---------------

``Pyrat Bay`` enables posterior sampling via the Markov-chain Monte
Carlo (MCMC) or Nested Sampling technique.  Use the ``sampler`` key to set
the sampling algorithm.   The following table list the available
options and references of their implementation:

=========== =================================== ================
``sampler`` Algorithm                           References
=========== =================================== ================
snooker     Snooker Differential-Evolution MCMC [Cubillos2017a]_
=========== =================================== ================

.. dynesty     Dynamic Nested Sampling             [Speagle2019]_

The '*snooker*' option implements the DEMC-z algorithm with
snooker proposals, described in [terBraak2008]_.

..  The '*dynesty*' option implements Dynamic Nested-sampling algorithm described in [Speagle2019]_.

MCMC Retrieval
^^^^^^^^^^^^^^

The ``nsamples`` key sets the total number of MCMC samples.

The ``nchains`` key sets the number of parallel MCMC chains to use.

The ``burnin`` key sets the number of iterations to remove (burn) at
the beginning of each chain.

.. The ``thinning`` optional key the thinning factor.

.. Gelman-Rubin stuff.


.. Nested-sampling Retrieval
   ^^^^^^^^^^^^^^^^^^^^^^^^^
   TBD.

----------------------------------------------------------------------

Examples
--------

Since this is an eclipse retrieval, the code requires a stellar
spectrum model to compute the planet-to-star flux ratios.  You will
also to download the configuration file shown above and the observing
filter files, e.g., with these shell commands:

.. code-block:: shell

   # Download Kurucz stellar model:
   wget http://kurucz.harvard.edu/grids/gridp00odfnew/fp00k2odfnew.pck

   tutorial_path=https://raw.githubusercontent.com/pcubillos/pyratbay/master/examples/tutorial
   # Download the configuration file:
   wget $tutorial_path/retrieval_eclipse.cfg

   # Download the filter files:
   for i in {0..9}
   do
   wget $tutorial_path/filters/tutorial_band0${i}.dat
   done
   mkdir filters/
   mv filters/tutorial_band0?.dat filters/


Also, be sure to have create the opacity file as described in
:ref:`opactutorial`.


As in a spectrum run, an MCMC run returns a '*pyrat*' object in an
interactive run.  The following Python script computes an opacity file
using the configuration file found at the top of this tutorial.
Just like before, to run the MCMC modeling, simply execute this command:

.. code-block:: python

    import pyratbay as pb

    pyrat = pb.run("retrieval_eclipse.cfg")
