.. _pyrat:

Pyrat Module
============

This module generates a spectrum by solving the radiative-transfer equation.


Initialization
--------------

The Initialization section reads and processes the input data, and calculates the varables that wont change.

Init
^^^^

Initialize the pyrat object.

Parse
^^^^^

Parse the input files.

Check Inputs
^^^^^^^^^^^^

Check that the user-input files make sense.

Generate Wavenumber
^^^^^^^^^^^^^^^^^^^

Calculate the spectral wavenumber array.

Atmospheric Model
^^^^^^^^^^^^^^^^^

Read the atmospheric model that contains the atmospheric species and the
radius (optional), pressure, temperature, and abundances profiles.
If the radius profile is not provided, ``pyrat`` will calculate it using
the hydrostatic-equilibrium equation.

Line-by-line Opacity
^^^^^^^^^^^^^^^^^^^^

Read the line-by-line opacity data stored in the TLI files.

Radius Profile
^^^^^^^^^^^^^^

Use the hydrostatic-equilibrium equation to calculate the layers radii
 if necessary.

Voigt Profiles
^^^^^^^^^^^^^^

Calculate a library of Voigt profiles, for a set of Doppler and Lorentz widths.
[Add Equations].

Read Cross-section
^^^^^^^^^^^^^^^^^^

Read the opacity files that contains the absorption as a function of
wavenumber and temperature.

Extinction-coefficient Table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If requested, calculate the extinction coefficient for a grid of
temperature, pressure, and wavenumber, for each species in the atmosphere.


Run
---

This section computes the spectrum.

Re-load Atmosphere
^^^^^^^^^^^^^^^^^^

If requested, update the atmospheric temperature, radius, and/or
abundances profiles.

Interpolate Cross-section Opacity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Interpolate the cross-section to the temperature of each layer.
Scale the extinction to the atmospheric number density:

.. math::
   e_{\rm CS} = \sum_j (e_{\rm interp} \prod_i \frac{n_{i}}{n_{0}})

to obtain the cross-section extinction coefficient, :math:`e_{\rm CS}`.

Optical Depth
^^^^^^^^^^^^^

First the code computes the extinction coefficient at each layer, considering the contribution fror all sources.

.. math::
  e =  e_{\rm LBL} + e_{\rm CS} + e_{\rm scat}
  :label: extintion-sum

Here is where we start making a difference between transit and eclipse geometry.

For transit geometry the code calculates the light path along
parallel rays ...
The distance traveled by a ray with impact parameter :math:`r_j`
between layers :math:`r_i` and :math:`r_{i+1}` is:

.. math::
   s_i (r_j) = \sqrt{r_i^2 - r_j^2} - \sqrt{r_{i+1}^2-r_j^2},

For eclipse geometry, the code calculates the light path for a ray diving
vertically into the planet:

.. math::
   s_i = r_{i} - r_{i+1}

Then, ``pyrat`` calculates the optical depth, :math:`\tau`.

Along impact-parameter ray, from the top of the atmospere,
:math:`R_{\rm top}`, to the closest-approach point:

.. math::
  \tau_j = \int_{R_{\rm top}}^{r_j} e {\rm d}s

Or between the top of the atmosphere and each atmospheric layer:

.. math::
  \tau_i = \int_{R_{\rm top}}^{r_i} e {\rm d}s
  :label: optical depth


Oputut Spectrum
^^^^^^^^^^^^^^^

Integrate the optical depth to compute the transmission or emission spectrum.

Modulation spectrum for transmission spectroscopy:

.. math::
   M_{\nu} = \frac{1}{R_{\rm s}^{2}} \left(R_{\rm top}^{2} +
       2 \int_{R_{\rm top}}^{R_{\rm deep}} b e^{-\tau(b)} {\rm d}b \right)
   :label: modulation

Intensity and flux-emission spectrum for eclipse spectroscopy:

.. math::
   I_{\nu} = \int_{0}^{\infty} B_{\nu} e^{-\tau} {\rm d}\tau
   :label: intensity

.. math::
   F_{\nu} = \pi \sum I_{\nu} \Delta\alpha
   :label: flux
