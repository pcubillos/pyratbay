// Copyright (c) 2016-2019 Patricio Cubillos and contributors.
// Pyrat Bay is currently proprietary software (see LICENSE).

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "ind.h"
#include "expn.h"

double
xi(double gamma, double tau){
  return 2.0/3.0 * ( (1/gamma) * (1 + (0.5*gamma*tau-1) * exp(-gamma*tau))   +
                     gamma * (1 - 0.5*pow(tau,2)) * expn(2, gamma*tau) + 1.0 );
}


PyDoc_STRVAR(TCEA__doc__,
"Generate a temperature profile based on the three-channel Eddington \n\
approximation model (Parmentier and Gillot 2014, AA 562), following  \n\
the parameterization of Line et al. (2013, ApJ 775).                 \n\
                                                                     \n\
Inputs                                                               \n\
------                                                               \n\
params: 1D float ndarray                                             \n\
  Array of free parameters:                                          \n\
   log10(kappa):  Planck thermal IR opacity in units cm^2/gr         \n\
   log10(gamma1): Visible-to-thermal stream Planck mean opacity ratio\n\
   log10(gamma2): Visible-to-thermal stream Planck mean opacity ratio\n\
   alpha:         Visible-stream partition (0.0--1.0)                \n\
   beta:          'catch-all' for albedo, emissivity, and day-night  \n\
                  redistribution (on the order of unity)             \n\
                                                                     \n\
pressure: 1D float ndarray                                           \n\
   Array of pressure values (in barye).                              \n\
R_star: Float                                                        \n\
   Stellar radius (in cm).                                           \n\
T_star: Float                                                        \n\
   Stellar effective temperature (in Kelvin degrees).                \n\
T_int:  Float                                                        \n\
   Planetary internal heat flux (in Kelvin degrees).                 \n\
sma:    Float                                                        \n\
   Semi-major axis (in cm).                                          \n\
grav:   Float                                                        \n\
   Planetary surface gravity (in cm s-2).                            \n\
                                                            \n\
Returns                                                     \n\
-------                                                     \n\
T: temperature array                                        \n\
                                                            \n\
Example                                                     \n\
-------                                                     \n\
>>> import pt as pt                                         \n\
>>> import scipy.constants as sc                            \n\
>>> import matplotlib.pyplot as plt                         \n\
>>> import numpy as np                                      \n\
                                                            \n\
>>> Rsun = 6.995e8 # Sun radius in meters                   \n\
                                                            \n\
>>> # Pressure array (barye):                               \n\
>>> press = np.logspace(2, -5, 100) *1e6                    \n\
                                                            \n\
>>> # Physical (fixed for each planet) parameters:          \n\
>>> Ts = 5040.0  # K                                        \n\
>>> Ti =  100.0  # K                                        \n\
>>> a  = 100.0 * 0.031 * sc.au # cm                         \n\
>>> Rs = 100.0 * 0.756 * Rsun  # cm                         \n\
>>> g  = 2192.8                # cm s-2                     \n\
                                                            \n\
>>> # Fitting parameters:                                   \n\
>>> kappa  = -1.5   # log10(3e-2)                           \n\
>>> gamma1 = -0.8   # log10(0.158)                          \n\
>>> gamma2 = -0.8   # log10(0.158)                          \n\
>>> alpha  = 0.5                                            \n\
>>> beta   = 1.0                                            \n\
>>> params = np.array([kappa, gamma1, gamma2, alpha, beta]) \n\
>>> T0 = pt.TCEA(params, press, Rs, Ts, Ti, a, g)           \n\
                                                            \n\
>>> plt.figure(1)                                           \n\
>>> plt.clf()                                               \n\
>>> plt.semilogy(T0, p, lw=2, color='b')                    \n\
>>> plt.ylim(press[0], press[-1])                           \n\
>>> plt.xlim(1000, 1700)                                    \n\
>>> plt.xlabel('Temperature  (K)')                          \n\
>>> plt.ylabel('Pressure  (bars)')                          \n\
                                                            \n\
Developers                                                  \n\
----------                                                  \n\
Madison Stemm      astromaddie@gmail.com                    \n\
Patricio Cubillos  pcubillos@fulbrightmail.org");

static PyObject *TCEA(PyObject *self, PyObject *args){
  PyArrayObject *freepars, *pressure, *temperature;
  double Rstar, Tstar, Tint, sma, grav,
         kappa, gamma1, gamma2, alpha, beta, tau, Tirr, xi1, xi2;
  int i, nlayers;     /* Auxilliary for-loop indices                        */
  npy_intp size[1];

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOddddd", &freepars, &pressure,
                                         &Rstar, &Tstar, &Tint, &sma, &grav))
    return NULL;

  /* Get array size:                                                        */
  size[0] = nlayers = (int)PyArray_DIM(pressure, 0);

  /* Allocate output:                                                       */
  temperature = (PyArrayObject *) PyArray_SimpleNew(1, size, NPY_DOUBLE);

  /* Unpack the model parameters:                                           */
  kappa  = pow(10.0, INDd(freepars, 0));
  gamma1 = pow(10.0, INDd(freepars, 1));
  gamma2 = pow(10.0, INDd(freepars, 2));
  alpha  = INDd(freepars, 3);
  beta   = INDd(freepars, 4);

  /* Stellar input temperature (at top of atmosphere):                      */
  Tirr = beta * sqrt(Rstar / (2.0*sma)) * Tstar;

  /* Gray IR optical depth:                                                 */
  for (i=0; i<nlayers; i++){
    tau = kappa * INDd(pressure,i) / grav;
    xi1 = xi(gamma1, tau);
    xi2 = xi(gamma2, tau);

    INDd(temperature,i) = pow(0.75 * (pow(Tint,4) * (2.0/3.0 + tau) +
                                      pow(Tirr,4) * (1-alpha) * xi1 +
                                      pow(Tirr,4) *    alpha  * xi2 ), 0.25);
  }

  return Py_BuildValue("N", temperature);
}


PyDoc_STRVAR(isothermal__doc__,
"Generate an isothermal temperature profile.\n\
                                            \n\
Inputs                                      \n\
------                                      \n\
T0: 1D float array                          \n\
   Atmospheric temperature (in Kelvin).     \n\
nlayers: integer                            \n\
   Number of atmospheric layers.            \n\
                                            \n\
Returns                                     \n\
-------                                     \n\
T: 1D float ndarray                         \n\
   Temperature profile.");

static PyObject *isothermal(PyObject *self, PyObject *args){
  PyArrayObject *T0, *temperature;
  int i, nlayers;     /* Auxilliary for-loop indices                        */
  npy_intp size[1];

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "Oi", &T0, &nlayers))
    return NULL;

  /* Get array size:                                                        */
  size[0] = nlayers;

  /* Allocate output:                                                       */
  temperature = (PyArrayObject *) PyArray_SimpleNew(1, size, NPY_DOUBLE);
  /* Set isothermal temperature:                                            */
  for (i=0; i<nlayers; i++){
    INDd(temperature,i) = INDd(T0,0);
  }
  return Py_BuildValue("N", temperature);
}


/* The module doc string                                                    */
PyDoc_STRVAR(pt__doc__,
  "Python wrapper for the temperature-profile models.");

/* A list of all the methods defined by this module.                        */
static PyMethodDef pt_methods[] = {
    {"TCEA",       TCEA,       METH_VARARGS, TCEA__doc__},
    {"isothermal", isothermal, METH_VARARGS, isothermal__doc__},
    {NULL,         NULL,       0,            NULL}
};


#if PY_MAJOR_VERSION >= 3
/* Module definition for Python 3.                                          */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "pt",
    pt__doc__,
    -1,
    pt_methods
};

/* When Python 3 imports a C module named 'X' it loads the module           */
/* then looks for a method named "PyInit_"+X and calls it.                  */
PyObject *PyInit_pt (void) {
  PyObject *module = PyModule_Create(&moduledef);
  import_array();
  return module;
}

#else
/* When Python 2 imports a C module named 'X' it loads the module           */
/* then looks for a method named "init"+X and calls it.                     */
void initpt(void){
  Py_InitModule3("pt", pt_methods, pt__doc__);
  import_array();
}
#endif
