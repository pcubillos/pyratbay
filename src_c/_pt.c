// Copyright (c) 2021 Patricio Cubillos
// Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "ind.h"
#include "expn.h"

double
xi(double gamma, double tau){
    return 2.0/3.0 * ( (1/gamma) * (1 + (0.5*gamma*tau-1) * exp(-gamma*tau)) +
                     gamma * (1 - 0.5*pow(tau,2)) * expn(2, gamma*tau) + 1.0 );
}


PyDoc_STRVAR(guillot__doc__,
"Generate a temperature profile based on the three-channel Eddington \n\
approximation model (Gillot 2010, AA 520), following the             \n\
the parameterization of Line et al. (2013, ApJ 775).                 \n\
                                                                     \n\
Inputs                                                               \n\
------                                                               \n\
params: 1D float ndarray                                             \n\
  Array of free parameters:                                          \n\
   log10(kappa):  Planck thermal IR opacity in units cm^2/gr         \n\
   log10(gamma1): Visible-to-thermal stream Planck mean opacity ratio\n\
   log10(gamma2): Visible-to-thermal stream Planck mean opacity ratio\n\
   alpha: Weight partition between visible-streams (0.0--1.0)        \n\
   t_irr: Stellar irradiation temperature (Kelvin degree)            \n\
          A good guess is the planet's equilibrium temperature.      \n\
   t_int: Planetary internal heat flux (in Kelvin degrees)           \n\
                                                                     \n\
pressure: 1D float ndarray                                           \n\
    Atmospheric pressure array (barye).                              \n\
grav: 1D flaot ndarray (optional)                                    \n\
    Atmospheric gravity at each pressure layer (cm s-2).             \n\
                                                            \n\
Returns                                                     \n\
-------                                                     \n\
temperature: 1D float ndarray                               \n\
    Temperature at each pressure layer (Kelvin degree).     \n\
                                                            \n\
Example                                                     \n\
-------                                                     \n\
>>> import _pt as pt                                        \n\
>>> import matplotlib.pyplot as plt                         \n\
>>> import numpy as np                                      \n\
                                                            \n\
>>> # Pressure (barye) and gravity (cm s-2) arrays:         \n\
>>> nlayers = 20                                            \n\
>>> press = np.logspace(-8, 2, nlayers) * 1e6               \n\
>>> grav = np.tile(2200.0, nlayers)                         \n\
                                                            \n\
>>> # Fitting parameters:                                   \n\
>>> kappa  = -1.5                                           \n\
>>> gamma1 = -0.8                                           \n\
>>> gamma2 = -0.8                                           \n\
>>> alpha  = 0.5                                            \n\
>>> t_irr = 1200.0                                          \n\
>>> t_int =  100.0                                          \n\
>>> params = np.array([kappa, gamma1, gamma2, alpha, t_irr, t_int])\n\
>>> temp = pt.guillot(params, press, grav)                  \n\
>>> print(temp)                                             \n\
[1046.89057361 1046.8906572  1046.89094586 1046.89194204 1046.89537746\n\
 1046.9072167  1046.94798848 1047.08827618 1047.57026844 1049.22033981\n\
 1054.80921126 1073.11749958 1127.5360275  1256.04683354 1458.34379995\n\
 1623.82740006 1659.07947584 1659.7176149  1660.94856336 1665.06440703]");

static PyObject *guillot(PyObject *self, PyObject *args){
    PyArrayObject *freepars, *pressure, *temperature, *gravity=NULL;
    double kappa, gamma1, gamma2, alpha, t_irr, t_int, tau, xi1, xi2, g;
    int i, nlayers;
    npy_intp size[1];

    /* Load inputs: */
    if (!PyArg_ParseTuple(args, "OO|O", &freepars, &pressure, &gravity))
        return NULL;

    /* Get array size: */
    size[0] = nlayers = (int)PyArray_DIM(pressure, 0);

    /* Allocate output: */
    temperature = (PyArrayObject *) PyArray_SimpleNew(1, size, NPY_DOUBLE);

    /* Unpack the model parameters: */
    kappa  = pow(10.0, INDd(freepars, 0));
    gamma1 = pow(10.0, INDd(freepars, 1));
    gamma2 = pow(10.0, INDd(freepars, 2));
    alpha  = INDd(freepars, 3);
    t_irr  = INDd(freepars, 4);
    t_int  = INDd(freepars, 5);

    /* Gray IR optical depth: */
    for (i=0; i<nlayers; i++){
        if (!gravity)
            g = 1.0;
        else
            g = INDd(gravity,i);
        tau = kappa * INDd(pressure,i) / g;
        xi1 = xi(gamma1, tau);
        xi2 = xi(gamma2, tau);

        INDd(temperature,i) = pow(
            0.75 * (pow(t_int,4) * (2.0/3.0 + tau) +
                    pow(t_irr,4) * (1-alpha) * xi1 +
                    pow(t_irr,4) *    alpha  * xi2 ), 0.25);
    }

    return Py_BuildValue("N", temperature);
}


PyDoc_STRVAR(isothermal__doc__,
"Generate an isothermal temperature profile.\n\
                                            \n\
Inputs                                      \n\
------                                      \n\
T0: 1D float array                          \n\
    Atmospheric temperature (in Kelvin).    \n\
nlayers: integer                            \n\
    Number of atmospheric layers.           \n\
                                            \n\
Returns                                     \n\
-------                                     \n\
T: 1D float ndarray                         \n\
    Temperature profile.");

static PyObject *isothermal(PyObject *self, PyObject *args){
    PyArrayObject *T0, *temperature;
    int i, nlayers;
    npy_intp size[1];

    /* Load inputs: */
    if (!PyArg_ParseTuple(args, "Oi", &T0, &nlayers))
      return NULL;

    size[0] = nlayers;
    temperature = (PyArrayObject *) PyArray_SimpleNew(1, size, NPY_DOUBLE);
    /* Set isothermal temperature: */
    for (i=0; i<nlayers; i++){
        INDd(temperature,i) = INDd(T0,0);
    }
    return Py_BuildValue("N", temperature);
}


/* The module doc string */
PyDoc_STRVAR(
    pt__doc__,
    "Python wrapper for the temperature-profile models."
);

/* A list of all the methods defined by this module. */
static PyMethodDef pt_methods[] = {
    {"guillot", guillot, METH_VARARGS, guillot__doc__},
    {"isothermal", isothermal, METH_VARARGS, isothermal__doc__},
    {NULL, NULL, 0, NULL}
};


/* Module definition for Python 3. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_pt",
    pt__doc__,
    -1,
    pt_methods
};

/* When Python 3 imports a C module named 'X' it loads the module */
/* then looks for a method named "PyInit_"+X and calls it.        */
PyObject *PyInit__pt (void) {
    PyObject *module = PyModule_Create(&moduledef);
    import_array();
    return module;
}
