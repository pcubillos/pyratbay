// Copyright (c) 2021 Patricio Cubillos
// Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "ind.h"
#include "simpson.h"


PyDoc_STRVAR(geth__doc__,
"Calculate the differentials for a Simpson-rule integration.\n\
                                                            \n\
Parameters                                                  \n\
----------                                                  \n\
h: 1D double ndarray                                        \n\
    Intervals between the X-axis samples.                   \n\
                                                            \n\
Returns                                                     \n\
-------                                                     \n\
hsum: 1D double ndarray                                     \n\
    Sums of interval pairs:                                 \n\
    hsum = [h0+h1, h2+h3, h4+h5, ...]                       \n\
hratio: 1D double ndarray                                   \n\
    Ratio of consecutive intervals:                         \n\
    hratio = [h0/h1, h2/h3, h4/h5, ...]                     \n\
hfactor: 1D double ndarray                                  \n\
    Factor interval:                                        \n\
    hfactor = [hsum0*hsum0/h0*h1, hsum1*hsum1/h2*h3, ...]   \n\
                                                            \n\
Notes                                                       \n\
-----                                                       \n\
If there's an even number of intervals, skip the first one.");

static PyObject *geth(PyObject *self, PyObject *args){
    PyArrayObject *h, *hsum, *hratio, *hfactor;
    npy_intp size[1];
    int i, j, even=0, n;

    /* Load inputs: */
    if (!PyArg_ParseTuple(args, "O", &h))
        return NULL;

    /* Get the number of intervals: */
    n = (int)PyArray_DIM(h, 0);

    /* Empty array case: */
    if (n==0){
        return Py_BuildValue("[i,i,i]", 0, 0, 0);
    }
      
    /* Check for even number of samples (odd number of intervals): */
    even = n%2;

    /* Allocate outputs: */ 
    size[0] = n/2;
    hsum    = (PyArrayObject *) PyArray_SimpleNew(1, size, NPY_DOUBLE);
    hratio  = (PyArrayObject *) PyArray_SimpleNew(1, size, NPY_DOUBLE);
    hfactor = (PyArrayObject *) PyArray_SimpleNew(1, size, NPY_DOUBLE);

    for (i=0; i<size[0]; i++){
        j = 2*i + even;
        INDd(hsum,   i) = INDd(h,(j  )) + INDd(h,(j+1));
        INDd(hratio, i) = INDd(h,(j  )) / INDd(h,(j+1));
        INDd(hfactor,i) = INDd(hsum,i)*INDd(hsum,i) / (INDd(h,j)*INDd(h,(j+1)));
    }
    return Py_BuildValue("[N,N,N]", hsum, hratio, hfactor);
}


PyDoc_STRVAR(simps__doc__,
"Wrapper for Simpson-rule integration.                        \n\
Based on the Scipy implementation of Simpson's integration,   \n\
Licensed under a BSD 3-Clause License.                        \n\
https://github.com/scipy/scipy/blob/master/LICENSE.txt        \n\
https://github.com/scipy/scipy/blob/v0.15.1/scipy/integrate/quadrature.py\n\
                                                              \n\
Parameters                                                    \n\
----------                                                    \n\
y: 1D double ndarray                                          \n\
    Function to integrate.                                    \n\
h: 1D double ndarray                                          \n\
    Intervals between function evaluations.                   \n\
hsum: 1D double ndarray                                       \n\
    Sums of interval pairs.                                   \n\
hratio: 1D double ndarray                                     \n\
    Ratio of consecutive intervals.                           \n\
hfactor: 1D double ndarray                                    \n\
    Factor interval.                                          \n\
                                                              \n\
Returns                                                       \n\
-------                                                       \n\
integ: Float                                                  \n\
    Integral of y over intervals h using the Simpson rule.    \n\
                                                              \n\
Notes                                                         \n\
-----                                                         \n\
- If there are even samples, use a trapezoidal integration for\n\
  the last interval.                                         \n\
- See geth for formula for hsum, hratio, and hfactor");


static PyObject *simps(PyObject *self, PyObject *args){
    PyArrayObject *y, *h, *hsum, *hratio, *hfactor;
    int n, even;
    double integ=0;

    if (!PyArg_ParseTuple(args, "OOOOO", &y, &h, &hsum, &hratio, &hfactor))
        return NULL;

    n = (int)PyArray_DIM(y, 0);
    /* Check if I have an even number of samples: */
    even = n%2 == 0;

    /* Simple case, nothing to integrate: */ 
    if (n < 2)
        return Py_BuildValue("d", 0.0);
    /* Simple case, do a trapezoidal integration: */
    if (n == 2)
        return Py_BuildValue("d", INDd(h,0) * 0.5*(INDd(y,0) + INDd(y,1)));

    /* Do Simpson integration: */
    integ = simpson(y, hsum, hratio, hfactor, n);
    /* Add trapezoidal rule for last interval if n is even: */
    if (even){
        integ += INDd(h,(n-2)) * 0.5*(INDd(y,(n-2)) + INDd(y,(n-1)));
    }

    return Py_BuildValue("d", integ);
}


PyDoc_STRVAR(simps2D__doc__,
"Wrapper for Simpson-rule integration on a 2D input.\n\
Replica of scipy.simps(y, x, axis=1, even='first'),           \n\
but with a custom number of integrating point at each wavelength.\n\
                                                              \n\
Parameters                                                    \n\
----------                                                    \n\
y: 2D double ndarray                                          \n\
    Function to integrate along second dimension, array of    \n\
    shape [ny0,ny1].                                          \n\
h: 1D double ndarray                                          \n\
    Intervals between function evaluations, array of size ny1-1.\n\
nint: 1D integer ndarray                                      \n\
    Number of values to integrate for each index in first dimension,\n\
    array of size ny0.                                        \n\
hsum: 1D double ndarray                                       \n\
    Sums of interval pairs.                                   \n\
hratio: 1D double ndarray                                     \n\
    Ratio of consecutive intervals.                           \n\
hfactor: 1D double ndarray                                    \n\
    Factor interval.                                          \n\
                                                              \n\
Returns                                                       \n\
-------                                                       \n\
integ: 1D float ndarray                                       \n\
    Integral of y over intervals h using the Simpson rule.    \n\
                                                              \n\
Notes                                                         \n\
-----                                                         \n\
- If there are even samples, use a trapezoidal integration for\n\
  the last interval.                                          \n\
- See geth for formula for hsum, hratio, and hfactor");

static PyObject *simps2D(PyObject *self, PyObject *args){
    PyArrayObject *y, *h, *hsum, *hratio, *hfactor, *nint, *integ;
    int i, n, nwave, even;
    npy_intp size[1];

    /* Load inputs: */
    if (!PyArg_ParseTuple(args, "OOOOOO",
            &y, &h, &nint, &hsum, &hratio, &hfactor))
        return NULL;

    /* Length of integrand: */
    nwave = size[0] = (int)PyArray_DIM(y, 1);
    integ = (PyArrayObject *) PyArray_SimpleNew(1, size, NPY_DOUBLE);

    for (i=0; i<nwave; i++){
        n = INDi(nint, i);
        /* Check if I have an even number of samples: */
        even = n%2 == 0;

        /* Simple case, nothing to integrate: */
        if (n < 2)
            INDd(integ,i) = 0.0;
        /* Simple case, do a trapezoidal integration: */
        else if (n == 2)
            INDd(integ,i) = INDd(h,0) * 0.5*(IND2d(y,0,i) + IND2d(y,1,i));
        /* Simpson integration: */
        else{
            INDd(integ,i) = simpson2(y, hsum, hratio, hfactor, n, i);
            if (even){
                /* Add trapezoidal rule for last interval if n is even: */
                INDd(integ,i) += INDd(h,(n-2))
                                 * 0.5*(IND2d(y,(n-2),i) + IND2d(y,(n-1),i));
            }
        }
    }
    return Py_BuildValue("N", integ);
}


/* The module doc string */
PyDoc_STRVAR(simpson__doc__, "Python wrapper for Simpson integration.");


/* A list of all the methods defined by this module. */
static PyMethodDef simpson_methods[] = {
    {"geth",      geth,       METH_VARARGS, geth__doc__},
    {"simps",     simps,      METH_VARARGS, simps__doc__},
    {"simps2D",   simps2D,    METH_VARARGS, simps2D__doc__},
    {NULL,        NULL,       0,            NULL}    /* sentinel */
};


/* Module definition for Python 3. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_simpson",
    simpson__doc__,
    -1,
    simpson_methods
};

/* When Python 3 imports a C module named 'X' it loads the module */
/* then looks for a method named "PyInit_"+X and calls it.        */
PyObject *PyInit__simpson (void) {
    PyObject *module = PyModule_Create(&moduledef);
    import_array();
    return module;
}
