// Copyright (c) 2021 Patricio Cubillos
// Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "ind.h"
#include "constants.h"
#include "utils.h"


PyDoc_STRVAR(blackbody_wn_2D__doc__,
"Compute the Planck emission function in wavenumber space:       \n\
   Bnu(T) = 2 h c**2 nu**3 / (exp(hc*nu/kT)-1),                  \n\
with units of erg s-1 sr-1 cm-2 cm.                              \n\
                                                                 \n\
Parameters                                                       \n\
----------                                                       \n\
wn: 1D float ndarray                                             \n\
    Wavenumber spectrum (cm-1).                                  \n\
temp: 1D float ndarray                                           \n\
    Temperature at each layer (K).                               \n\
B: 2D float ndarray [optional]                                   \n\
    Array to store the Planck emission of shape [nlayers, nwave].\n\
last: 1D integer ndarray [optional]                              \n\
    Indices of last layer to evaluate at each wavenumber.        \n\
                                                                 \n\
Returns                                                          \n\
-------                                                          \n\
(If B was not provided as input:)                                \n\
B: 2D float ndarray                                              \n\
    Planck emission function at wn (erg s-1 sr-1 cm-2 cm).");

static PyObject *blackbody_wn_2D(PyObject *self, PyObject *args){
    PyArrayObject *wn, *temp, *B=NULL, *last=NULL;
    int i, j, ilast, nwave, nlayers;
    npy_intp dims[2];
    double factor;

    /* Use dims as flag to detect whether B was passed as input argument: */
    dims[0] = 0;

    /* Load inputs: */
    if (!PyArg_ParseTuple(args, "OO|OO", &wn, &temp, &B, &last))
        return NULL;

    /* Get the spectrum size: */
    nwave   = (int)PyArray_DIM(wn, 0);
    nlayers = (int)PyArray_DIM(temp, 0);

    /* Initialize output array if necessary: */
    if (!B){
        dims[0] = nlayers;
        dims[1] = nwave;
        B = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    }

    /* Evaluate the Planck function: */
    for (i=0; i<nwave; i++){
        if (!last)
            ilast = nlayers - 1;
        else
            ilast = INDi(last,i);
        factor = 2 * H * LS*LS * pow(INDd(wn,i),3);
        for (j=0; j <= ilast; j++){
            IND2d(B,j,i) = factor
                           / (exp(H*LS*INDd(wn,i)/(KB*INDd(temp,j))) - 1.0);
        }
    }

    if (dims[0]==0)
        return Py_BuildValue("i", 1);
    return Py_BuildValue("N", B);
}


PyDoc_STRVAR(blackbody_wn__doc__,
"Calculate the Planck emission function in wavenumber space:\n\
   Bnu(T) = 2 h c**2 nu**3 / (exp(hc*nu/kT)-1),             \n\
with units of erg s-1 sr-1 cm-2 cm.                         \n\
                                                            \n\
Parameters                                                  \n\
----------                                                  \n\
wn: 1D float ndarray                                        \n\
    Wavenumber spectrum (cm-1).                             \n\
temp: Float                                                 \n\
    Temperature (Kelvin).                                   \n\
B: 1D float ndarray [optional]                              \n\
    Array to store the Planck emission.                     \n\
                                                            \n\
Returns                                                     \n\
-------                                                     \n\
(If B was not provided as input:)                           \n\
B: 1D float ndarray                                         \n\
    Planck emission function at wn (erg s-1 sr-1 cm-2 cm).");

static PyObject *blackbody_wn(PyObject *self, PyObject *args){
    PyArrayObject *B=NULL, *wn;
    int i;
    long nwave;
    npy_intp dims[1];
    double temp, factor;

    /* Use dims as flag to detect whether B was passed as input argument: */
    dims[0] = 0;

    /* Load inputs: */
    if (!PyArg_ParseTuple(args, "Od|O", &wn, &temp, &B))
        return NULL;

    /* Get the spectrum size: */
    nwave = PyArray_DIM(wn, 0);

    /* Initialize output array if necessary: */
    if (!B){
        dims[0] = nwave;
        B = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    }

    /* Evaluate the Planck function: */
    for (i=0; i<nwave; i++){
        factor = 2 * H * LS*LS * pow(INDd(wn,i),3);
        INDd(B,i) = factor / (exp(H*LS*INDd(wn,i)/(KB*temp)) - 1.0);
    }

    if (dims[0]==0)
        return Py_BuildValue("i", 1);
  return Py_BuildValue("N", B);
}


/* The module doc string    */
PyDoc_STRVAR(blackbody__doc__, "Wrapper for the Planck emission calculation.");

/* A list of all the methods defined by this module.                        */
static PyMethodDef blackbody_methods[] = {
    {"blackbody_wn_2D", blackbody_wn_2D, METH_VARARGS, blackbody_wn_2D__doc__},
    {"blackbody_wn",    blackbody_wn,    METH_VARARGS, blackbody_wn__doc__},
    {NULL,    NULL,  0,            NULL}  /* sentinel */
};


/* Module definition for Python 3.                                          */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_blackbody",
    blackbody__doc__,
    -1,
    blackbody_methods
};

/* When Python 3 imports a C module named 'X' it loads the module           */
/* then looks for a method named "PyInit_"+X and calls it.                  */
PyObject *PyInit__blackbody (void) {
    PyObject *module = PyModule_Create(&moduledef);
    import_array();
    return module;
}
