// Copyright (c) 2021 Patricio Cubillos
// Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "ind.h"
#include "utils.h"


PyDoc_STRVAR(
    ediff__doc__,
"Calculate the differences between consecutive elements of an array.\n\
                                              \n\
Parameters                                    \n\
----------                                    \n\
arr: 2D float ndarray                         \n\
   Input array to get the differences from.   \n\
                                              \n\
Returns                                       \n\
-------                                       \n\
diff: 1D float ndarray                        \n\
   Array wth differences.");


static PyObject *ediff(PyObject *self, PyObject *args){
    PyArrayObject *arr, *diff;
    int n, i;
    npy_intp size[1];

    /* Load inputs:                                                           */
    if (!PyArg_ParseTuple(args, "O", &arr))
      return NULL;

    size[0] = n = (int)PyArray_DIM(arr,0) - 1;
    diff = (PyArrayObject *) PyArray_SimpleNew(1, size, NPY_DOUBLE);

    for (i=0; i<n; i++)
        INDd(diff, i) = INDd(arr,(i+1)) - INDd(arr,i);
    return Py_BuildValue("N", diff);
}


PyDoc_STRVAR(
    arrbinsearch__doc__,
"Binary search of the indices in array for the closest element to\n\
each value.                                   \n\
                                              \n\
Parameters                                    \n\
----------                                    \n\
values: 1D float ndarray                      \n\
   Input values to search into array.         \n\
array: 1D float ndarray                       \n\
   Sorted array of values.differences from.   \n\
                                              \n\
Returns                                       \n\
-------                                       \n\
indices: 1D integer ndarray                   \n\
  Array indices of closest value.");


static PyObject *arrbinsearch(PyObject *self, PyObject *args){
  PyArrayObject *values, *array, *indices;
  int n, i;
  npy_intp size[1];

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OO", &values, &array))
    return NULL;

  size[0] = (int)PyArray_DIM(values, 0);
  n       = (int)PyArray_DIM(array,  0);
  indices = (PyArrayObject *) PyArray_SimpleNew(1, size, NPY_INT);

  for (i=0; i<size[0]; i++)
    INDi(indices, i) = binsearchapprox(array, INDd(values,i), 0, n);

  return Py_BuildValue("N", indices);
}


/* The module doc string    */
PyDoc_STRVAR(cutils__doc__, "Wrapper for the Planck emission calculation.");

/* A list of all the methods defined by this module.                        */
static PyMethodDef cutils_methods[] = {
    {"ediff", ediff, METH_VARARGS, ediff__doc__},
    {"arrbinsearch", arrbinsearch, METH_VARARGS, arrbinsearch__doc__},
    {NULL, NULL, 0, NULL}  // sentinel
};


/* Module definition for Python 3.                                          */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "cutils",
    cutils__doc__,
    -1,
    cutils_methods
};

/* When Python 3 imports a C module named 'X' it loads the module           */
/* then looks for a method named "PyInit_"+X and calls it.                  */
PyObject *PyInit_cutils (void) {
  PyObject *module = PyModule_Create(&moduledef);
  import_array();
  return module;
}

