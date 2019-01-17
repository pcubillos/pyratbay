// Copyright (c) 2016-2019 Patricio Cubillos and contributors.
// Pyrat Bay is currently proprietary software (see LICENSE).

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "ind.h"
#include "utils.h"

PyDoc_STRVAR(ediff__doc__,
"Calculate the differences between consecutive elements of an array.\n\
                                              \n\
Parameters                                    \n\
----------                                    \n\
arr: 2D float ndarray                         \n\
   Input array to get the differences from.   \n\
diff: 1D float ndarray                        \n\
   Array to store the differences.            \n\
n: Integer                                    \n\
   Number if elements in arr.");


static PyObject *ediff(PyObject *self, PyObject *args){
  PyArrayObject *arr, *diff;
  int n, i;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOi", &arr, &diff, &n))
    return NULL;

  for (i=0; i<n-1; i++)
    INDd(diff, i) = INDd(arr, (i+1)) - INDd(arr,i);
  return Py_BuildValue("i", 1);
}


PyDoc_STRVAR(arrbinsearch__doc__,
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
    {"ediff",        ediff,        METH_VARARGS, ediff__doc__},
    {"arrbinsearch", arrbinsearch, METH_VARARGS, arrbinsearch__doc__},
    {NULL,           NULL,         0,            NULL}          /* sentinel */
};

#if PY_MAJOR_VERSION >= 3
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

#else
/* When Python 2 imports a C module named 'X' it loads the module           */
/* then looks for a method named "init"+X and calls it.                     */
void initcutils(void){
  Py_InitModule3("cutils", cutils_methods, cutils__doc__);
  import_array();
}
#endif
