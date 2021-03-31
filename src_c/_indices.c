// Copyright (c) 2021 Patricio Cubillos
// Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "ind.h"
#include "utils.h"


PyDoc_STRVAR(ifirst__doc__,
"Get the first index where data is 1.                         \n\
                                                              \n\
Parameters                                                    \n\
----------                                                    \n\
data: 1D integer ndarray                                      \n\
    An array of (int) bools.                                  \n\
default_ret: Integer                                          \n\
    Default returned value when no value in data is 1.        \n\
                                                              \n\
Returns                                                       \n\
-------                                                       \n\
first: integer                                                \n\
   First index where data == 1.  Return default_ret otherwise.\n\
                                                              \n\
Examples                                                      \n\
--------                                                      \n\
>>> import numpy as np                                        \n\
>>> print(indices.ifirst(np.array([1,0,0])))                  \n\
0                                                             \n\
>>> print(indices.ifirst(np.array([0,1,0])))                  \n\
1                                                             \n\
>>> print(indices.ifirst(np.array([0,1,1])))                  \n\
1                                                             \n\
>>> print(indices.ifirst(np.array([0,0,0])))                  \n\
-1                                                            \n\
>>> default_ret = 0                                           \n\
>>> print(indices.ifirst(np.array([0,0,0]), default_ret))     \n\
0");

static PyObject *ifirst(PyObject *self, PyObject *args){
  PyArrayObject *data;
  int i, n, default_ret=-1;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "O|i", &data, &default_ret))
      return NULL;

  /* Get the number of intervals:                                           */
  n = (int)PyArray_DIM(data, 0);

  /* Check for even number of samples (odd number of intervals):            */
  for(i=0; i<n; i++){
      if (INDi(data,i) == 1)
          return Py_BuildValue("i", i);
  }
  return Py_BuildValue("i", default_ret);
}


PyDoc_STRVAR(ilast__doc__,
"Get the last index where data is 1.                         \n\
                                                             \n\
Parameters                                                   \n\
----------                                                   \n\
data: 1D integer ndarray                                     \n\
    An array of (int) bools.                                 \n\
default_ret: Integer                                         \n\
    Default returned value when no value in data is 1.       \n\
                                                             \n\
Returns                                                      \n\
-------                                                      \n\
last: integer                                                \n\
   Last index where data == 1.  Return default_ret otherwise.\n\
                                                             \n\
Examples                                                     \n\
--------                                                     \n\
>>> import numpy as np                                       \n\
>>> print(indices.ilast(np.array([1,0,0])))                  \n\
0                                                            \n\
>>> print(indices.ilast(np.array([0,1,0])))                  \n\
1                                                            \n\
>>> print(indices.ilast(np.array([0,1,1])))                  \n\
2                                                            \n\
>>> print(indices.ilast(np.array([0,0,0])))                  \n\
-1                                                           \n\
>>> default_ret = 0                                          \n\
>>> print(indices.ilast(np.array([0,0,0]), default_ret))     \n\
0");

static PyObject *ilast(PyObject *self, PyObject *args){
  PyArrayObject *data;
  int i, n, default_ret=-1;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "O|i", &data, &default_ret))
      return NULL;

  /* Get the number of intervals:                                           */
  n = (int)PyArray_DIM(data, 0);

  /* Check for even number of samples (odd number of intervals):            */
  for(i=n-1; i>=0; i--){
      if (INDi(data,i) == 1)
          return Py_BuildValue("i", i);
  }
  return Py_BuildValue("i", default_ret);
}


/* The module doc string                                                    */
PyDoc_STRVAR(indicesmod__doc__,
    "Efficient search for indices where a condition is met.");

/* A list of all the methods defined by this module.                        */
static PyMethodDef indices_methods[] = {
    {"ifirst", ifirst, METH_VARARGS, ifirst__doc__},
    {"ilast",  ilast,  METH_VARARGS, ilast__doc__},
    {NULL,     NULL,   0,            NULL}    /* sentinel            */
};


#if PY_MAJOR_VERSION >= 3
/* Module definition for Python 3.                                          */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_indices",
    indicesmod__doc__,
    -1,
    indices_methods
};

/* When Python 3 imports a C module named 'X' it loads the module           */
/* then looks for a method named "PyInit_"+X and calls it.                  */
PyObject *PyInit__indices (void) {
    PyObject *module = PyModule_Create(&moduledef);
    import_array();
    return module;
}

#else
/* When Python 2 imports a C module named 'X' it loads the module           */
/* then looks for a method named "init"+X and calls it.                     */
void init_indices(void){
    Py_InitModule3("_indices", indices_methods, indicesmod__doc__);
    import_array();
}
#endif
