// Copyright (c) 2016-2019 Patricio Cubillos and contributors.
// Pyrat Bay is currently proprietary software (see LICENSE).

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "ind.h"
#include "utils.h"


PyDoc_STRVAR(ifirst__doc__,
"Get the first index where data is 1.                       \n\
                                                            \n\
Parameters:                                                 \n\
-----------                                                 \n\
data: 1D integer ndarray                                    \n\
   An array of (int) bools.                                 \n\
                                                            \n\
Returns:                                                    \n\
--------                                                    \n\
first: integer                                              \n\
   First index where data == 1.  Return -1 otherwise.       \n\
");

static PyObject *ifirst(PyObject *self, PyObject *args){
  PyArrayObject *data;
  int i, n;  /* Auxilliary for-loop indices                              */

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "O", &data))
    return NULL;

  /* Get the number of intervals:                                           */
  n = (int)PyArray_DIM(data, 0);

  /* Check for even number of samples (odd number of intervals):            */
  for(i=0; i<n; i++){
      if (INDi(data,i) == 1)
          return Py_BuildValue("i", i);
  }
  return Py_BuildValue("i", -1);
}


PyDoc_STRVAR(ilast__doc__,
"Get the last index where data is 1.                        \n\
                                                            \n\
Parameters:                                                 \n\
-----------                                                 \n\
data: 1D integer ndarray                                    \n\
   An array of (int) bools.                                 \n\
                                                            \n\
Returns:                                                    \n\
--------                                                    \n\
last: integer                                               \n\
   Last index where data == 1.  Return -2 otherwise.        \n\
");

static PyObject *ilast(PyObject *self, PyObject *args){
  PyArrayObject *data;
  int i, n;  /* Auxilliary for-loop indices                              */

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "O", &data))
    return NULL;

  /* Get the number of intervals:                                           */
  n = (int)PyArray_DIM(data, 0);

  /* Check for even number of samples (odd number of intervals):            */
  for(i=n-1; i>=0; i--){
      if (INDi(data,i) == 1)
          return Py_BuildValue("i", i);
  }
  return Py_BuildValue("i", -2);
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
    "trapz",
    indicesmod__doc__,
    -1,
    indices_methods
};

/* When Python 3 imports a C module named 'X' it loads the module           */
/* then looks for a method named "PyInit_"+X and calls it.                  */
PyObject *PyInit_indices (void) {
  PyObject *module = PyModule_Create(&moduledef);
  import_array();
  return module;
}

#else
/* When Python 2 imports a C module named 'X' it loads the module           */
/* then looks for a method named "init"+X and calls it.                     */
void initindices(void){
  Py_InitModule3("indices", indices_methods, indicesmod__doc__);
  import_array();
}
#endif
