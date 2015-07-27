#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_8_API_VERSION
#include <numpy/arrayobject.h>

#include "ind.h"


PyDoc_STRVAR(trapz__doc__,
"Integrate the data array using the trapezoidal rule.       \n\
                                                            \n\
Parameters:                                                 \n\
-----------                                                 \n\
data: 1D double ndarray                                     \n\
   Sampled function (Y-axis) to integrate.                  \n\
intervals: 1D double ndarray                                \n\
   Intervals between the data samples (X-axis).             \n\
                                                            \n\
Returns:                                                    \n\
--------                                                    \n\
res: double                                                 \n\
   The integral of data over the given intervals.           \n\
");

static PyObject *trapz(PyObject *self, PyObject *args){
  PyArrayObject *data, *intervals;
  int i, nint;  /* Auxilliary for-loop indices                       */
  double res=0;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OO", &data, &intervals))
    return NULL;

  /* Get the number of intervals:                                           */
  nint = PyArray_DIM(intervals, 0);

  /* Empty array case:                                                      */
  if (nint < 1){
    return Py_BuildValue("d", 0.0);
  }

  /* Check for even number of samples (odd number of intervals):            */
  for(i=0; i<nint; i++){
    res += INDd(intervals,i) * (INDd(data,(i+1)) + INDd(data,i));
  }
  return Py_BuildValue("d", 0.5*res);
}


PyDoc_STRVAR(cumtrapz__doc__,
"Calculate the cumulative integral of a data array using the    \n\
trapezoidal-rule.  Stop calculation if it reaches the threshold.\n\
                                                                \n\
Parameters:                                                     \n\
-----------                                                     \n\
output: 1D double ndarray                                       \n\
   The output array to store the results.                       \n\
data: 1D double ndarray                                         \n\
   Sampled function (Y-axis) to integrate.                      \n\
intervals: 1D double ndarray                                    \n\
   Intervals between the data samples (X-axis).                 \n\
threshold: double                                               \n\
   Maximum threshold value to integrate.                        \n\
                                                                \n\
Returns:                                                        \n\
--------                                                        \n\
ind: integer                                                    \n\
   The index where the integration stopped.                     \n\
");

static PyObject *cumtrapz(PyObject *self, PyObject *args){
  PyArrayObject *output, *data, *intervals;
  int i, nint;       /* Auxilliary for-loop indices                         */
  double threshold;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOd", &output, &data, &intervals,
                                      &threshold))
    return NULL;

  /* Get the number of intervals:                                           */
  nint = PyArray_DIM(intervals, 0);

  /* First value is zero (zero-length interval):                            */
  INDd(output,0) = 0.0;
  /* Empty array case:                                                      */
  if (nint < 1){
    return Py_BuildValue("i", 0);
  }

  for (i=0; i<nint; i++){
    /* Integrate each interval:                                             */
    INDd(output, (i+1)) = INDd(output,i) +
            0.5*INDd(intervals,i) * (INDd(data,(i+1)) + INDd(data,i));
    /* If it reached the threshold, stop:                                   */
    if (INDd(output,(i+1)) >= threshold){
      return Py_BuildValue("i", (i+1));
    }
  }
  return Py_BuildValue("i", nint);
}


/* The module doc string                                                    */
PyDoc_STRVAR(trapzmod__doc__,
   "Python wrapper for trapezoidal-rule integration.");


/* A list of all the methods defined by this module.                        */
static PyMethodDef trapz_methods[] = {
    {"trapz",     trapz,      METH_VARARGS, trapz__doc__},
    {"cumtrapz",  cumtrapz,   METH_VARARGS, cumtrapz__doc__},
    {NULL,        NULL,       0,            NULL}    /* sentinel            */
};


#if PY_MAJOR_VERSION >= 3
/* Module definition for Python 3.                                          */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "trapz",
    trapzmod__doc__,
    -1,
    trapz_methods
};

/* When Python 3 imports a C module named 'X' it loads the module           */
/* then looks for a method named "PyInit_"+X and calls it.                  */
PyObject *PyInit_trapz (void) {
  PyObject *module = PyModule_Create(&moduledef);
  import_array();
  return module;
}

#else
/* When Python 2 imports a C module named 'X' it loads the module           */
/* then looks for a method named "init"+X and calls it.                     */
void inittrapz(void){
  Py_InitModule3("trapz", trapz_methods, trapzmod__doc__);
  import_array();
}
#endif
