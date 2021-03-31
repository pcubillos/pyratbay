// Copyright (c) 2021 Patricio Cubillos
// Pyrat Bay is open-source software under the GNU GPL-2.0 license (see LICENSE)

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "ind.h"
#include "utils.h"


PyDoc_STRVAR(trapz__doc__,
"Integrate the data array using the trapezoidal rule.       \n\
                                                            \n\
Parameters                                                  \n\
----------                                                  \n\
data: 1D double ndarray                                     \n\
    Sampled function (Y-axis) to integrate.                 \n\
intervals: 1D double ndarray                                \n\
    Intervals between the data samples (X-axis).            \n\
                                                            \n\
Returns                                                     \n\
-------                                                     \n\
res: double                                                 \n\
    The integral of data over the given intervals.          \n\
");

static PyObject *trapz(PyObject *self, PyObject *args){
    PyArrayObject *data, *intervals;
    int i, nint;
    double res=0;

    if (!PyArg_ParseTuple(args, "OO", &data, &intervals))
        return NULL;

    /* Get the number of intervals: */
    nint = (int)PyArray_DIM(intervals, 0);

    /* Empty array case: */
    if (nint < 1){
        return Py_BuildValue("d", 0.0);
    }

    for(i=0; i<nint; i++){
        res += INDd(intervals,i) * (INDd(data,(i+1)) + INDd(data,i));
    }
    return Py_BuildValue("d", 0.5*res);
}


PyDoc_STRVAR(trapz2D__doc__,
"Integrate the data array using the trapezoidal rule.       \n\
                                                            \n\
Parameters                                                  \n\
----------                                                  \n\
data: 2D double ndarray                                     \n\
    Sampled function (Y-axis) to integrate.                 \n\
intervals: 1D double ndarray                                \n\
    Intervals between the data samples (X-axis).            \n\
nint: 1D integer ndarray                                    \n\
    Number of intervals to integrate.                       \n\
                                                            \n\
Returns                                                     \n\
-------                                                     \n\
integ: 1D double ndarray                                    \n\
    The integral of data over the given intervals.          \n\
");


static PyObject *trapz2D(PyObject *self, PyObject *args){
    PyArrayObject *data, *intervals, *nint, *integ;
    int i, j, nwave;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "OOO", &data, &intervals, &nint))
        return NULL;

    nwave = dims[0] = (int)PyArray_DIM(data, 1);
    integ = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    for (j=0; j<nwave; j++){
        INDd(integ,j) = 0.0;
        for(i=0; i<INDi(nint,j); i++){
            INDd(integ,j) += INDd(intervals,i)
                           * (IND2d(data,i,j) + IND2d(data,(i+1),j));
        }
        INDd(integ,j) *= 0.5;
    }
    return Py_BuildValue("N", integ);
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
  nint = (int)PyArray_DIM(intervals, 0);

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


PyDoc_STRVAR(optdepth__doc__,
"Integrate optical depth using the trapezoidal rule.        \n\
                                                            \n\
Parameters                                                  \n\
----------                                                  \n\
data: 2D double ndarray                                     \n\
    Sampled function (Y-axis) to integrate.                 \n\
intervals: 1D double ndarray                                \n\
    Intervals between the data samples (X-axis).            \n\
taumax: Float                                               \n\
    Maximum optical depth to compute.                       \n\
ideep: 1D integer ndarray                                   \n\
    Flag of layer that reached taumax.                      \n\
ilay: Integer                                               \n\
    Current layer index.                                    \n\
                                                            \n\
Returns                                                     \n\
-------                                                     \n\
tau: 1D double ndarray                                      \n\
    2x the integral of data over the given intervals.       \n\
");

static PyObject *optdepth(PyObject *self, PyObject *args){
    PyArrayObject *data, *intervals, *ideep, *tau;
    int i, j, ilay, nint, nwave;
    double taumax;
    npy_intp dims[1];

    /* Load inputs: */
    if (!PyArg_ParseTuple(args, "OOdOi",
                          &data, &intervals, &taumax, &ideep, &ilay))
        return NULL;

    /* Get the number of intervals: */
    nint  = (int)PyArray_DIM(intervals, 0);
    nwave = (int)PyArray_DIM(data, 1);

    dims[0] = nwave;
    tau = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    for (j=0; j<nwave; j++){
        INDd(tau,j) = 0.0;

        if (INDi(ideep,j) < 0){
            for(i=0; i<nint; i++){
                INDd(tau,j) += INDd(intervals,i)
                               * (IND2d(data,(i+1),j) + IND2d(data,i,j));
            }
            /* I should now divide by two, but optical depth in   */
            /* transmission is twice the integral just calculated */
            if (INDd(tau,j) > taumax){
                INDi(ideep,j) = (int)ilay;
            }
        }
    }

    return Py_BuildValue("N", tau);
}


PyDoc_STRVAR(intensity__doc__,
"Intensity radiative-transfer integration under plane-parallel, LTE    \n\
approximation:                                                         \n\
  I = integ{B * exp(-tau/mu)} dtau/mu                                  \n\
                                                                       \n\
Parameters:                                                            \n\
-----------                                                            \n\
tau: 2D double ndarray                                                 \n\
   Optical depth as a function of altitude and wavelength.             \n\
ideep: 1D integer ndarray                                              \n\
   Bottom-layer index of the atmosphere as a function of wavelength.   \n\
planck: 2D dloble ndarray                                              \n\
   Planck blackbody emission as a funcrion of altitude and wavelength. \n\
mu: 1D float ndarray                                                   \n\
   Cosine of angles between normal and the the ray paths.              \n\
rtop: Integer                                                          \n\
   Top-layer index of the atmosphere.                                  \n\
                                                                       \n\
Returns:                                                               \n\
--------                                                               \n\
intensity: 2D double ndarray                                           \n\
   Intensity at the top of an atmosphere as a function of mu and       \n\
   wavelength.                                                         \n\
");

static PyObject *intensity(PyObject *self, PyObject *args){
  PyArrayObject *tau, *ideep, *bbody, *mu, *intensity, *dtau;
  int j, k, nwave, ntheta, rtop, last;
  double taumax;
  npy_intp idims[2], tdims[1];

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOOi", &tau, &ideep, &bbody, &mu, &rtop))
    return NULL;

  tdims[0] =          (int)PyArray_DIM(tau, 0);
  idims[1] = nwave  = (int)PyArray_DIM(tau, 1);
  idims[0] = ntheta = (int)PyArray_DIM(mu,  0);

  intensity = (PyArrayObject *) PyArray_SimpleNew(2, idims, NPY_DOUBLE);
  dtau      = (PyArrayObject *) PyArray_SimpleNew(1, tdims, NPY_DOUBLE);

  for (j=0; j<nwave; j++){
    last   = INDi(ideep, j);
    taumax = IND2d(tau,last,j);

    for (k=0; k<ntheta; k++){
      if (last-rtop == 1){
        IND2d(intensity,k,j) = IND2d(bbody,last,j);
      }else{
        /* Integral step: dtau = delta exp(-tau/mu)      */
        tdiff(dtau, tau, INDd(mu,k), rtop, last, j);
        /* Intensity trapezoidal integration:            */
        IND2d(intensity,k,j) = IND2d(bbody,last,j)*exp(-taumax/INDd(mu,k))
                               - itrapz(bbody, dtau, rtop, last, j);
      }
    }
  }

  Py_DECREF(dtau);
  return Py_BuildValue("N", intensity);
}


/* The module doc string                                                    */
PyDoc_STRVAR(trapzmod__doc__,
   "Python wrapper for trapezoidal-rule integration.");


/* A list of all the methods defined by this module.                        */
static PyMethodDef trapz_methods[] = {
    {"trapz",     trapz,      METH_VARARGS, trapz__doc__},
    {"trapz2D",   trapz2D,    METH_VARARGS, trapz2D__doc__},
    {"cumtrapz",  cumtrapz,   METH_VARARGS, cumtrapz__doc__},
    {"optdepth",  optdepth,   METH_VARARGS, optdepth__doc__},
    {"intensity", intensity,  METH_VARARGS, intensity__doc__},
    {NULL,        NULL,       0,            NULL}    /* sentinel            */
};


/* Module definition for Python 3.                                          */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_trapz",
    trapzmod__doc__,
    -1,
    trapz_methods
};

/* When Python 3 imports a C module named 'X' it loads the module           */
/* then looks for a method named "PyInit_"+X and calls it.                  */
PyObject *PyInit__trapz (void) {
  PyObject *module = PyModule_Create(&moduledef);
  import_array();
  return module;
}

