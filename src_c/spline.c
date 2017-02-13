// Copyright (c) 2016-2017 Patricio Cubillos and contributors.
// Pyrat Bay is currently proprietary software (see LICENSE).

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_8_API_VERSION
#include <numpy/arrayobject.h>
#include <math.h>

#include "ind.h"
#include "utils.h"
#include "spline.h"

PyDoc_STRVAR(splinterp__doc__,
"Cubic spline interpolation.\n\
                                                           \n\
Parameters                                                 \n\
----------                                                 \n\
yin: 1D double ndarray                                     \n\
   Input y values.                                         \n\
xin: 1D double ndarray                                     \n\
   Input x array.                                          \n\
xout: 1D double ndarray                                    \n\
   X array for inpterpolated values.                       \n\
extrap: double                                             \n\
   Value assigned to extrapolated points.                  \n\
                                                           \n\
Return                                                     \n\
------                                                     \n\
yout: 1D float ndarray                                     \n\
   Spline interpolated values.                             \n\
                                                           \n\
Notes                                                      \n\
-----                                                      \n\
This code is based on BART spline code written by          \n\
  Ryan Challener and Sarah Blumenthal.");

static PyObject *splinterp(PyObject *self, PyObject *args){
  PyArrayObject *xin, *yin, *xout, *yout, *dx, *z;
  npy_intp nin[1], nout[1];  /* Size of arrays                              */
  int i,                     /* For-loop counter                            */
      epleft=0, epright;     /* Number of extrapolating point on left/right */
  double extrap;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOd", &yin, &xin, &xout, &extrap))
    return NULL;

  /* Get the number of datapoints:                                          */
  nin [0] = PyArray_DIM(xin,  0)-1;
  nout[0] = PyArray_DIM(xout, 0);

  /* Allocate arrays:                                                       */ 
  dx   = (PyArrayObject *) PyArray_SimpleNew(1, nin,  NPY_DOUBLE);
  nin[0] += 1;
  z    = (PyArrayObject *) PyArray_SimpleNew(1, nin,  NPY_DOUBLE);
  yout = (PyArrayObject *) PyArray_SimpleNew(1, nout, NPY_DOUBLE);

  for (i=0; i<nin[0]-1; i++){
    INDd(dx, i) = INDd(xin,(i+1)) - INDd(xin,i);
  }

  /* Calculate second derivatives:                                          */
  tri(z, yin, dx, (int)nin[0]);

  /* Extrapolated values:                                                   */
  while(INDd(xout, epleft)  < INDd(xin,0)){
    INDd(yout, epleft) = extrap;
    epleft += 1;
  }
  epright = (int)nout[0]-1;
  while(INDd(xout, epright) > INDd(xin,(nin[0]-1))){
    INDd(yout,epright) = extrap;
    epright -= 1;
  }

  /* Number of points to interpolate:                                       */
  nout[0] = epright - epleft + 1;
  /* Interpolate:                                                           */
  spline3(xin, yin, (int)nin[0], xout, yout, (int)nout[0], z, dx, epleft);

  Py_DECREF(dx);
  Py_DECREF(z);

  return Py_BuildValue("N", yout);
}


PyDoc_STRVAR(splinterp_2D__doc__,
"Cubic spline interpolation for 2D array along the first axis,\n\
Repeating across second axis.                              \n\
                                                           \n\
Parameters                                                 \n\
----------                                                 \n\
yin: 2D double ndarray                                     \n\
   Input y values of shape [nin, nin2].                    \n\
xin: 1D double ndarray                                     \n\
   Input x array of length nin.                            \n\
z: 2D double ndarray                                       \n\
   Second derivatives of yin at xin (of shape [nin, nin2]).\n\
xout: 1D double ndarray                                    \n\
   X valuer where to interpolate (of length nout).         \n\
yout: 2D double ndarray                                    \n\
   X valuer where to interpolate (of shape [nout,nin2]).   \n\
lo: Integer                                                \n\
   Lower index along second axis of yin to interpolate.    \n\
hi: Integer                                                \n\
   Upper index along second axis of yin to interpolate.    \n\
                                                           \n\
Notes                                                      \n\
-----                                                      \n\
This code is based on BART spline code written by          \n\
  Ryan Challener and Sarah Blumenthal.");

static PyObject *splinterp_2D(PyObject *self, PyObject *args){
  PyArrayObject *xin, *yin, *z, *xout, *yout;
  int i, j, nin, nout, index, lo, hi;
  double dx, dy, deltax, a, b, c;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOOOii",
                        &yin, &xin, &z, &xout, &yout, &lo, &hi))
    return NULL;

  nin  = (int)PyArray_DIM(xin, 0);
  nout = (int)PyArray_DIM(xout,0);

  for (i=0; i<nout; i++){
    if (INDd(xout,i) < INDd(xin,0) || INDd(xout,i) > INDd(xin,(nin-1)))
      return Py_BuildValue("d", NAN);

    /* Binary search to find index:                                         */
    index = binsearchapprox(xin, INDd(xout,i), 0, nin-1);
    /* Enforce: x[i] <= xout (except if x[N-1] == xout):                    */
    if (index == nin-1 || INDd(xout,i) < INDd(xin,index)){
      index--;
    }

    /* If the requested x falls on a given x, no interpolation needed:      */
    if (INDd(xin,index) == INDd(xout,i)){
      for (j=lo; j<hi; j++){
        IND2d(yout,i,j) = IND2d(yin,index,j);
      }
      continue;
    }
    /* Else, spline-interpolate:                                            */
    dx = INDd(xin,(index+1)) - INDd(xin,index);
    deltax = INDd(xout,i) - INDd(xin,index);
    for (j=lo; j<hi; j++){
      dy = IND2d(yin,(index+1),j) - IND2d(yin,index,j);
      a = (IND2d(z,(index+1),j) - IND2d(z,index,j)) / (6*dx);
      b = 0.5 * IND2d(z,index,j);
      c = dy/dx - dx/6 * (IND2d(z,(index+1),j) + 2*IND2d(z,index,j));
      IND2d(yout,i,j) = IND2d(yin,index,j) + deltax*(c + deltax*(b + deltax*a));
    }
  }

  return Py_BuildValue("d", 0.0);
}


PyDoc_STRVAR(splinterp_pt__doc__,
"Cubic spline interpolation for a single point.            \n\
                                                           \n\
Parameters                                                 \n\
----------                                                 \n\
yin: 1D double ndarray                                     \n\
   Input y values.                                         \n\
xin: 1D double ndarray                                     \n\
   Input x array.                                          \n\
z: 1D double ndarray                                       \n\
   Second derivatives of yin at xin.                       \n\
nin: Integer                                               \n\
   Length of yin.                                          \n\
xout: Float                                                \n\
   X valuer where to interpolate.                          \n\
                                                           \n\
Return                                                     \n\
------                                                     \n\
yout: Float                                                \n\
   Spline interpolated value.                              \n\
                                                           \n\
Notes                                                      \n\
-----                                                      \n\
This code is based on BART spline code written by          \n\
  Ryan Challener and Sarah Blumenthal.");

static PyObject *splinterp_pt(PyObject *self, PyObject *args){
  PyArrayObject *xin, *yin, *z;
  int nin, index;
  double xout, yout, dx, dy, deltax, a, b, c;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOid", &yin, &xin, &z, &nin, &xout))
    return NULL;

  if (xout < INDd(xin,0) || xout > INDd(xin,(nin-1)))
    return Py_BuildValue("d", NAN);

  /* Binary search to find index:                                           */
  index = binsearchapprox(xin, xout, 0, nin-1);
  /* Enforce: x[i] <= xout (except if x[N-1] == xout):                      */
  if (index == nin-1 || xout < INDd(xin,index)){
    index--;
  }

  /* If the requested x falls on a given x, no interpolation is necessary:  */
  if (INDd(xin,index) == xout)
    return Py_BuildValue("d", INDd(yin,index));

  /* Calculate range of area in question:                                   */
  dx = INDd(xin,(index+1)) - INDd(xin,index);
  dy = INDd(yin,(index+1)) - INDd(yin,index);

  /* Else, spline-interpolate:                                              */
  deltax = xout - INDd(xin,index);
  a = (INDd(z,(index+1)) - INDd(z,index)) / (6*dx);
  b = 0.5 * INDd(z,index);
  c = dy/dx - dx/6 * (INDd(z,(index+1)) + 2*INDd(z,index));
  yout = INDd(yin,index) + deltax*(c + deltax*(b + deltax*a));

  return Py_BuildValue("d", yout);
}


PyDoc_STRVAR(spline_init__doc__,
"Cubic spline interpolation for a single point.            \n\
                                                           \n\
Parameters                                                 \n\
----------                                                 \n\
yin: 1D double ndarray                                     \n\
   Input y values.                                         \n\
xin: 1D double ndarray                                     \n\
   Input x array.                                          \n\
                                                           \n\
Return                                                     \n\
------                                                     \n\
z: 1D double ndarray                                       \n\
   Second derivatives of yin at xin.                       \n\
                                                           \n\
Notes                                                      \n\
-----                                                      \n\
This code is based on BART spline code written by          \n\
  Ryan Challener and Sarah Blumenthal.");

static PyObject *spline_init(PyObject *self, PyObject *args){
  PyArrayObject *xin, *yin, *dx, *z;
  npy_intp nin[1];  /* Size of arrays                                       */
  int i;            /* For-loop counter                                     */

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OO", &yin, &xin))
    return NULL;

  /* Array size - 1:                                                        */
  nin[0] = PyArray_DIM(yin,0) - 1;

  /* Allocate arrays:                                                       */
  dx   = (PyArrayObject *) PyArray_SimpleNew(1, nin,  NPY_DOUBLE);
  nin[0] += 1;
  z    = (PyArrayObject *) PyArray_SimpleNew(1, nin,  NPY_DOUBLE);

  for (i=0; i<nin[0]-1; i++){
    INDd(dx, i) = INDd(xin,(i+1)) - INDd(xin,i);
  }
  /* Calculate second derivatives:                                          */
  tri(z, yin, dx, (int)nin[0]);

  /* Free arrays:                                                           */
  Py_DECREF(dx);
  return Py_BuildValue("N", z);
}


/* The module doc string                                                    */
PyDoc_STRVAR(spline__doc__, "Python wrapper for Cubic Spline Interpolation.");


/* A list of all the methods defined by this module.                        */
static PyMethodDef spline_methods[] = {
    {"splinterp",    splinterp,    METH_VARARGS, splinterp__doc__},
    {"splinterp_pt", splinterp_pt, METH_VARARGS, splinterp_pt__doc__},
    {"splinterp_2D", splinterp_2D, METH_VARARGS, splinterp_2D__doc__},
    {"spline_init",  spline_init,  METH_VARARGS, spline_init__doc__},
    {NULL,           NULL,         0,            NULL}    /* sentinel       */
};


#if PY_MAJOR_VERSION >= 3
/* Module definition for Python 3.                                          */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spline",
    spline__doc__,
    -1,
    spline_methods
};

/* When Python 3 imports a C module named 'X' it loads the module           */
/* then looks for a method named "PyInit_"+X and calls it.                  */
PyObject *PyInit_splinterp (void) {
  PyObject *module = PyModule_Create(&moduledef);
  import_array();
  return module;
}

#else
/* When Python 2 imports a C module named 'X' it loads the module           */
/* then looks for a method named "init"+X and calls it.                     */
void initspline(void){
  Py_InitModule3("spline", spline_methods, spline__doc__);
  import_array();
}
#endif
