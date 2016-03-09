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
  nin [0] = (int)PyArray_DIM(xin,  0)-1;
  nout[0] = (int)PyArray_DIM(xout, 0);

  /* Allocate arrays:                                                       */ 
  dx   = (PyArrayObject *) PyArray_SimpleNew(1, nin,  NPY_DOUBLE);
  nin[0] += 1;
  z    = (PyArrayObject *) PyArray_SimpleNew(1, nin,  NPY_DOUBLE);
  yout = (PyArrayObject *) PyArray_SimpleNew(1, nout, NPY_DOUBLE);

  for (i=0; i<nin[0]-1; i++){
    INDd(dx, i) = INDd(xin,(i+1)) - INDd(xin,i);
  }

  /* Calculate second derivatives:                                          */
  tri(z, yin, dx, nin[0]);

  /* Extrapolated values:                                                   */
  while(INDd(xout, epleft)  < INDd(xin,0)){
    INDd(yout, epleft) = extrap;
    epleft += 1;
  }
  epright = nout[0]-1;
  while(INDd(xout, epright) > INDd(xin,(nin[0]-1))){
    INDd(yout,epright) = extrap;
    epright -= 1;
  }

  /* Number of points to interpolate:                                       */
  nout[0] = epright - epleft + 1;
  /* Interpolate:                                                           */
  spline3(xin, yin, nin[0], xout, yout, nout[0], z, dx, epleft);

  free(dx);
  free(z);

  return Py_BuildValue("N", yout);
}


/* The module doc string                                                    */
PyDoc_STRVAR(spline__doc__, "Python wrapper for Cubic Spline Interpolation.");


/* A list of all the methods defined by this module.                        */
static PyMethodDef spline_methods[] = {
    {"splinterp", splinterp,  METH_VARARGS, splinterp__doc__},
    {NULL,        NULL,       0,            NULL}    /* sentinel            */
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
