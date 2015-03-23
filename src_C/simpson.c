#include <Python.h>
#include <numpy/arrayobject.h>

/* Access to i-th value of array a:                                         */
#define INDd(a,i) *((double *)(a->data + i*a->strides[0]))

#include "simpson.h"

PyDoc_STRVAR(geth__doc__,
"Calculate the differentials for a Simpson-rule integration.\n\
                                                            \n\
Parameters:                                                 \n\
-----------                                                 \n\
h: 1D double ndarray                                        \n\
   Intervals between the X-axis samples.                    \n\
                                                            \n\
Returns:                                                    \n\
--------                                                    \n\
hsum: 1D double ndarray                                     \n\
   Sums of interval pairs.                                  \n\
hratio: 1D double ndarray                                   \n\
   Ratio of consecutive intervals.                          \n\
hfactor: 1D double ndarray                                  \n\
   Factor interval.                                         \n\
                                                            \n\
Notes:                                                      \n\
------                                                      \n\
- If there are even samples, skip the first interval.       \n\
- hsum    = [h0+h1, h2+h3, h4+h5, ...]                      \n\
- hratio  = [h1/h0, h3/h2, h5/h4, ...]                      \n\
- hfactor = [hsum0*hsum0/h0*h1, hsum1*hsum1/h2*h3, ...]     \n\
");

static PyObject *geth(PyObject *self, PyObject *args){
  PyArrayObject *h, *hsum, *hratio, *hfactor;
  npy_intp size[1];    /* Size of output numpy array                        */
  int i, j, even=0, n; /* Auxilliary for-loop indices                       */

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "O", &h))
    return NULL;

  /* Get the number of intervals:                                           */
  n = h->dimensions[0];

  /* Check for even number of samples (odd number of intervals):            */
  even = n%2;

  /* Allocate outputs:                                                      */ 
  size[0] = n/2;
  hsum    = (PyArrayObject *) PyArray_SimpleNew(1, size, PyArray_DOUBLE);
  hratio  = (PyArrayObject *) PyArray_SimpleNew(1, size, PyArray_DOUBLE);
  hfactor = (PyArrayObject *) PyArray_SimpleNew(1, size, PyArray_DOUBLE);

  for (i=0; i<size[0]; i++){
    j = 2*i + even;
    INDd(hsum,   i) = INDd(h,(j  )) + INDd(h,(j+1));
    INDd(hratio, i) = INDd(h,(j+1)) / INDd(h,(j  ));
    INDd(hfactor,i) = INDd(hsum,i)*INDd(hsum,i) / (INDd(h,(j))*INDd(h,(j+1)));
  }

  Py_XDECREF(size);
  return Py_BuildValue("[N,N,N]", hsum, hratio, hfactor);
}


PyDoc_STRVAR(simps__doc__,
"Wrapper for Simpson-rule integration.\n\
                                                              \n\
Parameters:                                                   \n\
-----------                                                   \n\
y: 1D double ndarray                                          \n\
   Function to integrate.                                     \n\
h: 1D double ndarray                                          \n\
   Intervals between function evaluations.                    \n\
hsum: 1D double ndarray                                       \n\
   Sums of interval pairs.                                    \n\
hratio: 1D double ndarray                                     \n\
   Ratio of consecutive intervals.                            \n\
hfactor: 1D double ndarray                                    \n\
   Factor interval.                                           \n\
                                                              \n\
Returns:                                                      \n\
--------                                                      \n\
integ: Float                                                  \n\
   Integral of y over intervals h using the Simpson rule.     \n\
                                                              \n\
Notes:                                                        \n\
------                                                        \n\
- If there are even samples, use a trapezoidal integration for\n\
  the first interval.                                         \n\
- See geth for formula for hsum, hratio, and hfactor");


static PyObject *simps(PyObject *self, PyObject *args){
  PyArrayObject *y, *h, *hsum, *hratio, *hfactor;
  int n, even;
  double integ=0;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOOO", &y, &h, &hsum, &hratio, &hfactor))
    return NULL;

  /* Length of integrand:                                                   */
  n = y->dimensions[0];
  /* Check if I have an even number of samples:                             */
  even = n%2 == 0;

  /* Simple case, nothing to integrate:                                     */ 
  if (n == 1)
    return Py_BuildValue("d", 0.0);
  /* Simple case, do a trapezoidal integration:                             */
  if (n == 2)
    return Py_BuildValue("d", INDd(h,0) * (INDd(y,0) + INDd(y,1))/2);

  /* Do Simpson integration (skip first if even):                           */
  integ = simpson(y, hsum, hratio, hfactor, (n-1)/2);

  /* Add trapezoidal rule for first interval if n is even:                  */
  if (even){
    integ += INDd(h,0) * (INDd(y,0) + INDd(y,1))/2;
  }

  return Py_BuildValue("d", integ);
}


/* The module doc string    */
PyDoc_STRVAR(simpson__doc__, "Python wrapper for Simpson integration.");


/* A list of all the methods defined by this module.                        */
static PyMethodDef simpson_methods[] = {
    {"geth",      geth,       METH_VARARGS, geth__doc__},
    {"simps",     simps,      METH_VARARGS, simps__doc__},
    {NULL,        NULL,       0,            NULL}    /* sentinel */
};


/* When Python imports a C module named 'X' it loads the module */
/* then looks for a method named "init"+X and calls it.         */
void initsimpson(void){
  Py_InitModule3("simpson", simpson_methods, simpson__doc__);
  import_array();
}
