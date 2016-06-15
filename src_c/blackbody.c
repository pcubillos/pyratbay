#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_8_API_VERSION
#include <numpy/arrayobject.h>

#include "ind.h"
#include "constants.h"
#include "utils.h"

PyDoc_STRVAR(planck__doc__,
"Calculate the Planck emission function in wavenumber space:\n\
   Bnu(T) = 2 h c**2 nu**3 / (exp(hc*nu/kT)-1),             \n\
all variables in cgs units.                                 \n\
                                                            \n\
Parameters:                                                 \n\
-----------                                                 \n\
B: 2D float ndarray                                         \n\
   Array to store the Planck emission [nlayers, nwave].     \n\
wn: 1D float ndarray                                        \n\
   Wavenumber spectrum (cm-1).                              \n\
temp: 1D float ndarray                                      \n\
   Temperature at each layer (K).                           \n\
last: 1D integer ndarray                                    \n\
   Indices of last layer to evaluate.");

static PyObject *planck(PyObject *self, PyObject *args){
  PyArrayObject *B,   /* Extinction coefficient [nwave, nlayers]            */
           *wn,       /* Wavenumber array [nwave]                           */
           *temp,     /* Temperature of layers [nlayers]                    */
           *last;     /* Index of last layer to evaluate [nwave]            */
  int i, j;    /* Auxilliary for-loop indices                               */
  long nwave;  /* Number of wavenumber spectral samples                     */
  double factor;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOO", &B, &wn, &temp, &last))
    return NULL;

  /* Get the spectrum size:                                                 */
  nwave = PyArray_DIM(wn, 0);

  /* Evaluate the Planck function:                                          */
  for (i=0; i<nwave; i++){
    factor = 2 * H * LS*LS * pow(INDd(wn,i),3);
    for (j=0; j < INDi(last, i); j++){
      IND2d(B,j,i) = factor / (exp(H*LS*INDd(wn,i)/(KB*INDd(temp,j))) - 1.0);
    }
  }

  return Py_BuildValue("i", 1);
}

/* The module doc string    */
PyDoc_STRVAR(blackbody__doc__, "Wrapper for the Planck emission calculation.");

/* A list of all the methods defined by this module.                        */
static PyMethodDef blackbody_methods[] = {
    {"planck", planck, METH_VARARGS, planck__doc__},
    {NULL,         NULL,       0,            NULL}              /* sentinel */
};


#if PY_MAJOR_VERSION >= 3
/* Module definition for Python 3.                                          */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "blackbody",
    blackbody__doc__,
    -1,
    blackbody_methods
};

/* When Python 3 imports a C module named 'X' it loads the module           */
/* then looks for a method named "PyInit_"+X and calls it.                  */
PyObject *PyInit_blackbody (void) {
  PyObject *module = PyModule_Create(&moduledef);
  import_array();
  return module;
}

#else
/* When Python 2 imports a C module named 'X' it loads the module           */
/* then looks for a method named "init"+X and calls it.                     */
void initblackbody(void){
  Py_InitModule3("blackbody", blackbody_methods, blackbody__doc__);
  import_array();
}
#endif
