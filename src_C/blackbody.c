#include <Python.h>
#include <numpy/arrayobject.h>

/* Access to i-th value of array a:                                         */
#define INDd(a,i) *((double *)(a->data + i*a->strides[0]))
#define INDi(a,i) *((int    *)(a->data + i*a->strides[0]))
#define IND2d(a,i,j) *((double *)(a->data + i*a->strides[0] + j*a->strides[1]))

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
   Array to store the Planck emission.                      \n\
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
  int i, j,   /* Auxilliary for-loop indices                                */
      nwave;  /* Number of wavenumber spectral samples                      */
  double factor;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOO", &B, &wn, &temp, &last))
    return NULL;

  /* Get the spectrum size:                                                 */
  nwave = wn->dimensions[0];

  /* Evaluate the Planck function:                                          */
  for (i=0; i<nwave; i++){
    factor = 2 * H * LS*LS * pow(INDd(wn,i),3);
    for (j=0; j < INDi(last, i); j++){
      IND2d(B,i,j) = factor / (exp(H*LS*INDd(wn,i)/(KB*INDd(temp,j))) - 1.0);
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


/* When Python imports a C module named 'X' it loads the module */
/* then looks for a method named "init"+X and calls it.         */
void initblackbody(void){
  Py_InitModule3("blackbody", blackbody_methods, blackbody__doc__);
  import_array();
}
