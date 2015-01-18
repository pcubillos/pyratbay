#include <Python.h>
#include <numpy/arrayobject.h>
#include "voigt.h"

/* Access to i-th value of array a:                                         */
#define INDd(a,i) *((double *)(a->data + i*a->strides[0]))
#define IND2i(a,i,j) *((int *)(a->data+i*a->strides[0]+j*a->strides[1]))
#define IND2d(a,i,j) *((double *)(a->data+i*a->strides[0]+j*a->strides[1]))


PyDoc_STRVAR(voigt__doc__,
"Calculate an oversampled Voigt profile for the given Lorentz and Doppler  \n\
widths.                                                                    \n\
                                                                           \n\
Parameters:                                                                \n\
-----------                                                                \n\
profile: 2D double ndarray                                                 \n\
   Array (osamp, nspec) where to put the calculated Voigt profiles.        \n\
lorentz: 1D double ndarray                                                 \n\
   Array of Lorentz widths in cm-1.                                        \n\
doppler: 1D double ndarray                                                 \n\
   Array of Doppler widths in cm-1.                                        \n\
psize: 2D integer ndarray                                                  \n\
   Array (nLorentz, nDoppler) with the half-size of spectral points of the \n\
   profiles.  Profiles not to be calculated have psize == 0.               \n\
osamp: Integer                                                             \n\
   Oversampling factor to calculate the Voigt profiles.                    \n\
dwn: Float                                                                 \n\
   Wavenumber step size in cm-1.                                           \n\
verb: Integer                                                              \n\
   Verbosity flag to print info to screen.                                 \n\
                                                                           \n\
Notes:                                                                     \n\
------                                                                     \n\
The profiles are stacked one next to the other in the profile array.       \n\
                                                                           \n\
Developers:                                                                \n\
-----------                                                                \n\
Patricio Rojo      U Cornell  pato@astro.cornell.edu (pato@oan.cl)         \n\
Patricio Cubillos  UCF        pcubillos@fulbrightmail.org                  \n\
                                                                           \n\
Modification History:                                                      \n\
---------------------                                                      \n\
2006        p. rojo      Writen as a component of the transit package.     \n\
2014-08-24  p. cubillos  Modified for use with the pyrat project.");

static PyObject *voigt(PyObject *self, PyObject *args){
  PyArrayObject *profile, *doppler, *lorentz, *psize;
  float **vprofile; /* Voigt profile for each (Dop,Lor) width               */
  double dwn;       /* Wavenumber sample step size                          */
  int nDop, nLor,  /* Number of Lorentz and Doppler width samples           */
      osamp,   /* Oversampling factor                                       */
      nwave,   /* Number of wavenumber samples of Voigt profile             */
      index=0, /* Profile index position                                    */
      verb,    /* Verbosity flag                                            */
      n, m, i, j; /* Auxilliary for-loop indices                            */

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOOidi", &profile, &lorentz, &doppler, &psize,
                                         &osamp, &dwn, &verb))
    return NULL;

  /* Get array sizes:                                                       */
  nLor = lorentz->dimensions[0];
  nDop = doppler->dimensions[0];

  /* Allocate Voigt profile per (Lorentz, Doppler) pair:                    */
  vprofile = (float **)calloc(osamp, sizeof(float *));

  for (m=0; m<nLor; m++){
    for   (n=0; n<nDop; n++){
      if (verb > 20)
        printf("DL Ratio: %.3f\n", INDd(doppler,n)/INDd(lorentz,m));
      /* If the profile size is > 0, calculate it:                          */
      if (IND2i(psize, m, n) != 0){
        /* Number of spectral samples:                                      */
        nwave = 2*IND2i(psize, m, n) + 1;
        /* Allocate profile array:                                          */
        *vprofile = (float *)calloc(nwave*osamp, sizeof(float));
        for(i=0; i<osamp; i++)
          vprofile[i] = vprofile[0] + i*nwave;

        if (verb > 10)
          printf("Calculating profile[%d, %d] = %d\n", m, n, IND2i(psize,m,n));
        /* Calculate voigt using a width that gives an integer number
           of 'dwn' spaced bins:                                            */
        j = voigtn(osamp, nwave, dwn*(long)(nwave/2),
                   INDd(lorentz,m), INDd(doppler,n), vprofile, -1,
                   nwave > _voigt_maxelements?VOIGT_QUICK:0);
        if (j != 1){
          printf("voigtn() returned error code %i.\n", j);
          return 0;
        }
        /* Store values in python-object profile:                           */
        for (i=0; i<osamp; i++){
          for (j=0; j<nwave; j++){
            //printf("%.5f,  ", vprofile[i][j]);
            IND2d(profile, i, (index+j)) = vprofile[i][j];
          }
          //printf("\n");
        }
        /* Free memory:                                                     */
        free(vprofile[0]);

        /* Update index of profile                                          */
        index += 2*IND2i(psize,m,n) + 1;
      }
      else{
        if (verb > 10)
          printf("Skip profile[%d, %d] calculation.\n", m, n);
      }
    }
    if (verb > 5 && verb <=10)
      printf("Calculated Voigt for Lorentz width %3d/%d.\n", m+1, nLor);
  }
  /* Free the allocated Voigt profile array and return:                     */
  free(vprofile);
  return Py_BuildValue("i", 1);
}


/* The module doc string    */
PyDoc_STRVAR(vprofile__doc__, "Python wrapper for Voigt profile calculation.");

/* A list of all the methods defined by this module.                        */
static PyMethodDef vprofile_methods[] = {
    {"voigt",      voigt,      METH_VARARGS, voigt__doc__},
    {NULL,         NULL,       0,            NULL}    /* sentinel */
};


/* When Python imports a C module named 'X' it loads the module */
/* then looks for a method named "init"+X and calls it.         */
void initvprofile(void){
  Py_InitModule3("vprofile", vprofile_methods, vprofile__doc__);
  import_array();
}
